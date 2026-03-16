"""Command-line interface for PeptoMatch."""

import argparse
import json
import logging
import sys
from pathlib import Path

import pandas as pd

from .utils import load_config, setup_logging, ensure_output_dir
from .io_loaders import load_composition_data, load_strain_table
from .composition_features import CompositionFeatureExtractor
from .genome_prior import GenomePriorBuilder
from .scoring import PeptoneRecommender
from .explain import RecommendationExplainer
from .taxonomy_priors import list_supported_genera
from .strain_db import StrainDB
from .ncbi_client import NCBIClient


def _get_config() -> dict:
    try:
        return load_config()
    except Exception:
        return load_config(Path(__file__).parent.parent.parent / "config" / "config.yaml")


def cmd_recommend(args):
    config = _get_config()
    comp_path = Path(config["data"]["composition_file"])
    comp_df = load_composition_data(comp_path, sheet_name=config["data"].get("composition_sheet", "data"))

    # Load strains from DB or Excel
    db = StrainDB(Path("data/strains.db"))
    if db.count() == 0:
        strain_path = Path(config["data"]["strain_file"])
        db.load_from_excel(strain_path)

    strain_df = db.get_strain_df()

    # Find strain
    if args.strain.isdigit():
        strain_id = int(args.strain)
    else:
        mask = (
            strain_df["full_name"].str.contains(args.strain, case=False, na=False) |
            strain_df["species"].str.contains(args.strain, case=False, na=False)
        )
        matches = strain_df[mask]
        if len(matches) == 0:
            print(f"\nNo strain found matching '{args.strain}'")
            print("\nAvailable strains (first 10):")
            for _, row in strain_df.head(10).iterrows():
                print(f"  {row['strain_id']}: {row['full_name']}")
            return 1
        strain_id = matches.iloc[0]["strain_id"]
        if len(matches) > 1:
            print(f"Multiple matches, using: {matches.iloc[0]['full_name']}")

    recommender = PeptoneRecommender(comp_df, strain_df, config)
    pf = None if args.all_peptones else config.get("peptone_filter")
    recommendations = recommender.recommend(strain_id, top_k=args.topk, peptone_filter=pf)

    explainer = RecommendationExplainer(comp_df, strain_df, config, language=args.language)
    recommendations = explainer.explain_batch(strain_id, recommendations, top_n_reasons=3)

    if args.format == "table":
        strain_info = strain_df[strain_df["strain_id"] == strain_id].iloc[0]
        print(f"\n{'='*60}")
        print(f"Peptone Recommendations for: {strain_info['full_name']}")
        print(f"{'='*60}")
        for _, row in recommendations.iterrows():
            print(f"\n{row['rank']}. {row['peptone']} (Score: {row['score']:.1f})")
            print(f"   {row['explanation']}")
        print(f"\n{'='*60}")
    elif args.format == "json":
        print(json.dumps(recommendations.to_dict(orient="records"), indent=2, ensure_ascii=False))
    elif args.format == "csv":
        out = ensure_output_dir() / f"recommendations_strain{strain_id}.csv"
        recommendations.to_csv(out, index=False, encoding="utf-8-sig")
        print(f"Saved to {out}")

    return 0


def cmd_list_strains(args):
    config = _get_config()
    db = StrainDB(Path("data/strains.db"))
    if db.count() == 0:
        strain_path = Path(config["data"]["strain_file"])
        db.load_from_excel(strain_path)

    strains = db.get_all()
    print(f"\n{'='*85}")
    print(f"{'ID':<5} {'Genus':<25} {'Species':<20} {'Strain':<15} {'GCF':<20}")
    print(f"{'-'*85}")
    for s in strains:
        gcf = s.get("gcf_accession") or "N/A"
        print(f"{s['strain_id']:<5} {(s.get('genus') or '')[:24]:<25} "
              f"{(s.get('species') or '')[:19]:<20} "
              f"{(s.get('strain_name') or '')[:14]:<15} {gcf[:19]:<20}")
    print(f"{'='*85}")
    print(f"Total: {len(strains)} strains")
    return 0


def cmd_list_peptones(args):
    config = _get_config()
    comp_path = Path(config["data"]["composition_file"])
    comp_df = load_composition_data(comp_path)
    extractor = CompositionFeatureExtractor(comp_df, config)
    features = extractor.compute_all_features()

    print(f"\n{'='*60}")
    print(f"{'Name':<15} {'FAA Total':<12} {'TAA Total':<12} {'MW Avg':<10}")
    print(f"{'-'*60}")
    for name in features.index:
        faa = features.loc[name].get("faa_total", 0)
        taa = features.loc[name].get("taa_total", 0)
        mw = features.loc[name].get("mw_avg", 0)
        print(f"{name:<15} {faa:<12.2f} {taa:<12.2f} {mw:<10.1f}")
    print(f"{'='*60}")
    print(f"Total: {len(features)} peptones")
    return 0


def cmd_search_ncbi(args):
    config = _get_config()
    ncbi_cfg = config.get("ncbi", {})
    client = NCBIClient(email=ncbi_cfg.get("email", "user@example.com"))

    print(f"\nSearching NCBI for '{args.taxon}'...")
    results = client.search_assemblies(
        taxon=args.taxon,
        assembly_level=args.level,
        limit=args.limit,
    )

    if not results:
        print("No assemblies found.")
        return 1

    print(f"\n{'='*90}")
    print(f"{'Accession':<18} {'Organism':<35} {'Level':<18} {'Size (Mb)':<10}")
    print(f"{'-'*90}")
    for r in results:
        size_mb = r.get("genome_size", 0) / 1_000_000
        print(f"{r['accession']:<18} {r['organism_name'][:34]:<35} "
              f"{r['assembly_level']:<18} {size_mb:<10.2f}")
    print(f"{'='*90}")

    if args.add:
        db = StrainDB(Path("data/strains.db"))
        if db.count() == 0:
            strain_path = Path(config["data"]["strain_file"])
            db.load_from_excel(strain_path)

        added = db.expand_from_ncbi(args.taxon, client,
                                     assembly_level=args.level,
                                     max_results=args.limit)
        print(f"\nAdded {len(added)} new strains to local DB")

    return 0


def cmd_build_priors(args):
    config = _get_config()
    comp_path = Path(config["data"]["composition_file"])

    db = StrainDB(Path("data/strains.db"))
    if db.count() == 0:
        strain_path = Path(config["data"]["strain_file"])
        db.load_from_excel(strain_path)

    strain_df = db.get_strain_df()
    builder = GenomePriorBuilder(strain_df, config)
    priors = builder.build_all_priors(use_gcf=True, force_rebuild=args.force)

    print(f"\nBuilt priors for {len(priors)} strains")
    sources = {}
    for p in priors.values():
        src = p.get("source", "unknown")
        sources[src] = sources.get(src, 0) + 1
    for src, cnt in sorted(sources.items()):
        print(f"  {src}: {cnt}")
    return 0


def main():
    parser = argparse.ArgumentParser(
        description="PeptoMatch - Genome-driven Peptone Recommendation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  peptomatch recommend --strain "plantarum" --topk 5
  peptomatch list-strains
  peptomatch list-peptones
  peptomatch search-ncbi --taxon "Lactobacillaceae" --limit 10 --add
  peptomatch build-priors
        """,
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    sub = parser.add_subparsers(dest="command")

    # recommend
    p = sub.add_parser("recommend", help="Get peptone recommendations")
    p.add_argument("--strain", "-s", required=True)
    p.add_argument("--topk", "-k", type=int, default=10)
    p.add_argument("--format", "-f", choices=["table", "json", "csv"], default="table")
    p.add_argument("--language", "-l", choices=["ko", "en"], default="ko")
    p.add_argument("--all-peptones", action="store_true", help="Include all peptones (ignore Sempio filter)")

    # list-strains
    sub.add_parser("list-strains", help="List available strains")

    # list-peptones
    sub.add_parser("list-peptones", help="List available peptones")

    # search-ncbi
    p = sub.add_parser("search-ncbi", help="Search NCBI for genomes")
    p.add_argument("--taxon", "-t", required=True, help="Taxon name or ID")
    p.add_argument("--level", default="complete_genome",
                   choices=["complete_genome", "chromosome", "scaffold", "contig"])
    p.add_argument("--limit", type=int, default=20)
    p.add_argument("--add", action="store_true", help="Add results to local strain DB")

    # build-priors
    p = sub.add_parser("build-priors", help="Build genome priors")
    p.add_argument("--force", action="store_true")

    args = parser.parse_args()
    setup_logging("DEBUG" if args.verbose else "INFO")

    if args.command is None:
        parser.print_help()
        return 0

    handlers = {
        "recommend": cmd_recommend,
        "list-strains": cmd_list_strains,
        "list-peptones": cmd_list_peptones,
        "search-ncbi": cmd_search_ncbi,
        "build-priors": cmd_build_priors,
    }

    handler = handlers.get(args.command)
    if handler:
        try:
            return handler(args)
        except Exception as e:
            print(f"Error: {e}")
            if args.verbose:
                import traceback
                traceback.print_exc()
            return 1
    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    sys.exit(main())

"""Genome-based prior features for bacterial strains.

Uses NCBIClient + KOAnnotator instead of eggNOG-mapper.
Falls back to expanded taxonomy priors.
"""

import json
import logging
from pathlib import Path
from typing import Any, Optional

import pandas as pd

from .utils import ensure_output_dir
from .kegg_pathway import analyze_ko_annotations
from .taxonomy_priors import get_taxonomy_prior
from .ncbi_client import NCBIClient
from .ko_annotator import KOAnnotator
from .kegg_client import KEGGClient

logger = logging.getLogger("peptomatch")


class GenomePriorBuilder:
    """Build genome-based prior features for strains."""

    def __init__(self, strain_df: pd.DataFrame, config: Optional[dict] = None):
        self.strain_df = strain_df
        self.config = config or {}
        self.priors: dict[int, dict] = {}

        self.cache_dir = ensure_output_dir("genome_cache")
        self.prior_file = ensure_output_dir() / "genome_prior_features.json"

        # Initialize NCBI client and KO annotator
        ncbi_cfg = self.config.get("ncbi", {})
        self.ncbi_client = NCBIClient(
            email=ncbi_cfg.get("email", "user@example.com"),
            api_key=ncbi_cfg.get("api_key"),
            cache_dir=Path(ncbi_cfg.get("cache_dir", "data/ncbi_cache")),
        )

        ann_cfg = self.config.get("annotation", {})
        self.ko_annotator = KOAnnotator(
            kofamscan_path=ann_cfg.get("kofamscan_path"),
            kofamscan_profiles=ann_cfg.get("kofamscan_profiles"),
        )
        self.annotation_strategy = ann_cfg.get("strategy", "auto")

        # KEGG REST API client (works on Streamlit Cloud)
        self.kegg_client = KEGGClient()

        self._load_cached_priors()

    def _load_cached_priors(self) -> None:
        if self.prior_file.exists():
            try:
                with open(self.prior_file, "r", encoding="utf-8") as f:
                    cached = json.load(f)
                    self.priors = {int(k): v for k, v in cached.items()}
                logger.info(f"Loaded {len(self.priors)} cached genome priors")
            except Exception as e:
                logger.warning(f"Could not load cached priors: {e}")

    def save_priors(self) -> Path:
        with open(self.prior_file, "w", encoding="utf-8") as f:
            json.dump(self.priors, f, indent=2, ensure_ascii=False)
        logger.info(f"Saved {len(self.priors)} genome priors to {self.prior_file}")
        return self.prior_file

    def get_prior(self, strain_id: int) -> dict[str, Any]:
        if strain_id in self.priors:
            return self.priors[strain_id]
        return self.build_prior(strain_id)

    def build_prior(self, strain_id: int, force_rebuild: bool = False) -> dict[str, Any]:
        if strain_id in self.priors and not force_rebuild:
            return self.priors[strain_id]

        strain_row = self.strain_df[self.strain_df["strain_id"] == strain_id]

        if strain_row.empty:
            logger.warning(f"Strain ID {strain_id} not found in strain table")
            prior = get_taxonomy_prior("")
            prior["source"] = "generic_default"
            prior["strain_id"] = strain_id
            self.priors[strain_id] = prior
            return prior

        strain_info = strain_row.iloc[0]
        genus = strain_info.get("genus", "")
        gcf = strain_info.get("GCF", "")

        logger.info(f"Building prior for strain {strain_id}: {genus} (GCF: {gcf or 'N/A'})")

        prior = None

        species = strain_info.get("species", "")

        # Try GCF-based prior first (NCBI download + KofamScan/GFF3)
        if gcf and pd.notna(gcf):
            prior = self._build_gcf_prior(gcf, strain_id, genus)

        # Try KEGG REST API (works on Streamlit Cloud)
        if prior is None and genus:
            prior = self._build_kegg_prior(genus, species, strain_id)

        # Fall back to taxonomy prior
        if prior is None:
            prior = get_taxonomy_prior(genus)
            prior["source"] = "taxonomy_fallback"
            prior["genus"] = genus

        prior["strain_id"] = strain_id
        self.priors[strain_id] = prior
        return prior

    def _build_gcf_prior(self, gcf: str, strain_id: int, genus: str = "") -> Optional[dict[str, Any]]:
        """Build prior from GCF accession using NCBI download + KO annotation."""
        # Check per-GCF cache
        cache_file = self.cache_dir / f"{gcf}_prior.json"
        if cache_file.exists():
            try:
                with open(cache_file, "r", encoding="utf-8") as f:
                    logger.info(f"Loaded cached prior for {gcf}")
                    return json.load(f)
            except Exception:
                pass

        try:
            # Download protein FASTA and/or GFF3
            protein_fasta = self.ncbi_client.download_protein_fasta(gcf)
            gff3_path = self.ncbi_client.download_gff3(gcf)

            if protein_fasta is None and gff3_path is None:
                logger.warning(f"Could not download any annotation for {gcf}")
                return None

            # Run KO annotation
            ko_list, source = self.ko_annotator.annotate(
                protein_fasta=protein_fasta,
                gff3_path=gff3_path,
                strategy=self.annotation_strategy,
            )

            if not ko_list:
                logger.warning(f"No KO annotations for {gcf}, falling back to taxonomy")
                return None

            # Calculate pathway completeness
            prior = analyze_ko_annotations(ko_list)
            prior["gcf"] = gcf
            prior["genus"] = genus
            prior["ko_count"] = len(ko_list)
            prior["annotation_source"] = source

            # Save to per-GCF cache
            with open(cache_file, "w", encoding="utf-8") as f:
                json.dump(prior, f, indent=2)

            logger.info(f"Built GCF prior for {gcf}: {len(ko_list)} KOs via {source}")
            return prior

        except Exception as e:
            logger.error(f"Failed to build GCF prior for {gcf}: {e}")
            return None

    def _build_kegg_prior(
        self, genus: str, species: str, strain_id: int
    ) -> Optional[dict[str, Any]]:
        """Build prior using KEGG REST API (online, no local tools needed)."""
        # Check cache first
        cache_file = self.cache_dir / f"kegg_{genus}_{species}_prior.json"
        if cache_file.exists():
            try:
                with open(cache_file, "r", encoding="utf-8") as f:
                    prior = json.load(f)
                    logger.info(f"Loaded cached KEGG prior for {genus} {species}")
                    return prior
            except Exception:
                pass

        try:
            ko_list, source, org_code = self.kegg_client.annotate_strain(genus, species)

            if not ko_list:
                logger.info(f"No KEGG data for {genus} {species}")
                return None

            prior = analyze_ko_annotations(ko_list)
            prior["source"] = "kegg_api"
            prior["genus"] = genus
            prior["species"] = species
            prior["kegg_org_code"] = org_code
            prior["ko_count"] = len(ko_list)

            # Save KO list for flowchart visualization
            ko_list_file = self.cache_dir / f"kegg_{org_code}_ko_list.json"
            with open(ko_list_file, "w") as f:
                json.dump(ko_list, f)

            # Save prior cache
            with open(cache_file, "w", encoding="utf-8") as f:
                json.dump(prior, f, indent=2)

            logger.info(
                f"Built KEGG prior for {genus} {species}: "
                f"{len(ko_list)} KOs via {org_code}"
            )
            return prior

        except Exception as e:
            logger.error(f"KEGG prior failed for {genus} {species}: {e}")
            return None

    def run_kegg_annotation(self, strain_id: int) -> dict[str, Any]:
        """Explicitly run KEGG API annotation for a strain (called from UI).

        Returns the updated prior dict with source info.
        """
        strain_row = self.strain_df[self.strain_df["strain_id"] == strain_id]
        if strain_row.empty:
            return {"error": "Strain not found"}

        info = strain_row.iloc[0]
        genus = info.get("genus", "")
        species = info.get("species", "")

        if not genus:
            return {"error": "No genus information"}

        # Force KEGG annotation
        prior = self._build_kegg_prior(genus, species, strain_id)

        if prior is None:
            return {
                "error": f"KEGG에 {genus} {species}에 해당하는 organism이 없습니다.",
                "suggestion": "NCBI GFF3 annotation 또는 taxonomy prior를 사용합니다.",
            }

        # Update cached priors
        prior["strain_id"] = strain_id
        self.priors[strain_id] = prior
        self.save_priors()

        return prior

    def build_all_priors(self, use_gcf: bool = True, force_rebuild: bool = False) -> dict[int, dict]:
        total = len(self.strain_df)
        gcf_success = 0
        taxonomy_fallback = 0

        for idx, row in self.strain_df.iterrows():
            strain_id = row["strain_id"]
            logger.info(f"Processing strain {strain_id} ({idx + 1}/{total})")

            if force_rebuild and strain_id in self.priors:
                del self.priors[strain_id]

            prior = self.build_prior(strain_id)

            if prior.get("source") == "gcf_annotation":
                gcf_success += 1
            else:
                taxonomy_fallback += 1

        logger.info(f"Built priors: {gcf_success} from GCF, {taxonomy_fallback} from taxonomy")
        self.save_priors()
        return self.priors

    def get_demand_scores(self, strain_id: int) -> dict[str, float]:
        """Convert prior to demand scores. Demand = 1 - biosynthesis_completeness."""
        prior = self.get_prior(strain_id)
        demand = {}

        aa_synth = prior.get("aa_biosynthesis", {})
        for aa, completeness in aa_synth.items():
            demand[f"demand_{aa}"] = 1.0 - completeness

        vit_synth = prior.get("vitamin_biosynthesis", {})
        for vit, completeness in vit_synth.items():
            demand[f"demand_vitamin_{vit}"] = 1.0 - completeness

        nuc_synth = prior.get("nucleotide_biosynthesis", 0.6)
        demand["demand_nucleotide"] = 1.0 - nuc_synth

        demand["transporter_bonus"] = prior.get("transporter_score", 0.5)
        return demand

    def get_prior_summary(self, strain_id: int) -> dict[str, Any]:
        prior = self.get_prior(strain_id)

        aa_synth = prior.get("aa_biosynthesis", {})
        sorted_aa = sorted(aa_synth.items(), key=lambda x: x[1])
        deficient_aa = [aa for aa, val in sorted_aa[:5] if val < 0.4]

        vit_synth = prior.get("vitamin_biosynthesis", {})
        deficient_vit = [vit for vit, val in vit_synth.items() if val < 0.4]

        return {
            "strain_id": strain_id,
            "source": prior.get("source", "unknown"),
            "genus": prior.get("genus", "unknown"),
            "gcf": prior.get("gcf", "N/A"),
            "ko_count": prior.get("ko_count", 0),
            "most_deficient_aa": deficient_aa,
            "deficient_vitamins": deficient_vit,
            "nucleotide_synthesis": prior.get("nucleotide_biosynthesis", 0.5),
            "peptide_utilization": prior.get("transporter_score", 0.5),
            "has_gcf_data": prior.get("source") == "gcf_annotation"
        }

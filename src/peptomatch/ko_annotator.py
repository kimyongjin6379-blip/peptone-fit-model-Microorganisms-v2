"""Lightweight KO annotation without eggNOG-mapper.

Strategies (in priority order):
1. KofamScan - HMMER-based, most accurate (~4GB profiles)
2. GFF3 product description matching - zero-install fallback
3. Taxonomy prior - when no genome annotation is available
"""

import logging
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Optional

logger = logging.getLogger("peptomatch")


# Build gene-name-to-KO reverse map from kegg_pathway.py definitions
def _build_gene_ko_map() -> dict[str, str]:
    """Build gene_name -> KO_ID map from kegg_pathway.py inline comments.

    Parses comments like:
        "K00765",  # hisG; ATP phosphoribosyltransferase
    to create: {"hisg": "K00765", "atp phosphoribosyltransferase": "K00765"}
    """
    from .kegg_pathway import (
        AA_BIOSYNTHESIS_KOS,
        VITAMIN_BIOSYNTHESIS_KOS,
        NUCLEOTIDE_BIOSYNTHESIS_KOS,
        TRANSPORTER_KOS,
    )

    gene_map: dict[str, str] = {}

    # Manually define gene name -> KO based on the comments in kegg_pathway.py
    _GENE_NAMES = {
        # Amino acid biosynthesis
        "K00765": ["hisG"], "K02501": ["hisI"], "K01814": ["hisA"],
        "K02500": ["hisH"], "K01693": ["hisB"], "K00013": ["hisD"],
        "K01089": ["hisC"],
        "K01754": ["ilvA", "threonine dehydratase"],
        "K01652": ["ilvB", "acetolactate synthase"],
        "K01653": ["ilvH"], "K00053": ["ilvC"], "K01687": ["ilvD"],
        "K00826": ["ilvE", "branched-chain amino acid aminotransferase"],
        "K01649": ["leuA", "2-isopropylmalate synthase"],
        "K01703": ["leuC"], "K01704": ["leuD"],
        "K00052": ["leuB", "3-isopropylmalate dehydrogenase"],
        "K00928": ["lysC", "aspartate kinase"],
        "K00133": ["asd", "aspartate-semialdehyde dehydrogenase"],
        "K01714": ["dapA", "dihydrodipicolinate synthase"],
        "K00215": ["dapB"], "K01439": ["dapE"], "K01778": ["dapF"],
        "K01586": ["lysA", "diaminopimelate decarboxylase"],
        "K00003": ["hom", "homoserine dehydrogenase"],
        "K00651": ["metA", "homoserine O-succinyltransferase"],
        "K01739": ["metB", "cystathionine gamma-synthase"],
        "K01760": ["metC", "cystathionine beta-lyase"],
        "K00548": ["metH", "methionine synthase"],
        "K00549": ["metE"],
        "K00800": ["aroA", "EPSP synthase"],
        "K01736": ["aroC", "chorismate synthase"],
        "K04518": ["pheA2", "prephenate dehydratase"],
        "K00832": ["tyrB", "aromatic-amino-acid transaminase"],
        "K01850": ["pheA", "chorismate mutase"],
        "K00872": ["thrB", "homoserine kinase"],
        "K01733": ["thrC", "threonine synthase"],
        "K01657": ["trpE", "anthranilate synthase"],
        "K01658": ["trpG"], "K00766": ["trpD"],
        "K01817": ["trpF"], "K01609": ["trpC"],
        "K01695": ["trpA", "tryptophan synthase alpha"],
        "K01696": ["trpB", "tryptophan synthase beta"],
        "K00814": ["GPT", "alanine transaminase"],
        "K00259": ["ald", "alanine dehydrogenase"],
        "K00611": ["argF", "ornithine carbamoyltransferase"],
        "K01940": ["argG", "argininosuccinate synthase"],
        "K01755": ["argH", "argininosuccinate lyase"],
        "K00930": ["argB"], "K00145": ["argC"], "K00821": ["argD"],
        "K01953": ["asnB", "asparagine synthase"],
        "K01914": ["asnA"],
        "K00813": ["aspC", "aspartate aminotransferase"],
        "K00812": ["aspB"],
        "K00640": ["cysE", "serine O-acetyltransferase"],
        "K01738": ["cysK", "cysteine synthase"],
        "K12339": ["cysM"],
        "K00262": ["gdhA", "glutamate dehydrogenase"],
        "K00265": ["gltB"], "K00266": ["gltD"],
        "K01915": ["glnA", "glutamine synthetase"],
        "K00600": ["glyA", "serine hydroxymethyltransferase"],
        "K00281": ["gcvP", "glycine dehydrogenase"],
        "K00931": ["proB", "glutamate 5-kinase"],
        "K00147": ["proA"], "K00286": ["proC"],
        "K00058": ["serA"], "K00831": ["serC"], "K01079": ["serB"],
        "K00220": ["tyrA", "prephenate dehydrogenase"],
        # Vitamins
        "K00941": ["thiD"], "K00788": ["thiE"], "K00946": ["thiL"],
        "K03147": ["thiC"],
        "K00794": ["ribH"], "K00793": ["ribE"], "K02858": ["ribD"],
        "K00082": ["ribA"],
        "K00767": ["nadC"], "K00969": ["nadD"],
        "K01916": ["nadE", "NAD+ synthase"], "K03517": ["nadB"],
        "K00077": ["panE"], "K00867": ["panK"], "K01918": ["panC"],
        "K00275": ["pdxH"], "K03472": ["pdxJ"], "K03473": ["pdxA"],
        "K00652": ["bioF"], "K00833": ["bioA"],
        "K01012": ["bioB", "biotin synthase"], "K00796": ["bioD", "folK"],
        "K00950": ["folC"], "K01930": ["folE"], "K00287": ["folA"],
        "K02232": ["cobA"], "K02224": ["cobI"], "K02229": ["cobO"],
        # Nucleotides
        "K00764": ["purF"], "K01945": ["purD"], "K01587": ["purN"],
        "K01952": ["purL"], "K01933": ["purM"], "K00602": ["purH"],
        "K01756": ["purB"],
        "K01955": ["carB"], "K00609": ["pyrB"], "K01465": ["pyrC"],
        "K00254": ["pyrD"], "K01591": ["pyrF"],
        # Transporters
        "K15580": ["oppA"], "K15581": ["oppB"], "K15582": ["oppC"],
        "K15583": ["oppD"], "K15584": ["oppF"],
        "K02035": ["dppA"], "K02036": ["dppB"], "K02037": ["dppC"],
        "K02038": ["dppD"], "K02039": ["dppF"],
        "K09969": ["aapJ"], "K09970": ["aapQ"], "K09971": ["aapM"],
        "K09972": ["aapP"],
        "K01999": ["livK"], "K01997": ["livH"], "K01998": ["livM"],
        "K01995": ["livG"], "K01996": ["livF"],
    }

    for ko_id, names in _GENE_NAMES.items():
        for name in names:
            gene_map[name.lower()] = ko_id

    return gene_map


# Module-level cache
_GENE_KO_MAP: Optional[dict[str, str]] = None


def _get_gene_ko_map() -> dict[str, str]:
    global _GENE_KO_MAP
    if _GENE_KO_MAP is None:
        _GENE_KO_MAP = _build_gene_ko_map()
    return _GENE_KO_MAP


class KOAnnotator:
    """Multi-strategy KO annotation."""

    def __init__(self, kofamscan_path: Optional[str] = None,
                 kofamscan_profiles: Optional[str] = None):
        # Priority: explicit arg > env var > PATH lookup
        self.kofamscan_path = (
            kofamscan_path
            or os.environ.get("KOFAMSCAN_PATH")
            or shutil.which("exec_annotation")
        )
        self.kofamscan_profiles = (
            kofamscan_profiles
            or os.environ.get("KOFAMSCAN_PROFILES")
        )
        # Verify binary actually exists
        self.has_kofamscan = (
            self.kofamscan_path is not None
            and Path(self.kofamscan_path).exists()
        )
        if self.has_kofamscan:
            logger.info(f"KofamScan found: {self.kofamscan_path}")

    def annotate(self, protein_fasta: Optional[Path] = None,
                 gff3_path: Optional[Path] = None,
                 strategy: str = "auto") -> tuple[list[str], str]:
        """Get KO IDs for a genome.

        Args:
            protein_fasta: Path to protein FASTA file
            gff3_path: Path to GFF3 annotation file
            strategy: 'auto', 'kofamscan', 'gff3'

        Returns:
            Tuple of (ko_list, source_label)
        """
        if strategy == "auto":
            # Try KofamScan first, then GFF3
            if self.has_kofamscan and protein_fasta:
                kos = self._annotate_kofamscan(protein_fasta)
                if kos:
                    return kos, "kofamscan"

            if gff3_path:
                kos = self._annotate_from_gff3(gff3_path)
                if kos:
                    return kos, "gff3_product"

        elif strategy == "kofamscan" and protein_fasta:
            kos = self._annotate_kofamscan(protein_fasta)
            if kos:
                return kos, "kofamscan"

        elif strategy == "gff3" and gff3_path:
            kos = self._annotate_from_gff3(gff3_path)
            if kos:
                return kos, "gff3_product"

        return [], "none"

    def _annotate_from_gff3(self, gff3_path: Path) -> list[str]:
        """Parse GFF3 product annotations and map to KOs via gene names."""
        gene_map = _get_gene_ko_map()
        found_kos: set[str] = set()

        try:
            with open(gff3_path, "r", encoding="utf-8", errors="replace") as f:
                for line in f:
                    if line.startswith("#"):
                        continue

                    parts = line.strip().split("\t")
                    if len(parts) < 9:
                        continue

                    feature_type = parts[2]
                    if feature_type not in ("CDS", "gene"):
                        continue

                    attrs = parts[8]

                    # Extract gene name from attributes
                    gene_match = re.search(r'gene=([^;]+)', attrs)
                    if gene_match:
                        gene_name = gene_match.group(1).strip().lower()
                        if gene_name in gene_map:
                            found_kos.add(gene_map[gene_name])

                    # Extract product description
                    product_match = re.search(r'product=([^;]+)', attrs)
                    if product_match:
                        product = product_match.group(1).strip().lower()
                        product = requests_unquote(product)
                        # Match against known gene/product names
                        for name, ko_id in gene_map.items():
                            if len(name) > 3 and name in product:
                                found_kos.add(ko_id)

            logger.info(f"GFF3 annotation matched {len(found_kos)} KOs from {gff3_path.name}")
            return list(found_kos)

        except Exception as e:
            logger.error(f"GFF3 parsing failed: {e}")
            return []

    def _annotate_kofamscan(self, protein_fasta: Path) -> list[str]:
        """Run KofamScan for KO annotation."""
        if not self.kofamscan_path:
            return []

        output_file = protein_fasta.parent / "kofamscan_result.txt"

        if output_file.exists():
            return self._parse_kofamscan_output(output_file)

        try:
            cmd = [self.kofamscan_path, "-o", str(output_file), "--cpu=2"]
            if self.kofamscan_profiles:
                cmd.extend(["-p", self.kofamscan_profiles])
            # Use config.yml if it exists alongside exec_annotation
            config_yml = Path(self.kofamscan_path).parent / "config.yml"
            if config_yml.exists():
                cmd.extend(["--config", str(config_yml)])
            cmd.append(str(protein_fasta))

            logger.info(f"Running KofamScan on {protein_fasta.name}...")
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)

            if result.returncode != 0:
                logger.error(f"KofamScan failed: {result.stderr[:500]}")
                return []

            return self._parse_kofamscan_output(output_file)

        except subprocess.TimeoutExpired:
            logger.error("KofamScan timed out")
            return []
        except FileNotFoundError:
            logger.error("KofamScan not found")
            return []

    def _parse_kofamscan_output(self, output_file: Path) -> list[str]:
        """Parse KofamScan output for significant KO hits."""
        kos: set[str] = set()

        try:
            with open(output_file, "r") as f:
                for line in f:
                    if line.startswith("#") or line.startswith("*"):
                        continue
                    # Format: * gene_name KO_ID ...  (* = significant hit)
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split()
                    if len(parts) >= 2:
                        ko = parts[1] if parts[0] == "*" else parts[0]
                        if ko.startswith("K") and len(ko) == 6:
                            kos.add(ko)

            logger.info(f"KofamScan found {len(kos)} KOs")
            return list(kos)

        except Exception as e:
            logger.error(f"Failed to parse KofamScan output: {e}")
            return []


def requests_unquote(s: str) -> str:
    """URL-decode a string (GFF3 uses URL encoding for special chars)."""
    try:
        from urllib.parse import unquote
        return unquote(s)
    except ImportError:
        return s

"""NCBI Datasets REST API v2 and Entrez client.

Handles assembly search, protein/GFF3 download, and taxonomy lookup
using only the `requests` library (no Biopython dependency).
"""

import json
import logging
import time
import zipfile
from pathlib import Path
from typing import Any, Optional
from xml.etree import ElementTree

import requests

logger = logging.getLogger("peptomatch")

DATASETS_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


class NCBIClient:
    """NCBI API client for genome assembly search and download."""

    def __init__(self, email: str = "user@example.com", api_key: Optional[str] = None,
                 cache_dir: Optional[Path] = None):
        self.email = email
        self.api_key = api_key
        self.rate_delay = 0.12 if api_key else 0.34
        self._last_request = 0.0
        self.cache_dir = cache_dir or Path("data/ncbi_cache")
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.session = requests.Session()
        self.session.headers.update({"Accept": "application/json"})

    def _rate_limit(self):
        elapsed = time.time() - self._last_request
        if elapsed < self.rate_delay:
            time.sleep(self.rate_delay - elapsed)
        self._last_request = time.time()

    # ── Assembly Search ──────────────────────────────────────────────────

    def search_assemblies(
        self,
        taxon: str,
        assembly_level: Optional[str] = None,
        refseq_only: bool = True,
        limit: int = 20,
    ) -> list[dict[str, Any]]:
        """Search NCBI Assembly DB by taxon name or ID.

        Args:
            taxon: Taxonomy name (e.g., 'Lactobacillaceae') or NCBI tax ID
            assembly_level: Filter: 'complete_genome', 'chromosome', 'scaffold', 'contig'
            refseq_only: Only return RefSeq assemblies
            limit: Max results

        Returns:
            List of assembly info dicts
        """
        self._rate_limit()
        url = f"{DATASETS_BASE}/genome/taxon/{requests.utils.quote(taxon)}/dataset_report"

        # Map user-friendly level names to API filter values
        level_map = {
            "complete_genome": "complete_genome",
            "chromosome": "chromosome",
            "scaffold": "scaffold",
            "contig": "contig",
        }

        params: dict[str, Any] = {
            "page_size": min(limit, 100),
        }
        if assembly_level:
            params["filters.assembly_level"] = level_map.get(assembly_level, assembly_level)
        if refseq_only:
            params["filters.assembly_source"] = "refseq"

        try:
            resp = self.session.get(url, params=params, timeout=30)
            resp.raise_for_status()
            data = resp.json()
        except requests.RequestException as e:
            logger.error(f"Assembly search failed for '{taxon}': {e}")
            return []

        reports = data.get("reports", [])
        results = []
        for r in reports[:limit]:
            info = r.get("assembly_info", {})
            org = r.get("organism", {})
            stats = r.get("assembly_stats", {})
            results.append({
                "accession": r.get("accession", ""),
                "organism_name": org.get("organism_name", ""),
                "tax_id": org.get("tax_id", ""),
                "assembly_level": info.get("assembly_level", ""),
                "assembly_name": info.get("assembly_name", ""),
                "refseq_category": info.get("refseq_category", ""),
                "genome_size": int(stats.get("total_sequence_length", 0) or 0),
                "gc_percent": stats.get("gc_percent", 0),
                "gene_count": stats.get("total_number_of_chromosomes", 0),
            })

        logger.info(f"Found {len(results)} assemblies for '{taxon}'")
        return results

    def get_assembly_info(self, accession: str) -> Optional[dict[str, Any]]:
        """Get details for a specific assembly accession."""
        self._rate_limit()
        url = f"{DATASETS_BASE}/genome/accession/{accession}/dataset_report"

        try:
            resp = self.session.get(url, timeout=30)
            resp.raise_for_status()
            data = resp.json()
            reports = data.get("reports", [])
            if reports:
                return reports[0]
        except requests.RequestException as e:
            logger.error(f"Failed to get assembly info for {accession}: {e}")

        return None

    # ── Downloads ────────────────────────────────────────────────────────

    def download_protein_fasta(self, accession: str) -> Optional[Path]:
        """Download protein FASTA for a genome accession.

        Returns:
            Path to protein.faa file, or None on failure
        """
        out_dir = self.cache_dir / accession
        out_dir.mkdir(parents=True, exist_ok=True)
        protein_file = out_dir / "protein.faa"

        if protein_file.exists() and protein_file.stat().st_size > 0:
            logger.info(f"Using cached proteins for {accession}")
            return protein_file

        self._rate_limit()
        url = f"{DATASETS_BASE}/genome/accession/{accession}/download"
        params = {"include_annotation_type": "PROT_FASTA"}

        zip_path = out_dir / "download.zip"

        try:
            logger.info(f"Downloading proteins for {accession}...")
            resp = self.session.get(url, params=params, timeout=300, stream=True)
            resp.raise_for_status()

            with open(zip_path, "wb") as f:
                for chunk in resp.iter_content(chunk_size=8192):
                    f.write(chunk)

            # Extract protein.faa from zip
            with zipfile.ZipFile(zip_path, "r") as zf:
                prot_files = [n for n in zf.namelist() if n.endswith("protein.faa")]
                if not prot_files:
                    logger.error(f"No protein.faa in archive for {accession}")
                    return None
                with zf.open(prot_files[0]) as src, open(protein_file, "wb") as dst:
                    dst.write(src.read())

            zip_path.unlink(missing_ok=True)
            logger.info(f"Downloaded {protein_file.stat().st_size / 1024:.1f} KB for {accession}")
            return protein_file

        except Exception as e:
            logger.error(f"Protein download failed for {accession}: {e}")
            zip_path.unlink(missing_ok=True)
            return None

    def download_gff3(self, accession: str) -> Optional[Path]:
        """Download GFF3 annotation for a genome accession.

        Returns:
            Path to genomic.gff file, or None on failure
        """
        out_dir = self.cache_dir / accession
        out_dir.mkdir(parents=True, exist_ok=True)
        gff_file = out_dir / "genomic.gff"

        if gff_file.exists() and gff_file.stat().st_size > 0:
            logger.info(f"Using cached GFF3 for {accession}")
            return gff_file

        self._rate_limit()
        url = f"{DATASETS_BASE}/genome/accession/{accession}/download"
        params = {"include_annotation_type": "GENOME_GFF"}

        zip_path = out_dir / "download_gff.zip"

        try:
            logger.info(f"Downloading GFF3 for {accession}...")
            resp = self.session.get(url, params=params, timeout=300, stream=True)
            resp.raise_for_status()

            with open(zip_path, "wb") as f:
                for chunk in resp.iter_content(chunk_size=8192):
                    f.write(chunk)

            with zipfile.ZipFile(zip_path, "r") as zf:
                gff_files = [n for n in zf.namelist() if n.endswith(".gff")]
                if not gff_files:
                    logger.error(f"No .gff in archive for {accession}")
                    return None
                with zf.open(gff_files[0]) as src, open(gff_file, "wb") as dst:
                    dst.write(src.read())

            zip_path.unlink(missing_ok=True)
            logger.info(f"Downloaded GFF3 for {accession}")
            return gff_file

        except Exception as e:
            logger.error(f"GFF3 download failed for {accession}: {e}")
            zip_path.unlink(missing_ok=True)
            return None

    # ── Taxonomy ─────────────────────────────────────────────────────────

    def search_taxonomy(self, query: str) -> Optional[str]:
        """Search NCBI Taxonomy for a tax ID."""
        self._rate_limit()
        params = {
            "db": "taxonomy",
            "term": query,
            "retmax": 1,
            "email": self.email,
            "retmode": "json",
        }
        if self.api_key:
            params["api_key"] = self.api_key

        try:
            resp = self.session.get(f"{ENTREZ_BASE}/esearch.fcgi", params=params, timeout=15)
            resp.raise_for_status()
            data = resp.json()
            ids = data.get("esearchresult", {}).get("idlist", [])
            return ids[0] if ids else None
        except Exception as e:
            logger.error(f"Taxonomy search failed for '{query}': {e}")
            return None

    def get_taxonomy_info(self, tax_id: str) -> Optional[dict[str, Any]]:
        """Get taxonomy details by ID."""
        self._rate_limit()
        params = {
            "db": "taxonomy",
            "id": tax_id,
            "retmode": "xml",
            "email": self.email,
        }
        if self.api_key:
            params["api_key"] = self.api_key

        try:
            resp = self.session.get(f"{ENTREZ_BASE}/efetch.fcgi", params=params, timeout=15)
            resp.raise_for_status()

            root = ElementTree.fromstring(resp.content)
            taxon = root.find(".//Taxon")
            if taxon is None:
                return None

            sci_name = taxon.findtext("ScientificName", "")
            parts = sci_name.split()

            lineage_items = taxon.findall(".//LineageEx/Taxon")
            lineage = [t.findtext("ScientificName", "") for t in lineage_items]

            return {
                "tax_id": tax_id,
                "scientific_name": sci_name,
                "rank": taxon.findtext("Rank", ""),
                "lineage": lineage,
                "genus": parts[0] if parts else "",
                "species": parts[1] if len(parts) > 1 else "",
            }

        except Exception as e:
            logger.error(f"Taxonomy fetch failed for {tax_id}: {e}")
            return None

    @staticmethod
    def resolve_genus(organism_name: str) -> str:
        """Extract genus from an organism name string."""
        if not organism_name:
            return ""
        return organism_name.strip().split()[0]

"""KEGG REST API client for online KO annotation.

Uses the free KEGG REST API (https://rest.kegg.jp/) to:
1. Search for organisms by species name
2. Retrieve KO (KEGG Orthology) assignments for all genes
3. Cache results locally for reuse

This works on Streamlit Cloud — no local tools required.
"""

import json
import logging
import time
from pathlib import Path
from typing import Optional

import requests

from .utils import ensure_output_dir

logger = logging.getLogger("peptomatch")

KEGG_REST_BASE = "https://rest.kegg.jp"


class KEGGClient:
    """KEGG REST API client for organism search and KO retrieval."""

    def __init__(self, cache_dir: Optional[Path] = None):
        self.cache_dir = cache_dir or ensure_output_dir("kegg_cache")
        self.session = requests.Session()
        self._last_request = 0.0
        # KEGG asks for max 3 requests/second for free use
        self.rate_delay = 0.35

    def _rate_limit(self):
        elapsed = time.time() - self._last_request
        if elapsed < self.rate_delay:
            time.sleep(self.rate_delay - elapsed)
        self._last_request = time.time()

    # ── Organism Search ──────────────────────────────────────────────

    def search_organism(self, query: str) -> list[dict[str, str]]:
        """Search KEGG organisms by name.

        Args:
            query: Species name (e.g., "Lactobacillus plantarum")

        Returns:
            List of dicts with 'org_code', 'tax_id', 'name'
        """
        self._rate_limit()
        url = f"{KEGG_REST_BASE}/find/genome/{requests.utils.quote(query)}"

        try:
            resp = self.session.get(url, timeout=15)
            resp.raise_for_status()
        except requests.RequestException as e:
            logger.error(f"KEGG organism search failed for '{query}': {e}")
            return []

        results = []
        for line in resp.text.strip().split("\n"):
            if not line.strip():
                continue
            # Format: genome:T00123\tOrganism name; description
            parts = line.split("\t", 1)
            if len(parts) < 2:
                continue
            genome_id = parts[0].replace("genome:", "").strip()
            name = parts[1].strip()
            results.append({
                "genome_id": genome_id,
                "name": name,
            })

        logger.info(f"KEGG search '{query}': {len(results)} results")
        return results

    def find_organism_code(self, genus: str, species: str = "") -> Optional[str]:
        """Find KEGG 3/4-letter organism code for a species.

        Args:
            genus: Genus name (e.g., "Lactobacillus")
            species: Species name (e.g., "plantarum")

        Returns:
            Organism code (e.g., "lpl") or None
        """
        query = f"{genus} {species}".strip()

        # Check cache
        cache_file = self.cache_dir / f"org_{genus}_{species}.json"
        if cache_file.exists():
            try:
                with open(cache_file, "r") as f:
                    data = json.load(f)
                    if data.get("org_code"):
                        return data["org_code"]
            except Exception:
                pass

        # Search KEGG for organism list
        self._rate_limit()
        url = f"{KEGG_REST_BASE}/list/organism"

        try:
            resp = self.session.get(url, timeout=30)
            resp.raise_for_status()
        except requests.RequestException as e:
            logger.error(f"KEGG organism list failed: {e}")
            return None

        # Parse: T00123\torg_code\tOrganism name\tLineage
        best_match = None
        genus_only_match = None
        genus_lower = genus.lower()
        species_lower = species.lower() if species else ""

        for line in resp.text.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            org_code = parts[1].strip()
            org_name = parts[2].strip().lower()

            # Exact genus + species match (highest priority)
            if species_lower and genus_lower in org_name and species_lower in org_name:
                best_match = org_code
                break
            # Genus-only match (keep first hit as fallback)
            elif genus_only_match is None and genus_lower in org_name:
                genus_only_match = org_code

        if not best_match:
            best_match = genus_only_match

        if best_match:
            # Cache it
            cache_file.parent.mkdir(parents=True, exist_ok=True)
            with open(cache_file, "w") as f:
                json.dump({"org_code": best_match, "query": query}, f)
            logger.info(f"Found KEGG organism code: {best_match} for {query}")

        return best_match

    # ── KO Retrieval ─────────────────────────────────────────────────

    def get_ko_list(self, org_code: str) -> list[str]:
        """Get all KO IDs assigned to genes in an organism.

        Args:
            org_code: KEGG organism code (e.g., "lpl", "eco")

        Returns:
            List of KO IDs (e.g., ["K00765", "K00013", ...])
        """
        # Check cache
        cache_file = self.cache_dir / f"{org_code}_ko_list.json"
        if cache_file.exists():
            try:
                with open(cache_file, "r") as f:
                    data = json.load(f)
                    logger.info(f"Loaded cached KO list for {org_code}: {len(data)} KOs")
                    return data
            except Exception:
                pass

        self._rate_limit()
        url = f"{KEGG_REST_BASE}/link/ko/{org_code}"

        try:
            resp = self.session.get(url, timeout=60)
            resp.raise_for_status()
        except requests.RequestException as e:
            logger.error(f"KEGG KO retrieval failed for {org_code}: {e}")
            return []

        ko_set: set[str] = set()
        for line in resp.text.strip().split("\n"):
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                ko_raw = parts[1].strip()
                # Format: ko:K00765
                ko_id = ko_raw.replace("ko:", "")
                if ko_id.startswith("K") and len(ko_id) == 6:
                    ko_set.add(ko_id)

        ko_list = sorted(ko_set)

        # Cache
        if ko_list:
            cache_file.parent.mkdir(parents=True, exist_ok=True)
            with open(cache_file, "w") as f:
                json.dump(ko_list, f)
            logger.info(f"Retrieved {len(ko_list)} KOs for {org_code} from KEGG")

        return ko_list

    def annotate_strain(
        self,
        genus: str,
        species: str = "",
    ) -> tuple[list[str], str, Optional[str]]:
        """Full pipeline: find organism → get KO list.

        Args:
            genus: Genus name
            species: Species name

        Returns:
            Tuple of (ko_list, source_label, org_code)
        """
        org_code = self.find_organism_code(genus, species)
        if not org_code:
            logger.warning(f"No KEGG organism found for {genus} {species}")
            return [], "not_found", None

        ko_list = self.get_ko_list(org_code)
        if not ko_list:
            return [], "no_kos", org_code

        return ko_list, "kegg_api", org_code

    def clear_cache(self, org_code: Optional[str] = None):
        """Clear cached data."""
        if org_code:
            cache_file = self.cache_dir / f"{org_code}_ko_list.json"
            cache_file.unlink(missing_ok=True)
        else:
            import shutil
            if self.cache_dir.exists():
                shutil.rmtree(self.cache_dir)
                self.cache_dir.mkdir(parents=True, exist_ok=True)

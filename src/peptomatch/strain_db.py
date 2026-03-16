"""SQLite-based strain database with NCBI expansion capability."""

import logging
import sqlite3
from pathlib import Path
from typing import Any, Optional

import pandas as pd

from .ncbi_client import NCBIClient

logger = logging.getLogger("peptomatch")

CREATE_TABLE_SQL = """
CREATE TABLE IF NOT EXISTS strains (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    strain_id INTEGER UNIQUE NOT NULL,
    genus TEXT,
    species TEXT,
    strain_name TEXT,
    gcf_accession TEXT,
    assembly_level TEXT,
    source TEXT DEFAULT 'manual',
    taxonomy_id TEXT,
    domain TEXT,
    notes TEXT,
    full_name TEXT
)
"""


class StrainDB:
    """Local strain database with NCBI expansion."""

    def __init__(self, db_path: Optional[Path] = None):
        self.db_path = db_path or Path("data/strains.db")
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self.conn = sqlite3.connect(str(self.db_path), check_same_thread=False)
        self.conn.row_factory = sqlite3.Row
        self.conn.execute(CREATE_TABLE_SQL)
        self.conn.commit()

    def close(self):
        self.conn.close()

    def _next_strain_id(self) -> int:
        cur = self.conn.execute("SELECT MAX(strain_id) FROM strains")
        row = cur.fetchone()
        return (row[0] or 0) + 1

    def add_strain(self, genus: str, species: str, strain_name: str = "",
                   gcf: str = "", assembly_level: str = "",
                   source: str = "manual", domain: str = "Bacteria",
                   taxonomy_id: str = "", notes: str = "") -> int:
        """Add a strain. Returns strain_id."""
        strain_id = self._next_strain_id()
        full_name = f"{genus} {species} {strain_name}".strip()

        try:
            self.conn.execute(
                """INSERT INTO strains
                   (strain_id, genus, species, strain_name, gcf_accession,
                    assembly_level, source, taxonomy_id, domain, notes, full_name)
                   VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                (strain_id, genus, species, strain_name, gcf or None,
                 assembly_level, source, taxonomy_id, domain, notes, full_name)
            )
            self.conn.commit()
            return strain_id
        except sqlite3.IntegrityError:
            logger.warning(f"Strain {full_name} already exists")
            return -1

    def load_from_excel(self, excel_path: Path) -> int:
        """Import strains from the existing strain list Excel file."""
        from .io_loaders import load_strain_table

        df = load_strain_table(excel_path)
        count = 0

        for _, row in df.iterrows():
            genus = row.get("genus", "") or ""
            species = row.get("species", "") or ""
            strain_name = row.get("strain_name", "") or ""
            gcf = row.get("GCF", "") or ""

            if not genus and not species:
                continue

            sid = self.add_strain(
                genus=genus, species=species, strain_name=strain_name,
                gcf=gcf, source="excel_import"
            )
            if sid > 0:
                count += 1

        logger.info(f"Imported {count} strains from {excel_path.name}")
        return count

    def search_local(self, query: str) -> list[dict]:
        """Search local strains by name/genus/species."""
        cur = self.conn.execute(
            """SELECT * FROM strains
               WHERE genus LIKE ? OR species LIKE ? OR strain_name LIKE ?
                  OR full_name LIKE ? OR gcf_accession LIKE ?
               ORDER BY strain_id""",
            tuple(f"%{query}%" for _ in range(5))
        )
        return [dict(r) for r in cur.fetchall()]

    def get_all(self) -> list[dict]:
        cur = self.conn.execute("SELECT * FROM strains ORDER BY strain_id")
        return [dict(r) for r in cur.fetchall()]

    def get_strain(self, strain_id: int) -> Optional[dict]:
        cur = self.conn.execute("SELECT * FROM strains WHERE strain_id = ?", (strain_id,))
        row = cur.fetchone()
        return dict(row) if row else None

    def count(self) -> int:
        cur = self.conn.execute("SELECT COUNT(*) FROM strains")
        return cur.fetchone()[0]

    def list_genera(self) -> list[str]:
        cur = self.conn.execute(
            "SELECT DISTINCT genus FROM strains WHERE genus IS NOT NULL ORDER BY genus"
        )
        return [r[0] for r in cur.fetchall()]

    def expand_from_ncbi(self, taxon: str, ncbi_client: NCBIClient,
                         assembly_level: str = "complete_genome",
                         max_results: int = 20) -> list[dict]:
        """Search NCBI for assemblies and add to local DB.

        Returns list of newly added strain dicts.
        """
        assemblies = ncbi_client.search_assemblies(
            taxon=taxon,
            assembly_level=assembly_level,
            refseq_only=True,
            limit=max_results,
        )

        added = []
        for asm in assemblies:
            accession = asm.get("accession", "")
            organism = asm.get("organism_name", "")

            # Check if already exists
            cur = self.conn.execute(
                "SELECT strain_id FROM strains WHERE gcf_accession = ?",
                (accession,)
            )
            if cur.fetchone():
                continue

            parts = organism.split()
            genus = parts[0] if parts else ""
            species = parts[1] if len(parts) > 1 else ""
            strain_name = " ".join(parts[2:]) if len(parts) > 2 else ""

            sid = self.add_strain(
                genus=genus, species=species, strain_name=strain_name,
                gcf=accession,
                assembly_level=asm.get("assembly_level", ""),
                source="ncbi_search",
                taxonomy_id=str(asm.get("tax_id", "")),
                notes=f"NCBI search: {taxon}"
            )

            if sid > 0:
                added.append({
                    "strain_id": sid, "genus": genus, "species": species,
                    "strain_name": strain_name, "GCF": accession,
                    "organism_name": organism,
                })

        logger.info(f"Added {len(added)} strains from NCBI search '{taxon}'")
        return added

    def get_strain_df(self) -> pd.DataFrame:
        """Export as DataFrame compatible with existing scoring code."""
        rows = self.get_all()
        if not rows:
            return pd.DataFrame(columns=[
                "strain_id", "genus", "species", "strain_name", "GCF", "notes", "full_name"
            ])

        df = pd.DataFrame(rows)
        # Rename columns for compatibility
        rename_map = {"gcf_accession": "GCF"}
        df = df.rename(columns=rename_map)

        # Ensure required columns
        for col in ["strain_id", "genus", "species", "strain_name", "GCF", "notes", "full_name"]:
            if col not in df.columns:
                df[col] = ""

        return df

    def delete_strain(self, strain_id: int) -> bool:
        cur = self.conn.execute("DELETE FROM strains WHERE strain_id = ?", (strain_id,))
        self.conn.commit()
        return cur.rowcount > 0

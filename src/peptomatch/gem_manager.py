"""GEM (Genome-scale Metabolic Model) management for PeptoMatch.

Responsibilities
----------------
1. Manage local cache of GEM files (SBML .xml) per strain (keyed by GCF accession)
2. Wrap gapseq invocation (find -> draft -> fill) — runs via WSL on Windows
3. Validate / list / fetch GEMs for FBA consumption

gapseq runs ONLY on Linux/WSL — this module assumes that environment when
generating, but caching/loading works on any OS once the SBML file exists.

Typical flow
------------
    mgr = GEMManager()
    if not mgr.has_gem("GCF_000014425.1"):
        mgr.generate_gem(gcf="GCF_000014425.1",
                         genome_fasta=Path("/mnt/d/.../genome.fna"))
    gem_path = mgr.get_gem_path("GCF_000014425.1")
    # → pass to FBASimulator
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

logger = logging.getLogger("peptomatch.gem")

# ── Paths ──────────────────────────────────────────────────────
DEFAULT_GEM_CACHE_DIR = Path("outputs/gem_cache")
DEFAULT_GAPSEQ_WORK_DIR = Path("outputs/gapseq_work")


@dataclass
class GEMInfo:
    """Metadata for a single GEM file."""
    gcf: str
    path: Path
    size_bytes: int
    exists: bool

    def to_dict(self) -> dict:
        return {
            "gcf": self.gcf,
            "path": str(self.path),
            "size_bytes": self.size_bytes,
            "exists": self.exists,
        }


class GEMManager:
    """Local cache manager for genome-scale metabolic models."""

    def __init__(
        self,
        cache_dir: Optional[Path] = None,
        gapseq_work_dir: Optional[Path] = None,
        wsl_distro: Optional[str] = None,
    ):
        """
        Parameters
        ----------
        cache_dir : Path
            Where final SBML files live (e.g. outputs/gem_cache/{GCF}.xml)
        gapseq_work_dir : Path
            Scratch space for gapseq find/draft/fill intermediates
        wsl_distro : str | None
            WSL distribution name (e.g. "Ubuntu-22.04"); None = default distro.
            On Linux this is unused.
        """
        self.cache_dir = Path(cache_dir or DEFAULT_GEM_CACHE_DIR)
        self.work_dir = Path(gapseq_work_dir or DEFAULT_GAPSEQ_WORK_DIR)
        self.wsl_distro = wsl_distro
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.work_dir.mkdir(parents=True, exist_ok=True)

    # ── cache lookup ─────────────────────────────────────────
    def gem_path_for(self, gcf: str) -> Path:
        """Canonical path: {cache_dir}/{GCF}.xml"""
        gcf = gcf.strip()
        return self.cache_dir / f"{gcf}.xml"

    def has_gem(self, gcf: str) -> bool:
        p = self.gem_path_for(gcf)
        return p.exists() and p.stat().st_size > 0

    def get_gem_path(self, gcf: str) -> Optional[Path]:
        """Return Path if present, else None."""
        return self.gem_path_for(gcf) if self.has_gem(gcf) else None

    def list_gems(self) -> list[GEMInfo]:
        """All cached GEMs."""
        out: list[GEMInfo] = []
        for p in sorted(self.cache_dir.glob("*.xml")):
            out.append(GEMInfo(
                gcf=p.stem,
                path=p,
                size_bytes=p.stat().st_size,
                exists=True,
            ))
        return out

    # ── gapseq invocation ────────────────────────────────────
    def _wsl_prefix(self) -> list[str]:
        """Return the wsl.exe prefix used to call linux commands from Windows."""
        if self.wsl_distro:
            return ["wsl.exe", "-d", self.wsl_distro, "--"]
        return ["wsl.exe", "--"]

    def _is_windows(self) -> bool:
        import platform
        return platform.system().lower() == "windows"

    def _gapseq_available(self) -> bool:
        """Best-effort check that gapseq is callable."""
        cmd = self._wsl_prefix() + ["which", "gapseq"] if self._is_windows() else ["which", "gapseq"]
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            return r.returncode == 0 and r.stdout.strip() != ""
        except Exception:
            return False

    def generate_gem(
        self,
        gcf: str,
        genome_fasta: Path,
        bitscore_cutoff: int = 200,
        timeout_seconds: int = 7200,  # 2h default
    ) -> Path:
        """Run gapseq find → draft → fill for a strain genome.

        Parameters
        ----------
        gcf : str
            GCF accession (used as cache key + filename)
        genome_fasta : Path
            Genome FASTA file (.fna or .fasta) — accessible from WSL/Linux
        bitscore_cutoff : int
            gapseq find -b parameter (default 200)
        timeout_seconds : int
            Hard upper limit for the full pipeline

        Returns
        -------
        Path to the generated SBML file in cache_dir.

        Notes
        -----
        - Long-running (typically 30 min ~ several hours per strain)
        - Requires gapseq installed in WSL (or Linux)
        - Should typically be run once per strain, then cached
        """
        if not self._gapseq_available():
            raise RuntimeError(
                "gapseq is not callable. On Windows, ensure WSL is installed "
                "and gapseq is on PATH inside the WSL distro."
            )

        work = self.work_dir / gcf
        work.mkdir(parents=True, exist_ok=True)

        prefix = self._wsl_prefix() if self._is_windows() else []

        # Convert Windows path → /mnt/x/... if running through WSL
        fasta_arg = self._to_wsl_path(genome_fasta) if self._is_windows() else str(genome_fasta)
        work_arg = self._to_wsl_path(work) if self._is_windows() else str(work)

        # Step 1: find pathways
        cmd_find = prefix + [
            "bash", "-lc",
            f"cd {work_arg} && gapseq find -p all -b {bitscore_cutoff} {fasta_arg}",
        ]
        logger.info(f"[gapseq find] {' '.join(cmd_find)}")
        subprocess.run(cmd_find, check=True, timeout=timeout_seconds)

        # Step 2: draft
        # gapseq find produces {basename}-Pathways.tbl etc; draft consumes them
        basename = genome_fasta.stem
        cmd_draft = prefix + [
            "bash", "-lc",
            f"cd {work_arg} && gapseq draft -r {basename}-all-Reactions.tbl "
            f"-t {basename}-Transporter.tbl -p {basename}-all-Pathways.tbl "
            f"-c {fasta_arg}",
        ]
        logger.info(f"[gapseq draft] {' '.join(cmd_draft)}")
        subprocess.run(cmd_draft, check=True, timeout=timeout_seconds)

        # Step 3: fill
        cmd_fill = prefix + [
            "bash", "-lc",
            f"cd {work_arg} && gapseq fill -m {basename}-draft.RDS "
            f"-c {basename}-rxnWeights.RDS -g {basename}-rxnXgenes.RDS "
            f"-n MRS",  # default medium for gap-filling; can be parameterized
        ]
        logger.info(f"[gapseq fill] {' '.join(cmd_fill)}")
        subprocess.run(cmd_fill, check=True, timeout=timeout_seconds)

        # Locate produced SBML and copy to cache
        produced = next(work.glob(f"{basename}*.xml"), None)
        if produced is None:
            raise RuntimeError(
                f"gapseq finished but no SBML found under {work}. "
                "Check gapseq output for errors."
            )
        target = self.gem_path_for(gcf)
        shutil.copy2(produced, target)
        logger.info(f"GEM cached: {target}")
        return target

    @staticmethod
    def _to_wsl_path(p: Path) -> str:
        """Convert C:\\foo\\bar → /mnt/c/foo/bar for WSL invocation."""
        s = str(p.resolve())
        if len(s) >= 2 and s[1] == ":":
            drive = s[0].lower()
            rest = s[2:].replace("\\", "/")
            return f"/mnt/{drive}{rest}"
        return s.replace("\\", "/")

    # ── validation ───────────────────────────────────────────
    def validate_gem(self, gcf_or_path) -> dict:
        """Light-weight SBML validation via cobra.

        Returns dict: {"valid": bool, "n_reactions": int, "n_metabolites": int,
                       "objective": str, "error": str|None}
        """
        if isinstance(gcf_or_path, (str, bytes)):
            path = self.get_gem_path(gcf_or_path)  # treat as GCF
            if path is None:
                path = Path(gcf_or_path)  # fall back: treat as path
        else:
            path = Path(gcf_or_path)

        if not path.exists():
            return {"valid": False, "error": f"file not found: {path}"}

        try:
            import cobra  # local import — keep module importable without cobra
            model = cobra.io.read_sbml_model(str(path))
            return {
                "valid": True,
                "n_reactions": len(model.reactions),
                "n_metabolites": len(model.metabolites),
                "n_genes": len(model.genes),
                "objective": str(model.objective.expression),
                "error": None,
            }
        except ImportError:
            return {"valid": False, "error": "cobra not installed (pip install cobra)"}
        except Exception as e:
            return {"valid": False, "error": str(e)}


# ── helpers used by API layer ──────────────────────────────────
_DEFAULT_MGR: Optional[GEMManager] = None


def get_default_manager() -> GEMManager:
    """Shared singleton for the FastAPI process."""
    global _DEFAULT_MGR
    if _DEFAULT_MGR is None:
        _DEFAULT_MGR = GEMManager()
    return _DEFAULT_MGR

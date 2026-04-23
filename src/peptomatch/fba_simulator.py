"""FBA (Flux Balance Analysis) simulator for PeptoMatch.

Wraps COBRApy to:
1. Load a GEM (SBML) for a strain
2. Convert basal medium + peptone composition → exchange-reaction bounds
3. Run FBA to predict specific growth rate (μ, h⁻¹)
4. Extract shadow prices (bottleneck nutrients)
5. (Optional) Optimize medium composition under constraints

Conversion convention
---------------------
For each KEGG-mapped compound we set:

    EX_{cpd_id}_e.lower_bound = -mmol/L   (uptake is negative flux in COBRApy)

mmol/L is computed from:
  - basal medium: read directly from basal_media_composition.xlsx
  - peptone:       per-group unit conversion (see UNIT_FACTORS below)

Units cheat sheet (composition_template.xlsx — confirmed with analytical lab)
----------------------------------------------------------------------------
  faa_*, taa_*       : %w/w  (g AA per 100 g powder)
  mineral_*          : mg/kg (ppm)
  orgacid_*          : mg/kg
  sugar_*            : mg/kg
  nucleotide_*       : mg/kg (assumed)
  vitB_*             : mg/kg (assumed)

Unit conversion to mmol/L:
  %w/w (AA):    val × peptone_conc_g_per_L × 10 / MW
  mg/kg:        val × peptone_conc_g_per_L / 1000 / MW
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Optional

import pandas as pd

logger = logging.getLogger("peptomatch.fba")

# Per-column-group unit → multiplier that takes (raw_value × peptone_conc_g_per_L)
# and produces mg/L (before dividing by MW(g/mol) to get mmol/L).
#
# Derivations:
#   %w/w  = g per 100 g powder
#         = (val/100) g AA per g powder × peptone_conc g powder/L
#         = val × peptone_conc / 100 g AA/L
#         = val × peptone_conc × 10 mg/L
#   mg/kg = mg per 1000 g powder
#         = val × peptone_conc / 1000 mg/L
UNIT_FACTORS = {
    "faa_":        10.0,     # %w/w → mg/L = val × conc × 10
    "taa_":        10.0,     # %w/w
    "mineral_":    1e-3,     # mg/kg → mg/L = val × conc / 1000
    "orgacid_":    1e-3,     # mg/kg
    "sugar_":      1e-3,     # mg/kg
    "nucleotide_": 1e-3,     # mg/kg (assumed)
    "vitB_":       1e-3,     # mg/kg (assumed)
}

# Default amino-acid molecular weights (g/mol) for common AAs.
# Used when peptone composition gives mass values but no MW column.
AA_MW = {
    "Ala": 89.09, "Arg": 174.20, "Asn": 132.12, "Asp": 133.10,
    "Cys": 121.16, "Cys2": 240.30, "Gln": 146.15, "Glu": 147.13,
    "Gly": 75.07, "His": 155.16, "Hyp": 131.13, "Ile": 131.17,
    "Leu": 131.17, "Lys": 146.19, "Met": 149.21, "Orn": 132.16,
    "Phe": 165.19, "Pro": 115.13, "Ser": 105.09, "Thr": 119.12,
    "Trp": 204.23, "Tyr": 181.19, "Val": 117.15, "Cit": 175.19,
    "GABA": 103.12,
}

# Mapping from composition column name → KEGG compound id (e.g. cpd00033 = Glycine)
# Add entries here as more KEGG mappings are confirmed.
# Currently 45/89 columns are KEGG-mapped per the project status doc.
KEGG_COMPOUND_MAP = {
    # Free amino acids → KEGG / ModelSEED cpd ids
    "faa_Alanine":        "cpd00035",
    "faa_Arginine":       "cpd00051",
    "faa_Asparagine":     "cpd00132",
    "faa_Aspartic acid":  "cpd00041",
    "faa_Cysteine":       "cpd00084",
    "faa_Glutamic acid":  "cpd00023",
    "faa_Glutamine":      "cpd00053",
    "faa_Glycine":        "cpd00033",
    "faa_Histidine":      "cpd00119",
    "faa_Isoleucine":     "cpd00322",
    "faa_Leucine":        "cpd00107",
    "faa_Lysine":         "cpd00039",
    "faa_Methionine":     "cpd00060",
    "faa_Phenylalanine":  "cpd00066",
    "faa_Proline":        "cpd00129",
    "faa_Serine":         "cpd00054",
    "faa_Threonine":      "cpd00161",
    "faa_Tryptophan":     "cpd00065",
    "faa_Tyrosine":       "cpd00069",
    "faa_Valine":         "cpd00156",
    # Sugars
    "sugar_Glucose":      "cpd00027",
    "sugar_Fructose":     "cpd00082",
    "sugar_Sucrose":      "cpd00076",
    "sugar_Lactose":      "cpd00208",
    "sugar_Maltose":      "cpd00179",
    # Minerals
    "mineral_Na":         "cpd00971",
    "mineral_K":          "cpd00205",
    "mineral_Mg":         "cpd00254",
    "mineral_Ca":         "cpd00063",
    # Nucleotides (released bases)
    "nucleotide_AMP":     "cpd00018",
    "nucleotide_GMP":     "cpd00126",
    "nucleotide_UMP":     "cpd00091",
    "nucleotide_IMP":     "cpd00246",
    "nucleotide_CMP":     "cpd00046",
    "nucleotide_Hypoxanthine": "cpd00226",
    # Organic acids
    "orgacid_Citric":     "cpd00137",
    "orgacid_Malic":      "cpd00130",
    "orgacid_Succinic":   "cpd00036",
    "orgacid_Lactic":     "cpd00159",
    "orgacid_Acetic":     "cpd00029",
    # B vitamins
    "vitB_B1":            "cpd00305",  # Thiamine
    "vitB_B2":            "cpd00220",  # Riboflavin
    "vitB_B3":            "cpd00218",  # Niacin
    "vitB_B6":            "cpd00263",  # Pyridoxine
    "vitB_B9":            "cpd00393",  # Folate
}


# ── Loaders ────────────────────────────────────────────────────

def load_basal_media(path: Path) -> pd.DataFrame:
    """Load basal_media_composition.xlsx and clean it up."""
    df = pd.read_excel(path)
    df = df.dropna(subset=["media_id", "component"])
    df["g_per_L"] = pd.to_numeric(df["g_per_L"], errors="coerce")
    df["MW (g/mol)"] = pd.to_numeric(df["MW (g/mol)"], errors="coerce")
    return df


def load_peptone_composition(path: Path) -> pd.DataFrame:
    """Load composition_template.xlsx (one row per peptone)."""
    return pd.read_excel(path)


# ── Composition → mmol/L ───────────────────────────────────────

def basal_medium_to_mmol_L(basal_df: pd.DataFrame, media_id: str) -> dict[str, float]:
    """Returns {kegg_cpd: mmol/L} for one basal medium (e.g. MRS).

    Skips rows with missing kegg_cpd / MW (e.g. complex extracts replaced by peptone).
    """
    sub = basal_df[basal_df["media_id"] == media_id]
    out: dict[str, float] = {}
    for _, r in sub.iterrows():
        cpd = str(r.get("kegg_cpd", "")).strip()
        mw = r.get("MW (g/mol)")
        g_l = r.get("g_per_L")
        if not cpd or cpd in ("-", "nan", "NaN") or pd.isna(mw) or pd.isna(g_l):
            continue
        try:
            out[cpd] = float(g_l) / float(mw) * 1000.0  # mmol/L
        except (ValueError, ZeroDivisionError):
            continue
    return out


def peptone_to_mmol_L(
    composition_row: pd.Series,
    peptone_conc_g_per_L: float,
) -> dict[str, float]:
    """Convert one peptone's composition row → {kegg_cpd: mmol/L}.

    Composition values are interpreted as mg / 100 g powder (default).
    """
    out: dict[str, float] = {}
    for col, cpd in KEGG_COMPOUND_MAP.items():
        if col not in composition_row.index:
            continue
        val = composition_row[col]
        if pd.isna(val):
            continue
        # Robust numeric coercion: handles "< LOQ", "N/D", "-", blank strings, etc.
        try:
            val_f = float(val)
        except (TypeError, ValueError):
            # Non-numeric sentinel → treat as undetected = 0
            continue
        if val_f == 0:
            continue

        # MW resolution: prefer AA_MW for amino acid columns, else skip when unknown
        mw = _resolve_mw(col)
        if mw is None:
            continue

        # Pick unit-conversion factor based on column prefix
        factor = None
        for prefix, f in UNIT_FACTORS.items():
            if col.startswith(prefix):
                factor = f
                break
        if factor is None:
            # Unknown group — skip to avoid silent wrong scaling
            logger.warning(f"No unit factor for column '{col}'; skipping")
            continue

        # factor × val × peptone_conc = mg/L; / MW(g/mol) = mmol/L
        mg_per_L = val_f * peptone_conc_g_per_L * factor
        mmol_per_L = mg_per_L / mw
        out[cpd] = out.get(cpd, 0.0) + mmol_per_L
    return out


def _resolve_mw(col_name: str) -> Optional[float]:
    """Return MW (g/mol) for a composition column, or None if unknown."""
    if col_name.startswith("faa_") or col_name.startswith("taa_"):
        aa_full = col_name.split("_", 1)[1].strip().lower()
        # Match by full name → 3-letter code
        from .composition_features import CompositionFeatureExtractor as _CFE
        code = _CFE.AA_NAME_MAP.get(aa_full)
        return AA_MW.get(code) if code else None
    # Sugars
    sugar_mw = {
        "Glucose": 180.16, "Fructose": 180.16, "Sucrose": 342.30,
        "Lactose": 342.30, "Maltose": 342.30,
    }
    if col_name.startswith("sugar_"):
        return sugar_mw.get(col_name.split("_", 1)[1])
    # Minerals (atomic weights for the cation)
    mineral_mw = {"Na": 22.99, "K": 39.10, "Mg": 24.31, "Ca": 40.08}
    if col_name.startswith("mineral_"):
        return mineral_mw.get(col_name.split("_", 1)[1])
    # Organic acids (free acid form)
    orgacid_mw = {
        "Citric": 192.12, "Malic": 134.09, "Succinic": 118.09,
        "Lactic": 90.08, "Acetic": 60.05,
    }
    if col_name.startswith("orgacid_"):
        return orgacid_mw.get(col_name.split("_", 1)[1])
    # Nucleotides
    nuc_mw = {
        "AMP": 347.22, "GMP": 363.22, "UMP": 324.18,
        "IMP": 348.21, "CMP": 323.20, "Hypoxanthine": 136.11,
    }
    if col_name.startswith("nucleotide_"):
        return nuc_mw.get(col_name.split("_", 1)[1])
    # B vitamins
    vit_mw = {"B1": 265.36, "B2": 376.36, "B3": 123.11, "B6": 169.18, "B9": 441.40}
    if col_name.startswith("vitB_"):
        return vit_mw.get(col_name.split("_", 1)[1])
    return None


# ── FBA Simulator ──────────────────────────────────────────────

class FBASimulator:
    """Per-strain FBA wrapper.

    Construction is cheap, but `optimize()` calls cobra.Model.optimize() which
    does the linear-programming work. Reuse one simulator per strain.
    """

    def __init__(self, gem_path: Path):
        self.gem_path = Path(gem_path)
        self._model = None  # lazily loaded

    @property
    def model(self):
        """Lazily load the SBML model so the constructor is cheap."""
        if self._model is None:
            try:
                import cobra
            except ImportError as e:
                raise ImportError(
                    "cobra is required for FBA. Install with: pip install cobra"
                ) from e
            self._model = cobra.io.read_sbml_model(str(self.gem_path))
            logger.info(
                f"Loaded GEM {self.gem_path.name}: "
                f"{len(self._model.reactions)} reactions, "
                f"{len(self._model.metabolites)} metabolites"
            )
        return self._model

    # ── medium setup ──────────────────────────────────────

    # Trace essentials that every realistic aqueous medium provides.
    # These stay open during reset_medium() so set_medium() doesn't have to
    # duplicate them in every basal/peptone spec. ModelSEED cpd ids.
    #
    # Philosophy: gapseq models gap-filled with ALLmed assume trace AA/vitamin
    # availability. When restricting to MRS+peptone bounds, we'd lose those
    # implicit assumptions and get μ=0 for auxotrophs (e.g. Lactobacillus
    # needs Pro/Cys/Trp often below detection in peptone analyses).
    #
    # Bound levels:
    #   Large (100+): ubiquitous inorganics (water, H+, CO2, etc.)
    #   Medium (1~10): major ions
    #   Trace AA (0.1): mimics scavenging from peptide bonds + carryover
    #   Trace vit (0.001): vitamins carry over in trace from YE remnants
    #   Nucleobases (0.05): salvage pathway substrates
    _TRACE_ESSENTIALS = {
        # ── Inorganics (very large) ──────────────────────────
        "cpd00001": 1000.0,   # H2O
        "cpd00007": 20.0,     # O2
        "cpd00011": 1000.0,   # CO2
        "cpd00067": 1000.0,   # H+
        "cpd00009": 10.0,     # Phosphate
        "cpd00013": 10.0,     # NH3
        "cpd00048": 10.0,     # Sulfate
        # ── Major ions ───────────────────────────────────────
        "cpd00063": 1.0,      # Ca2+
        "cpd00099": 10.0,     # Cl-
        "cpd00205": 10.0,     # K+
        "cpd00254": 1.0,      # Mg2+
        "cpd00971": 10.0,     # Na+
        # ── Trace metals ─────────────────────────────────────
        "cpd00030": 0.1,      # Mn2+  (KEY for lactobacilli)
        "cpd00034": 0.1,      # Zn2+
        "cpd00058": 0.01,     # Cu2+
        "cpd00149": 0.01,     # Co2+
        "cpd00244": 0.01,     # Ni2+
        "cpd10515": 0.1,      # Fe2+
        "cpd10516": 0.1,      # Fe3+
        # ── Amino acid scavenging (trace, auxotroph rescue) ──
        # Peptones provide free AAs, but analyses often miss <LOQ items.
        # This 0.1 mmol/L baseline simulates peptide-bound / trace carryover.
        "cpd00023": 0.1,      # Glu
        "cpd00033": 0.1,      # Gly
        "cpd00035": 0.1,      # Ala
        "cpd00039": 0.1,      # Lys
        "cpd00041": 0.1,      # Asp
        "cpd00051": 0.1,      # Arg
        "cpd00053": 0.1,      # Gln
        "cpd00054": 0.1,      # Ser
        "cpd00060": 0.1,      # Met
        "cpd00065": 0.1,      # Trp  (often missing in analyses)
        "cpd00066": 0.1,      # Phe
        "cpd00069": 0.1,      # Tyr
        "cpd00084": 0.1,      # Cys  (often < LOQ)
        "cpd00107": 0.1,      # Leu
        "cpd00119": 0.1,      # His
        "cpd00129": 0.1,      # Pro  (bottleneck for LR/LA/LP!)
        "cpd00132": 0.1,      # Asn
        "cpd00156": 0.1,      # Val
        "cpd00161": 0.1,      # Thr
        "cpd00322": 0.1,      # Ile
        # ── B vitamins (trace from YE remnants) ──────────────
        # Note: raised to 0.01 (was 0.001) because gapseq ALLmed supplies
        # these at -1 and lactobacilli have very high FAD/FMN demand.
        "cpd00218": 0.01,     # Niacin (B3)
        "cpd00220": 0.01,     # Riboflavin (B2)
        "cpd00263": 0.01,     # Pyridoxine (B6)
        "cpd00305": 0.01,     # Thiamine (B1)
        "cpd00393": 0.01,     # Folate (B9)
        "cpd00644": 0.01,     # Pantothenate (B5)
        "cpd00104": 0.01,     # Biotin (B7)
        "cpd00166": 0.01,     # Calomel/Calcium-pantothenate salt
        # Flavin cofactors — direct biomass precursors
        "cpd00015": 0.01,     # FAD
        "cpd00050": 0.01,     # FMN
        # B6 (multiple forms — LP needs specifically pyridoxal, not pyridoxine)
        "cpd00215": 0.01,     # Pyridoxal
        "cpd00016": 0.01,     # Pyridoxal-5-phosphate (PLP)
        # Cell-wall precursor (some LAB can't synthesize DAP — e.g. LA)
        "cpd00516": 0.05,     # meso-2,6-Diaminopimelate (peptidoglycan cross-link)
        # ── Nucleobases / nucleosides (salvage pathway) ──────
        "cpd00226": 0.05,     # Hypoxanthine
        "cpd00092": 0.05,     # Uracil
        "cpd00182": 0.05,     # Adenine
        "cpd00307": 0.05,     # Cytosine
        "cpd00309": 0.05,     # Thymine
        "cpd00311": 0.05,     # Guanosine
        "cpd00018": 0.05,     # AMP
        "cpd00038": 0.01,     # GTP
        # ── Lipid-synthesis precursors (gapseq ALLmed items) ─
        # gapseq gap-fills biomass against ALLmed which provides these for free.
        # Without them, cpd15xxx lipoteichoic-acid / cardiolipin / PG synthesis
        # is blocked and μ collapses to 0 for Gram-positive strains.
        "cpd01080": 0.5,      # ocdca (stearic acid, C18:0) — key lipid C source
        "cpd03847": 0.1,      # Myristic acid (C14:0)
        "cpd00214": 0.1,      # Hexadecanoic / palmitic acid (C16:0)
        "cpd00080": 0.1,      # Glycerol-3-phosphate (phospholipid backbone)
        "cpd00098": 0.1,      # Choline
        # ── Cofactor precursors (often "leaked" at low level in vivo) ──
        "cpd00010": 0.01,     # CoA
        "cpd00006": 0.01,     # NADP
        "cpd00003": 0.01,     # NAD
        "cpd00793": 0.001,    # Thiamine phosphate (TPP)
        "cpd00028": 0.1,      # Heme (raised: gapseq ALLmed had -1, LAB biomass demand)
        # ── Polyamines (lactobacilli often require) ──────────
        "cpd00118": 0.01,     # Putrescine
        "cpd00264": 0.01,     # Spermidine
        # ── Dipeptides (gapseq ALLmed carryover, mimic peptide utilisation) ─
        "cpd11588": 0.05,     # gly-pro
        "cpd11590": 0.05,     # met-ala
        "cpd11592": 0.05,     # gly-glu
        "cpd01017": 0.05,     # cys-gly
        "cpd15605": 0.05,     # gly-phe
        "cpd15606": 0.05,     # gly-tyr
        # ── Fermentation/central-metabolism intermediates ────
        "cpd00281": 0.05,     # GABA (γ-aminobutyrate)
        "cpd00276": 0.05,     # Glu-Gln (GLUM)
        "cpd00106": 0.05,     # Fumarate
        "cpd00036": 0.05,     # Succinate
        "cpd00794": 0.05,     # Trehalose
    }

    def reset_medium(self, keep_trace_essentials: bool = True) -> None:
        """Close all uptake exchanges to 0, optionally preserving trace essentials.

        Parameters
        ----------
        keep_trace_essentials : bool
            If True (default), leave small uptake bounds open for H2O, H+, CO2,
            O2, phosphate, NH3, sulfate, and common trace metals (Mg, Mn, Fe,
            Zn, etc.). These are present in every realistic medium and their
            absence forces μ=0 even with full amino-acid/sugar supply.
        """
        if not hasattr(self, "_ex_index"):
            self._ex_index = self._build_exchange_index()

        for rxn in self._exchange_reactions():
            rxn.lower_bound = 0

        if keep_trace_essentials:
            for cpd, mmol in self._TRACE_ESSENTIALS.items():
                rxn_id = self._ex_index.get(cpd)
                if rxn_id is None:
                    continue
                self.model.reactions.get_by_id(rxn_id).lower_bound = -abs(mmol)

    def _exchange_reactions(self):
        return [r for r in self.model.reactions if r.id.startswith("EX_")]

    def _build_exchange_index(self) -> dict[str, str]:
        """Build {cpd_id: full_exchange_rxn_id} map.

        Handles different SBML naming conventions:
          - gapseq/ModelSEED:  EX_cpd00027_e0
          - BiGG:              EX_glc__D_e
          - COBRA legacy:      EX_cpd00027(e)
          - Plain:             EX_cpd00027_e
        """
        import re
        idx: dict[str, str] = {}
        # Match EX_<something>_<compartment> where compartment is e, e0, etc.
        pattern = re.compile(r"^EX_(.+?)(?:_e\d*|\(e\d*\))$")
        for rxn in self._exchange_reactions():
            m = pattern.match(rxn.id)
            if m:
                core = m.group(1)
                # Keep first occurrence; gapseq typically has one per compound
                idx.setdefault(core, rxn.id)
        return idx

    def _set_uptake(self, kegg_cpd: str, mmol_per_L: float) -> bool:
        """Set lower_bound on exchange reaction for this compound.

        Returns True if an exchange was found and updated.
        """
        # Lazy-build the index once
        if not hasattr(self, "_ex_index"):
            self._ex_index = self._build_exchange_index()

        rxn_id = self._ex_index.get(kegg_cpd)
        if rxn_id is None:
            # Fallback: try direct ID variants
            for candidate in (f"EX_{kegg_cpd}_e0", f"EX_{kegg_cpd}_e",
                              f"EX_{kegg_cpd}(e)"):
                if candidate in self.model.reactions:
                    rxn_id = candidate
                    self._ex_index[kegg_cpd] = rxn_id
                    break
        if rxn_id is None:
            return False
        # Use max-allowed uptake semantics: if the reaction is already open
        # wider than the requested supply (e.g. trace-essential or a prior
        # basal/peptone contribution), keep the wider bound.  Setting a
        # SMALLER |lb| would artificially tighten uptake and violate LP
        # monotonicity when sources are combined (basal + peptone + supp).
        rxn = self.model.reactions.get_by_id(rxn_id)
        new_lb = -abs(mmol_per_L)
        if rxn.lower_bound > new_lb:  # current is less-open (closer to 0)
            rxn.lower_bound = new_lb
        return True

    def set_medium(
        self,
        basal_mmol: dict[str, float],
        peptone_mmol: Optional[dict[str, float]] = None,
        supplements_mmol: Optional[list[dict[str, float]]] = None,
        reset_first: bool = True,
    ) -> dict[str, Any]:
        """Apply combined medium to model exchanges.

        Parameters
        ----------
        basal_mmol       : dict {cpd_id: mmol/L}  — base salts/glucose/etc.
        peptone_mmol     : dict {cpd_id: mmol/L}  — primary N source (one peptone)
        supplements_mmol : list of dicts          — additional boosters
                           e.g. [{"cpd00129": 3.2, ...}] from Gistex LS Ferm 5 g/L
                           Multiple supplements are summed.
        reset_first      : reset all exchanges to closed first (default True)

        Returns
        -------
        dict with keys:
          applied: list of (cpd_id, mmol_per_L, source) actually set
          missing: list of cpd_ids that have no matching exchange in this GEM
        """
        if reset_first:
            self.reset_medium()

        peptone_mmol = peptone_mmol or {}
        supplements_mmol = supplements_mmol or []

        # Combine: peptone + supplements sum on top of basal
        combined: dict[str, tuple[float, str]] = {}
        for cpd, mm in basal_mmol.items():
            combined[cpd] = (mm, "basal")
        for cpd, mm in peptone_mmol.items():
            if cpd in combined:
                combined[cpd] = (combined[cpd][0] + mm, combined[cpd][1] + "+peptone")
            else:
                combined[cpd] = (mm, "peptone")
        for sup in supplements_mmol:
            for cpd, mm in sup.items():
                if cpd in combined:
                    combined[cpd] = (combined[cpd][0] + mm, combined[cpd][1] + "+supp")
                else:
                    combined[cpd] = (mm, "supplement")

        applied: list[tuple[str, float, str]] = []
        missing: list[str] = []
        for cpd, (mm, src) in combined.items():
            if self._set_uptake(cpd, mm):
                applied.append((cpd, mm, src))
            else:
                missing.append(cpd)
        return {"applied": applied, "missing": missing}

    # ── optimization ──────────────────────────────────────
    def predict_growth(self) -> float:
        """Run FBA → return predicted specific growth rate (h⁻¹)."""
        sol = self.model.optimize()
        if sol.status != "optimal":
            logger.warning(f"FBA non-optimal: status={sol.status}")
            return 0.0
        return float(sol.objective_value or 0.0)

    def get_shadow_prices(self, top_n: int = 10) -> list[dict]:
        """Return top-N bottleneck nutrients by |shadow price|.

        Shadow price = ∂μ/∂(metabolite supply) — high |sp| means small extra
        supply yields a large growth gain (= bottleneck).
        """
        sol = self.model.optimize()
        if sol.status != "optimal" or sol.shadow_prices is None:
            return []
        sp = sol.shadow_prices
        # Filter to extracellular metabolites only (those that appear in EX_*)
        rows = [
            {"metabolite": met, "shadow_price": float(val)}
            for met, val in sp.items()
            if abs(float(val)) > 1e-9
        ]
        rows.sort(key=lambda r: abs(r["shadow_price"]), reverse=True)
        return rows[:top_n]

    # ── medium optimization (LP-based) ────────────────────
    def optimize_medium(
        self,
        candidate_supplements: dict[str, tuple[float, float]],
        target_growth: Optional[float] = None,
    ) -> dict[str, Any]:
        """Find smallest set of supplements maximizing growth.

        Skeleton — full implementation pending.

        Parameters
        ----------
        candidate_supplements : {cpd_id: (min_mmol_L, max_mmol_L)}
        target_growth : if set, find minimum total supplement to reach μ ≥ target

        Returns
        -------
        {"recommended": [{"cpd": ..., "mmol_L": ..., "delta_growth": ...}, ...],
         "predicted_growth": float}
        """
        # Baseline growth
        baseline = self.predict_growth()
        results: list[dict] = []
        for cpd, (lo, hi) in candidate_supplements.items():
            rxn_id = f"EX_{cpd}_e"
            if rxn_id not in self.model.reactions:
                continue
            rxn = self.model.reactions.get_by_id(rxn_id)
            old = rxn.lower_bound
            try:
                rxn.lower_bound = -hi  # try max supplementation
                mu = self.predict_growth()
                results.append({
                    "cpd": cpd,
                    "mmol_L": hi,
                    "delta_growth": mu - baseline,
                })
            finally:
                rxn.lower_bound = old
        results.sort(key=lambda r: r["delta_growth"], reverse=True)
        return {
            "baseline_growth": baseline,
            "recommended": results[:10],
        }


# ── High-level helper ──────────────────────────────────────────

def predict_growth_for_recommendation(
    gem_path: Path,
    basal_df: pd.DataFrame,
    media_id: str,
    composition_df: pd.DataFrame,
    peptone_name: str,
    peptone_conc_g_per_L: float,
) -> dict[str, Any]:
    """One-shot FBA prediction for a (strain, media, peptone) triple.

    Returns
    -------
    {
        "predicted_mu": float (h⁻¹),
        "n_exchanges_set": int,
        "missing_compounds": list[str],
        "shadow_prices_top": [{"metabolite": ..., "shadow_price": ...}, ...]
    }
    """
    sim = FBASimulator(gem_path)

    # Build bounds
    basal_mmol = basal_medium_to_mmol_L(basal_df, media_id)

    pep_row = composition_df[composition_df["Sample_name"] == peptone_name]
    if pep_row.empty:
        peptone_mmol: dict[str, float] = {}
    else:
        peptone_mmol = peptone_to_mmol_L(pep_row.iloc[0], peptone_conc_g_per_L)

    info = sim.set_medium(basal_mmol, peptone_mmol)
    mu = sim.predict_growth()
    sp = sim.get_shadow_prices(top_n=10)

    return {
        "predicted_mu": mu,
        "n_exchanges_set": len(info["applied"]),
        "missing_compounds": info["missing"],
        "shadow_prices_top": sp,
    }

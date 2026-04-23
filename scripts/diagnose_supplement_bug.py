"""Diagnose why Mode B (with Gistex LS Ferm supplement) μ DECREASES vs Mode A.

LP monotonicity says: adding more substrate cannot decrease optimum μ.
If μ_B < μ_A, there MUST be a bound that gets STRICTER in Mode B.
This script finds that bound.

Usage:
    python diagnose_supplement_bug.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from peptomatch.fba_simulator import (
    FBASimulator,
    basal_medium_to_mmol_L,
    peptone_to_mmol_L,
)

STRAIN = "LR"
GEM = ROOT / "outputs" / "gem_cache" / f"{STRAIN}.xml"
BASAL_XLSX = ROOT / "data" / "basal_media_composition.xlsx"
COMP_XLSX = ROOT / "data" / "composition_template.xlsx"

PEPTONE = "SOY-1"
PEPTONE_G = 25.0
YE = "Gistex LS Ferm"
YE_G = 5.0


def dump_bounds(sim: FBASimulator) -> dict[str, tuple[float, float]]:
    """Snapshot every EX reaction bound."""
    return {
        r.id: (r.lower_bound, r.upper_bound)
        for r in sim.model.reactions
        if r.id.startswith("EX_")
    }


def main() -> None:
    print(f"=== {STRAIN}: {PEPTONE} {PEPTONE_G} g/L, ±{YE} {YE_G} g/L ===\n")

    sim = FBASimulator(str(GEM))
    basal = basal_medium_to_mmol_L(pd.read_excel(BASAL_XLSX), media_id="MRS")

    comp = pd.read_excel(COMP_XLSX)
    pep_row = comp[comp["Sample_name"] == PEPTONE].iloc[0]
    ye_row = comp[comp["Sample_name"] == YE].iloc[0]

    pep_mmol = peptone_to_mmol_L(pep_row, PEPTONE_G)
    ye_mmol = peptone_to_mmol_L(ye_row, YE_G)

    # ── Mode A: no supplement ───────────────────────────────
    sim.set_medium(basal, pep_mmol)
    mu_a = sim.predict_growth()
    bounds_a = dump_bounds(sim)
    print(f"Mode A (sole peptone):       μ = {mu_a:.4f}")

    # ── Mode B: + Gistex ────────────────────────────────────
    sim.set_medium(basal, pep_mmol, supplements_mmol=[ye_mmol])
    mu_b = sim.predict_growth()
    bounds_b = dump_bounds(sim)
    print(f"Mode B (+ Gistex {YE_G} g/L): μ = {mu_b:.4f}")
    print(f"Δμ = {mu_b - mu_a:+.4f}\n")

    # ── Find STRICTER bounds in B vs A ──────────────────────
    # A bound is stricter when lower_bound is LESS NEGATIVE (closer to 0)
    # or upper_bound is SMALLER (for outflow).
    stricter: list[tuple[str, tuple, tuple]] = []
    for rxn_id, (lb_a, ub_a) in bounds_a.items():
        lb_b, ub_b = bounds_b.get(rxn_id, (lb_a, ub_a))
        # Lower bound getting stricter = less uptake allowed
        # In COBRA, more negative lower_bound = more uptake allowed,
        # so stricter = lb_b > lb_a (closer to 0)
        if lb_b > lb_a + 1e-12:
            stricter.append(
                (rxn_id, (lb_a, ub_a), (lb_b, ub_b))
            )

    if stricter:
        print(f"!! {len(stricter)} exchange(s) became STRICTER in Mode B:")
        print(f"   (lower_bound moved toward 0, i.e. less uptake allowed)\n")
        for rxn_id, a, b in stricter[:30]:
            print(f"    {rxn_id:30s}  A: lb={a[0]:>12.4g}  →  B: lb={b[0]:>12.4g}")
        print()
    else:
        print("✓ No bounds became stricter in Mode B.\n")
        print("  Then μ_B < μ_A is NOT caused by set_medium.")
        print("  Check solver / duplicate solver calls / stale shadow_prices.\n")

    # ── Find LOOSER bounds too, for sanity ──────────────────
    looser = sum(
        1 for rxn_id, (lb_a, _) in bounds_a.items()
        if bounds_b.get(rxn_id, (lb_a, 0))[0] < lb_a - 1e-12
    )
    print(f"(Info: {looser} bounds became LOOSER in Mode B — expected.)\n")

    # ── Check Gistex composition for anomalies ──────────────
    print("=== Gistex LS Ferm compounds being applied ===")
    applied_cpds = sorted(ye_mmol.items(), key=lambda kv: -kv[1])
    for cpd, mm in applied_cpds[:15]:
        print(f"    {cpd}: {mm:.6g} mmol/L")
    print(f"    ... ({len(applied_cpds)} total)")


if __name__ == "__main__":
    main()

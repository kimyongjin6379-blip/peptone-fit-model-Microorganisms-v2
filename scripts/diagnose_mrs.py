"""
MRS + SOY-1 조건에서 왜 μ=0 나오는지 진단.

확인 사항:
1. MRS basal이 실제로 어떤 cpd_id에 얼마의 mmol/L을 주는가
2. SOY-1 펩톤이 어떤 AA에 얼마의 mmol/L을 주는가
3. LR GEM의 growth에 필수적인 exchange 중 빠진 것은?
4. 펩톤 농도를 10배 올리면 μ가 변하는가?
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import pandas as pd
from peptomatch.fba_simulator import (
    FBASimulator,
    load_basal_media,
    load_peptone_composition,
    basal_medium_to_mmol_L,
    peptone_to_mmol_L,
)

BASAL_XLSX = ROOT / "data" / "basal_media_composition.xlsx"
COMPO_XLSX = ROOT / "data" / "composition_template.xlsx"
GEM_DIR = ROOT / "outputs" / "gem_cache"

# ModelSEED cpd -> human-readable name (small helper)
CPD_NAMES = {
    "cpd00001": "H2O", "cpd00007": "O2", "cpd00009": "Phosphate",
    "cpd00011": "CO2", "cpd00013": "NH3", "cpd00023": "L-Glu",
    "cpd00027": "D-Glucose", "cpd00033": "Gly", "cpd00035": "Ala",
    "cpd00039": "Lys", "cpd00041": "Asp", "cpd00048": "Sulfate",
    "cpd00051": "Arg", "cpd00053": "Gln", "cpd00054": "Ser",
    "cpd00060": "Met", "cpd00063": "Ca", "cpd00065": "Trp",
    "cpd00066": "Phe", "cpd00067": "H+", "cpd00069": "Tyr",
    "cpd00082": "Fructose", "cpd00084": "Cys", "cpd00099": "Cl-",
    "cpd00107": "Leu", "cpd00119": "His", "cpd00129": "Pro",
    "cpd00132": "Asn", "cpd00137": "Citrate", "cpd00156": "Val",
    "cpd00159": "Lactate", "cpd00161": "Thr", "cpd00205": "K+",
    "cpd00218": "Niacin(B3)", "cpd00220": "Riboflavin(B2)",
    "cpd00226": "Hypoxanthine", "cpd00254": "Mg", "cpd00263": "B6",
    "cpd00305": "Thiamine(B1)", "cpd00322": "Ile", "cpd00393": "Folate(B9)",
    "cpd00971": "Na+", "cpd10515": "Fe2+", "cpd00030": "Mn",
    "cpd00034": "Zn", "cpd00058": "Cu",
}


def name(cpd):
    return CPD_NAMES.get(cpd, cpd)


def main():
    basal_df = load_basal_media(BASAL_XLSX)
    compo_df = load_peptone_composition(COMPO_XLSX)

    print("=" * 70)
    print("1) MRS basal → mmol/L (KEGG-mapped only)")
    print("=" * 70)
    mrs = basal_medium_to_mmol_L(basal_df, "MRS")
    for cpd, mm in sorted(mrs.items(), key=lambda kv: -kv[1]):
        print(f"  {cpd:<12} {name(cpd):<18} {mm:>10.4f} mmol/L")
    print(f"  TOTAL: {len(mrs)} components")

    print()
    print("=" * 70)
    print("2) MRS basal xlsx — RAW rows (including non-KEGG-mapped)")
    print("=" * 70)
    mrs_raw = basal_df[basal_df["media_id"] == "MRS"]
    for _, r in mrs_raw.iterrows():
        cpd = str(r.get("kegg_cpd", "-"))
        comp = str(r.get("component", "-"))
        gl = r.get("g_per_L", "-")
        mw = r.get("MW (g/mol)", "-")
        marker = "✓" if cpd and cpd not in ("-", "nan", "NaN") and pd.notna(r.get("MW (g/mol)")) else "✗"
        print(f"  {marker} {comp:<25} g/L={gl:<8} MW={mw:<8} kegg={cpd}")

    print()
    print("=" * 70)
    print("3) SOY-1 peptone → mmol/L @ 20 g/L")
    print("=" * 70)
    pep_row = compo_df[compo_df["Sample_name"] == "SOY-1"]
    if pep_row.empty:
        first_col = compo_df.columns[0]
        print(f"  Sample_name not found; first column = '{first_col}'")
        pep_row = compo_df[compo_df[first_col] == "SOY-1"]
    if pep_row.empty:
        print("  [FATAL] SOY-1 row missing entirely")
        return
    pep20 = peptone_to_mmol_L(pep_row.iloc[0], 20.0)
    for cpd, mm in sorted(pep20.items(), key=lambda kv: -kv[1]):
        print(f"  {cpd:<12} {name(cpd):<18} {mm:>10.4f} mmol/L")
    print(f"  TOTAL: {len(pep20)} compounds")

    print()
    print("=" * 70)
    print("4) Peptone 200 g/L (10x) — does more AA supply help?")
    print("=" * 70)
    pep200 = peptone_to_mmol_L(pep_row.iloc[0], 200.0)
    for cpd, mm in sorted(pep200.items(), key=lambda kv: -kv[1])[:10]:
        print(f"  {cpd:<12} {name(cpd):<18} {mm:>10.4f} mmol/L")

    print()
    print("=" * 70)
    print("5) Growth sweep: LR on MRS + SOY-1 at varying peptone conc")
    print("=" * 70)
    gem = GEM_DIR / "LR.xml"
    if not gem.exists():
        print("  LR.xml not found")
        return
    sim = FBASimulator(gem)
    for conc in [20, 50, 100, 200, 500, 1000]:
        pep_mmol = peptone_to_mmol_L(pep_row.iloc[0], conc)
        sim.set_medium(mrs, pep_mmol)
        mu = sim.predict_growth()
        print(f"  peptone = {conc:>5} g/L   applied_bounds={len(mrs)+len(pep_mmol)}  μ = {mu:.4f}")

    print()
    print("=" * 70)
    print("6) LR GEM's biomass reaction — which metabolites does it demand?")
    print("=" * 70)
    sim2 = FBASimulator(gem)
    bio = sim2.model.objective.expression
    print(f"  Objective: {str(bio)[:200]}...")

    # Find the biomass reaction directly
    biomass_rxns = [r for r in sim2.model.reactions if "biomass" in r.id.lower() or "bio" in r.id.lower()[:4]]
    print(f"  Biomass-like reactions: {[r.id for r in biomass_rxns]}")
    if biomass_rxns:
        br = biomass_rxns[0]
        reactants = [(m.id, coef) for m, coef in br.metabolites.items() if coef < 0]
        print(f"  Top biomass precursors (first 15):")
        for mid, c in sorted(reactants, key=lambda x: x[1])[:15]:
            print(f"    {mid:<20} coeff={c:+.4f}")


if __name__ == "__main__":
    main()

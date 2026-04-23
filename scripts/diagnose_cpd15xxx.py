"""
LR GEM에서 cpd15xxx shadow-price dead-end의 정체 파악.

- biomass에 얼마나 쓰이는가 (coeff)
- 생산 reaction 개수 / blocked 여부
- exchange / sink / demand 유무
- 만약 sink만 열어주면 growth 살아나는지 테스트
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import cobra
from peptomatch.fba_simulator import (
    FBASimulator, load_basal_media, load_peptone_composition,
    basal_medium_to_mmol_L, peptone_to_mmol_L,
)

BASAL_XLSX = ROOT / "data" / "basal_media_composition.xlsx"
COMPO_XLSX = ROOT / "data" / "composition_template.xlsx"
GEM = ROOT / "outputs" / "gem_cache" / "LR.xml"

TARGETS = [
    "cpd15311", "cpd15421", "cpd15540", "cpd15748",
    "cpd15775", "cpd15547", "cpd15739", "cpd15526",
    "cpd15757", "cpd15793",
]

def main():
    model = cobra.io.read_sbml_model(str(GEM))
    print(f"LR GEM: {len(model.reactions)} rxn, {len(model.metabolites)} met")

    # biomass reaction
    biomass = next((r for r in model.reactions if "bio1" in r.id.lower() or "biomass" in r.id.lower()), None)
    print(f"Biomass: {biomass.id if biomass else 'NOT FOUND'}")

    print()
    print("=" * 78)
    print(f"{'cpd':<10} {'name':<30} {'biomass_coeff':<14} {'producers':<10} {'consumers':<10}")
    print("=" * 78)

    # each target: find metabolite + inspect
    for cpd in TARGETS:
        # match all compartments
        mets = [m for m in model.metabolites if m.id.startswith(cpd + "_")]
        if not mets:
            print(f"{cpd:<10} NOT FOUND")
            continue
        for m in mets:
            # biomass coeff
            bio_coef = biomass.metabolites.get(m, 0) if biomass else 0
            # producers / consumers (excluding biomass)
            producers = [r for r in m.reactions if r.metabolites[m] > 0 and r != biomass]
            consumers = [r for r in m.reactions if r.metabolites[m] < 0 and r != biomass]
            name = (m.name or "")[:28]
            print(f"{m.id:<10} {name:<30} {bio_coef:<+14.4f} {len(producers):<10} {len(consumers):<10}")

    # For each cpd found in cytosol, try adding a sink reaction and see if growth recovers
    print()
    print("=" * 78)
    print("Test: add sink reaction for each dead-end & check growth")
    print("=" * 78)

    # load media first
    basal_df = load_basal_media(BASAL_XLSX)
    compo_df = load_peptone_composition(COMPO_XLSX)
    mrs = basal_medium_to_mmol_L(basal_df, "MRS")
    pep = peptone_to_mmol_L(
        compo_df[compo_df["Sample_name"] == "SOY-1"].iloc[0], 25.0
    )

    sim = FBASimulator(GEM)
    sim.set_medium(mrs, pep)
    mu_base = sim.predict_growth()
    print(f"Baseline mu (no sinks added): {mu_base:.4f}")

    # re-load fresh model for mutation
    m2 = cobra.io.read_sbml_model(str(GEM))
    # re-apply same bounds
    sim2 = FBASimulator(GEM)
    sim2.set_medium(mrs, pep)
    # inject sinks on sim2.model
    added = []
    for cpd in TARGETS:
        for met in sim2.model.metabolites:
            if met.id.startswith(cpd + "_c"):  # cytosol
                sink_id = f"SK_{met.id}"
                if sink_id in sim2.model.reactions:
                    continue
                try:
                    sim2.model.add_boundary(met, type="sink", lb=-10, ub=10)
                    added.append(met.id)
                except Exception as e:
                    print(f"  [skip] {met.id}: {e}")
    mu_sink = sim2.model.optimize().objective_value
    print(f"After adding sinks (lb=-10) for {len(added)} metabolites: mu = {mu_sink:.4f}")

    # sink (lb=-10) lets metabolite be DRAINED, but if we want it SUPPLIED we need
    # a demand/source where lb is more negative (or ub positive on reverse direction).
    # cobra's add_boundary(type="sink", lb=-1000, ub=1000) allows both directions.
    print()
    print("-" * 78)
    print("Retry: allow bidirectional flow (supply-capable sink)")
    print("-" * 78)
    sim3 = FBASimulator(GEM)
    sim3.set_medium(mrs, pep)
    added2 = []
    for cpd in TARGETS:
        for met in sim3.model.metabolites:
            if met.id.startswith(cpd + "_c"):
                try:
                    sim3.model.add_boundary(met, type="sink", lb=-1000, ub=1000)
                    added2.append(met.id)
                except Exception:
                    pass
    mu_bi = sim3.model.optimize().objective_value
    print(f"Bidirectional sinks for {len(added2)} mets: mu = {mu_bi:.4f}")

    if mu_bi and mu_bi > 1e-6:
        sol = sim3.model.optimize()
        print("Sink fluxes (nonzero):")
        for met_id in added2:
            flux = sol.fluxes.get(f"SK_{met_id}", 0)
            if abs(flux) > 1e-6:
                direction = "SUPPLIED from outside" if flux > 0 else "drained"
                print(f"  SK_{met_id:<18} flux = {flux:+.4f}   ({direction})")


if __name__ == "__main__":
    main()

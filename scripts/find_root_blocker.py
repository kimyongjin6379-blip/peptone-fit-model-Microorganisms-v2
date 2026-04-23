"""
LR GEM의 진짜 root blocker 찾기.

전략:
1. 모든 exchange/sink을 활짝 열어서(-1000) growth가 되는지 확인 → GEM 자체 문제인지
2. 그래도 안 되면 biomass-infeasible (GEM 재생성 필요)
3. 된다면, 어떤 cpd가 반드시 필요한지 blockedreactions + shadow price 계층적으로 추적
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import cobra

GEM = ROOT / "outputs" / "gem_cache" / "LR.xml"

def main():
    model = cobra.io.read_sbml_model(str(GEM))
    n_ex = sum(1 for r in model.reactions if r.id.startswith("EX_"))
    print(f"LR GEM: {len(model.reactions)} rxn, {n_ex} exchange rxn")

    # Test 0: default (what gapseq left it in)
    sol = model.optimize()
    print(f"[T0] Default bounds mu = {sol.objective_value:.4f}")

    # Test 1: open ALL exchanges wide → is GEM feasible at all?
    m1 = cobra.io.read_sbml_model(str(GEM))
    for r in m1.reactions:
        if r.id.startswith("EX_"):
            r.lower_bound = -1000
            r.upper_bound = 1000
    sol1 = m1.optimize()
    print(f"[T1] All exchanges OPEN (-1000,+1000) mu = {sol1.objective_value:.4f}")
    if sol1.objective_value is None or sol1.objective_value < 1e-6:
        print("    → GEM biomass infeasible even with unlimited input. Dead-end in network.")

    # Test 2: close all, apply only MRS-like minimal, see shadow prices of biomass precursors
    m2 = cobra.io.read_sbml_model(str(GEM))
    for r in m2.reactions:
        if r.id.startswith("EX_"):
            r.lower_bound = 0
    # supply only essentials (glucose + O2/CO2/H2O/H+/Pi/NH3/SO4 + all 20 AAs + B vitamins + nucleobases + metals)
    essentials = {
        "cpd00001": 1000, "cpd00007": 20, "cpd00009": 10, "cpd00011": 1000,
        "cpd00013": 10, "cpd00048": 10, "cpd00067": 1000, "cpd00027": 20,
        "cpd00023": 1, "cpd00033": 1, "cpd00035": 1, "cpd00039": 1, "cpd00041": 1,
        "cpd00051": 1, "cpd00053": 1, "cpd00054": 1, "cpd00060": 1, "cpd00065": 1,
        "cpd00066": 1, "cpd00069": 1, "cpd00084": 1, "cpd00107": 1, "cpd00119": 1,
        "cpd00129": 1, "cpd00132": 1, "cpd00156": 1, "cpd00161": 1, "cpd00322": 1,
        "cpd00218": 0.01, "cpd00220": 0.01, "cpd00263": 0.01, "cpd00305": 0.01,
        "cpd00393": 0.01, "cpd00644": 0.01, "cpd00104": 0.01,
        "cpd00226": 0.5, "cpd00092": 0.5, "cpd00018": 0.5, "cpd00182": 0.5,
        "cpd00030": 0.1, "cpd00034": 0.1, "cpd00058": 0.1, "cpd10515": 0.1,
        "cpd00063": 1, "cpd00099": 10, "cpd00205": 10, "cpd00254": 1,
        "cpd00971": 10,
    }
    applied = 0
    for cpd, val in essentials.items():
        for cand in (f"EX_{cpd}_e0", f"EX_{cpd}_e", f"EX_{cpd}(e)"):
            if cand in [r.id for r in m2.reactions]:
                m2.reactions.get_by_id(cand).lower_bound = -val
                applied += 1
                break
    sol2 = m2.optimize()
    print(f"[T2] Essentials-only (~{applied} exchanges) mu = {sol2.objective_value:.4f}")

    # Test 3: also add ALL lipids/fatty acids as exchanges (if GEM has them)
    m3 = m2.copy()
    lipid_cpds = ["cpd11493", "cpd15237", "cpd15268", "cpd15269", "cpd15270", "cpd15271",
                  "cpd15272", "cpd15273", "cpd15274", "cpd15275", "cpd15276", "cpd15277",
                  "cpd15278", "cpd15279", "cpd15280", "cpd15281", "cpd15282", "cpd15283",
                  "cpd15284", "cpd15285", "cpd15286", "cpd15287", "cpd15288", "cpd15289",
                  # actual dead-ends from previous diag
                  "cpd15311", "cpd15421", "cpd15540", "cpd15748", "cpd15775",
                  "cpd15547", "cpd15739", "cpd15526", "cpd15757", "cpd15793",
                  # fatty acids
                  "cpd00214", "cpd15269", "cpd15270"]
    added = 0
    for cpd in lipid_cpds:
        for met in m3.metabolites:
            if met.id.startswith(cpd + "_c"):
                try:
                    m3.add_boundary(met, type="sink", lb=-100, ub=1000)
                    added += 1
                except Exception:
                    pass
                break
    sol3 = m3.optimize()
    print(f"[T3] + bidirectional sinks for {added} lipid mets: mu = {sol3.objective_value:.4f}")

    # If T3 still zero, print top shadow prices
    if sol3.objective_value is None or sol3.objective_value < 1e-6:
        print()
        print("Top 20 negative shadow prices after T3:")
        try:
            sp = sol3.shadow_prices.sort_values()
            for mid, v in sp.head(20).items():
                met = m3.metabolites.get_by_id(mid) if mid in [mm.id for mm in m3.metabolites] else None
                nm = met.name[:30] if met else ""
                print(f"  {mid:<20} {nm:<32} sp={v:+.3f}")
        except Exception as e:
            print(f"  (shadow price unavailable: {e})")

    # Test 4: scan for blocked reactions among biomass precursors
    print()
    print("=" * 70)
    print("[T4] Biomass precursors that CANNOT be produced (blocked)")
    print("=" * 70)
    biomass = next((r for r in model.reactions if "bio1" in r.id.lower()), None)
    if biomass is None:
        return
    # For each reactant, check if its producing reactions have flux capacity
    from cobra.flux_analysis import find_blocked_reactions
    # this is slow so only on precursors
    precursors = [m for m, c in biomass.metabolites.items() if c < 0]
    print(f"  {len(precursors)} biomass precursors")

    # A quick proxy: for each precursor check if any producer has nonzero upper_bound path
    # Simpler check: with all exchanges wide open, find metabolites whose shadow price is highly negative
    # OR whose producing reactions are all at zero flux.
    m_open = cobra.io.read_sbml_model(str(GEM))
    for r in m_open.reactions:
        if r.id.startswith("EX_"):
            r.lower_bound = -1000
    sol_open = m_open.optimize()
    if sol_open.objective_value and sol_open.objective_value > 1e-6:
        print(f"  With all EX open, mu = {sol_open.objective_value:.4f} → GEM is alive; issue is medium")
    else:
        print(f"  With all EX open, still mu = {sol_open.objective_value}. Hunting blocked precursors...")
        for m in precursors:
            try:
                mm = m_open.metabolites.get_by_id(m.id)
            except KeyError:
                continue
            # max production: temporarily make this metabolite the objective via a demand
            try:
                with m_open:
                    demand = m_open.add_boundary(mm, type="demand")
                    m_open.objective = demand
                    sol_d = m_open.optimize()
                    prod_max = sol_d.objective_value or 0
                if prod_max < 1e-6:
                    print(f"  BLOCKED: {m.id:<25} {(m.name or '')[:30]:<32} max_prod = {prod_max:.4f}")
            except Exception as e:
                print(f"  [err {m.id}] {e}")


if __name__ == "__main__":
    main()

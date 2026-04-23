"""Check whether LR GEM has Co2+ (cpd00149) exchange, and whether the
trace-essentials pass actually sets its bound."""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import cobra
from peptomatch.fba_simulator import (
    FBASimulator, load_basal_media, load_peptone_composition,
    basal_medium_to_mmol_L, peptone_to_mmol_L,
)

GEM = ROOT / "outputs" / "gem_cache" / "LR.xml"
BASAL = ROOT / "data" / "basal_media_composition.xlsx"
COMPO = ROOT / "data" / "composition_template.xlsx"

sim = FBASimulator(GEM)

# 1. Does EX_cpd00149_* exist?
co_rxns = [r for r in sim.model.reactions if "cpd00149" in r.id.lower()]
print("Reactions mentioning cpd00149:")
for r in co_rxns:
    print(f"  {r.id:<25} lb={r.lower_bound} ub={r.upper_bound}  {r.reaction}")

# 2. After set_medium, what is Co2+ exchange bound?
basal_df = load_basal_media(BASAL)
compo_df = load_peptone_composition(COMPO)
mrs = basal_medium_to_mmol_L(basal_df, "MRS")
pep = peptone_to_mmol_L(compo_df[compo_df["Sample_name"] == "SOY-1"].iloc[0], 25.0)

info = sim.set_medium(mrs, pep)
# Find Co2+ exchange
for rxn in sim.model.reactions:
    if rxn.id.startswith("EX_") and "cpd00149" in rxn.id:
        print(f"After set_medium: {rxn.id}  lb={rxn.lower_bound}")
mu = sim.predict_growth()
print(f"\nmu = {mu:.4f}")

# 3. Raise cobalt heavily
sim2 = FBASimulator(GEM)
sim2.set_medium(mrs, pep)
for rxn in sim2.model.reactions:
    if rxn.id.startswith("EX_cpd00149"):
        rxn.lower_bound = -10  # flood with cobalt
mu2 = sim2.predict_growth()
print(f"mu with cobalt -10: {mu2:.4f}")

# 4. Top shadow prices again
sp = sim2.get_shadow_prices(top_n=20)
print("\nTop 20 shadow prices after cobalt flood:")
for row in sp:
    print(f"  {row['metabolite']:<25} sp={row['shadow_price']:+.3f}")

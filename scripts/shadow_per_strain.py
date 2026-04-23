"""각 균주별 (LR/LA/LP/EF/BS) 현재 상태 shadow price top 15."""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from peptomatch.fba_simulator import (
    FBASimulator, load_basal_media, load_peptone_composition,
    basal_medium_to_mmol_L, peptone_to_mmol_L,
)

BASAL = ROOT / "data" / "basal_media_composition.xlsx"
COMPO = ROOT / "data" / "composition_template.xlsx"
GEM_DIR = ROOT / "outputs" / "gem_cache"

basal_df = load_basal_media(BASAL)
compo_df = load_peptone_composition(COMPO)
pep_row = compo_df[compo_df["Sample_name"] == "SOY-1"].iloc[0]

CASES = [("LR", "MRS"), ("LA", "MRS"), ("LP", "MRS"), ("EF", "BHI"), ("BS", "LB")]
for ab, media in CASES:
    gem = GEM_DIR / f"{ab}.xml"
    if not gem.exists():
        continue
    sim = FBASimulator(gem)
    basal = basal_medium_to_mmol_L(basal_df, media)
    pep = peptone_to_mmol_L(pep_row, 25.0)
    sim.set_medium(basal, pep)
    mu = sim.predict_growth()
    print(f"\n=== {ab} on {media} + SOY-1 25g/L  ==>  mu = {mu:.4f} ===")
    sp = sim.get_shadow_prices(top_n=15)
    for row in sp:
        mid = row["metabolite"]
        try:
            met = sim.model.metabolites.get_by_id(mid)
            nm = (met.name or "")[:35]
        except KeyError:
            nm = ""
        print(f"  {mid:<20} {nm:<37} sp={row['shadow_price']:+.3f}")

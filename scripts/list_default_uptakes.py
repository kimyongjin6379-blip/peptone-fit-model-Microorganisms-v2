"""LR GEM default에서 열려 있는(lb<0) exchange 전체 출력.
이것이 ALLmed gap-fill 당시 의도된 minimal open set."""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
import cobra
GEM = ROOT / "outputs" / "gem_cache" / "LR.xml"
model = cobra.io.read_sbml_model(str(GEM))

open_ex = [r for r in model.reactions if r.id.startswith("EX_") and r.lower_bound < 0]
print(f"Default open uptakes: {len(open_ex)}")
print()
# try to resolve readable name via metabolite
for r in sorted(open_ex, key=lambda x: x.lower_bound):
    met = list(r.metabolites.keys())[0]
    print(f"  {r.id:<25} lb={r.lower_bound:>8.2f}  {met.name[:40]}")

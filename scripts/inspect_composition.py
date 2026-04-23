"""
composition_template.xlsx의 원본 값 덤프.

SOY-1 행의 아미노산/미네랄/유기산 원본 숫자를 보여서
단위(mg/100g vs mg/kg vs g/100g)를 추정.
"""
import sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
COMPO_XLSX = ROOT / "data" / "composition_template.xlsx"

df = pd.read_excel(COMPO_XLSX)

print(f"Columns ({len(df.columns)}): ")
for c in df.columns:
    print(f"  {c}")
print()

# Find SOY-1 row
name_col = "Sample_name" if "Sample_name" in df.columns else df.columns[0]
print(f"Name column: {name_col}")
print()

row = df[df[name_col] == "SOY-1"]
if row.empty:
    print("SOY-1 row not found. First 3 rows:")
    print(df.head(3).T)
    sys.exit(1)

r = row.iloc[0]

print("=" * 70)
print("SOY-1 row — grouped by prefix")
print("=" * 70)

groups = {
    "faa_":       "Free amino acids",
    "taa_":       "Total amino acids",
    "mineral_":   "Minerals",
    "sugar_":     "Sugars",
    "orgacid_":   "Organic acids",
    "nucleotide_": "Nucleotides",
    "vitB_":      "B vitamins",
}

for prefix, label in groups.items():
    cols = [c for c in df.columns if c.startswith(prefix)]
    if not cols:
        continue
    print(f"\n[{label}] ({len(cols)} cols)")
    for c in cols:
        v = r[c]
        print(f"  {c:<30} = {v}")

# Also print non-grouped columns
print("\n[Other columns]")
all_grouped = set()
for prefix in groups:
    all_grouped.update(c for c in df.columns if c.startswith(prefix))
for c in df.columns:
    if c in all_grouped or c == name_col:
        continue
    v = r[c]
    print(f"  {c:<30} = {v}")

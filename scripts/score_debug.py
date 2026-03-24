"""Debug scoring breakdown for L. acidophilus."""
import sys
sys.path.insert(0, "src")

import numpy as np
from peptomatch.io_loaders import load_composition_data
from peptomatch.composition_features import CompositionFeatureExtractor
from peptomatch.kegg_client import KEGGClient
from peptomatch.kegg_pathway import analyze_ko_annotations

# Load composition
comp_df = load_composition_data("data/composition_template.xlsx")

SEMPIO = [
    "SOY-1", "SOY-N+", "SOY-L", "SOY-P", "WHEAT-1", "PEA-1", "RICE-1",
    "SOY-BIO", "SOY-BIO N50", "WHEAT-BIO", "PEA-BIO", "RICE-BIO",
    "Fish Collagen", "PPR Type2", "PPR Type3", "PPR Type4",
]

extractor = CompositionFeatureExtractor(comp_df, {"peptone_filter": SEMPIO})
supply_scores = extractor.compute_supply_scores(norm_samples=SEMPIO)

# Get L. acidophilus KOs via KEGG
client = KEGGClient()
ko_list, source, org_code = client.annotate_strain("Lactobacillus", "acidophilus")
print(f"KEGG: {len(ko_list)} KOs via {org_code}")

prior = analyze_ko_annotations(ko_list)

# Build demand dict
demand = {}
for aa, completeness in prior.get("aa_biosynthesis", {}).items():
    demand[f"demand_{aa}"] = 1.0 - completeness
for vit, completeness in prior.get("vitamin_biosynthesis", {}).items():
    demand[f"demand_vitamin_{vit}"] = 1.0 - completeness
nuc_bio = prior.get("nucleotide_biosynthesis", 0.5)
if isinstance(nuc_bio, dict):
    demand["demand_nucleotide"] = 1.0 - nuc_bio.get("purine", 0.5)
else:
    demand["demand_nucleotide"] = 1.0 - float(nuc_bio)
demand["transporter_bonus"] = prior.get("transporter_score", 0.5)

for sugar, data in prior.get("sugar_metabolism", {}).items():
    val = data if isinstance(data, (int, float)) else data.get("completeness", 0)
    demand[f"sugar_utilization_{sugar}"] = val

for mineral, data in prior.get("mineral_transport", {}).items():
    val = data if isinstance(data, (int, float)) else data.get("completeness", 0)
    demand[f"mineral_demand_{mineral}"] = val

for acid, data in prior.get("organic_acid_metabolism", {}).items():
    val = data if isinstance(data, (int, float)) else data.get("completeness", 0)
    demand[f"orgacid_utilization_{acid}"] = val

# Print demand
print("\n=== L. acidophilus Demand Profile ===")
print("\nAA demand (1.0 = fully needed from outside):")
for k in sorted(demand.keys()):
    if k.startswith("demand_") and not k.startswith("demand_vitamin") and not k.startswith("demand_nuc"):
        print(f"  {k:30s} = {demand[k]:.3f}")

print("\nSugar utilization:")
for k in sorted(demand.keys()):
    if "sugar" in k:
        print(f"  {k:30s} = {demand[k]:.3f}")

print("\nMineral demand:")
for k in sorted(demand.keys()):
    if "mineral" in k:
        print(f"  {k:30s} = {demand[k]:.3f}")

print("\nOrganic acid utilization:")
for k in sorted(demand.keys()):
    if "orgacid" in k:
        print(f"  {k:30s} = {demand[k]:.3f}")

print(f"\nTransporter bonus: {demand.get('transporter_bonus', 0):.3f}")

from peptomatch.scoring import STRAIN_TYPE_PRESETS

# Use LAB preset for L. acidophilus
lab_preset = STRAIN_TYPE_PRESETS["LAB"]["weights"]
default_preset = STRAIN_TYPE_PRESETS["default"]["weights"]

print("\n=== Weight Comparison: LAB vs Default ===")
for key in sorted(lab_preset.keys()):
    lab_v = lab_preset[key]
    def_v = default_preset[key]
    marker = " <<<" if lab_v != def_v else ""
    print(f"  {key:25s}  LAB={lab_v:4.1f}  Default={def_v:4.1f}{marker}")

w = {
    "faa_abundance": lab_preset["faa_abundance"],
    "taa_abundance": lab_preset["taa_abundance"],
    "mw_low": lab_preset["mw_low"],
    "mw_medium": lab_preset["mw_medium"],
    "mw_high": lab_preset["mw_high"],
    "vitamin_b": lab_preset["vitamin_b"],
    "nucleotides": lab_preset["nucleotides"],
    "aa_match": lab_preset["aa_biosynthesis_gap"],
    "transporter": lab_preset["transporter_bonus"],
    "sugar": lab_preset["sugar"],
    "mineral": lab_preset["mineral"],
    "orgacid": lab_preset["orgacid"],
    "nitrogen_quality": lab_preset["nitrogen_quality"],
}

# Score each peptone
all_results = {}
for pep in SEMPIO:
    if pep not in supply_scores.index:
        continue
    s = supply_scores.loc[pep]
    tb = demand.get("transporter_bonus", 0.5)

    c = {}
    c["faa"] = s.get("supply_faa", 0) * w["faa_abundance"]
    c["taa"] = s.get("supply_taa", 0) * w["taa_abundance"]
    c["low_mw"] = s.get("supply_low_mw", 0) * w["mw_low"] * (1 + tb)
    c["med_mw"] = s.get("supply_medium_mw", 0) * w["mw_medium"] * tb
    c["high_mw"] = s.get("supply_high_mw", 0) * w["mw_high"] * (1 - tb * 0.3)

    vit_dem = np.mean([demand.get(f"demand_vitamin_{v}", 0.5) for v in ["B1","B2","B3","B6","B9"]])
    c["vitamin"] = s.get("supply_vitamin", 0) * vit_dem * w["vitamin_b"]

    nuc_dem = demand.get("demand_nucleotide", 0.5)
    c["nucleotide"] = s.get("supply_nucleotide", 0) * nuc_dem * w["nucleotides"]

    all_aa = ["His","Ile","Leu","Lys","Met","Phe","Thr","Trp","Val",
              "Ala","Arg","Asn","Asp","Cys","Glu","Gln","Gly","Pro","Ser","Tyr"]
    aa_score = 0
    aa_cnt = 0
    for aa in all_aa:
        sk = f"supply_{aa}"
        dk = f"demand_{aa}"
        if sk in s.index and dk in demand:
            aa_score += s[sk] * demand[dk]
            aa_cnt += 1
    c["aa_match"] = (aa_score / aa_cnt * w["aa_match"]) if aa_cnt > 0 else 0

    c["transporter"] = tb * (s.get("supply_medium_mw", 0) + s.get("supply_low_mw", 0) * 0.5) * w["transporter"]

    # Sugar
    sugar_score = 0
    sugar_cnt = 0
    for sug in ["glucose", "sucrose", "lactose", "maltose"]:
        sv = s.get(f"supply_sugar_{sug}", 0)
        dv = demand.get(f"sugar_utilization_{sug}", 0)
        if sv > 0 or dv > 0:
            sugar_score += sv * dv
            sugar_cnt += 1
    avg_util = np.mean([demand.get(f"sugar_utilization_{su}", 0) for su in ["glucose","sucrose","lactose","maltose"]])
    sugar_score += s.get("supply_sugar_total", 0) * avg_util * 0.5
    c["sugar"] = sugar_score * w["sugar"] if sugar_cnt > 0 else 0

    # Mineral
    min_score = 0
    for sm, dm in [("K","K"),("Mg","Mg"),("Ca","Ca"),("Na","Na")]:
        sv = s.get(f"supply_mineral_{sm}", 0)
        dv = demand.get(f"mineral_demand_{dm}", 0.5)
        min_score += sv * dv
    fe_d = demand.get("mineral_demand_Fe", 0)
    mn_d = demand.get("mineral_demand_Mn", 0)
    min_score += s.get("supply_mineral_Mg", 0) * (fe_d + mn_d) * 0.2
    c["mineral"] = min_score * w["mineral"]

    # Orgacid
    oa_score = 0
    for sa, da in [("lactate","lactate"),("citrate","citrate"),("acetate","acetate"),("succinate","succinate"),("malate","malate")]:
        sv = s.get(f"supply_orgacid_{sa}", 0)
        dv = demand.get(f"orgacid_utilization_{da}", 0)
        if dv >= 0.5:
            oa_score += sv * dv * 0.3
        else:
            oa_score -= sv * (1 - dv) * 0.3
    c["orgacid"] = oa_score * w["orgacid"]

    c["nq"] = s.get("supply_nitrogen_quality", 0) * w["nitrogen_quality"]

    c["TOTAL"] = sum(c.values())
    all_results[pep] = c

# Sort and print
sorted_r = sorted(all_results.items(), key=lambda x: x[1]["TOTAL"], reverse=True)

print("\n" + "=" * 120)
print("=== SCORE BREAKDOWN (L. acidophilus) ===")
print("=" * 120)

cols = ["faa", "taa", "aa_match", "low_mw", "med_mw", "high_mw", "vitamin", "nucleotide", "transporter", "sugar", "mineral", "orgacid", "nq"]
header = f"{'Product':18s} {'TOTAL':>6s} |"
for col in cols:
    header += f" {col[:5]:>5s}"
print(header)
print("-" * len(header))

for name, c in sorted_r:
    line = f"{name:18s} {c['TOTAL']:6.3f} |"
    for col in cols:
        line += f" {c[col]:5.3f}"
    print(line)

# Old vs New summary
print("\n=== OLD vs NEW Component Contribution ===")
print(f"{'Product':18s}  {'Old(AA+MW+Vit)':>16s}  {'New(Sug+Min+OA+NQ)':>20s}  {'Total':>8s}")
print("-" * 70)
for name, c in sorted_r:
    old = c["faa"] + c["taa"] + c["aa_match"] + c["low_mw"] + c["med_mw"] + c["high_mw"] + c["vitamin"] + c["nucleotide"] + c["transporter"]
    new = c["sugar"] + c["mineral"] + c["orgacid"] + c["nq"]
    total = c["TOTAL"]
    old_pct = old / total * 100 if total > 0 else 0
    new_pct = new / total * 100 if total > 0 else 0
    print(f"{name:18s}  {old:6.3f} ({old_pct:4.0f}%)      {new:6.3f} ({new_pct:4.0f}%)          {total:6.3f}")

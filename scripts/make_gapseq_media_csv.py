"""
gapseq gap-fill용 MRS.csv / BHI.csv 생성.

출력 형식 (gapseq 표준):
    compounds,name,maxFlux
    cpd00027,D-Glucose,111.0
    ...

전략:
  1. basal salts/sugars    : basal_media_composition.xlsx에서 KEGG 매핑된 것
  2. YE/BE 유래 영양소      : composition_template.xlsx의 실측 TAA% 데이터
                              (MRS 레시피의 YE 5g/L + BE 10g/L 기준 자동 계산)
  3. 나머지 cofactor/지방산 : hardcoded supplement (nucleobase/B vit/지방산/DAP 등)
  4. Tween 80 대체          : 지방산 (oleate, palmitate, stearate)
  5. 범용 inorganic         : H2O, O2, CO2, H+, NH3, SO4, Cl-, Na, K, Ca, Fe, Mn, Zn, Cu, Co

두 버전 생성:
  - MRS_strict.csv  : 기본 레시피만 (엄격, gap-fill 실패 가능)
  - MRS_rich.csv    : YE/BE 유래 AA 포함 (권장) ← 실측 조성 사용

같은 방식으로 BHI_strict.csv / BHI_rich.csv 생성.
"""
import sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from peptomatch.fba_simulator import (
    KEGG_COMPOUND_MAP, AA_MW, UNIT_FACTORS, _resolve_mw,
)

XLSX = ROOT / "data" / "basal_media_composition.xlsx"
COMPO_XLSX = ROOT / "data" / "composition_template.xlsx"
OUT_DIR = ROOT / "data" / "gapseq_media"
OUT_DIR.mkdir(exist_ok=True)


# ── YE/BE 실측 조성 → mmol/L 변환 ───────────────────────────────
# composition_template.xlsx 의 'Yeast extract'/'Beef extract' 행에서
# TAA_* 컬럼 (%w/w) 을 읽어, 해당 배지 레시피의 g/L 에 맞게 mmol/L 계산.
#
# MRS 레시피 표준:
#   - Peptone (from animal tissue)  10 g/L   ← 본 모델에선 사용자 펩톤 25 g/L
#   - Beef extract                  10 g/L
#   - Yeast extract                  5 g/L
#   - Glucose                       20 g/L
#   - Tween 80                       1 g/L
#   - K2HPO4                         2 g/L
#   - Sodium acetate                 5 g/L
#   - Ammonium citrate               2 g/L
#   - MgSO4·7H2O                   0.1 g/L
#   - MnSO4·4H2O                  0.05 g/L
#
# BHI 레시피 (approximate):
#   - Brain heart infusion          17.5 g/L   ← 동물 추출물, BE로 대체
#   - Peptone                       10 g/L    ← 사용자 펩톤
#   - NaCl                           5 g/L
#   - Na2HPO4                      2.5 g/L
#   - Glucose                        2 g/L
#
# 변환:
#   mg/L = %w/w × g_per_L × 10   (fba_simulator.UNIT_FACTORS["taa_"] = 10)
#   mmol/L = mg/L / MW

MEDIA_EXTRACT_RECIPE = {
    # media_id → {ingredient_sample_name: g_per_L}
    "MRS": {"Yeast extract": 5.0, "Beef extract": 10.0},
    "BHI": {"Beef extract": 17.5},   # BHI의 'brain heart infusion'을 BE로 근사
    "TSB": {},                        # TSB는 gapseq 내장 사용 예정이라 비워둠
    # LB recipe: Tryptone 10 g/L + Yeast extract 5 g/L + NaCl 10 g/L (실측 조성 사용)
    "LB":  {"Yeast extract": 5.0, "LPS Tryptone": 10.0},
}


def load_extract_composition() -> pd.DataFrame:
    """composition_template.xlsx 로드."""
    return pd.read_excel(COMPO_XLSX)


def extract_to_mmol_L(compo_df: pd.DataFrame, sample_name: str, g_per_L: float) -> dict[str, tuple[str, float]]:
    """YE 또는 BE 1종 → {cpd: (name, mmol/L)}.
    TAA 우선 사용 (전체 아미노산 풀, 단백분해 후 가용).
    """
    sub = compo_df[compo_df["Sample_name"] == sample_name]
    if sub.empty:
        print(f"  [WARN] '{sample_name}' not found in composition_template.xlsx")
        return {}
    row = sub.iloc[0]
    out: dict[str, tuple[str, float]] = {}

    # TAA 컬럼만 사용 (FAA는 TAA의 부분집합 — 중복 방지)
    for col, cpd in KEGG_COMPOUND_MAP.items():
        if not col.startswith("faa_"):
            continue
        taa_col = "taa_" + col.split("_", 1)[1]
        if taa_col not in row.index:
            continue
        val = row[taa_col]
        if pd.isna(val):
            continue
        try:
            val_f = float(val)
        except (TypeError, ValueError):
            continue
        if val_f <= 0:
            continue

        mw = _resolve_mw(col)  # faa_ 와 taa_ 는 같은 MW
        if mw is None:
            continue

        factor = UNIT_FACTORS["taa_"]  # 10.0 (%w/w → mg/L multiplier per g/L)
        mg_per_L = val_f * g_per_L * factor
        mmol_per_L = mg_per_L / mw

        aa_name = col.split("_", 1)[1]
        prev = out.get(cpd, (aa_name, 0.0))[1]
        out[cpd] = (aa_name, prev + mmol_per_L)
    return out


def compute_extract_supplement(compo_df: pd.DataFrame, media_id: str) -> dict[str, tuple[str, float]]:
    """배지별 YE/BE 레시피 → 총 AA 공급량 (mmol/L)."""
    recipe = MEDIA_EXTRACT_RECIPE.get(media_id, {})
    total: dict[str, tuple[str, float]] = {}
    for sample_name, g_per_L in recipe.items():
        partial = extract_to_mmol_L(compo_df, sample_name, g_per_L)
        for cpd, (name, mmol) in partial.items():
            prev = total.get(cpd, (name, 0.0))[1]
            total[cpd] = (name, prev + mmol)
    return total


# ── 나머지 supplement (AA 외: vit/nucleobase/lipid/cofactor) ────
# YE/BE에 AA 외에 B-vit, 뉴클레오사이드, 미량 cofactor가 풍부하지만
# composition_template.xlsx의 YE/BE 행에 해당 실측 데이터가 없음.
# → 문헌값 기반 근사로 보충 (gap-fill 알고리즘이 필요 시 사용).
NON_AA_SUPPLEMENT = {
    # B vitamins (YE 5 g/L 기준 대략적 공급량)
    "cpd00218": ("Niacin", 0.1),
    "cpd00220": ("Riboflavin", 0.1),
    "cpd00263": ("Pyridoxine", 0.1),
    "cpd00215": ("Pyridoxal", 0.1),
    "cpd00305": ("Thiamine", 0.1),
    "cpd00393": ("Folate", 0.1),
    "cpd00644": ("Pantothenate", 0.1),
    "cpd00104": ("Biotin", 0.1),
    "cpd00016": ("Pyridoxal-5-P", 0.1),
    # Nucleobases / nucleosides (salvage)
    "cpd00226": ("Hypoxanthine", 0.5),
    "cpd00092": ("Uracil", 0.5),
    "cpd00182": ("Adenine", 0.5),
    "cpd00307": ("Cytosine", 0.5),
    "cpd00309": ("Thymine", 0.5),
    "cpd00311": ("Guanosine", 0.5),
    "cpd00018": ("AMP", 0.1),
    "cpd00038": ("GTP", 0.1),
    # Cofactors (trace)
    "cpd00010": ("CoA", 0.05),
    "cpd00006": ("NADP", 0.05),
    "cpd00003": ("NAD", 0.05),
    "cpd00028": ("Heme", 0.05),
    # Polyamines (LAB often need)
    "cpd00118": ("Putrescine", 0.1),
    "cpd00264": ("Spermidine", 0.1),
    # Dipeptides (단백분해 산물 — gapseq draft에 있을 수 있음)
    "cpd11588": ("Gly-Pro", 0.5),
    "cpd11590": ("Met-Ala", 0.5),
    "cpd11592": ("Gly-Glu", 0.5),
    # Cell-wall precursor
    "cpd00516": ("meso-DAP", 0.5),
    "cpd00098": ("Choline", 0.5),
    "cpd00080": ("Glycerol-3-P", 0.5),
}

# Tween 80 대체 — fatty acids (C16, C18, oleate)
TWEEN80_SUBSTITUTE = {
    "cpd01080": ("Stearate (C18:0)", 1.0),
    "cpd00214": ("Palmitate (C16:0)", 1.0),
    "cpd03847": ("Myristate (C14:0)", 0.5),
    "cpd15237": ("Oleate (C18:1)", 1.0),
}

# 범용 inorganics (모든 배지 공통)
INORGANICS = {
    "cpd00001": ("H2O", 100.0),
    "cpd00007": ("O2", 20.0),
    "cpd00011": ("CO2", 100.0),
    "cpd00067": ("H+", 100.0),
    "cpd00013": ("NH3", 10.0),
    "cpd00048": ("Sulfate", 10.0),
    "cpd00099": ("Cl-", 10.0),
    "cpd00205": ("K+", 10.0),
    "cpd00971": ("Na+", 10.0),
    "cpd00063": ("Ca2+", 1.0),
    "cpd00254": ("Mg2+", 1.0),   # 필수 (ATP·리보솜·DNA pol)
    "cpd00030": ("Mn2+", 0.5),   # LAB 필수 (SOD·SPX 대체 효소)
    "cpd00149": ("Co2+", 0.1),
    "cpd00058": ("Cu2+", 0.01),
    "cpd10515": ("Fe2+", 0.1),
    "cpd10516": ("Fe3+", 0.1),
    "cpd00034": ("Zn2+", 0.1),
    "cpd11574": ("Molybdate", 0.01),
}

# BHI 기본 레시피 (g/L) — 단순 화학성분만 (BE는 extract로 따로 처리)
BHI_BASAL = [
    # (cpd, name, g_per_L, MW)
    ("cpd00027", "D-Glucose", 2.0, 180.16),
    ("cpd00971", "Na+",        5.0 / 58.44, 1.0),
    ("cpd00099", "Cl-",        5.0 / 58.44, 1.0),
    ("cpd00009", "Phosphate",  2.5 / 141.96, 1.0),
]


def load_basal(media_id: str) -> list[tuple[str, str, float]]:
    """basal xlsx에서 특정 media의 KEGG-매핑된 성분 → (cpd, name, mmol/L) 리스트."""
    df = pd.read_excel(XLSX)
    df["media_id"] = df["media_id"].astype(str).str.strip()
    rows = df[df["media_id"] == media_id]

    out = []
    for _, r in rows.iterrows():
        cpd = str(r.get("kegg_cpd", "")).strip()
        if not cpd or cpd in ("-", "nan", "NaN"):
            continue
        mw_raw = r.get("MW (g/mol)")
        try:
            mw = float(mw_raw)
        except (TypeError, ValueError):
            continue
        gl = float(r.get("g_per_L", 0))
        if gl <= 0 or mw <= 0:
            continue
        mmol = gl * 1000.0 / mw
        name = str(r.get("ion/compound", r.get("component", cpd))).strip()
        out.append((cpd, name, round(mmol, 4)))
    return out


def write_gapseq_csv(path: Path, entries: dict[str, tuple[str, float]]) -> None:
    """dict[cpd → (name, maxFlux)] → gapseq CSV."""
    lines = ["compounds,name,maxFlux"]
    for cpd, (name, flux) in entries.items():
        safe_name = name.replace(",", " ")
        lines.append(f"{cpd},{safe_name},{round(flux, 4)}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"  wrote {path.name}  ({len(entries)} compounds)")


def make_media(
    media_id: str,
    basal_rows: list[tuple[str, str, float]],
    compo_df: pd.DataFrame,
    include_supplement: bool = True,
) -> dict:
    """entries = inorganics + basal + (YE/BE 실측 AA + non-AA supplement if rich)"""
    entries: dict[str, tuple[str, float]] = {}
    # 1) inorganics
    for cpd, (name, mmol) in INORGANICS.items():
        entries[cpd] = (name, mmol)
    # 2) basal
    for cpd, name, mmol in basal_rows:
        prev = entries.get(cpd, (name, 0.0))[1]
        entries[cpd] = (name, max(mmol, prev))
    # 3) rich supplement
    if include_supplement:
        # 3a) 실측 YE/BE AA (max, 중복 누적 않고 최댓값)
        extract_aa = compute_extract_supplement(compo_df, media_id)
        for cpd, (name, mmol) in extract_aa.items():
            prev = entries.get(cpd, (name, 0.0))[1]
            entries[cpd] = (name, max(mmol, prev))
        # 3b) non-AA (vit/nuc/cofactor) — hardcoded, 실측 AA를 덮어쓰지 않음
        for cpd, (name, mmol) in NON_AA_SUPPLEMENT.items():
            entries.setdefault(cpd, (name, mmol))
        # 3c) Tween 80 대체 지방산
        for cpd, (name, mmol) in TWEEN80_SUBSTITUTE.items():
            entries.setdefault(cpd, (name, mmol))
    return entries


def print_extract_summary(compo_df: pd.DataFrame, media_id: str) -> None:
    """배지별 YE/BE AA 공급량 요약 출력."""
    recipe = MEDIA_EXTRACT_RECIPE.get(media_id, {})
    if not recipe:
        return
    print(f"  extract recipe: {recipe}")
    aa = compute_extract_supplement(compo_df, media_id)
    total_mmol = sum(v[1] for v in aa.values())
    print(f"  → {len(aa)} AA types, total {total_mmol:.1f} mmol/L")
    # top 5
    top = sorted(aa.items(), key=lambda kv: -kv[1][1])[:5]
    for cpd, (name, mmol) in top:
        print(f"    {cpd}  {name:<20} {mmol:>7.3f} mmol/L")


def main():
    compo_df = load_extract_composition()
    print(f"Loaded composition_template.xlsx: {len(compo_df)} samples")
    ye = compo_df[compo_df["Sample_name"] == "Yeast extract"]
    be = compo_df[compo_df["Sample_name"] == "Beef extract"]
    print(f"  Yeast extract row: {'OK' if len(ye) else 'MISSING'}")
    print(f"  Beef extract row : {'OK' if len(be) else 'MISSING'}")
    print()

    # ── MRS ──────────────────────────────────────
    print("▶ MRS")
    mrs_basal = load_basal("MRS")
    print(f"  basal (xlsx): {len(mrs_basal)} KEGG-mapped rows")
    print_extract_summary(compo_df, "MRS")

    mrs_strict = make_media("MRS", mrs_basal, compo_df, include_supplement=False)
    write_gapseq_csv(OUT_DIR / "MRS_strict.csv", mrs_strict)

    mrs_rich = make_media("MRS", mrs_basal, compo_df, include_supplement=True)
    write_gapseq_csv(OUT_DIR / "MRS_rich.csv", mrs_rich)

    # ── MRS_rich_cys (for L. gasseri: MRS + 0.01% cysteine) ──
    # 0.01% w/v = 0.1 g/L cysteine. MW=121.16 g/mol → 0.826 mmol/L.
    # 안전 여유로 1.0 mmol/L 설정. cpd00084 = L-Cysteine.
    mrs_rich_cys = dict(mrs_rich)
    prev_cys = mrs_rich_cys.get("cpd00084", ("L-Cysteine", 0.0))[1]
    mrs_rich_cys["cpd00084"] = ("L-Cysteine", max(prev_cys, 1.0))
    write_gapseq_csv(OUT_DIR / "MRS_rich_cys.csv", mrs_rich_cys)

    # ── BHI ──────────────────────────────────────
    print("\n▶ BHI")
    bhi_basal = []
    for cpd, name, gl, mw in BHI_BASAL:
        bhi_basal.append((cpd, name, round(gl * 1000.0 / mw, 4) if mw > 1 else round(gl, 4)))
    print(f"  basal (hardcoded): {len(bhi_basal)} rows")
    print_extract_summary(compo_df, "BHI")

    bhi_strict = make_media("BHI", bhi_basal, compo_df, include_supplement=False)
    write_gapseq_csv(OUT_DIR / "BHI_strict.csv", bhi_strict)

    bhi_rich = make_media("BHI", bhi_basal, compo_df, include_supplement=True)
    write_gapseq_csv(OUT_DIR / "BHI_rich.csv", bhi_rich)

    # ── TSB (참고용, gapseq 내장과 비교 가능) ────
    print("\n▶ TSB (reference)")
    tsb_basal = load_basal("TSB")
    print(f"  basal (xlsx): {len(tsb_basal)} KEGG-mapped rows")
    tsb_rich = make_media("TSB", tsb_basal, compo_df, include_supplement=True)
    write_gapseq_csv(OUT_DIR / "TSB_rich.csv", tsb_rich)

    # ── LB (E. coli용 — Tryptone 10 + YE 5 + NaCl 10 g/L) ──
    print("\n▶ LB")
    lb_basal = load_basal("LB")
    # LB 표준 NaCl 10 g/L = 171 mmol/L (INORGANICS 기본 10 보다 보강)
    lb_basal.append(("cpd00971", "Na+", round(10.0 * 1000.0 / 58.44, 4)))
    lb_basal.append(("cpd00099", "Cl-", round(10.0 * 1000.0 / 58.44, 4)))
    print(f"  basal (xlsx+NaCl): {len(lb_basal)} rows")
    print_extract_summary(compo_df, "LB")
    lb_rich = make_media("LB", lb_basal, compo_df, include_supplement=True)
    write_gapseq_csv(OUT_DIR / "LB_rich.csv", lb_rich)

    print(f"\nDone. All files saved to: {OUT_DIR}")


if __name__ == "__main__":
    main()

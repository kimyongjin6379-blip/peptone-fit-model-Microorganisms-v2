"""
End-to-end FBA 테스트: LR + MRS + SOY-1 조건으로 μ 예측

gap-fill 배지(ALLmed)의 이론치 μ가 아닌, 실제 MRS+펩톤 조성으로
bounds를 설정해 진짜 실험 조건의 μ를 뽑아냅니다.
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
    predict_growth_for_recommendation,
)

# ── 경로 설정 ────────────────────────────────────────
BASAL_XLSX = ROOT / "data" / "basal_media_composition.xlsx"
COMPO_XLSX = ROOT / "data" / "composition_template.xlsx"
GEM_DIR = ROOT / "outputs" / "gem_cache"

# ── 테스트 설정 ──────────────────────────────────────
TEST_CASES = [
    # (strain_abbrev, basal_media_id, peptone_name, peptone_conc_g_per_L)
    # 조건 A: 전량 대체 → peptone 25 g/L (YE, BE 제거)
    ("LR", "MRS", "SOY-1", 25.0),
    ("LA", "MRS", "SOY-1", 25.0),
    ("LP", "MRS", "SOY-1", 25.0),
    ("EF", "BHI", "SOY-1", 25.0),
    ("BS", "LB",  "SOY-1", 25.0),
]


def main():
    print("=" * 70)
    print("FBA end-to-end test: MRS/BHI/LB basal + real peptone composition")
    print("=" * 70)

    # 데이터 로드
    if not BASAL_XLSX.exists():
        print(f"[FATAL] basal_media_composition.xlsx not found at {BASAL_XLSX}")
        sys.exit(1)

    basal_df = load_basal_media(BASAL_XLSX)
    print(f"Loaded basal media: {sorted(basal_df['media_id'].dropna().unique())}")

    compo_df = load_peptone_composition(COMPO_XLSX)
    pep_names = compo_df.get("Sample_name", compo_df.iloc[:, 0]).tolist()
    print(f"Loaded peptones ({len(pep_names)} rows): first few = {pep_names[:5]}")
    print()

    # 케이스별 테스트
    print(f"{'Strain':<6} {'Basal':<5} {'Peptone':<10} {'mu (h^-1)':<12} "
          f"{'Bounds':<10} {'Missing':<8}")
    print("-" * 70)

    for abbrev, media, peptone, conc in TEST_CASES:
        gem = GEM_DIR / f"{abbrev}.xml"
        if not gem.exists():
            print(f"{abbrev:<6} SKIP (GEM not found at {gem})")
            continue

        try:
            result = predict_growth_for_recommendation(
                gem_path=gem,
                basal_df=basal_df,
                media_id=media,
                composition_df=compo_df,
                peptone_name=peptone,
                peptone_conc_g_per_L=conc,
            )
            mu = result["predicted_mu"]
            n_applied = result["n_exchanges_set"]
            n_missing = len(result["missing_compounds"])
            print(f"{abbrev:<6} {media:<5} {peptone:<10} {mu:<12.3f} "
                  f"{n_applied:<10} {n_missing:<8}")
        except Exception as e:
            print(f"{abbrev:<6} ERROR: {type(e).__name__}: {e}")

    # 세부 진단: LR 한 케이스 깊이 파보기
    print()
    print("=" * 70)
    print("DETAIL: LR + MRS + SOY-1 shadow prices (top bottlenecks)")
    print("=" * 70)

    gem = GEM_DIR / "LR.xml"
    if gem.exists():
        sim = FBASimulator(gem)
        basal_mmol = basal_medium_to_mmol_L(basal_df, "MRS")
        pep_row = compo_df[compo_df["Sample_name"] == "SOY-1"]
        if pep_row.empty:
            # fallback: first column might be the name column
            first_col = compo_df.columns[0]
            pep_row = compo_df[compo_df[first_col] == "SOY-1"]
        if not pep_row.empty:
            pep_mmol = peptone_to_mmol_L(pep_row.iloc[0], 20.0)
        else:
            print("[WARN] SOY-1 not found in composition template — using basal only")
            pep_mmol = {}

        info = sim.set_medium(basal_mmol, pep_mmol)
        print(f"Applied exchanges: {len(info['applied'])}")
        print(f"Missing from GEM : {len(info['missing'])} "
              f"(e.g. {info['missing'][:5]})")
        mu = sim.predict_growth()
        print(f"Predicted mu     : {mu:.4f} h^-1")

        sp = sim.get_shadow_prices(top_n=10)
        print("Top 10 bottleneck metabolites (shadow price):")
        for row in sp:
            print(f"  {row['metabolite']:<20} sp = {row['shadow_price']:+.3f}")


if __name__ == "__main__":
    main()

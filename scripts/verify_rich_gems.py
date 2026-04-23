"""
새로 gap-fill한 10종 GEM 검증 스크립트.

목적:
  gapseq fill 단계의 μ는 "rich medium에서 feasibility test".
  실제로 쓸 조건 (basal + SOY-1 25 g/L peptone) 에서 μ를 내서
  GEM이 진짜 쓸 수 있는지 확인.

검증 기준 (PASS/FAIL):
  1. μ > 0.1 (생리적으로 의미 있는 성장)
  2. 10종 μ 표준편차 ≥ 0.05 (균주 차이가 있어야 함)
  3. Shadow price 상위 5개가 AA/vitamin/mineral/lipid 계열 (의심 cpd 없어야)
  4. Infeasible 없음

실행 (Windows VS Code):
  cd D:\\folder1\\peptomatch
  python scripts/verify_rich_gems.py
"""
import sys
from pathlib import Path
from statistics import mean, pstdev

# Windows cp949 콘솔에서 유니코드 심볼(μ, ≥, ⁻¹, ✓) 출력 지원
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from peptomatch.fba_simulator import (
    FBASimulator, load_basal_media, load_peptone_composition,
    basal_medium_to_mmol_L, peptone_to_mmol_L,
)

BASAL_XLSX = ROOT / "data" / "basal_media_composition.xlsx"
COMPO_XLSX = ROOT / "data" / "composition_template.xlsx"
GEM_DIR = ROOT / "outputs" / "gem_cache"

# 균주 → basal medium 매핑 (실제 실험 환경 기준)
STRAIN_MEDIA = {
    # ── 기존 10종 ──
    "LP": "MRS", "LR": "MRS", "LA": "MRS", "LC": "MRS",
    "LPC": "MRS", "LS": "MRS", "STT": "MRS",
    "EF": "BHI",
    "BS": "TSB", "BC": "TSB",
    # ── Tier 1 신규 6종 (2026-04) ──
    "LRE": "MRS",   # L. reuteri
    "LB":  "MRS",   # L. delbrueckii subsp. bulgaricus
    "LF":  "MRS",   # L. fermentum
    "LH":  "MRS",   # L. helveticus
    "LG":  "MRS",   # L. gasseri (gap-fill은 MRS_rich_cys, 검증은 일반 MRS + 펩톤으로도 자체 cys 공급 가능)
    "EC":  "LB",    # E. coli BL21 — LB basal
}
PEPTONE = "SOY-1"
PEPTONE_CONC = 25.0   # g/L

# μ threshold
MU_MIN_PASS = 0.1
MU_STDDEV_MIN = 0.05

# 의심 metabolite 키워드 (shadow price에 뜨면 WARNING)
SUSPICIOUS_KEYWORDS = [
    "sulfite", "dimethyl", "propanediol", "xylitol", "arabitol",
    "mannitol", "sorbitol", "formate",  # 과발현된 side path 의심
]

# 정상 bottleneck 키워드 (있으면 OK)
NORMAL_KEYWORDS = [
    "amino", "acid", "thiamine", "riboflavin", "niacin", "folate",
    "pantothenate", "biotin", "pyridox", "cobalamin", "heme",
    "fe", "mg", "mn", "zn", "ca", "cu", "co",
    "stearate", "palmitate", "oleate", "myristate",
    "putrescine", "spermidine", "glyc", "lact", "acetate",
    "dap", "diaminopimelate", "cysteine", "methionine",
    "nucleotide", "adenine", "guanine", "cytosine", "uracil", "thymine",
    "phosphate", "sulfate", "nitrogen", "ammoni",
]


def is_normal_bottleneck(name: str) -> bool:
    nm = (name or "").lower()
    return any(kw in nm for kw in NORMAL_KEYWORDS)


def is_suspicious(name: str) -> bool:
    nm = (name or "").lower()
    return any(kw in nm for kw in SUSPICIOUS_KEYWORDS)


def _extract_sp(sim, top_n: int = 5) -> list[dict]:
    sp_top = sim.get_shadow_prices(top_n=top_n)
    out = []
    for row in sp_top:
        mid = row["metabolite"]
        try:
            met = sim.model.metabolites.get_by_id(mid)
            nm = met.name or mid
        except KeyError:
            nm = mid
        out.append({
            "id": mid, "name": nm, "sp": row["shadow_price"],
            "normal": is_normal_bottleneck(nm),
            "suspicious": is_suspicious(nm),
        })
    return out


def verify_one(strain: str, media_id: str, basal_df, pep_row,
               ye_row=None, ye_g_per_L: float = 0.0) -> dict:
    """Run FBA in two modes and report both.

    Mode A (sole peptone) : basal + peptone 25 g/L  (matches real experiment)
    Mode B (peptone + YE) : basal + peptone 25 g/L + Gistex LS Ferm 5 g/L
                            (recommendation if vitamin bottleneck)
    """
    gem = GEM_DIR / f"{strain}.xml"
    result = {"strain": strain, "media": media_id, "gem_exists": gem.exists()}

    if not gem.exists():
        result["status"] = "MISSING"
        return result

    try:
        sim = FBASimulator(gem)
        basal = basal_medium_to_mmol_L(basal_df, media_id)
        pep = peptone_to_mmol_L(pep_row, PEPTONE_CONC)

        # Mode A: sole peptone (실험 조건)
        sim.set_medium(basal, pep)
        result["mu_sole"] = sim.predict_growth()
        result["sp_sole"] = _extract_sp(sim, top_n=5)

        # Mode B: + YE 5 g/L (Gistex LS Ferm)
        if ye_row is not None and ye_g_per_L > 0:
            ye_mmol = peptone_to_mmol_L(ye_row, ye_g_per_L)
            sim.set_medium(basal, pep, supplements_mmol=[ye_mmol])
            result["mu_ye"] = sim.predict_growth()
            result["sp_ye"] = _extract_sp(sim, top_n=5)
            result["delta_mu"] = result["mu_ye"] - result["mu_sole"]
            result["improvement_pct"] = (
                (result["delta_mu"] / result["mu_sole"] * 100)
                if result["mu_sole"] > 0 else 0.0
            )

        result["status"] = "OK"
    except Exception as e:
        result["status"] = "ERROR"
        result["error"] = str(e)
    return result


def _print_sp(sp_list):
    for sp in sp_list:
        marker = "!" if sp["suspicious"] else ("✓" if sp["normal"] else "?")
        print(f"    [{marker}] {sp['id']:<18} {sp['name'][:35]:<37} sp={sp['sp']:+.3f}")


def print_strain_result(r: dict):
    print(f"\n{'='*70}")
    print(f"  {r['strain']}  on  {r['media']} + {PEPTONE} {PEPTONE_CONC} g/L")
    print('='*70)
    if r["status"] != "OK":
        print(f"  STATUS: {r['status']}")
        if "error" in r:
            print(f"  ERROR: {r['error']}")
        return

    mu_s = r["mu_sole"]
    flag = "✓" if mu_s >= MU_MIN_PASS else "✗"
    print(f"  [A] Peptone only        {flag} μ = {mu_s:.4f} h⁻¹")
    print(f"      Shadow price top 5:")
    _print_sp(r["sp_sole"])

    if "mu_ye" in r:
        mu_y = r["mu_ye"]
        dmu = r["delta_mu"]
        pct = r["improvement_pct"]
        flag2 = "✓" if mu_y >= MU_MIN_PASS else "✗"
        print()
        print(f"  [B] + YE 5 g/L          {flag2} μ = {mu_y:.4f} h⁻¹  "
              f"(Δμ={dmu:+.4f}, {pct:+.1f}%)")
        print(f"      Shadow price top 5:")
        _print_sp(r["sp_ye"])


YE_SAMPLE = "Gistex LS Ferm"   # Sempio 실험실 실제 YE
YE_G_PER_L = 5.0


def main():
    print(f"\n{'#'*70}")
    print(f"#  GEM verification - {PEPTONE} {PEPTONE_CONC} g/L on strain-specific basal")
    print(f"#  Boost compare: [B] + {YE_SAMPLE} {YE_G_PER_L} g/L")
    print(f"#  GEM dir: {GEM_DIR}")
    print('#' * 70)

    basal_df = load_basal_media(BASAL_XLSX)
    compo_df = load_peptone_composition(COMPO_XLSX)
    pep_row = compo_df[compo_df["Sample_name"] == PEPTONE].iloc[0]

    ye_sub = compo_df[compo_df["Sample_name"] == YE_SAMPLE]
    if ye_sub.empty:
        print(f"  ⚠ YE sample '{YE_SAMPLE}' not found — Mode B skipped")
        ye_row = None
    else:
        ye_row = ye_sub.iloc[0]

    results = []
    for strain, media in STRAIN_MEDIA.items():
        r = verify_one(strain, media, basal_df, pep_row, ye_row, YE_G_PER_L)
        results.append(r)
        print_strain_result(r)

    # ── Summary ─────────────────────────────────────────
    print(f"\n\n{'#'*70}")
    print("#  SUMMARY")
    print('#'*70)

    ok = [r for r in results if r["status"] == "OK"]
    mus_sole = [r["mu_sole"] for r in ok]
    mus_ye = [r["mu_ye"] for r in ok if "mu_ye" in r]
    passing = [r for r in ok if r["mu_sole"] >= MU_MIN_PASS]

    print(f"\n  GEMs loaded: {len(ok)}/{len(results)}")
    print(f"  μ_sole ≥ {MU_MIN_PASS}: {len(passing)}/{len(results)}")

    if mus_sole:
        print(f"  μ_sole range: {min(mus_sole):.4f} ~ {max(mus_sole):.4f}  "
              f"(stddev={pstdev(mus_sole):.4f})")
    if mus_ye:
        print(f"  μ_+YE range:  {min(mus_ye):.4f} ~ {max(mus_ye):.4f}  "
              f"(stddev={pstdev(mus_ye):.4f})")

    # 검증 기준 체크 (Mode A 기준)
    print(f"\n  Criterion checks (Mode A - sole peptone):")
    c1 = len(passing) == len(results)
    c2 = bool(mus_sole) and pstdev(mus_sole) >= MU_STDDEV_MIN
    suspicious_found = []
    for r in ok:
        for sp in r.get("sp_sole", []):
            if sp["suspicious"]:
                suspicious_found.append((r["strain"], sp["name"]))
    c3 = len(suspicious_found) == 0
    c4 = all(r["status"] == "OK" for r in results)

    print(f"    1. 모든 균주 μ_sole ≥ {MU_MIN_PASS}         : {'PASS' if c1 else 'FAIL'}")
    print(f"    2. μ_sole 균주간 stddev ≥ {MU_STDDEV_MIN}     : {'PASS' if c2 else 'FAIL'}")
    print(f"    3. 의심 shadow price 없음                 : {'PASS' if c3 else 'FAIL'}")
    print(f"    4. 모든 GEM 로드·최적화 성공              : {'PASS' if c4 else 'FAIL'}")

    if suspicious_found:
        print(f"\n  ⚠ Suspicious bottlenecks:")
        for st, nm in suspicious_found:
            print(f"    {st}: {nm}")

    # 균주별 비교 한눈에
    print(f"\n  μ per strain  (sole → +YE 5 g/L):")
    header = f"    {'Strain':<5} {'Base':<4}  {'μ_sole':>8}  {'μ_+YE':>8}  {'Δμ':>7}  {'%':>6}  Recommend"
    print(header)
    print(f"    {'-'*5} {'-'*4}  {'-'*8}  {'-'*8}  {'-'*7}  {'-'*6}  ---------")
    for r in results:
        if r["status"] != "OK":
            print(f"    {r['strain']:<5} {r['media']:<4}  {r['status']}")
            continue
        ms = r["mu_sole"]
        if "mu_ye" not in r:
            print(f"    {r['strain']:<5} {r['media']:<4}  {ms:>8.4f}")
            continue
        my = r["mu_ye"]
        dmu = r["delta_mu"]
        pct = r["improvement_pct"]
        # Recommendation: YE 추가 시 개선 10%↑ 이면 권장
        rec = "★ +YE 권장" if pct >= 10.0 else ("= 펩톤만" if pct < 1.0 else "○ 선택적")
        print(f"    {r['strain']:<5} {r['media']:<4}  {ms:>8.4f}  {my:>8.4f}  "
              f"{dmu:>+7.4f}  {pct:>+5.1f}%  {rec}")

    all_pass = c1 and c2 and c3 and c4
    print(f"\n  {'='*40}")
    print(f"  OVERALL: {'✅ READY FOR TIER 1' if all_pass else '❌ FIX GEMs FIRST'}")
    print(f"  {'='*40}\n")


if __name__ == "__main__":
    main()

"""16 균주 × 7 펩톤 FBA 예측 매트릭스 + 인터랙티브 헤트맵.

실행:
    python scripts/fba_heatmap.py

출력:
    outputs/fba_heatmap_data.csv           — 전체 결과 (224행)
    outputs/fba_heatmap_sole.html          — 펩톤만 (인터랙티브)
    outputs/fba_heatmap_with_ye.html       — +YE 5g/L
    outputs/fba_heatmap_delta.html         — YE 효과 (Δμ)

CSV 컬럼:
    strain, media, peptone, peptone_g_per_L,
    mu_sole, mu_with_ye, delta_mu, pct_improvement,
    top_bottleneck_sole, top_bottleneck_with_ye
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from peptomatch.fba_simulator import (
    FBASimulator,
    basal_medium_to_mmol_L,
    peptone_to_mmol_L,
)

# ── 설정 ──────────────────────────────────────────────────
GEM_DIR = ROOT / "outputs" / "gem_cache"
BASAL_XLSX = ROOT / "data" / "basal_media_composition.xlsx"
COMP_XLSX = ROOT / "data" / "composition_template.xlsx"
OUT_DIR = ROOT / "outputs"
OUT_DIR.mkdir(exist_ok=True)

PEPTONE_G_PER_L = 25.0   # 실제 실험: 총 질소원 25 g/L를 단일 펩톤으로 교체
YE_NAME = "Gistex LS Ferm"
YE_G_PER_L = 5.0

# 균주 → 배지 매핑
STRAIN_MEDIA = {
    "LP":  "MRS", "LR":  "MRS", "LA":  "MRS", "LC":  "MRS",
    "LPC": "MRS", "LS":  "MRS", "STT": "MRS", "LRE": "MRS",
    "LB":  "MRS", "LF":  "MRS", "LH":  "MRS", "LG":  "MRS",
    "EF":  "BHI",
    "BS":  "TSB", "BC":  "TSB",
    "EC":  "LB",
}

PEPTONES = ["SOY-1", "SOY-N+", "SOY-L", "SOY-B", "RICE-1", "WHEAT-1", "PEA-1"]


# ── 실행 ──────────────────────────────────────────────────
def top_bottleneck(sim: FBASimulator) -> str:
    sp = sim.get_shadow_prices(top_n=1)
    if not sp:
        return "-"
    return f"{sp[0]['metabolite']} ({sp[0]['shadow_price']:.2f})"


def run_one(strain: str, media: str, peptone: str,
            basal_df: pd.DataFrame, comp_df: pd.DataFrame,
            ye_row) -> dict:
    gem = GEM_DIR / f"{strain}.xml"
    if not gem.exists():
        return {"strain": strain, "media": media, "peptone": peptone,
                "mu_sole": None, "mu_with_ye": None,
                "error": "GEM missing"}

    sim = FBASimulator(str(gem))
    basal = basal_medium_to_mmol_L(basal_df, media_id=media)

    pep_row = comp_df[comp_df["Sample_name"] == peptone]
    if pep_row.empty:
        return {"strain": strain, "media": media, "peptone": peptone,
                "mu_sole": None, "mu_with_ye": None,
                "error": "peptone missing"}
    pep_mmol = peptone_to_mmol_L(pep_row.iloc[0], PEPTONE_G_PER_L)

    # Mode A: 펩톤만
    sim.set_medium(basal, pep_mmol)
    mu_sole = sim.predict_growth()
    bn_sole = top_bottleneck(sim)

    # Mode B: +Gistex
    ye_mmol = peptone_to_mmol_L(ye_row, YE_G_PER_L)
    sim.set_medium(basal, pep_mmol, supplements_mmol=[ye_mmol])
    mu_ye = sim.predict_growth()
    bn_ye = top_bottleneck(sim)

    delta = mu_ye - mu_sole
    pct = (delta / mu_sole * 100) if mu_sole > 1e-6 else 0.0

    return {
        "strain": strain, "media": media, "peptone": peptone,
        "peptone_g_per_L": PEPTONE_G_PER_L,
        "mu_sole": round(mu_sole, 4),
        "mu_with_ye": round(mu_ye, 4),
        "delta_mu": round(delta, 4),
        "pct_improvement": round(pct, 1),
        "top_bottleneck_sole": bn_sole,
        "top_bottleneck_with_ye": bn_ye,
    }


def main() -> None:
    basal_df = pd.read_excel(BASAL_XLSX)
    comp_df = pd.read_excel(COMP_XLSX)

    ye_row_df = comp_df[comp_df["Sample_name"] == YE_NAME]
    if ye_row_df.empty:
        print(f"[FATAL] '{YE_NAME}' not found in composition_template.")
        sys.exit(1)
    ye_row = ye_row_df.iloc[0]

    rows = []
    total = len(STRAIN_MEDIA) * len(PEPTONES)
    n = 0
    for strain, media in STRAIN_MEDIA.items():
        for peptone in PEPTONES:
            n += 1
            r = run_one(strain, media, peptone, basal_df, comp_df, ye_row)
            rows.append(r)
            mu_s = r.get("mu_sole")
            mu_y = r.get("mu_with_ye")
            if mu_s is None:
                print(f"  [{n:3d}/{total}] {strain:4s} / {peptone:8s}  "
                      f"ERROR: {r.get('error')}")
            else:
                print(f"  [{n:3d}/{total}] {strain:4s} / {peptone:8s}  "
                      f"mu_sole={mu_s:.3f}  mu_YE={mu_y:.3f}  "
                      f"Δ={r['pct_improvement']:+.1f}%")

    df = pd.DataFrame(rows)
    csv_path = OUT_DIR / "fba_heatmap_data.csv"
    df.to_csv(csv_path, index=False, encoding="utf-8-sig")
    print(f"\n[OK] CSV: {csv_path}")

    # ── Heatmaps ─────────────────────────────────────────
    strains = list(STRAIN_MEDIA.keys())

    def make_matrix(col: str) -> list[list[float]]:
        mat = []
        for s in strains:
            row = []
            for p in PEPTONES:
                sub = df[(df["strain"] == s) & (df["peptone"] == p)]
                val = sub.iloc[0][col] if not sub.empty else None
                row.append(val)
            mat.append(row)
        return mat

    def make_hover(col_mu: str, col_bn: str) -> list[list[str]]:
        mat = []
        for s in strains:
            row = []
            for p in PEPTONES:
                sub = df[(df["strain"] == s) & (df["peptone"] == p)]
                if sub.empty:
                    row.append("N/A")
                else:
                    r = sub.iloc[0]
                    row.append(
                        f"<b>{s} × {p}</b><br>"
                        f"배지: {r['media']}<br>"
                        f"μ: {r[col_mu]:.4f} h⁻¹<br>"
                        f"Bottleneck: {r[col_bn]}"
                    )
            mat.append(row)
        return mat

    mat_sole = make_matrix("mu_sole")
    mat_ye = make_matrix("mu_with_ye")
    mat_delta = make_matrix("pct_improvement")

    hover_sole = make_hover("mu_sole", "top_bottleneck_sole")
    hover_ye = make_hover("mu_with_ye", "top_bottleneck_with_ye")

    y_labels = [f"{s} ({STRAIN_MEDIA[s]})" for s in strains]

    # (1) 펩톤만
    fig1 = go.Figure(data=go.Heatmap(
        z=mat_sole, x=PEPTONES, y=y_labels,
        text=[[f"{v:.3f}" if v is not None else "" for v in row] for row in mat_sole],
        texttemplate="%{text}",
        customdata=hover_sole,
        hovertemplate="%{customdata}<extra></extra>",
        colorscale="Viridis",
        colorbar=dict(title="μ (h⁻¹)"),
    ))
    fig1.update_layout(
        title=f"FBA 예측 μ — 펩톤만 ({PEPTONE_G_PER_L} g/L, single source)",
        xaxis_title="펩톤",
        yaxis_title="균주 (배지)",
        yaxis=dict(autorange="reversed"),
        height=700, width=1000,
    )
    p1 = OUT_DIR / "fba_heatmap_sole.html"
    fig1.write_html(p1)
    print(f"[OK] HTML: {p1}")

    # (2) +YE 5g/L
    fig2 = go.Figure(data=go.Heatmap(
        z=mat_ye, x=PEPTONES, y=y_labels,
        text=[[f"{v:.3f}" if v is not None else "" for v in row] for row in mat_ye],
        texttemplate="%{text}",
        customdata=hover_ye,
        hovertemplate="%{customdata}<extra></extra>",
        colorscale="Viridis",
        colorbar=dict(title="μ (h⁻¹)"),
    ))
    fig2.update_layout(
        title=f"FBA 예측 μ — 펩톤 + {YE_NAME} {YE_G_PER_L} g/L",
        xaxis_title="펩톤",
        yaxis_title="균주 (배지)",
        yaxis=dict(autorange="reversed"),
        height=700, width=1000,
    )
    p2 = OUT_DIR / "fba_heatmap_with_ye.html"
    fig2.write_html(p2)
    print(f"[OK] HTML: {p2}")

    # (3) YE 효과 (%)
    fig3 = go.Figure(data=go.Heatmap(
        z=mat_delta, x=PEPTONES, y=y_labels,
        text=[[f"{v:+.1f}%" if v is not None else "" for v in row] for row in mat_delta],
        texttemplate="%{text}",
        colorscale="RdYlGn",
        zmid=0,
        colorbar=dict(title="μ 개선 (%)"),
    ))
    fig3.update_layout(
        title=f"YE 보충 효과 (%) — Mode B vs Mode A",
        xaxis_title="펩톤",
        yaxis_title="균주 (배지)",
        yaxis=dict(autorange="reversed"),
        height=700, width=1000,
    )
    p3 = OUT_DIR / "fba_heatmap_delta.html"
    fig3.write_html(p3)
    print(f"[OK] HTML: {p3}")

    # ── 요약 출력 ────────────────────────────────────────
    print("\n" + "=" * 70)
    print("  SUMMARY — Mode A (펩톤만)")
    print("=" * 70)
    pivot = df.pivot(index="strain", columns="peptone", values="mu_sole")
    pivot = pivot.reindex(index=strains, columns=PEPTONES)
    print(pivot.to_string(float_format=lambda v: f"{v:.3f}" if pd.notna(v) else "  -  "))

    print("\n  펩톤별 최고 성능 균주:")
    for p in PEPTONES:
        best = df[df["peptone"] == p].nlargest(1, "mu_sole")
        if not best.empty:
            r = best.iloc[0]
            print(f"    {p:8s}: {r['strain']} ({r['mu_sole']:.3f})")

    print("\n  균주별 최적 펩톤:")
    for s in strains:
        best = df[df["strain"] == s].nlargest(1, "mu_sole")
        if not best.empty:
            r = best.iloc[0]
            print(f"    {s:4s}: {r['peptone']:8s} (μ={r['mu_sole']:.3f})")


if __name__ == "__main__":
    main()

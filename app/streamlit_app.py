"""PeptoMatch - Streamlit Web Application."""

import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from peptomatch.utils import load_config
from peptomatch.io_loaders import load_composition_data
from peptomatch.composition_features import CompositionFeatureExtractor
from peptomatch.genome_prior import GenomePriorBuilder
from peptomatch.scoring import PeptoneRecommender
from peptomatch.explain import RecommendationExplainer
from peptomatch.strain_db import StrainDB
from peptomatch.ncbi_client import NCBIClient
from peptomatch.taxonomy_priors import list_supported_genera

st.set_page_config(
    page_title="PeptoMatch",
    page_icon="🧬",
    layout="wide",
)


@st.cache_resource
def get_config():
    try:
        return load_config(PROJECT_ROOT / "config" / "config.yaml")
    except Exception:
        return load_config()


@st.cache_resource
def get_composition_df():
    config = get_config()
    return load_composition_data(
        Path(config["data"]["composition_file"]),
        sheet_name=config["data"].get("composition_sheet", "data"),
    )


@st.cache_resource
def get_strain_db():
    config = get_config()
    db = StrainDB(PROJECT_ROOT / "data" / "strains.db")
    if db.count() == 0:
        strain_path = Path(config["data"]["strain_file"])
        if strain_path.exists():
            db.load_from_excel(strain_path)
    return db


@st.cache_resource
def get_ncbi_client():
    config = get_config()
    ncbi_cfg = config.get("ncbi", {})
    return NCBIClient(
        email=ncbi_cfg.get("email", "user@example.com"),
        cache_dir=PROJECT_ROOT / "data" / "ncbi_cache",
    )


def main():
    st.title("🧬 PeptoMatch")
    st.caption("Genome-driven Peptone Recommendation System")

    tab1, tab2, tab3 = st.tabs(["🎯 펩톤 추천", "🔍 균주 브라우저", "📊 펩톤 탐색기"])

    # ── Tab 1: Peptone Recommendation ─────────────────────────────────
    with tab1:
        render_recommend_tab()

    # ── Tab 2: Strain Browser ─────────────────────────────────────────
    with tab2:
        render_strain_browser_tab()

    # ── Tab 3: Peptone Explorer ───────────────────────────────────────
    with tab3:
        render_peptone_explorer_tab()


def render_recommend_tab():
    config = get_config()
    comp_df = get_composition_df()
    db = get_strain_db()
    strain_df = db.get_strain_df()

    if strain_df.empty:
        st.warning("균주 데이터가 없습니다. '균주 브라우저' 탭에서 균주를 추가해주세요.")
        return

    col1, col2 = st.columns([2, 1])

    with col1:
        # Strain selection
        strain_options = {
            f"{r['strain_id']}: {r['full_name']}": r["strain_id"]
            for _, r in strain_df.iterrows()
        }
        selected = st.selectbox("균주 선택", options=list(strain_options.keys()))
        strain_id = strain_options[selected]

    with col2:
        top_k = st.slider("추천 수", 3, 16, 10)
        language = st.radio("언어", ["ko", "en"], horizontal=True)
        sempio_only = st.checkbox("Sempio 소재만", value=True)

    if st.button("🔬 분석 시작", type="primary", use_container_width=True):
        with st.spinner("유전체 기반 분석 중..."):
            recommender = PeptoneRecommender(comp_df, strain_df, config)
            pf = config.get("peptone_filter") if sempio_only else None
            recommendations = recommender.recommend(strain_id, top_k=top_k, peptone_filter=pf)

            explainer = RecommendationExplainer(comp_df, strain_df, config, language=language)
            recommendations = explainer.explain_batch(strain_id, recommendations, top_n_reasons=3)

            # Strain summary
            summary = explainer.get_strain_summary(strain_id)

        # Display results
        st.subheader("균주 영양 프로파일")
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("균주", summary.get("name", ""))
        c2.metric("정보 출처", summary.get("prior_source", ""))
        c3.metric("수송체 활성", summary.get("transporter_level", ""))
        c4.metric("뉴클레오타이드 합성", summary.get("nucleotide_synthesis", ""))

        if summary.get("key_aa_deficiencies"):
            st.info(f"**주요 아미노산 결핍**: {', '.join(summary['key_aa_deficiencies'])}")
        if summary.get("vitamin_deficiencies"):
            st.info(f"**비타민 결핍**: {', '.join(summary['vitamin_deficiencies'])}")

        st.subheader(f"Top-{top_k} 펩톤 추천")

        # Score bar chart
        fig = px.bar(
            recommendations,
            x="peptone", y="score",
            color="score",
            color_continuous_scale="Viridis",
            labels={"peptone": "펩톤", "score": "매칭 점수"},
        )
        fig.update_layout(height=400, showlegend=False)
        st.plotly_chart(fig, use_container_width=True)

        # Detailed table
        display_df = recommendations[["rank", "peptone", "score", "explanation"]].copy()
        display_df.columns = ["순위", "펩톤", "점수", "추천 이유"]
        st.dataframe(display_df, use_container_width=True, hide_index=True)

        # Download
        csv = recommendations.to_csv(index=False, encoding="utf-8-sig")
        st.download_button("📥 결과 다운로드 (CSV)", csv,
                           f"peptomatch_strain{strain_id}.csv", "text/csv")


def render_strain_browser_tab():
    config = get_config()
    db = get_strain_db()

    st.subheader("보유 균주 목록")
    strains = db.get_all()

    if strains:
        df = pd.DataFrame(strains)
        display_cols = ["strain_id", "genus", "species", "strain_name", "gcf_accession", "source"]
        available_cols = [c for c in display_cols if c in df.columns]
        st.dataframe(df[available_cols], use_container_width=True, hide_index=True)
        st.caption(f"총 {len(strains)}개 균주")
    else:
        st.info("등록된 균주가 없습니다.")

    st.divider()
    st.subheader("NCBI에서 균주 검색 및 추가")

    col1, col2, col3 = st.columns([2, 1, 1])
    with col1:
        taxon = st.text_input("분류군 검색", placeholder="예: Lactobacillaceae, Bifidobacterium, Bacillus subtilis")
    with col2:
        level = st.selectbox("Assembly Level", ["complete_genome", "chromosome", "scaffold"])
    with col3:
        limit = st.number_input("최대 결과", 5, 100, 20)

    if st.button("🔍 NCBI 검색") and taxon:
        client = get_ncbi_client()
        with st.spinner(f"NCBI에서 '{taxon}' 검색 중..."):
            results = client.search_assemblies(taxon=taxon, assembly_level=level, limit=limit)

        if results:
            results_df = pd.DataFrame(results)
            results_df["genome_size_mb"] = results_df["genome_size"] / 1_000_000
            st.dataframe(
                results_df[["accession", "organism_name", "assembly_level", "genome_size_mb"]],
                use_container_width=True, hide_index=True,
            )

            if st.button("📥 로컬 DB에 추가", type="primary"):
                added = db.expand_from_ncbi(taxon, client, assembly_level=level, max_results=limit)
                st.success(f"{len(added)}개 균주가 추가되었습니다!")
                st.rerun()
        else:
            st.warning("검색 결과가 없습니다.")

    st.divider()
    st.subheader("지원 Genus 목록 (Taxonomy Prior)")
    genera = list_supported_genera()
    st.write(", ".join(sorted(genera)))


def render_peptone_explorer_tab():
    comp_df = get_composition_df()
    config = get_config()
    extractor = CompositionFeatureExtractor(comp_df, config)
    features = extractor.compute_all_features()

    st.subheader("펩톤 성분 비교")

    # Filter option
    sempio_filter = config.get("peptone_filter", [])
    show_sempio_only = st.checkbox("Sempio 소재만 표시", value=True, key="explorer_filter")
    peptone_list = [p for p in features.index.tolist() if p in sempio_filter] if (show_sempio_only and sempio_filter) else features.index.tolist()

    selected_peptones = st.multiselect("비교할 펩톤 선택 (2-5개)", peptone_list, default=peptone_list[:3])

    if len(selected_peptones) >= 2:
        # FAA comparison
        faa_cols = [c for c in features.columns if c.startswith("faa_") and c != "faa_total" and c != "faa_mean"]
        if faa_cols:
            faa_data = features.loc[selected_peptones, faa_cols].T
            faa_data.index = [c.replace("faa_", "") for c in faa_data.index]

            melted = faa_data.reset_index().melt(id_vars="index")
            color_col = [c for c in melted.columns if c not in ("index", "value")][0]
            fig = px.bar(
                melted,
                x="index", y="value", color=color_col,
                barmode="group",
                labels={"index": "아미노산", "value": "함량", color_col: "펩톤"},
                title="유리 아미노산 (FAA) 비교",
            )
            fig.update_layout(height=500)
            st.plotly_chart(fig, use_container_width=True)

        # Summary table
        summary_cols = ["faa_total", "taa_total", "mw_avg", "mw_pct_low", "vitamin_total", "nucleotide_total"]
        available = [c for c in summary_cols if c in features.columns]
        if available:
            st.subheader("요약 비교")
            summary = features.loc[selected_peptones, available].copy()
            summary.columns = [c.replace("_", " ").title() for c in summary.columns]
            st.dataframe(summary, use_container_width=True)

    elif len(selected_peptones) == 1:
        st.info("비교를 위해 2개 이상의 펩톤을 선택해주세요.")

    # Full composition table
    st.divider()
    st.subheader("전체 펩톤 성분 데이터")

    cols_to_show = ["faa_total", "taa_total", "mw_avg", "mw_pct_low"]
    available = [c for c in cols_to_show if c in features.columns]
    if available:
        st.dataframe(features[available].sort_values(available[0], ascending=False),
                     use_container_width=True)


if __name__ == "__main__":
    main()

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
from peptomatch.compare import StrainComparator
from peptomatch.kegg_viz import KEGGVisualizer

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

    tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
        "🎯 펩톤 추천",
        "🔍 균주 브라우저",
        "📊 펩톤 탐색기",
        "⚖️ 균주 비교",
        "🧪 KEGG 경로",
        "📄 PDF 리포트",
    ])

    with tab1:
        render_recommend_tab()
    with tab2:
        render_strain_browser_tab()
    with tab3:
        render_peptone_explorer_tab()
    with tab4:
        render_compare_tab()
    with tab5:
        render_kegg_tab()
    with tab6:
        render_pdf_tab()


# ── Tab 1: 펩톤 추천 ───────────────────────────────────────────────

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

    if st.button("🔬 분석 시작", type="primary", width="stretch"):
        with st.spinner("유전체 기반 분석 중..."):
            recommender = PeptoneRecommender(comp_df, strain_df, config)
            pf = config.get("peptone_filter") if sempio_only else None
            recommendations = recommender.recommend(strain_id, top_k=top_k, peptone_filter=pf)

            explainer = RecommendationExplainer(comp_df, strain_df, config, language=language)
            recommendations = explainer.explain_batch(strain_id, recommendations, top_n_reasons=3)
            summary = explainer.get_strain_summary(strain_id)

        # Store for PDF tab reuse
        st.session_state["last_recommendations"] = recommendations
        st.session_state["last_summary"] = summary
        st.session_state["last_strain_id"] = strain_id

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

        fig = px.bar(
            recommendations,
            x="peptone", y="score",
            color="score",
            color_continuous_scale="Viridis",
            labels={"peptone": "펩톤", "score": "매칭 점수"},
        )
        fig.update_layout(height=400, showlegend=False)
        st.plotly_chart(fig, use_container_width=True)

        display_df = recommendations[["rank", "peptone", "score", "explanation"]].copy()
        display_df.columns = ["순위", "펩톤", "점수", "추천 이유"]
        st.dataframe(display_df, width="stretch", hide_index=True)

        csv = recommendations.to_csv(index=False, encoding="utf-8-sig")
        st.download_button("📥 결과 다운로드 (CSV)", csv,
                           f"peptomatch_strain{strain_id}.csv", "text/csv")


# ── Tab 2: 균주 브라우저 ───────────────────────────────────────────

def render_strain_browser_tab():
    config = get_config()
    db = get_strain_db()

    st.subheader("보유 균주 목록")
    strains = db.get_all()

    if strains:
        df = pd.DataFrame(strains)
        display_cols = ["strain_id", "genus", "species", "strain_name", "gcf_accession", "source"]
        available_cols = [c for c in display_cols if c in df.columns]
        st.dataframe(df[available_cols], width="stretch", hide_index=True)
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
        st.session_state["ncbi_results"] = results

    results = st.session_state.get("ncbi_results", [])

    if results:
        st.markdown(f"**{len(results)}건 검색됨** — 추가할 균주를 선택하세요:")

        selected_indices = []
        for i, r in enumerate(results):
            already = db.conn.execute(
                "SELECT 1 FROM strains WHERE gcf_accession = ?", (r["accession"],)
            ).fetchone()
            label = f"`{r['accession']}` — {r['organism_name']} ({r['genome_size']/1_000_000:.2f} Mb)"
            if already:
                st.checkbox(label, value=False, disabled=True, key=f"ncbi_{i}",
                            help="이미 DB에 존재")
            else:
                if st.checkbox(label, value=True, key=f"ncbi_{i}"):
                    selected_indices.append(i)

        if selected_indices:
            if st.button(f"📥 선택한 {len(selected_indices)}개 균주 추가", type="primary"):
                selected_assemblies = [results[i] for i in selected_indices]
                added = db.add_assemblies(selected_assemblies)
                st.success(f"{len(added)}개 균주가 추가되었습니다!")
                del st.session_state["ncbi_results"]
                st.rerun()
    elif "ncbi_results" in st.session_state:
        st.warning("검색 결과가 없습니다.")

    st.divider()
    st.subheader("지원 Genus 목록 (Taxonomy Prior)")
    genera = list_supported_genera()
    st.write(", ".join(sorted(genera)))


# ── Tab 3: 펩톤 탐색기 ────────────────────────────────────────────

def render_peptone_explorer_tab():
    comp_df = get_composition_df()
    config = get_config()
    extractor = CompositionFeatureExtractor(comp_df, config)
    features = extractor.compute_all_features()

    st.subheader("펩톤 성분 비교")

    sempio_filter = config.get("peptone_filter", [])
    show_sempio_only = st.checkbox("Sempio 소재만 표시", value=True, key="explorer_filter")
    peptone_list = (
        [p for p in features.index.tolist() if p in sempio_filter]
        if (show_sempio_only and sempio_filter)
        else features.index.tolist()
    )

    selected_peptones = st.multiselect("비교할 펩톤 선택 (2-5개)", peptone_list, default=peptone_list[:3])

    if len(selected_peptones) >= 2:
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

        summary_cols = ["faa_total", "taa_total", "mw_avg", "mw_pct_low", "vitamin_total", "nucleotide_total"]
        available = [c for c in summary_cols if c in features.columns]
        if available:
            st.subheader("요약 비교")
            summary = features.loc[selected_peptones, available].copy()
            summary.columns = [c.replace("_", " ").title() for c in summary.columns]
            st.dataframe(summary, width="stretch")

    elif len(selected_peptones) == 1:
        st.info("비교를 위해 2개 이상의 펩톤을 선택해주세요.")

    st.divider()
    st.subheader("전체 펩톤 성분 데이터")

    cols_to_show = ["faa_total", "taa_total", "mw_avg", "mw_pct_low"]
    available = [c for c in cols_to_show if c in features.columns]
    if available:
        st.dataframe(features[available].sort_values(available[0], ascending=False),
                     width="stretch")


# ── Tab 4: 균주 비교 ──────────────────────────────────────────────

def render_compare_tab():
    config = get_config()
    comp_df = get_composition_df()
    db = get_strain_db()
    strain_df = db.get_strain_df()

    if strain_df.empty:
        st.warning("균주 데이터가 없습니다.")
        return

    st.subheader("⚖️ 균주 간 비교 분석")
    st.caption("2~4개 균주를 선택하여 영양 요구도와 펩톤 매칭을 비교합니다.")

    strain_options = {
        f"{r['strain_id']}: {r['full_name']}": r["strain_id"]
        for _, r in strain_df.iterrows()
    }

    selected_labels = st.multiselect(
        "비교할 균주 선택 (2-4개)",
        options=list(strain_options.keys()),
        max_selections=4,
        key="compare_strains",
    )

    if len(selected_labels) < 2:
        st.info("비교를 위해 2개 이상의 균주를 선택하세요.")
        return

    strain_ids = [strain_options[s] for s in selected_labels]
    sempio_only = st.checkbox("Sempio 소재만", value=True, key="compare_sempio")

    if st.button("📊 비교 분석", type="primary", width="stretch"):
        with st.spinner("비교 분석 중..."):
            comparator = StrainComparator(comp_df, strain_df, config)

            # 1. Radar chart
            st.subheader("아미노산 요구도 비교")
            radar_fig = comparator.radar_chart(strain_ids)
            st.plotly_chart(radar_fig, use_container_width=True)

            # 2. Heatmap
            st.subheader("영양소 요구도 히트맵")
            heatmap_fig = comparator.heatmap_chart(strain_ids)
            st.plotly_chart(heatmap_fig, use_container_width=True)

            # 3. Recommendation comparison
            st.subheader("펩톤 추천 비교")
            pf = config.get("peptone_filter") if sempio_only else None
            recs_by_strain = comparator.compare_recommendations(
                strain_ids, top_k=5, peptone_filter=pf
            )

            score_fig = comparator.score_comparison_chart(recs_by_strain)
            st.plotly_chart(score_fig, use_container_width=True)

            # 4. Side-by-side tables
            cols = st.columns(len(strain_ids))
            for col, (label, recs) in zip(cols, recs_by_strain.items()):
                with col:
                    st.markdown(f"**{label}**")
                    display = recs[["rank", "peptone", "score"]].copy()
                    display.columns = ["순위", "펩톤", "점수"]
                    display["점수"] = display["점수"].round(1)
                    st.dataframe(display, width="stretch", hide_index=True)

            # 5. Demand profile table
            st.subheader("상세 요구도 비교")
            demand_df = comparator.compare_demand(strain_ids)
            # Show only AA demands
            aa_cols = [c for c in demand_df.columns if c.startswith("demand_") and "vitamin" not in c and "nucleotide" not in c]
            if aa_cols:
                display_demand = demand_df[aa_cols].copy()
                display_demand.columns = [c.replace("demand_", "") for c in display_demand.columns]
                st.dataframe(display_demand, width="stretch")


# ── Tab 5: KEGG 경로 ──────────────────────────────────────────────

def render_kegg_tab():
    config = get_config()
    db = get_strain_db()
    strain_df = db.get_strain_df()

    if strain_df.empty:
        st.warning("균주 데이터가 없습니다.")
        return

    st.subheader("🧪 KEGG 생합성 경로 시각화")
    st.caption("균주의 아미노산·비타민 생합성 경로 완성도를 시각화합니다.")

    strain_options = {
        f"{r['strain_id']}: {r['full_name']}": r["strain_id"]
        for _, r in strain_df.iterrows()
    }
    selected = st.selectbox("균주 선택", options=list(strain_options.keys()), key="kegg_strain")
    strain_id = strain_options[selected]

    if st.button("🔬 경로 분석", type="primary", key="kegg_btn", width="stretch"):
        with st.spinner("KEGG 경로 분석 중..."):
            viz = KEGGVisualizer(strain_df, config)

            # Overview heatmap
            st.subheader("경로 종합 Overview")
            overview_fig = viz.overview_chart(strain_id)
            st.plotly_chart(overview_fig, use_container_width=True)

            # AA biosynthesis bar chart
            st.subheader("아미노산 생합성 완성도")
            aa_fig = viz.aa_pathway_chart(strain_id)
            st.plotly_chart(aa_fig, use_container_width=True)

            # Vitamin chart
            st.subheader("비타민 생합성 완성도")
            vit_fig = viz.vitamin_chart(strain_id)
            st.plotly_chart(vit_fig, use_container_width=True)

            # Pathway detail flowcharts for most deficient AAs
            prior = viz.prior_builder.get_prior(strain_id)
            aa_synth = prior.get("aa_biosynthesis", {})
            if aa_synth:
                sorted_aa = sorted(aa_synth.items(), key=lambda x: x[1])
                worst_3 = [aa for aa, val in sorted_aa[:3] if val < 0.8]

                if worst_3:
                    st.subheader("🔬 주요 결핍 아미노산 경로 플로우차트")
                    st.caption("각 효소를 박스로 표시하고 화살표로 연결합니다. 초록=존재, 빨강=결핍, ★=핵심 효소")
                    for aa in worst_3:
                        try:
                            detail_fig = viz.pathway_detail_chart(strain_id, aa)
                            st.plotly_chart(detail_fig, use_container_width=True)
                        except Exception as e:
                            st.warning(f"{aa} 경로 상세를 표시할 수 없습니다: {e}")

            # Vitamin pathway flowcharts for deficient vitamins
            vit_synth = prior.get("vitamin_biosynthesis", {})
            if vit_synth:
                worst_vit = [v for v, val in sorted(vit_synth.items(), key=lambda x: x[1]) if val < 0.8][:3]
                if worst_vit:
                    st.subheader("🔬 주요 결핍 비타민 경로 플로우차트")
                    for vit in worst_vit:
                        try:
                            vit_detail = viz.pathway_detail_chart(strain_id, vit, pathway_source="vitamin")
                            st.plotly_chart(vit_detail, use_container_width=True)
                        except Exception as e:
                            st.warning(f"비타민 {vit} 경로 상세를 표시할 수 없습니다: {e}")

            # Interactive pathway explorer
            st.subheader("🗺️ 경로 탐색기")
            explore_col1, explore_col2 = st.columns(2)
            with explore_col1:
                aa_list = list(sorted(aa_synth.keys()))
                selected_aa = st.selectbox("아미노산 경로 선택", aa_list, key="explore_aa")
            with explore_col2:
                vit_list = list(sorted(vit_synth.keys())) if vit_synth else []
                selected_vit = st.selectbox("비타민 경로 선택", ["(선택 안 함)"] + vit_list, key="explore_vit")

            if selected_aa:
                try:
                    explore_fig = viz.pathway_detail_chart(strain_id, selected_aa)
                    st.plotly_chart(explore_fig, use_container_width=True)
                except Exception as e:
                    st.warning(f"경로를 표시할 수 없습니다: {e}")

            if selected_vit and selected_vit != "(선택 안 함)":
                try:
                    explore_vit_fig = viz.pathway_detail_chart(strain_id, selected_vit, pathway_source="vitamin")
                    st.plotly_chart(explore_vit_fig, use_container_width=True)
                except Exception as e:
                    st.warning(f"비타민 경로를 표시할 수 없습니다: {e}")


# ── Tab 6: PDF 리포트 ─────────────────────────────────────────────

def render_pdf_tab():
    config = get_config()
    comp_df = get_composition_df()
    db = get_strain_db()
    strain_df = db.get_strain_df()

    if strain_df.empty:
        st.warning("균주 데이터가 없습니다.")
        return

    st.subheader("📄 PDF 리포트 생성")
    st.caption("분석 결과를 PDF 문서로 다운로드합니다.")

    # Check if we have cached results from the recommend tab
    has_cached = "last_recommendations" in st.session_state

    strain_options = {
        f"{r['strain_id']}: {r['full_name']}": r["strain_id"]
        for _, r in strain_df.iterrows()
    }

    col1, col2 = st.columns([2, 1])
    with col1:
        selected = st.selectbox("균주 선택", options=list(strain_options.keys()), key="pdf_strain")
        strain_id = strain_options[selected]
    with col2:
        top_k = st.slider("추천 수", 3, 16, 10, key="pdf_topk")
        sempio_only = st.checkbox("Sempio 소재만", value=True, key="pdf_sempio")
        include_charts = st.checkbox("차트 이미지 포함", value=True, key="pdf_charts")

    if has_cached:
        cached_sid = st.session_state.get("last_strain_id")
        if cached_sid == strain_id:
            st.success("💡 '펩톤 추천' 탭의 분석 결과를 재사용합니다.")

    if st.button("📄 PDF 생성", type="primary", width="stretch"):
        with st.spinner("PDF 리포트 생성 중..."):
            # Generate recommendations (or reuse cached)
            if has_cached and st.session_state.get("last_strain_id") == strain_id:
                recommendations = st.session_state["last_recommendations"]
                summary = st.session_state["last_summary"]
            else:
                recommender = PeptoneRecommender(comp_df, strain_df, config)
                pf = config.get("peptone_filter") if sempio_only else None
                recommendations = recommender.recommend(strain_id, top_k=top_k, peptone_filter=pf)

                explainer = RecommendationExplainer(comp_df, strain_df, config, language="ko")
                recommendations = explainer.explain_batch(strain_id, recommendations, top_n_reasons=3)
                summary = explainer.get_strain_summary(strain_id)

            # Generate chart images
            charts = {}
            if include_charts:
                try:
                    # Score bar chart
                    fig = px.bar(
                        recommendations, x="peptone", y="score",
                        color="score", color_continuous_scale="Viridis",
                    )
                    fig.update_layout(height=350, showlegend=False, title="Peptone Match Scores")
                    charts["Peptone Match Scores"] = fig.to_image(format="png", width=900, height=350)
                except Exception:
                    pass  # kaleido might not be installed

                try:
                    viz = KEGGVisualizer(strain_df, config)
                    aa_fig = viz.aa_pathway_chart(strain_id)
                    aa_fig.update_layout(height=350)
                    charts["AA Biosynthesis Completeness"] = aa_fig.to_image(format="png", width=900, height=350)
                except Exception:
                    pass

            # Generate PDF
            from peptomatch.report_pdf import ReportGenerator
            generator = ReportGenerator(comp_df, strain_df, config)
            pdf_bytes = generator.generate(
                strain_id=strain_id,
                recommendations=recommendations,
                summary=summary,
                charts=charts,
            )

            strain_name = summary.get("name", f"strain{strain_id}").replace(" ", "_")
            filename = f"PeptoMatch_{strain_name}.pdf"

            st.download_button(
                label="📥 PDF 다운로드",
                data=pdf_bytes,
                file_name=filename,
                mime="application/pdf",
                type="primary",
                width="stretch",
            )

            st.success(f"PDF 리포트가 준비되었습니다! ({len(pdf_bytes)/1024:.1f} KB)")


if __name__ == "__main__":
    main()

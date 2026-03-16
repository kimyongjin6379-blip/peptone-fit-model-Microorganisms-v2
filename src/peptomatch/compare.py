"""Strain comparison module — side-by-side demand profile analysis."""

from typing import Any, Optional

import pandas as pd
import numpy as np
import plotly.graph_objects as go

from .genome_prior import GenomePriorBuilder
from .scoring import PeptoneRecommender


class StrainComparator:
    """Compare multiple strains on demand profiles and peptone recommendations."""

    def __init__(
        self,
        composition_df: pd.DataFrame,
        strain_df: pd.DataFrame,
        config: Optional[dict] = None,
    ):
        self.composition_df = composition_df
        self.strain_df = strain_df
        self.config = config or {}
        self.prior_builder = GenomePriorBuilder(strain_df, config)

    # ── public API ──────────────────────────────────────────────────

    def compare_demand(self, strain_ids: list[int]) -> pd.DataFrame:
        """Build a DataFrame of demand scores for given strains.

        Returns:
            DataFrame (strains × demand features)
        """
        rows = []
        for sid in strain_ids:
            demand = self.prior_builder.get_demand_scores(sid)
            row = self._strain_label(sid)
            row.update(demand)
            rows.append(row)
        return pd.DataFrame(rows).set_index("label")

    def compare_recommendations(
        self,
        strain_ids: list[int],
        top_k: int = 5,
        peptone_filter: Optional[list[str]] = None,
    ) -> dict[str, pd.DataFrame]:
        """Get top-K recommendations per strain for side-by-side comparison."""
        recommender = PeptoneRecommender(
            self.composition_df, self.strain_df, self.config
        )
        result = {}
        for sid in strain_ids:
            pf = peptone_filter
            if pf is None:
                pf = self.config.get("peptone_filter")
            recs = recommender.recommend(sid, top_k=top_k, peptone_filter=pf)
            label = self._strain_label(sid)["label"]
            result[label] = recs
        return result

    # ── plotly charts ───────────────────────────────────────────────

    def radar_chart(self, strain_ids: list[int]) -> go.Figure:
        """Create a radar chart comparing AA demand profiles."""
        aa_list = [
            "His", "Ile", "Leu", "Lys", "Met", "Phe", "Thr", "Trp", "Val",
            "Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "Pro", "Ser", "Tyr",
        ]
        fig = go.Figure()

        for sid in strain_ids:
            demand = self.prior_builder.get_demand_scores(sid)
            label = self._strain_label(sid)["label"]
            values = [demand.get(f"demand_{aa}", 0) for aa in aa_list]
            # close the polygon
            values_closed = values + [values[0]]
            aa_closed = aa_list + [aa_list[0]]

            fig.add_trace(go.Scatterpolar(
                r=values_closed,
                theta=aa_closed,
                name=label,
                fill="toself",
                opacity=0.55,
            ))

        fig.update_layout(
            polar=dict(radialaxis=dict(visible=True, range=[0, 1])),
            title="아미노산 요구도 (Demand) 비교",
            height=550,
        )
        return fig

    def heatmap_chart(self, strain_ids: list[int]) -> go.Figure:
        """Heatmap of demand scores (AA + vitamins)."""
        aa_list = [
            "His", "Ile", "Leu", "Lys", "Met", "Phe", "Thr", "Trp", "Val",
            "Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "Pro", "Ser", "Tyr",
        ]
        vit_list = ["B1", "B2", "B3", "B6", "B9"]
        features = aa_list + [f"Vit {v}" for v in vit_list] + ["Nucleotide"]

        labels = []
        z_values = []

        for sid in strain_ids:
            demand = self.prior_builder.get_demand_scores(sid)
            label = self._strain_label(sid)["label"]
            labels.append(label)

            row = [demand.get(f"demand_{aa}", 0) for aa in aa_list]
            row += [demand.get(f"demand_vitamin_{v}", 0) for v in vit_list]
            row += [demand.get("demand_nucleotide", 0)]
            z_values.append(row)

        fig = go.Figure(data=go.Heatmap(
            z=z_values,
            x=features,
            y=labels,
            colorscale="RdYlGn_r",  # red=high demand, green=low
            zmin=0, zmax=1,
            colorbar_title="Demand",
            text=[[f"{v:.0%}" for v in row] for row in z_values],
            texttemplate="%{text}",
        ))
        fig.update_layout(
            title="영양소 요구도 히트맵",
            height=max(300, 120 * len(strain_ids)),
            xaxis_title="영양소",
            yaxis_title="균주",
        )
        return fig

    def score_comparison_chart(
        self,
        recs_by_strain: dict[str, pd.DataFrame],
    ) -> go.Figure:
        """Grouped bar chart of peptone scores across strains."""
        # Collect all peptones that appear
        all_peptones = set()
        for df in recs_by_strain.values():
            all_peptones.update(df["peptone"].tolist())
        all_peptones = sorted(all_peptones)

        fig = go.Figure()
        for label, df in recs_by_strain.items():
            score_map = dict(zip(df["peptone"], df["score"]))
            scores = [score_map.get(p, 0) for p in all_peptones]
            fig.add_trace(go.Bar(name=label, x=all_peptones, y=scores))

        fig.update_layout(
            barmode="group",
            title="펩톤 매칭 점수 비교",
            xaxis_title="펩톤",
            yaxis_title="매칭 점수",
            height=450,
        )
        return fig

    # ── helpers ──────────────────────────────────────────────────────

    def _strain_label(self, strain_id: int) -> dict[str, Any]:
        row = self.strain_df[self.strain_df["strain_id"] == strain_id]
        if row.empty:
            return {"strain_id": strain_id, "label": f"Strain {strain_id}"}
        r = row.iloc[0]
        genus = r.get("genus", "")
        species = r.get("species", "")
        strain = r.get("strain_name", "")
        short = f"{genus} {species}"
        if strain:
            short += f" {strain}"
        return {"strain_id": strain_id, "label": short.strip()}

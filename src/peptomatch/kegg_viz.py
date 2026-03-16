"""KEGG pathway visualization — biosynthesis completeness maps."""

from typing import Any, Optional

import plotly.graph_objects as go
import plotly.express as px

from .kegg_pathway import (
    AA_BIOSYNTHESIS_KOS,
    VITAMIN_BIOSYNTHESIS_KOS,
    NUCLEOTIDE_BIOSYNTHESIS_KOS,
    TRANSPORTER_KOS,
)
from .genome_prior import GenomePriorBuilder

import pandas as pd


class KEGGVisualizer:
    """Visualize KEGG pathway completeness for strains."""

    def __init__(
        self,
        strain_df: pd.DataFrame,
        config: Optional[dict] = None,
    ):
        self.strain_df = strain_df
        self.config = config or {}
        self.prior_builder = GenomePriorBuilder(strain_df, config)

    # ── AA biosynthesis pathway map ─────────────────────────────────

    def aa_pathway_chart(self, strain_id: int) -> go.Figure:
        """Stepped bar chart of AA biosynthesis completeness per amino acid.

        Colors: green (≥ 80%), yellow (40-79%), red (<40%).
        """
        prior = self.prior_builder.get_prior(strain_id)
        aa_synth = prior.get("aa_biosynthesis", {})
        label = self._label(strain_id)

        # Sort: essential first, then non-essential
        essential = ["His", "Ile", "Leu", "Lys", "Met", "Phe", "Thr", "Trp", "Val"]
        non_essential = ["Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "Pro", "Ser", "Tyr"]
        ordered = essential + non_essential

        aa_names = []
        completeness_vals = []
        colors = []

        for aa in ordered:
            comp = aa_synth.get(aa, 0.5)
            aa_names.append(aa)
            completeness_vals.append(comp)
            if comp >= 0.8:
                colors.append("#2ecc71")   # green
            elif comp >= 0.4:
                colors.append("#f1c40f")   # yellow
            else:
                colors.append("#e74c3c")   # red

        fig = go.Figure(data=[
            go.Bar(
                x=aa_names,
                y=completeness_vals,
                marker_color=colors,
                text=[f"{v:.0%}" for v in completeness_vals],
                textposition="outside",
            )
        ])

        fig.add_shape(
            type="line", x0=-0.5, x1=len(essential) - 0.5,
            y0=-0.06, y1=-0.06, line=dict(color="royalblue", width=4),
        )
        fig.add_annotation(
            x=len(essential) / 2 - 0.5, y=-0.12,
            text="필수 아미노산", showarrow=False,
            font=dict(size=11, color="royalblue"),
        )
        fig.add_shape(
            type="line", x0=len(essential) - 0.5, x1=len(ordered) - 0.5,
            y0=-0.06, y1=-0.06, line=dict(color="grey", width=4),
        )
        fig.add_annotation(
            x=len(essential) + len(non_essential) / 2 - 0.5, y=-0.12,
            text="비필수 아미노산", showarrow=False,
            font=dict(size=11, color="grey"),
        )

        fig.update_layout(
            title=f"아미노산 생합성 경로 완성도 — {label}",
            yaxis=dict(title="완성도", range=[-0.2, 1.15], tickformat=".0%"),
            xaxis_title="아미노산",
            height=480,
            showlegend=False,
        )
        return fig

    # ── Vitamin / nucleotide / transporter summary ──────────────────

    def vitamin_chart(self, strain_id: int) -> go.Figure:
        """Horizontal bar chart for vitamin biosynthesis completeness."""
        prior = self.prior_builder.get_prior(strain_id)
        vit_synth = prior.get("vitamin_biosynthesis", {})
        label = self._label(strain_id)

        vitamins = sorted(vit_synth.keys())
        vals = [vit_synth.get(v, 0) for v in vitamins]
        colors = ["#2ecc71" if v >= 0.8 else "#f1c40f" if v >= 0.4 else "#e74c3c" for v in vals]
        vit_labels = [f"비타민 {v}" for v in vitamins]

        fig = go.Figure(data=[
            go.Bar(
                y=vit_labels, x=vals,
                orientation="h",
                marker_color=colors,
                text=[f"{v:.0%}" for v in vals],
                textposition="outside",
            )
        ])
        fig.update_layout(
            title=f"비타민 생합성 경로 — {label}",
            xaxis=dict(title="완성도", range=[0, 1.15], tickformat=".0%"),
            height=300,
            showlegend=False,
        )
        return fig

    # ── Detailed pathway step view ──────────────────────────────────

    def pathway_detail_chart(self, strain_id: int, aa: str) -> go.Figure:
        """Show individual KO steps for a specific AA biosynthesis pathway.

        Each step (KO) is shown as a colored box: green=found, red=missing.
        """
        prior = self.prior_builder.get_prior(strain_id)
        label = self._label(strain_id)

        if aa not in AA_BIOSYNTHESIS_KOS:
            raise ValueError(f"Unknown amino acid: {aa}")

        pathway = AA_BIOSYNTHESIS_KOS[aa]
        kos = pathway["kos"]
        essential = set(pathway.get("essential_kos", []))

        # Determine which KOs are found (from prior source)
        # We re-derive from the prior's source info
        found_kos = self._get_found_kos(strain_id)

        step_names = []
        found_vals = []
        colors = []
        annotations = []

        for ko in kos:
            # Extract enzyme short name from the AA_BIOSYNTHESIS_KOS comments
            name = ko
            is_found = ko in found_kos
            is_essential = ko in essential

            step_names.append(name)
            found_vals.append(1 if is_found else 0)
            if is_found:
                colors.append("#2ecc71")  # green
            else:
                colors.append("#e74c3c")  # red

            label_text = "✓" if is_found else "✗"
            if is_essential:
                label_text += " ★"
            annotations.append(label_text)

        fig = go.Figure(data=[
            go.Bar(
                x=step_names,
                y=[1] * len(step_names),
                marker_color=colors,
                text=annotations,
                textposition="inside",
                textfont=dict(size=16, color="white"),
            )
        ])

        completeness = prior.get("aa_biosynthesis", {}).get(aa, 0)

        fig.update_layout(
            title=f"{pathway['name']} (M{pathway['pathway'][1:]}) — 완성도: {completeness:.0%}",
            yaxis=dict(visible=False),
            xaxis_title="효소 (KO ID)",
            height=250,
            showlegend=False,
            annotations=[
                dict(x=0.01, y=1.15, xref="paper", yref="paper",
                     text="✓ = 존재 | ✗ = 결핍 | ★ = 핵심 효소",
                     showarrow=False, font=dict(size=11)),
            ],
        )
        return fig

    # ── Comprehensive overview ──────────────────────────────────────

    def overview_chart(self, strain_id: int) -> go.Figure:
        """Combined overview: AA + vitamins + nucleotide + transporter as a single heatmap row."""
        prior = self.prior_builder.get_prior(strain_id)
        aa_synth = prior.get("aa_biosynthesis", {})
        vit_synth = prior.get("vitamin_biosynthesis", {})
        nuc = prior.get("nucleotide_biosynthesis", 0.5)
        trans = prior.get("transporter_score", 0.5)

        categories = []
        values = []

        essential = ["His", "Ile", "Leu", "Lys", "Met", "Phe", "Thr", "Trp", "Val"]
        non_essential = ["Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "Pro", "Ser", "Tyr"]

        for aa in essential + non_essential:
            categories.append(aa)
            values.append(aa_synth.get(aa, 0.5))

        categories.append("")  # spacer
        values.append(None)

        for v in sorted(vit_synth.keys()):
            categories.append(f"B{v[1:]}" if v.startswith("B") else v)
            values.append(vit_synth.get(v, 0.5))

        categories.append("")
        values.append(None)
        categories.append("Nucleotide")
        values.append(nuc)
        categories.append("Transporter")
        values.append(trans)

        label = self._label(strain_id)

        fig = go.Figure(data=go.Heatmap(
            z=[values],
            x=categories,
            y=[label],
            colorscale="RdYlGn",
            zmin=0, zmax=1,
            colorbar_title="완성도",
            text=[[f"{v:.0%}" if v is not None else "" for v in values]],
            texttemplate="%{text}",
        ))
        fig.update_layout(
            title=f"생합성 경로 종합 — {label}",
            height=200,
            xaxis_title="경로",
        )
        return fig

    # ── helpers ──────────────────────────────────────────────────────

    def _label(self, strain_id: int) -> str:
        row = self.strain_df[self.strain_df["strain_id"] == strain_id]
        if row.empty:
            return f"Strain {strain_id}"
        r = row.iloc[0]
        return f"{r.get('genus', '')} {r.get('species', '')} {r.get('strain_name', '')}".strip()

    def _get_found_kos(self, strain_id: int) -> set[str]:
        """Try to retrieve the raw KO list for this strain from cache."""
        import json
        from pathlib import Path
        from .utils import ensure_output_dir

        prior = self.prior_builder.get_prior(strain_id)
        gcf = prior.get("gcf", "")
        if not gcf:
            return set()

        cache_dir = ensure_output_dir("genome_cache")
        # Check for KO list file
        ko_file = cache_dir / f"{gcf}_ko_list.json"
        if ko_file.exists():
            try:
                with open(ko_file, "r") as f:
                    return set(json.load(f))
            except Exception:
                pass

        # Fall back: reconstruct from pathway completeness
        # If a pathway is 100% complete, all KOs are present
        found = set()
        aa_synth = prior.get("aa_biosynthesis", {})
        for aa, comp in aa_synth.items():
            if aa in AA_BIOSYNTHESIS_KOS:
                pathway = AA_BIOSYNTHESIS_KOS[aa]
                essential_kos = pathway.get("essential_kos", [])
                # If completeness > 0, assume proportional KOs found
                n_found = int(round(comp * len(essential_kos)))
                for ko in essential_kos[:n_found]:
                    found.add(ko)
        return found

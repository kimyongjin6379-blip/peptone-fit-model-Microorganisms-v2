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

    # ── KO enzyme name lookup ────────────────────────────────────────

    # Maps KO IDs → short gene names (extracted from kegg_pathway.py comments)
    _KO_GENE_NAMES: dict[str, str] = {
        # Histidine
        "K00765": "hisG", "K02501": "hisI", "K01814": "hisA",
        "K02500": "hisH", "K01693": "hisB", "K00013": "hisD", "K01089": "hisC",
        # Isoleucine / Valine / Leucine
        "K01754": "ilvA", "K01652": "ilvB", "K01653": "ilvH",
        "K00053": "ilvC", "K01687": "ilvD", "K00826": "ilvE",
        "K01649": "leuA", "K01703": "leuC", "K01704": "leuD", "K00052": "leuB",
        # Lysine
        "K00928": "lysC", "K00133": "asd", "K01714": "dapA",
        "K00215": "dapB", "K01439": "dapE", "K01778": "dapF", "K01586": "lysA",
        # Methionine
        "K00003": "hom", "K00651": "metA", "K01739": "metB",
        "K01760": "metC", "K00548": "metH", "K00549": "metE",
        # Phenylalanine / Tyrosine
        "K01609": "trpC", "K00800": "aroA", "K01736": "aroC",
        "K04518": "pheA2", "K00832": "tyrB", "K01850": "pheA",
        "K00220": "tyrA",
        # Threonine
        "K00872": "thrB", "K01733": "thrC",
        # Tryptophan
        "K01657": "trpE", "K01658": "trpG", "K00766": "trpD",
        "K01817": "trpF", "K01695": "trpA", "K01696": "trpB",
        # Non-essential
        "K00814": "GPT", "K00259": "ald",
        "K00611": "argF", "K01940": "argG", "K01755": "argH",
        "K00930": "argB", "K00145": "argC", "K00821": "argD",
        "K01953": "asnB", "K01914": "asnA",
        "K00813": "aspC", "K00812": "aspB",
        "K00640": "cysE", "K01738": "cysK", "K12339": "cysM",
        "K00262": "gdhA", "K00265": "gltB", "K00266": "gltD", "K01915": "glnA",
        "K00600": "glyA", "K00281": "gcvP",
        "K00931": "proB", "K00147": "proA", "K00286": "proC",
        "K00058": "serA", "K00831": "serC", "K01079": "serB",
        # Vitamins
        "K00941": "thiD", "K00788": "thiE", "K00946": "thiL", "K03147": "thiC",
        "K00794": "ribH", "K00793": "ribE", "K02858": "ribD", "K00082": "ribA",
        "K00767": "nadC", "K00969": "nadD", "K01916": "nadE", "K03517": "nadB",
        "K00077": "panE", "K00867": "panK", "K01918": "panC",
        "K00275": "pdxH", "K03472": "pdxJ", "K03473": "pdxA",
        "K00652": "bioF", "K00833": "bioA", "K01012": "bioB", "K00796": "folK",
        "K00950": "folC", "K01930": "folE", "K00287": "folA",
        "K02232": "cobA", "K02224": "cobI", "K02229": "cobO",
    }

    # ── Detailed pathway step view (flowchart) ────────────────────────

    def pathway_detail_chart(
        self,
        strain_id: int,
        aa: str,
        pathway_source: str = "aa",
    ) -> go.Figure:
        """Draw a flowchart of individual KO enzyme steps for a pathway.

        Each enzyme is a colored box connected by arrows.
        Green = enzyme found, Red = enzyme missing.
        ★ marks essential (key) enzymes.

        Args:
            strain_id: Strain ID
            aa: Amino acid 3-letter code (e.g. "His") or vitamin key (e.g. "B1")
            pathway_source: "aa" for amino acid, "vitamin" for vitamin pathways
        """
        prior = self.prior_builder.get_prior(strain_id)

        if pathway_source == "vitamin":
            if aa not in VITAMIN_BIOSYNTHESIS_KOS:
                raise ValueError(f"Unknown vitamin: {aa}")
            pathway = VITAMIN_BIOSYNTHESIS_KOS[aa]
            completeness = prior.get("vitamin_biosynthesis", {}).get(aa, 0)
        else:
            if aa not in AA_BIOSYNTHESIS_KOS:
                raise ValueError(f"Unknown amino acid: {aa}")
            pathway = AA_BIOSYNTHESIS_KOS[aa]
            completeness = prior.get("aa_biosynthesis", {}).get(aa, 0)

        kos = pathway["kos"]
        essential = set(pathway.get("essential_kos", []))
        found_kos = self._get_found_kos(strain_id)

        return self._draw_flowchart(
            kos=kos,
            essential=essential,
            found_kos=found_kos,
            pathway_name=pathway["name"],
            module_id=pathway["pathway"],
            completeness=completeness,
        )

    def _draw_flowchart(
        self,
        kos: list[str],
        essential: set[str],
        found_kos: set[str],
        pathway_name: str,
        module_id: str,
        completeness: float,
    ) -> go.Figure:
        """Render a node-arrow flowchart using Plotly shapes.

        Layout: left-to-right, wrapping to a new row when exceeding
        ``max_per_row`` nodes.
        """
        n = len(kos)
        max_per_row = min(n, 5)
        n_rows = (n + max_per_row - 1) // max_per_row

        # ── sizing constants ──
        box_w = 1.6
        box_h = 0.9
        gap_x = 0.9          # horizontal gap between boxes
        gap_y = 1.6          # vertical gap between rows
        step_x = box_w + gap_x
        step_y = box_h + gap_y

        # Canvas bounds
        total_w = max_per_row * step_x - gap_x + 1.0
        total_h = n_rows * step_y + 0.5

        fig = go.Figure()

        shapes: list[dict] = []
        annotations: list[dict] = []
        arrow_x: list[float | None] = []
        arrow_y: list[float | None] = []

        node_centers: list[tuple[float, float]] = []

        for idx, ko in enumerate(kos):
            row = idx // max_per_row
            col = idx % max_per_row
            cx = col * step_x + box_w / 2 + 0.3
            cy = total_h - row * step_y - box_h / 2 - 0.3
            node_centers.append((cx, cy))

            is_found = ko in found_kos
            is_essential = ko in essential
            gene = self._KO_GENE_NAMES.get(ko, ko)

            fill = "#27ae60" if is_found else "#c0392b"
            border = "#1e8449" if is_found else "#922b21"

            # Rectangle shape
            shapes.append(dict(
                type="rect",
                x0=cx - box_w / 2, y0=cy - box_h / 2,
                x1=cx + box_w / 2, y1=cy + box_h / 2,
                fillcolor=fill,
                line=dict(color=border, width=2),
                layer="below",
            ))

            # Gene name label (top line)
            top_label = gene
            if is_essential:
                top_label = f"★ {gene}"
            annotations.append(dict(
                x=cx, y=cy + 0.12,
                text=f"<b>{top_label}</b>",
                showarrow=False,
                font=dict(size=12, color="white"),
                xanchor="center", yanchor="middle",
            ))
            # KO ID label (bottom line)
            status = "✓" if is_found else "✗"
            annotations.append(dict(
                x=cx, y=cy - 0.18,
                text=f"{ko}  {status}",
                showarrow=False,
                font=dict(size=9, color="rgba(255,255,255,0.85)"),
                xanchor="center", yanchor="middle",
            ))

        # ── arrows between consecutive nodes ──
        for i in range(len(node_centers) - 1):
            x0, y0 = node_centers[i]
            x1, y1 = node_centers[i + 1]

            same_row = (i // max_per_row) == ((i + 1) // max_per_row)

            if same_row:
                # Horizontal arrow: from right edge of box i to left edge of box i+1
                ax_start = x0 + box_w / 2
                ax_end = x1 - box_w / 2
                ay = y0
                # Draw line segments
                arrow_x.extend([ax_start, ax_end, None])
                arrow_y.extend([ay, ay, None])
                # Arrowhead annotation
                annotations.append(dict(
                    x=ax_end, y=ay,
                    ax=ax_end - 0.25, ay=ay,
                    xref="x", yref="y", axref="x", ayref="y",
                    showarrow=True,
                    arrowhead=3, arrowsize=1.5, arrowwidth=2,
                    arrowcolor="#555",
                ))
            else:
                # Row wrap: go down from bottom of box i, across, then into top of box i+1
                mid_y = (y0 - box_h / 2 + y1 + box_h / 2) / 2

                # Segment 1: down from box i
                arrow_x.extend([x0, x0, None])
                arrow_y.extend([y0 - box_h / 2, mid_y, None])
                # Segment 2: horizontal to x1
                arrow_x.extend([x0, x1, None])
                arrow_y.extend([mid_y, mid_y, None])
                # Segment 3: down to box i+1
                arrow_x.extend([x1, x1, None])
                arrow_y.extend([mid_y, y1 + box_h / 2, None])
                # Arrowhead
                annotations.append(dict(
                    x=x1, y=y1 + box_h / 2,
                    ax=x1, ay=y1 + box_h / 2 + 0.25,
                    xref="x", yref="y", axref="x", ayref="y",
                    showarrow=True,
                    arrowhead=3, arrowsize=1.5, arrowwidth=2,
                    arrowcolor="#555",
                ))

        # Arrow lines as a Scatter trace
        fig.add_trace(go.Scatter(
            x=arrow_x, y=arrow_y,
            mode="lines",
            line=dict(color="#555", width=2),
            hoverinfo="skip",
            showlegend=False,
        ))

        # Legend items
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode="markers",
            marker=dict(size=12, color="#27ae60", symbol="square"),
            name="효소 존재 (Found)",
        ))
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode="markers",
            marker=dict(size=12, color="#c0392b", symbol="square"),
            name="효소 결핍 (Missing)",
        ))

        # Title annotation
        title_text = (
            f"{pathway_name} ({module_id}) — "
            f"완성도: {completeness:.0%}"
        )
        annotations.append(dict(
            x=0.5, y=1.08, xref="paper", yref="paper",
            text=f"<b>{title_text}</b>",
            showarrow=False, font=dict(size=14),
            xanchor="center",
        ))
        # Legend annotation
        annotations.append(dict(
            x=0.01, y=-0.06, xref="paper", yref="paper",
            text="★ = 핵심 효소 (Essential)  |  ✓ = 존재  |  ✗ = 결핍",
            showarrow=False, font=dict(size=10, color="#666"),
            xanchor="left",
        ))

        fig.update_layout(
            shapes=shapes,
            annotations=annotations,
            xaxis=dict(visible=False, range=[0, total_w]),
            yaxis=dict(visible=False, range=[0, total_h + 0.3],
                       scaleanchor="x", scaleratio=1),
            height=max(280, n_rows * 160 + 80),
            margin=dict(l=10, r=10, t=50, b=40),
            showlegend=True,
            legend=dict(orientation="h", yanchor="bottom", y=-0.15, x=0.5, xanchor="center"),
            plot_bgcolor="white",
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
        group_labels = []  # for annotation

        essential = ["His", "Ile", "Leu", "Lys", "Met", "Phe", "Thr", "Trp", "Val"]
        non_essential = ["Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "Pro", "Ser", "Tyr"]

        for aa in essential:
            categories.append(aa)
            values.append(aa_synth.get(aa, 0.0))
            group_labels.append("Essential AA")

        for aa in non_essential:
            categories.append(aa)
            values.append(aa_synth.get(aa, 0.0))
            group_labels.append("Non-essential AA")

        # Spacer with unique label (not empty string)
        categories.append(" ")
        values.append(float("nan"))
        group_labels.append("")

        vit_keys = sorted(vit_synth.keys()) if vit_synth else ["B1", "B2", "B3", "B6", "B9"]
        for v in vit_keys:
            display = f"B{v[1:]}" if v.startswith("B") else v
            categories.append(display)
            values.append(vit_synth.get(v, 0.0))
            group_labels.append("Vitamin")

        # Second spacer with different unique label
        categories.append("  ")
        values.append(float("nan"))
        group_labels.append("")

        categories.append("Nucleotide")
        values.append(nuc if isinstance(nuc, (int, float)) else 0.5)
        group_labels.append("Other")

        categories.append("Transporter")
        values.append(trans if isinstance(trans, (int, float)) else 0.5)
        group_labels.append("Other")

        label = self._label(strain_id)

        # Build text labels (show % for valid values, empty for spacers)
        text_row = []
        for v in values:
            if v is not None and not (isinstance(v, float) and v != v):  # check for NaN
                text_row.append(f"{v:.0%}")
            else:
                text_row.append("")

        fig = go.Figure(data=go.Heatmap(
            z=[values],
            x=categories,
            y=[label],
            colorscale="RdYlGn",
            zmin=0, zmax=1,
            colorbar_title="완성도",
            text=[text_row],
            texttemplate="%{text}",
            hovertemplate="<b>%{x}</b>: %{text}<extra></extra>",
            xgap=1,
            ygap=1,
        ))

        # Add group separator annotations
        n_ess = len(essential)
        n_noness = len(non_essential)
        n_vit = len(vit_keys)

        fig.update_layout(
            title=f"생합성 경로 종합 — {label}",
            height=240,
            xaxis=dict(
                title="",
                tickangle=-45,
                side="bottom",
            ),
            annotations=[
                dict(x=n_ess / 2 - 0.5, y=1.15, text="<b>Essential AA</b>",
                     showarrow=False, xref="x", yref="paper", font=dict(size=10)),
                dict(x=n_ess + n_noness / 2 - 0.5, y=1.15, text="<b>Non-essential AA</b>",
                     showarrow=False, xref="x", yref="paper", font=dict(size=10)),
                dict(x=n_ess + n_noness + 1 + n_vit / 2 - 0.5, y=1.15, text="<b>Vitamin</b>",
                     showarrow=False, xref="x", yref="paper", font=dict(size=10)),
                dict(x=n_ess + n_noness + 1 + n_vit + 1 + 0.5, y=1.15, text="<b>Other</b>",
                     showarrow=False, xref="x", yref="paper", font=dict(size=10)),
            ],
            margin=dict(t=80, b=10),
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
        """Try to retrieve the raw KO list for this strain from cache.

        Search order:
        1. GCF-based KO list cache (from NCBI/KofamScan annotation)
        2. KEGG API KO list cache (from KEGG REST API)
        3. Fallback: reconstruct from pathway completeness (inaccurate)
        """
        import json
        from .utils import ensure_output_dir

        prior = self.prior_builder.get_prior(strain_id)
        cache_dir = ensure_output_dir("genome_cache")

        # 1. Check GCF-based KO list
        gcf = prior.get("gcf", "")
        if gcf:
            ko_file = cache_dir / f"{gcf}_ko_list.json"
            if ko_file.exists():
                try:
                    with open(ko_file, "r") as f:
                        return set(json.load(f))
                except Exception:
                    pass

        # 2. Check KEGG API KO list
        kegg_org = prior.get("kegg_org_code", "")
        if kegg_org:
            kegg_cache = ensure_output_dir("kegg_cache")
            kegg_file = kegg_cache / f"{kegg_org}_ko_list.json"
            if kegg_file.exists():
                try:
                    with open(kegg_file, "r") as f:
                        return set(json.load(f))
                except Exception:
                    pass

            # Also check genome_cache for kegg-sourced files
            kegg_file2 = cache_dir / f"kegg_{kegg_org}_ko_list.json"
            if kegg_file2.exists():
                try:
                    with open(kegg_file2, "r") as f:
                        return set(json.load(f))
                except Exception:
                    pass

        # 3. Fallback: reconstruct from pathway completeness (inaccurate)
        found = set()
        aa_synth = prior.get("aa_biosynthesis", {})
        for aa, comp in aa_synth.items():
            if aa in AA_BIOSYNTHESIS_KOS:
                pathway = AA_BIOSYNTHESIS_KOS[aa]
                all_kos = pathway.get("kos", [])
                # Use ALL kos, not just essential, for better reconstruction
                n_found = int(round(comp * len(all_kos)))
                for ko in all_kos[:n_found]:
                    found.add(ko)
        return found

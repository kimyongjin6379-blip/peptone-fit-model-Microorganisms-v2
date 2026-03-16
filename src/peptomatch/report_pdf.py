"""PDF report generation for PeptoMatch results.

Uses fpdf2 for lightweight PDF creation.
Korean text is transliterated to English for PDF compatibility (Helvetica font).
"""

import tempfile
from pathlib import Path
from typing import Any, Optional

import pandas as pd
from fpdf import FPDF

from .genome_prior import GenomePriorBuilder


def _safe(text: str) -> str:
    """Replace Korean / non-latin characters with ASCII equivalents for PDF."""
    replacements = {
        "높음": "High",
        "보통": "Moderate",
        "낮음": "Low",
        "없음": "None",
        "합성 능력 부족": "low biosynthesis",
        "→": "->",
        "의 높은": "'s high",
        "함량으로 보충 가능": "content can supplement",
        "유리 아미노산 함량": "free AA content",
        "저분자 펩타이드 비율": "low-MW peptide ratio",
        "이 높아 빠른 흡수에 유리": " favors rapid absorption",
        "균주의 펩타이드/아미노산 수송체 활성이 높아": "High transporter activity",
        "중분자 펩타이드 활용에 유리": "favors medium-MW peptides",
        "합성 능력 부족 → ": "deficiency -> supplemented by ",
        "에서 보충 가능": "",
        "핵산 요구도가 높은 균주에 적합한 뉴클레오타이드 함량": "Suitable nucleotide content for high-demand strain",
        "은 균형 잡힌 아미노산 프로파일을 제공": " provides balanced AA profile",
        "영양 요구를 충족": "meets nutritional needs",
        "총": "total:",
        "비타민": "Vitamin",
    }
    result = str(text)
    for ko, en in replacements.items():
        result = result.replace(ko, en)

    # Strip any remaining non-latin1 characters
    safe_chars = []
    for ch in result:
        try:
            ch.encode("latin-1")
            safe_chars.append(ch)
        except UnicodeEncodeError:
            safe_chars.append("?")
    return "".join(safe_chars)


class ReportGenerator:
    """Generate PDF reports for PeptoMatch analysis results."""

    def __init__(
        self,
        composition_df: pd.DataFrame,
        strain_df: pd.DataFrame,
        config: Optional[dict] = None,
        language: str = "ko",
    ):
        self.composition_df = composition_df
        self.strain_df = strain_df
        self.config = config or {}
        self.language = language
        self.prior_builder = GenomePriorBuilder(strain_df, config)

    def generate(
        self,
        strain_id: int,
        recommendations: pd.DataFrame,
        summary: dict[str, Any],
        charts: Optional[dict[str, bytes]] = None,
    ) -> bytes:
        """Generate a PDF report and return as bytes."""
        label = summary.get("name", f"Strain {strain_id}")
        pdf = self._create_pdf(_safe(label))

        self._add_strain_profile(pdf, summary)
        self._add_biosynthesis_table(pdf, strain_id)

        if charts:
            for chart_name, png_bytes in charts.items():
                self._add_chart_image(pdf, chart_name, png_bytes)

        pdf.add_page()
        self._add_recommendations(pdf, recommendations)

        return pdf.output()

    def _create_pdf(self, strain_label: str) -> FPDF:
        pdf = FPDF()
        pdf.alias_nb_pages()
        pdf.set_auto_page_break(auto=True, margin=20)
        pdf.add_page()

        pdf.set_font("Helvetica", "B", 18)
        pdf.cell(0, 12, "PeptoMatch Analysis Report", align="C",
                 new_x="LMARGIN", new_y="NEXT")
        pdf.set_font("Helvetica", "", 11)
        pdf.cell(0, 8, f"Strain: {strain_label}", align="C",
                 new_x="LMARGIN", new_y="NEXT")
        pdf.ln(6)
        pdf.line(10, pdf.get_y(), 200, pdf.get_y())
        pdf.ln(6)

        return pdf

    def _add_strain_profile(self, pdf: FPDF, summary: dict[str, Any]):
        pdf.set_font("Helvetica", "B", 13)
        pdf.cell(0, 8, "1. Strain Profile", new_x="LMARGIN", new_y="NEXT")
        pdf.ln(2)

        info_rows = [
            ("Strain Name", _safe(summary.get("name", ""))),
            ("Strain", _safe(summary.get("strain_name", ""))),
            ("GCF Accession", _safe(summary.get("gcf", "N/A"))),
            ("Data Source", _safe(summary.get("prior_source", ""))),
            ("Transporter Activity", _safe(summary.get("transporter_level", ""))),
            ("Nucleotide Synthesis", _safe(summary.get("nucleotide_synthesis", ""))),
        ]

        for label, value in info_rows:
            pdf.set_font("Helvetica", "B", 10)
            pdf.cell(55, 7, label, border=1)
            pdf.set_font("Helvetica", "", 10)
            pdf.cell(0, 7, value, border=1, new_x="LMARGIN", new_y="NEXT")

        pdf.ln(4)

        if summary.get("key_aa_deficiencies"):
            pdf.set_font("Helvetica", "B", 10)
            pdf.cell(55, 7, "Key AA Deficiencies", border=1)
            pdf.set_font("Helvetica", "", 10)
            pdf.cell(0, 7, _safe(", ".join(summary["key_aa_deficiencies"])),
                     border=1, new_x="LMARGIN", new_y="NEXT")

        if summary.get("vitamin_deficiencies"):
            pdf.set_font("Helvetica", "B", 10)
            pdf.cell(55, 7, "Vitamin Deficiencies", border=1)
            pdf.set_font("Helvetica", "", 10)
            pdf.cell(0, 7, _safe(", ".join(summary["vitamin_deficiencies"])),
                     border=1, new_x="LMARGIN", new_y="NEXT")

        pdf.ln(6)

    def _add_biosynthesis_table(self, pdf: FPDF, strain_id: int):
        prior = self.prior_builder.get_prior(strain_id)
        aa_synth = prior.get("aa_biosynthesis", {})
        if not aa_synth:
            return

        pdf.set_font("Helvetica", "B", 13)
        pdf.cell(0, 8, "2. Amino Acid Biosynthesis Completeness",
                 new_x="LMARGIN", new_y="NEXT")
        pdf.ln(2)

        essential = ["His", "Ile", "Leu", "Lys", "Met", "Phe", "Thr", "Trp", "Val"]
        non_essential = ["Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly",
                         "Pro", "Ser", "Tyr"]

        col_w = 19

        # Essential row
        pdf.set_font("Helvetica", "B", 8)
        pdf.cell(35, 7, "Category", border=1, align="C")
        for aa in essential:
            pdf.cell(col_w, 7, aa, border=1, align="C")
        pdf.ln()

        pdf.set_font("Helvetica", "", 8)
        pdf.cell(35, 7, "Essential", border=1, align="C")
        for aa in essential:
            val = aa_synth.get(aa, 0.5)
            pdf.set_fill_color(*self._completeness_color(val))
            pdf.cell(col_w, 7, f"{val:.0%}", border=1, align="C", fill=True)
        pdf.ln()
        pdf.ln(3)

        # Non-essential row
        ne_per_row = 10
        for start in range(0, len(non_essential), ne_per_row):
            chunk = non_essential[start:start + ne_per_row]
            pdf.set_font("Helvetica", "B", 8)
            pdf.cell(35, 7, "Category" if start == 0 else "", border=1, align="C")
            for aa in chunk:
                pdf.cell(col_w, 7, aa, border=1, align="C")
            pdf.ln()

            pdf.set_font("Helvetica", "", 8)
            pdf.cell(35, 7, "Non-Essential" if start == 0 else "", border=1, align="C")
            for aa in chunk:
                val = aa_synth.get(aa, 0.5)
                pdf.set_fill_color(*self._completeness_color(val))
                pdf.cell(col_w, 7, f"{val:.0%}", border=1, align="C", fill=True)
            pdf.ln()

        pdf.ln(6)

    def _add_recommendations(self, pdf: FPDF, recommendations: pd.DataFrame):
        pdf.set_font("Helvetica", "B", 13)
        pdf.cell(0, 8, "3. Peptone Recommendations", new_x="LMARGIN", new_y="NEXT")
        pdf.ln(2)

        col_widths = [12, 35, 20, 123]
        headers = ["Rank", "Peptone", "Score", "Recommendation Reason"]

        pdf.set_font("Helvetica", "B", 9)
        for w, h in zip(col_widths, headers):
            pdf.cell(w, 7, h, border=1, align="C")
        pdf.ln()

        pdf.set_font("Helvetica", "", 8)
        for _, row in recommendations.iterrows():
            rank = str(row.get("rank", ""))
            peptone = _safe(str(row.get("peptone", "")))
            score = f"{row.get('score', 0):.1f}"
            explanation = _safe(str(row.get("explanation", "")))

            max_chars = 85
            if len(explanation) > max_chars:
                explanation = explanation[:max_chars] + "..."

            pdf.cell(col_widths[0], 7, rank, border=1, align="C")
            pdf.cell(col_widths[1], 7, peptone, border=1)
            pdf.cell(col_widths[2], 7, score, border=1, align="C")
            pdf.cell(col_widths[3], 7, explanation, border=1)
            pdf.ln()

        pdf.ln(6)
        pdf.set_font("Helvetica", "I", 8)
        pdf.cell(0, 6,
                 "Generated by PeptoMatch - Genome-driven Peptone Recommendation System",
                 align="C")

    def _add_chart_image(self, pdf: FPDF, title: str, png_bytes: bytes):
        pdf.ln(4)
        pdf.set_font("Helvetica", "B", 11)
        pdf.cell(0, 8, _safe(title), new_x="LMARGIN", new_y="NEXT")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            f.write(png_bytes)
            tmp_path = f.name
        try:
            pdf.image(tmp_path, x=15, w=180)
        except Exception:
            pdf.set_font("Helvetica", "I", 9)
            pdf.cell(0, 7, "(Chart image could not be loaded)")
        finally:
            Path(tmp_path).unlink(missing_ok=True)

        pdf.ln(4)

    @staticmethod
    def _completeness_color(val: float) -> tuple[int, int, int]:
        if val >= 0.8:
            return (144, 238, 144)
        elif val >= 0.4:
            return (255, 255, 150)
        else:
            return (255, 180, 180)

"""Generate human-readable explanations for peptone recommendations."""

import logging
from typing import Any, Optional

import numpy as np
import pandas as pd

from .composition_features import CompositionFeatureExtractor
from .genome_prior import GenomePriorBuilder

logger = logging.getLogger("peptone_recommender")


# Templates for explanation generation
EXPLANATION_TEMPLATES = {
    # AA-related
    "aa_deficiency": "{aa} 합성 능력 부족 ({completeness:.0%}) → {peptone}의 높은 {aa} 함량으로 보충 가능",
    "aa_deficiency_en": "Low {aa} biosynthesis ({completeness:.0%}) → {peptone} provides high {aa} content",

    # MW-related
    "low_mw_benefit": "{peptone}의 저분자 펩타이드 비율 ({pct:.1f}%)이 높아 빠른 흡수에 유리",
    "low_mw_benefit_en": "{peptone} has high low-MW peptide content ({pct:.1f}%) for rapid absorption",

    "transporter_match": "균주의 펩타이드/아미노산 수송체 활성이 높아 ({score:.0%}) 중분자 펩타이드 활용에 유리",
    "transporter_match_en": "Strain has high transporter activity ({score:.0%}), favoring medium-MW peptides",

    # Vitamin-related
    "vitamin_deficiency": "{vitamin} 합성 능력 부족 → {peptone}에서 보충 가능",
    "vitamin_deficiency_en": "Low {vitamin} biosynthesis → {peptone} can supplement",

    # General
    "high_faa": "{peptone}의 높은 유리 아미노산 함량 (총 {value:.1f})이 영양 요구를 충족",
    "high_faa_en": "{peptone} has high free amino acid content (total: {value:.1f})",

    "balanced_profile": "{peptone}은 균형 잡힌 아미노산 프로파일을 제공",
    "balanced_profile_en": "{peptone} provides a balanced amino acid profile",

    "nucleotide_supply": "핵산 요구도가 높은 균주에 적합한 뉴클레오타이드 함량",
    "nucleotide_supply_en": "Suitable nucleotide content for strains with high nucleotide demand",
}


class RecommendationExplainer:
    """Generate explanations for peptone recommendations."""

    def __init__(
        self,
        composition_df: pd.DataFrame,
        strain_df: pd.DataFrame,
        config: Optional[dict] = None,
        language: str = "ko"
    ):
        """Initialize explainer.

        Args:
            composition_df: Peptone composition DataFrame
            strain_df: Strain table DataFrame
            config: Optional configuration
            language: Output language ('ko' for Korean, 'en' for English)
        """
        self.composition_df = composition_df
        self.strain_df = strain_df
        self.config = config or {}
        self.language = language

        self.composition_extractor = CompositionFeatureExtractor(composition_df, config)
        self.genome_prior_builder = GenomePriorBuilder(strain_df, config)

        # Precompute features
        self.all_features = self.composition_extractor.compute_all_features()
        self.supply_scores = self.composition_extractor.compute_supply_scores()

    def _get_template(self, key: str) -> str:
        """Get template for current language.

        Args:
            key: Template key

        Returns:
            Template string
        """
        if self.language == "en":
            return EXPLANATION_TEMPLATES.get(f"{key}_en", EXPLANATION_TEMPLATES.get(key, ""))
        return EXPLANATION_TEMPLATES.get(key, "")

    def explain_recommendation(
        self,
        strain_id: int,
        peptone_name: str,
        score_components: Optional[dict[str, float]] = None,
        top_n_reasons: int = 5
    ) -> list[str]:
        """Generate explanation for why a peptone is recommended.

        Args:
            strain_id: Strain ID
            peptone_name: Recommended peptone name
            score_components: Optional precomputed score components
            top_n_reasons: Number of top reasons to include

        Returns:
            List of explanation strings
        """
        explanations = []

        # Get strain prior
        prior = self.genome_prior_builder.get_prior(strain_id)
        demand = self.genome_prior_builder.get_demand_scores(strain_id)

        # Get peptone features
        if peptone_name not in self.all_features.index:
            return ["펩톤 정보를 찾을 수 없습니다." if self.language == "ko" else "Peptone information not found."]

        features = self.all_features.loc[peptone_name]
        supply = self.supply_scores.loc[peptone_name]

        # Collect potential reasons with scores
        reasons = []

        # 1. Check AA deficiencies and matching supplies
        aa_synth = prior.get("aa_biosynthesis", {})
        all_aa = ["His", "Ile", "Leu", "Lys", "Met", "Phe", "Thr", "Trp", "Val",
                  "Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "Pro", "Ser", "Tyr"]

        for aa in all_aa:
            completeness = aa_synth.get(aa, 0.5)
            supply_score = supply.get(f"supply_{aa}", 0)

            # If deficient (< 70%) and peptone supplies it
            if completeness < 0.7 and supply_score > 0.2:
                impact = (1 - completeness) * supply_score
                reasons.append({
                    "type": "aa_deficiency",
                    "impact": impact,
                    "text": self._get_template("aa_deficiency").format(
                        aa=aa,
                        completeness=completeness,
                        peptone=peptone_name
                    )
                })

        # 2. MW distribution match
        transporter_score = prior.get("transporter_score", 0.5)
        low_mw_pct = features.get("mw_pct_low", 0)

        if low_mw_pct > 20:
            reasons.append({
                "type": "low_mw",
                "impact": low_mw_pct / 100 * 0.8,
                "text": self._get_template("low_mw_benefit").format(
                    peptone=peptone_name,
                    pct=low_mw_pct
                )
            })

        if transporter_score > 0.6:
            reasons.append({
                "type": "transporter",
                "impact": transporter_score * 0.7,
                "text": self._get_template("transporter_match").format(
                    score=transporter_score
                )
            })

        # 3. Vitamin deficiency match
        vit_synth = prior.get("vitamin_biosynthesis", {})
        for vit, completeness in vit_synth.items():
            vit_supply = features.get(f"vitamin_{vit.lower()}", 0)
            if completeness < 0.3 and vit_supply > 0:
                reasons.append({
                    "type": "vitamin",
                    "impact": (1 - completeness) * 0.5,
                    "text": self._get_template("vitamin_deficiency").format(
                        vitamin=f"비타민 {vit}" if self.language == "ko" else f"Vitamin {vit}",
                        peptone=peptone_name
                    )
                })

        # 4. High FAA content
        faa_total = features.get("faa_total", 0)
        if faa_total > 5:  # Arbitrary threshold
            reasons.append({
                "type": "faa",
                "impact": min(faa_total / 20, 1.0) * 0.6,
                "text": self._get_template("high_faa").format(
                    peptone=peptone_name,
                    value=faa_total
                )
            })

        # 5. Nucleotide supply if needed
        nuc_demand = demand.get("demand_nucleotide", 0.5)
        nuc_supply = supply.get("supply_nucleotide", 0)
        if nuc_demand > 0.5 and nuc_supply > 0.3:
            reasons.append({
                "type": "nucleotide",
                "impact": nuc_demand * nuc_supply * 0.5,
                "text": self._get_template("nucleotide_supply")
            })

        # Sort by impact and take top N
        reasons.sort(key=lambda x: x["impact"], reverse=True)
        explanations = [r["text"] for r in reasons[:top_n_reasons]]

        # Add fallback if no specific reasons
        if not explanations:
            explanations.append(
                self._get_template("balanced_profile").format(peptone=peptone_name)
            )

        return explanations

    def explain_batch(
        self,
        strain_id: int,
        recommendations: pd.DataFrame,
        top_n_reasons: int = 3
    ) -> pd.DataFrame:
        """Add explanations to a recommendations DataFrame.

        Args:
            strain_id: Strain ID
            recommendations: Recommendations DataFrame with 'peptone' column
            top_n_reasons: Number of reasons per peptone

        Returns:
            DataFrame with added 'explanation' column
        """
        explanations = []

        for _, row in recommendations.iterrows():
            peptone = row["peptone"]
            components = row.get("components", None)

            explanation = self.explain_recommendation(
                strain_id,
                peptone,
                score_components=components,
                top_n_reasons=top_n_reasons
            )

            explanations.append(" | ".join(explanation))

        result = recommendations.copy()
        result["explanation"] = explanations

        return result

    def get_strain_summary(self, strain_id: int) -> dict[str, Any]:
        """Get a summary of strain characteristics for recommendation context.

        Args:
            strain_id: Strain ID

        Returns:
            Summary dictionary
        """
        strain_row = self.strain_df[self.strain_df["strain_id"] == strain_id]

        if strain_row.empty:
            return {"error": "Strain not found"}

        strain_info = strain_row.iloc[0]
        prior = self.genome_prior_builder.get_prior(strain_id)

        # Find key deficiencies
        aa_synth = prior.get("aa_biosynthesis", {})
        most_deficient = sorted(aa_synth.items(), key=lambda x: x[1])[:3]

        vit_synth = prior.get("vitamin_biosynthesis", {})
        deficient_vitamins = [v for v, c in vit_synth.items() if c < 0.4]

        summary = {
            "strain_id": strain_id,
            "name": f"{strain_info.get('genus', '')} {strain_info.get('species', '')}",
            "strain_name": strain_info.get("strain_name", ""),
            "gcf": strain_info.get("GCF", "N/A"),
            "prior_source": prior.get("source", "unknown"),
            "key_aa_deficiencies": [f"{aa} ({comp:.0%})" for aa, comp in most_deficient],
            "vitamin_deficiencies": deficient_vitamins,
            "transporter_level": "높음" if prior.get("transporter_score", 0.5) > 0.6 else "보통",
            "nucleotide_synthesis": f"{prior.get('nucleotide_biosynthesis', 0.5):.0%}"
        }

        if self.language == "en":
            summary["transporter_level"] = "high" if prior.get("transporter_score", 0.5) > 0.6 else "moderate"

        return summary

    def generate_report(
        self,
        strain_id: int,
        recommendations: pd.DataFrame
    ) -> str:
        """Generate a full text report for recommendations.

        Args:
            strain_id: Strain ID
            recommendations: Recommendations DataFrame

        Returns:
            Formatted report string
        """
        summary = self.get_strain_summary(strain_id)

        lines = []

        # Header
        if self.language == "ko":
            lines.append(f"# {summary['name']} 펩톤 추천 리포트")
            lines.append("")
            lines.append("## 균주 정보")
            lines.append(f"- 균주명: {summary['strain_name']}")
            lines.append(f"- GCF: {summary['gcf']}")
            lines.append(f"- 정보 출처: {summary['prior_source']}")
            lines.append("")
            lines.append("## 주요 영양 요구")
            lines.append(f"- 주요 아미노산 결핍: {', '.join(summary['key_aa_deficiencies'])}")
            lines.append(f"- 비타민 결핍: {', '.join(summary['vitamin_deficiencies']) or '없음'}")
            lines.append(f"- 펩타이드 수송체 활성: {summary['transporter_level']}")
            lines.append("")
            lines.append("## 추천 펩톤 (상위 순위)")
            lines.append("")
        else:
            lines.append(f"# Peptone Recommendation Report for {summary['name']}")
            lines.append("")
            lines.append("## Strain Information")
            lines.append(f"- Strain: {summary['strain_name']}")
            lines.append(f"- GCF: {summary['gcf']}")
            lines.append(f"- Data source: {summary['prior_source']}")
            lines.append("")
            lines.append("## Key Nutritional Requirements")
            lines.append(f"- Key AA deficiencies: {', '.join(summary['key_aa_deficiencies'])}")
            lines.append(f"- Vitamin deficiencies: {', '.join(summary['vitamin_deficiencies']) or 'None'}")
            lines.append(f"- Peptide transporter activity: {summary['transporter_level']}")
            lines.append("")
            lines.append("## Recommended Peptones (Top Ranked)")
            lines.append("")

        # Recommendations
        for _, row in recommendations.iterrows():
            rank = row.get("rank", "")
            peptone = row["peptone"]
            score = row["score"]

            lines.append(f"### {rank}. {peptone} (점수: {score:.1f})" if self.language == "ko"
                        else f"### {rank}. {peptone} (Score: {score:.1f})")

            explanation = self.explain_recommendation(strain_id, peptone, top_n_reasons=3)
            for exp in explanation:
                lines.append(f"  - {exp}")

            lines.append("")

        return "\n".join(lines)


def format_score_breakdown(components: dict[str, float], language: str = "ko") -> str:
    """Format score components for display.

    Args:
        components: Dictionary of score components
        language: Display language

    Returns:
        Formatted string
    """
    labels = {
        "ko": {
            "faa_abundance": "유리 아미노산",
            "taa_abundance": "총 아미노산",
            "low_mw": "저분자 펩타이드",
            "medium_mw": "중분자 펩타이드",
            "high_mw": "고분자 펩타이드",
            "vitamins": "비타민",
            "nucleotides": "뉴클레오타이드",
            "aa_match": "아미노산 매칭",
            "transporter_bonus": "수송체 보너스"
        },
        "en": {
            "faa_abundance": "Free AA",
            "taa_abundance": "Total AA",
            "low_mw": "Low MW",
            "medium_mw": "Medium MW",
            "high_mw": "High MW",
            "vitamins": "Vitamins",
            "nucleotides": "Nucleotides",
            "aa_match": "AA Match",
            "transporter_bonus": "Transporter"
        }
    }

    label_map = labels.get(language, labels["en"])

    parts = []
    for key, value in sorted(components.items(), key=lambda x: -x[1]):
        label = label_map.get(key, key)
        parts.append(f"{label}: {value:.2f}")

    return " | ".join(parts)

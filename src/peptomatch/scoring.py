"""Scoring and recommendation logic for peptone-strain matching."""

import logging
from typing import Any, Optional

import numpy as np
import pandas as pd

from .composition_features import CompositionFeatureExtractor
from .genome_prior import GenomePriorBuilder
from .media_config import get_default_media, get_media_responsibility, MEDIA_CONFIGS

logger = logging.getLogger("peptone_recommender")

# ============================================================================
# Strain-type-specific weight presets (literature-based)
#
# Rationale:
# - LAB (lactic acid bacteria): fermentative, carbon-source dependent,
#   highly auxotrophic but sugar is often the growth-limiting factor
# - Enterobacteriaceae (E. coli etc.): strong biosynthetic capability,
#   AA supplementation is the main benefit of peptone
# - Bacillus: versatile metabolism, balanced needs
# - Corynebacterium: AA-producing but often needs vitamins/minerals
# - Bifidobacterium: similar to LAB but stricter anaerobe, sugar-dependent
#
# Note: These are educated estimates. Will be calibrated with experimental
# data when available. Confidence level is indicated per preset.
# ============================================================================

STRAIN_TYPE_PRESETS: dict[str, dict] = {
    "LAB": {
        "description": "Lactic Acid Bacteria (fermentative, peptone as nitrogen source)",
        "confidence": "medium (literature-based, not yet experimentally validated)",
        # Design rationale:
        # Peptone is added at only 2.0-2.5% in media. Base medium (MRS etc.)
        # already provides ~20g/L glucose. Peptone's primary role is nitrogen/
        # peptide supply, NOT carbon source. Sugar in peptone is a minor bonus.
        # Minerals from peptone can supplement base medium salts.
        "weights": {
            "faa_abundance": 0.9,     # FAA is important — LAB peptide transporter active
            "taa_abundance": 0.7,     # TAA provides peptide-bound AA
            "mw_low": 1.2,            # di/tri-peptides well absorbed by LAB
            "mw_medium": 0.8,
            "mw_high": 0.4,
            "vitamin_b": 0.8,         # LAB often auxotrophic for B vitamins
            "nucleotides": 0.5,
            "aa_biosynthesis_gap": 1.3,  # AA matching still important for LAB
            "transporter_bonus": 0.6,
            "sugar": 0.8,             # minor bonus only — base medium supplies most sugar
            "mineral": 1.2,           # Mg/Mn important, supplements base medium
            "orgacid": 0.3,           # small bonus — LAB tolerant but base medium handles pH
            "nitrogen_quality": 0.8,  # AN/TN ratio matters for nitrogen utilization
        },
    },
    "Enterobacteriaceae": {
        "description": "E. coli, Salmonella, etc. (strong biosynthesis, AA-focused)",
        "confidence": "medium (literature-based)",
        "weights": {
            "faa_abundance": 1.2,
            "taa_abundance": 1.0,
            "mw_low": 1.2,
            "mw_medium": 1.0,
            "mw_high": 0.6,
            "vitamin_b": 0.3,         # E. coli synthesizes most vitamins
            "nucleotides": 0.3,
            "aa_biosynthesis_gap": 2.0,  # AA matching is most important
            "transporter_bonus": 0.5,
            "sugar": 0.3,             # base medium provides glucose; peptone sugar negligible
            "mineral": 0.4,           # base medium salts sufficient for E. coli
            "orgacid": 0.2,
            "nitrogen_quality": 0.8,
        },
    },
    "Bacillus": {
        "description": "Bacillus spp. (versatile, balanced needs)",
        "confidence": "medium (literature-based)",
        "weights": {
            "faa_abundance": 1.0,
            "taa_abundance": 0.8,
            "mw_low": 1.0,
            "mw_medium": 1.0,
            "mw_high": 0.6,
            "vitamin_b": 0.5,
            "nucleotides": 0.4,
            "aa_biosynthesis_gap": 1.5,
            "transporter_bonus": 0.5,
            "sugar": 0.5,             # base medium supplies carbon
            "mineral": 0.8,           # Bacillus needs Mn/Mg for sporulation enzymes
            "orgacid": 0.3,
            "nitrogen_quality": 0.7,
        },
    },
    "Corynebacterium": {
        "description": "Corynebacterium spp. (AA producer, vitamin-dependent)",
        "confidence": "low (limited literature)",
        "weights": {
            "faa_abundance": 0.8,
            "taa_abundance": 0.8,
            "mw_low": 1.0,
            "mw_medium": 1.0,
            "mw_high": 0.6,
            "vitamin_b": 1.2,         # often needs B vitamins
            "nucleotides": 0.5,
            "aa_biosynthesis_gap": 1.2,
            "transporter_bonus": 0.4,
            "sugar": 0.4,             # base medium provides carbon
            "mineral": 1.0,           # trace metals important for AA production
            "orgacid": 0.3,
            "nitrogen_quality": 0.6,
        },
    },
    "Bifidobacterium": {
        "description": "Bifidobacterium spp. (strict anaerobe, nitrogen-focused)",
        "confidence": "medium (literature-based)",
        "weights": {
            "faa_abundance": 0.8,
            "taa_abundance": 0.7,
            "mw_low": 1.0,
            "mw_medium": 0.8,
            "mw_high": 0.4,
            "vitamin_b": 0.6,
            "nucleotides": 0.6,
            "aa_biosynthesis_gap": 1.2,
            "transporter_bonus": 0.5,
            "sugar": 0.8,             # minor bonus — base medium supplies sugar
            "mineral": 1.0,
            "orgacid": 0.3,
            "nitrogen_quality": 0.7,
        },
    },
    "default": {
        "description": "Default weights — peptone as nitrogen/peptide source",
        "confidence": "low (generic)",
        # Peptone at 2-2.5% in media: primary role is nitrogen supply.
        # Sugar/mineral are secondary supplements to base medium.
        "weights": {
            "faa_abundance": 1.0,
            "taa_abundance": 0.8,
            "mw_low": 1.2,
            "mw_medium": 1.0,
            "mw_high": 0.6,
            "vitamin_b": 0.5,
            "nucleotides": 0.4,
            "aa_biosynthesis_gap": 1.5,
            "transporter_bonus": 0.5,
            "sugar": 0.5,             # minor — base medium provides most sugar
            "mineral": 0.6,           # supplementary to base medium salts
            "orgacid": 0.3,           # small factor
            "nitrogen_quality": 0.6,
        },
    },
}

# Genus → strain type mapping
GENUS_TO_TYPE: dict[str, str] = {
    # LAB
    "Lactobacillus": "LAB",
    "Lactiplantibacillus": "LAB",
    "Lacticaseibacillus": "LAB",
    "Levilactobacillus": "LAB",
    "Lentilactobacillus": "LAB",
    "Latilactobacillus": "LAB",
    "Limosilactobacillus": "LAB",
    "Ligilactobacillus": "LAB",
    "Lactococcus": "LAB",
    "Streptococcus": "LAB",
    "Leuconostoc": "LAB",
    "Pediococcus": "LAB",
    "Enterococcus": "LAB",
    "Weissella": "LAB",
    # Enterobacteriaceae
    "Escherichia": "Enterobacteriaceae",
    "Salmonella": "Enterobacteriaceae",
    "Klebsiella": "Enterobacteriaceae",
    "Enterobacter": "Enterobacteriaceae",
    "Serratia": "Enterobacteriaceae",
    # Bacillus
    "Bacillus": "Bacillus",
    "Priestia": "Bacillus",
    "Brevibacillus": "Bacillus",
    # Corynebacterium
    "Corynebacterium": "Corynebacterium",
    "Brevibacterium": "Corynebacterium",
    # Bifidobacterium
    "Bifidobacterium": "Bifidobacterium",
}


def get_strain_type(genus: str) -> str:
    """Get strain type preset name from genus."""
    return GENUS_TO_TYPE.get(genus, "default")


def get_weight_preset(genus: str) -> dict:
    """Get weight preset for a genus. Returns (preset_name, weights, confidence)."""
    strain_type = get_strain_type(genus)
    preset = STRAIN_TYPE_PRESETS[strain_type]
    return {
        "type": strain_type,
        "weights": preset["weights"],
        "confidence": preset["confidence"],
        "description": preset["description"],
    }


class PeptoneRecommender:
    """Rule-based peptone recommendation system."""

    def __init__(
        self,
        composition_df: pd.DataFrame,
        strain_df: pd.DataFrame,
        config: Optional[dict] = None
    ):
        """Initialize recommender.

        Args:
            composition_df: Peptone composition DataFrame
            strain_df: Strain table DataFrame
            config: Optional configuration dictionary
        """
        self.composition_df = composition_df
        self.strain_df = strain_df
        self.config = config or {}

        # Initialize feature extractors
        self.composition_extractor = CompositionFeatureExtractor(composition_df, config)
        self.genome_prior_builder = GenomePriorBuilder(strain_df, config)

        # Compute supply scores, normalized within the peptone filter subset
        norm_samples = self.config.get("peptone_filter", None)
        self.supply_scores = self.composition_extractor.compute_supply_scores(
            norm_samples=norm_samples
        )

        # Get weights from config
        self.weights = self.config.get("weights", {})

        logger.info("Initialized PeptoneRecommender")

    def recommend(
        self,
        strain_id: int,
        top_k: int = 10,
        min_score: float = 0.0,
        peptone_filter: Optional[list[str]] = None,
        media_key: Optional[str] = None,
    ) -> pd.DataFrame:
        """Get top-K peptone recommendations for a strain.

        Args:
            strain_id: Strain ID from strain table
            top_k: Number of recommendations to return
            min_score: Minimum score threshold
            peptone_filter: If provided, only consider these peptone names
            media_key: Media config key (e.g. "MRS_2.0"). Auto-detected if None.

        Returns:
            DataFrame with recommendations (peptone, score, components)
        """
        # Get strain info
        strain_info = self.strain_df[self.strain_df["strain_id"] == strain_id]

        if strain_info.empty:
            raise ValueError(f"Strain ID {strain_id} not found")

        strain_row = strain_info.iloc[0]
        genus = strain_row.get("genus", "")
        logger.info(
            f"Computing recommendations for strain {strain_id}: "
            f"{genus} {strain_row.get('species', '')}"
        )

        # Resolve media key and get responsibility weights
        if media_key is None:
            media_key = get_default_media(genus)
        responsibility = get_media_responsibility(media_key)
        media_cfg = MEDIA_CONFIGS.get(media_key, MEDIA_CONFIGS["MRS_2.5"])
        logger.info(
            f"Using media config: {media_key} "
            f"({media_cfg.get('display_name', media_key)})"
        )

        # Get demand scores for this strain
        demand_scores = self.genome_prior_builder.get_demand_scores(strain_id)

        # Store media info in demand for UI access
        demand_scores["_media_key"] = media_key
        demand_scores["_media_display_name"] = media_cfg.get("display_name", media_key)
        demand_scores["_media_description"] = media_cfg.get("description", "")
        demand_scores["_media_peptone_g_per_L"] = media_cfg.get("peptone_g_per_L", 25.0)
        demand_scores["_media_responsibility"] = responsibility

        # Also keep legacy preset info for backward compat
        preset = get_weight_preset(genus)
        demand_scores["_weight_preset"] = preset["type"]
        demand_scores["_weight_confidence"] = preset["confidence"]
        demand_scores["_weight_description"] = preset["description"]

        # Determine which peptones to score
        if peptone_filter is None:
            peptone_filter = self.config.get("peptone_filter", None)

        peptone_names = self.supply_scores.index
        if peptone_filter:
            peptone_names = [p for p in peptone_names if p in peptone_filter]

        # Compute match scores for each peptone
        results = []

        for peptone_name in peptone_names:
            supply = self.supply_scores.loc[peptone_name]
            score, components = self._compute_match_score(
                supply, demand_scores, responsibility=responsibility
            )

            results.append({
                "peptone": peptone_name,
                "score": score,
                "components": components
            })

        # Create DataFrame and sort
        df_results = pd.DataFrame(results)
        df_results = df_results.sort_values("score", ascending=False)
        df_results = df_results[df_results["score"] >= min_score]
        df_results = df_results.head(top_k)
        df_results = df_results.reset_index(drop=True)
        df_results["rank"] = df_results.index + 1

        return df_results

    def _compute_match_score(
        self,
        supply: pd.Series,
        demand: dict[str, float],
        strain_weights: Optional[dict[str, float]] = None,
        responsibility: Optional[dict[str, float]] = None,
    ) -> tuple[float, dict[str, float]]:
        """Compute match score between peptone supply and strain demand.

        Args:
            supply: Supply scores for a peptone
            demand: Demand scores for a strain
            strain_weights: Legacy strain-type-specific weights (deprecated, use responsibility)
            responsibility: Media-based peptone responsibility weights (preferred)

        Returns:
            Tuple of (total_score, component_scores)
        """
        components = {}

        if responsibility is not None:
            # NEW: Media-based responsibility scoring
            # Final weight = base_strength × peptone_responsibility
            base_strength = {
                "faa_abundance": 1.0, "taa_abundance": 0.8,
                "mw_low": 1.2, "mw_medium": 1.0, "mw_high": 0.6,
                "vitamin_b": 0.8, "nucleotides": 0.5,
                "aa_match": 1.5, "transporter": 0.5,
                "sugar": 1.5, "mineral": 1.0, "orgacid": 0.5, "nitrogen_quality": 0.7,
            }

            # Map weight keys to responsibility keys
            resp_key_map = {
                "faa_abundance": "faa_abundance",
                "taa_abundance": "taa_abundance",
                "mw_low": "mw_low",
                "mw_medium": "mw_medium",
                "mw_high": "mw_high",
                "vitamin_b": "vitamin_b",
                "nucleotides": "nucleotides",
                "aa_match": "aa_biosynthesis_gap",
                "transporter": "transporter_bonus",
                "sugar": "sugar",
                "mineral": "mineral",
                "orgacid": "orgacid",
                "nitrogen_quality": "nitrogen_quality",
            }

            resp = responsibility
            w = {}
            for wkey, bval in base_strength.items():
                rkey = resp_key_map.get(wkey, wkey)
                w[wkey] = bval * resp.get(rkey, 0.5)
        elif strain_weights is not None:
            # LEGACY: strain-type-specific weights (backward compat)
            sw = strain_weights
            defaults = STRAIN_TYPE_PRESETS["default"]["weights"]
            w = {
                "faa_abundance": sw.get("faa_abundance", defaults["faa_abundance"]),
                "taa_abundance": sw.get("taa_abundance", defaults["taa_abundance"]),
                "mw_low": sw.get("mw_low", defaults["mw_low"]),
                "mw_medium": sw.get("mw_medium", defaults["mw_medium"]),
                "mw_high": sw.get("mw_high", defaults["mw_high"]),
                "vitamin_b": sw.get("vitamin_b", defaults["vitamin_b"]),
                "nucleotides": sw.get("nucleotides", defaults["nucleotides"]),
                "aa_match": sw.get("aa_biosynthesis_gap", defaults["aa_biosynthesis_gap"]),
                "transporter": sw.get("transporter_bonus", defaults["transporter_bonus"]),
                "sugar": sw.get("sugar", defaults["sugar"]),
                "mineral": sw.get("mineral", defaults["mineral"]),
                "orgacid": sw.get("orgacid", defaults["orgacid"]),
                "nitrogen_quality": sw.get("nitrogen_quality", defaults["nitrogen_quality"]),
            }
        else:
            # Fallback: use default preset weights
            defaults = STRAIN_TYPE_PRESETS["default"]["weights"]
            w = {
                "faa_abundance": defaults["faa_abundance"],
                "taa_abundance": defaults["taa_abundance"],
                "mw_low": defaults["mw_low"],
                "mw_medium": defaults["mw_medium"],
                "mw_high": defaults["mw_high"],
                "vitamin_b": defaults["vitamin_b"],
                "nucleotides": defaults["nucleotides"],
                "aa_match": defaults["aa_biosynthesis_gap"],
                "transporter": defaults["transporter_bonus"],
                "sugar": defaults["sugar"],
                "mineral": defaults["mineral"],
                "orgacid": defaults["orgacid"],
                "nitrogen_quality": defaults["nitrogen_quality"],
            }

        # 1. Base composition scores (supply only)
        faa_score = supply.get("supply_faa", 0) * w["faa_abundance"]
        components["faa_abundance"] = faa_score

        taa_score = supply.get("supply_taa", 0) * w["taa_abundance"]
        components["taa_abundance"] = taa_score

        # 2. MW distribution score
        transporter_bonus = demand.get("transporter_bonus", 0.5)

        # Prefer low MW for high transporter strains
        low_mw_score = supply.get("supply_low_mw", 0) * w["mw_low"] * (1 + transporter_bonus)
        components["low_mw"] = low_mw_score

        medium_mw_score = supply.get("supply_medium_mw", 0) * w["mw_medium"] * transporter_bonus
        components["medium_mw"] = medium_mw_score

        high_mw_score = supply.get("supply_high_mw", 0) * w["mw_high"] * (1 - transporter_bonus * 0.3)
        components["high_mw"] = high_mw_score

        # 3. Vitamin and nucleotide scores weighted by demand
        vitamin_demand = np.mean([demand.get(f"demand_vitamin_{v}", 0.5) for v in ["B1", "B2", "B3", "B6", "B9"]])
        vitamin_score = supply.get("supply_vitamin", 0) * vitamin_demand * w["vitamin_b"]
        components["vitamins"] = vitamin_score

        nucleotide_demand = demand.get("demand_nucleotide", 0.5)
        nucleotide_score = supply.get("supply_nucleotide", 0) * nucleotide_demand * w["nucleotides"]
        components["nucleotides"] = nucleotide_score

        # 4. AA-specific matching (supply AAs the strain can't make)
        all_aa = ["His", "Ile", "Leu", "Lys", "Met", "Phe", "Thr", "Trp", "Val",
                  "Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "Pro", "Ser", "Tyr"]
        aa_match_score = 0
        aa_count = 0

        for aa in all_aa:
            supply_key = f"supply_{aa}"
            demand_key = f"demand_{aa}"

            if supply_key in supply.index and demand_key in demand:
                # High demand + high supply = good match
                aa_match = supply[supply_key] * demand[demand_key]
                aa_match_score += aa_match
                aa_count += 1

        if aa_count > 0:
            aa_match_score = (aa_match_score / aa_count) * w["aa_match"]
        components["aa_match"] = aa_match_score

        # 5. Transporter bonus for peptide-rich peptones
        transporter_score = transporter_bonus * (
            supply.get("supply_medium_mw", 0) + supply.get("supply_low_mw", 0) * 0.5
        ) * w["transporter"]
        components["transporter_bonus"] = transporter_score

        # ===== NEW SCORING COMPONENTS =====

        # 6. Sugar/carbon source matching
        # Score = Σ(sugar_supply × strain_can_use) for each sugar type
        sugar_types = ["glucose", "sucrose", "lactose", "maltose"]
        sugar_score = 0
        sugar_count = 0
        for sugar in sugar_types:
            supply_key = f"supply_sugar_{sugar}"
            demand_key = f"sugar_utilization_{sugar}"
            s_val = supply.get(supply_key, 0)
            d_val = demand.get(demand_key, 0)
            if s_val > 0 or d_val > 0:
                sugar_score += s_val * d_val
                sugar_count += 1
        # Also add total sugar weighted by average utilization
        avg_utilization = np.mean([demand.get(f"sugar_utilization_{s}", 0) for s in sugar_types])
        sugar_total_bonus = supply.get("supply_sugar_total", 0) * avg_utilization * 0.5
        sugar_score += sugar_total_bonus
        if sugar_count > 0:
            sugar_score = sugar_score * w["sugar"]
        components["sugar_carbon"] = sugar_score

        # 7. Mineral supply matching
        # Score = Σ(mineral_supply × strain_has_transporter) for each mineral
        mineral_types = [("K", "K"), ("Mg", "Mg"), ("Ca", "Ca"), ("Na", "Na")]
        mineral_score = 0
        for supply_mineral, demand_mineral in mineral_types:
            supply_key = f"supply_mineral_{supply_mineral}"
            demand_key = f"mineral_demand_{demand_mineral}"
            s_val = supply.get(supply_key, 0)
            d_val = demand.get(demand_key, 0.5)  # Default 0.5 = assume moderate need
            mineral_score += s_val * d_val
        # Add Fe and Mn demand even without direct supply data (bonus for Mg-rich)
        fe_demand = demand.get("mineral_demand_Fe", 0)
        mn_demand = demand.get("mineral_demand_Mn", 0)
        # Mg-rich peptones often co-supply trace metals
        mineral_score += supply.get("supply_mineral_Mg", 0) * (fe_demand + mn_demand) * 0.2
        mineral_score = mineral_score * w["mineral"]
        components["mineral_supply"] = mineral_score

        # 8. Organic acid utilization (replaces simple penalty)
        # If strain CAN use the acid → positive; if CANNOT → penalty
        orgacid_types = [("lactate", "lactate"), ("citrate", "citrate"),
                         ("acetate", "acetate"), ("succinate", "succinate"),
                         ("malate", "malate")]
        orgacid_score = 0
        for supply_acid, demand_acid in orgacid_types:
            supply_key = f"supply_orgacid_{supply_acid}"
            demand_key = f"orgacid_utilization_{demand_acid}"
            s_val = supply.get(supply_key, 0)
            d_val = demand.get(demand_key, 0)  # 0 = cannot use

            if d_val >= 0.5:
                # Strain CAN use this acid → positive contribution
                orgacid_score += s_val * d_val * 0.3
            else:
                # Strain CANNOT use this acid → penalty
                orgacid_score -= s_val * (1 - d_val) * 0.3
        orgacid_score = orgacid_score * w["orgacid"]
        components["orgacid_utilization"] = orgacid_score

        # 9. Nitrogen quality bonus (AN/TN ratio)
        nq_score = supply.get("supply_nitrogen_quality", 0) * w["nitrogen_quality"]
        components["nitrogen_quality"] = nq_score

        # Total score
        total_score = sum(components.values())

        # Normalize to 0-100 scale
        max_possible = sum(w.values()) * 2  # Rough estimate
        normalized_score = max(0, (total_score / max_possible) * 100)

        return normalized_score, components

    def recommend_all_strains(
        self,
        top_k: int = 5
    ) -> dict[int, pd.DataFrame]:
        """Get recommendations for all strains.

        Args:
            top_k: Number of recommendations per strain

        Returns:
            Dictionary of strain_id -> recommendations DataFrame
        """
        all_recommendations = {}

        for _, row in self.strain_df.iterrows():
            strain_id = row["strain_id"]
            try:
                recommendations = self.recommend(strain_id, top_k=top_k)
                all_recommendations[strain_id] = recommendations
            except Exception as e:
                logger.warning(f"Failed to get recommendations for strain {strain_id}: {e}")

        return all_recommendations

    def get_peptone_ranking_for_strain(
        self,
        strain_id: int
    ) -> pd.DataFrame:
        """Get full peptone ranking for a strain with detailed scores.

        Args:
            strain_id: Strain ID

        Returns:
            DataFrame with all peptones ranked
        """
        return self.recommend(strain_id, top_k=len(self.supply_scores))

    def compare_peptones(
        self,
        strain_id: int,
        peptone_names: list[str]
    ) -> pd.DataFrame:
        """Compare specific peptones for a strain.

        Args:
            strain_id: Strain ID
            peptone_names: List of peptone names to compare

        Returns:
            Comparison DataFrame
        """
        demand_scores = self.genome_prior_builder.get_demand_scores(strain_id)

        results = []
        for peptone in peptone_names:
            if peptone not in self.supply_scores.index:
                logger.warning(f"Peptone '{peptone}' not found")
                continue

            supply = self.supply_scores.loc[peptone]
            score, components = self._compute_match_score(supply, demand_scores)

            result = {"peptone": peptone, "total_score": score}
            result.update(components)
            results.append(result)

        return pd.DataFrame(results)

    def get_strain_info(self, strain_id: int) -> dict[str, Any]:
        """Get detailed information about a strain.

        Args:
            strain_id: Strain ID

        Returns:
            Strain information dictionary
        """
        strain_row = self.strain_df[self.strain_df["strain_id"] == strain_id]

        if strain_row.empty:
            raise ValueError(f"Strain ID {strain_id} not found")

        info = strain_row.iloc[0].to_dict()

        # Add prior summary
        prior_summary = self.genome_prior_builder.get_prior_summary(strain_id)
        info["prior_summary"] = prior_summary

        return info

    def get_peptone_info(self, peptone_name: str) -> dict[str, Any]:
        """Get detailed information about a peptone.

        Args:
            peptone_name: Peptone sample name

        Returns:
            Peptone information dictionary
        """
        return self.composition_extractor.get_peptone_profile(peptone_name)

    def recommend_blend(
        self,
        strain_id: int,
        max_components: int = 3,
        top_k: int = 5,
        peptone_filter: Optional[list[str]] = None,
        min_ratio: float = 0.1,
        max_ratio: float = 0.8,
        media_key: Optional[str] = None,
    ) -> list:
        """Get top blended peptone recommendations for a strain.

        Args:
            strain_id: Strain ID from strain table
            max_components: Maximum components per blend (2 or 3)
            top_k: Number of top blends to return
            peptone_filter: Optional list of peptone names to consider
            min_ratio: Minimum ratio per component
            max_ratio: Maximum ratio per component
            media_key: Media config key (e.g. "MRS_2.0"). Auto-detected if None.

        Returns:
            List of BlendResult objects
        """
        from .blend_optimizer import BlendOptimizer

        # Get strain info for media selection
        strain_info = self.strain_df[self.strain_df["strain_id"] == strain_id]
        genus = strain_info.iloc[0].get("genus", "") if not strain_info.empty else ""

        # Resolve media key and get responsibility weights
        if media_key is None:
            media_key = get_default_media(genus)
        responsibility = get_media_responsibility(media_key)

        # Get demand scores
        demand_scores = self.genome_prior_builder.get_demand_scores(strain_id)

        # Determine peptone filter
        if peptone_filter is None:
            peptone_filter = self.config.get("peptone_filter", None)

        # Create scoring function with media responsibility baked in
        def scoring_fn_with_weights(supply, demand):
            return self._compute_match_score(supply, demand, responsibility=responsibility)

        # Create optimizer
        optimizer = BlendOptimizer(
            supply_scores=self.supply_scores,
            min_ratio=min_ratio,
            max_ratio=max_ratio,
        )

        # Find best blends
        results = optimizer.find_best_blends(
            demand_scores=demand_scores,
            scoring_fn=scoring_fn_with_weights,
            peptone_filter=peptone_filter,
            max_components=max_components,
            top_k=top_k,
        )

        logger.info(
            f"Found {len(results)} blend recommendations for strain {strain_id}"
        )

        return results


class CompositionOnlyRecommender:
    """Simplified recommender for when genome priors are not available."""

    def __init__(
        self,
        composition_df: pd.DataFrame,
        config: Optional[dict] = None
    ):
        """Initialize composition-only recommender.

        Args:
            composition_df: Peptone composition DataFrame
            config: Optional configuration
        """
        self.composition_df = composition_df
        self.config = config or {}
        self.extractor = CompositionFeatureExtractor(composition_df, config)
        norm_samples = self.config.get("peptone_filter", None)
        self.supply_scores = self.extractor.compute_supply_scores(
            norm_samples=norm_samples
        )

        logger.info("Initialized CompositionOnlyRecommender (no genome priors)")

    def recommend_general(
        self,
        organism_type: str = "lactobacillus",
        top_k: int = 10
    ) -> pd.DataFrame:
        """Get general recommendations based on organism type.

        Args:
            organism_type: Type of organism (lactobacillus, bifidobacterium, generic)
            top_k: Number of recommendations

        Returns:
            Recommendations DataFrame
        """
        # Use generic demands based on organism type
        from .genome_prior import TAXONOMY_PRIORS

        if organism_type.lower() in ["lactobacillus", "lactic", "lab"]:
            # Use Lactiplantibacillus as reference
            prior = TAXONOMY_PRIORS.get("Lactiplantibacillus", TAXONOMY_PRIORS["_generic"])
        elif organism_type.lower() in ["bifidobacterium", "bifido"]:
            prior = TAXONOMY_PRIORS.get("Bifidobacterium", TAXONOMY_PRIORS["_generic"])
        else:
            prior = TAXONOMY_PRIORS["_generic"]

        # Convert to demand
        demand = {}
        for aa, val in prior.get("aa_biosynthesis", {}).items():
            demand[f"demand_{aa}"] = 1.0 - val
        for vit, val in prior.get("vitamin_biosynthesis", {}).items():
            demand[f"demand_vitamin_{vit}"] = 1.0 - val
        demand["demand_nucleotide"] = 1.0 - prior.get("nucleotide_biosynthesis", 0.5)
        demand["transporter_bonus"] = prior.get("transporter_score", 0.5)

        # Score peptones
        results = []
        for peptone_name in self.supply_scores.index:
            supply = self.supply_scores.loc[peptone_name]
            score = self._simple_score(supply, demand)
            results.append({
                "peptone": peptone_name,
                "score": score,
                "organism_type": organism_type
            })

        df = pd.DataFrame(results)
        df = df.sort_values("score", ascending=False).head(top_k)
        df = df.reset_index(drop=True)
        df["rank"] = df.index + 1

        return df

    def _simple_score(
        self,
        supply: pd.Series,
        demand: dict[str, float]
    ) -> float:
        """Compute simple match score.

        Args:
            supply: Supply scores
            demand: Demand scores

        Returns:
            Match score
        """
        score = 0

        # FAA and TAA base scores
        score += supply.get("supply_faa", 0) * 1.0
        score += supply.get("supply_taa", 0) * 0.8

        # MW preference
        transporter = demand.get("transporter_bonus", 0.5)
        score += supply.get("supply_low_mw", 0) * 1.2 * (1 + transporter)
        score += supply.get("supply_medium_mw", 0) * transporter

        # Vitamin and nucleotide
        score += supply.get("supply_vitamin", 0) * 0.5
        score += supply.get("supply_nucleotide", 0) * 0.4

        return score * 20  # Scale to ~0-100

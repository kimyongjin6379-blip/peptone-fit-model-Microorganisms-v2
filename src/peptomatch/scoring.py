"""Scoring and recommendation logic for peptone-strain matching."""

import logging
from typing import Any, Optional

import numpy as np
import pandas as pd

from .composition_features import CompositionFeatureExtractor
from .genome_prior import GenomePriorBuilder

logger = logging.getLogger("peptone_recommender")


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

        # Compute supply scores once
        self.supply_scores = self.composition_extractor.compute_supply_scores()

        # Get weights from config
        self.weights = self.config.get("weights", {})

        logger.info("Initialized PeptoneRecommender")

    def recommend(
        self,
        strain_id: int,
        top_k: int = 10,
        min_score: float = 0.0,
        peptone_filter: Optional[list[str]] = None
    ) -> pd.DataFrame:
        """Get top-K peptone recommendations for a strain.

        Args:
            strain_id: Strain ID from strain table
            top_k: Number of recommendations to return
            min_score: Minimum score threshold
            peptone_filter: If provided, only consider these peptone names

        Returns:
            DataFrame with recommendations (peptone, score, components)
        """
        # Get strain info
        strain_info = self.strain_df[self.strain_df["strain_id"] == strain_id]

        if strain_info.empty:
            raise ValueError(f"Strain ID {strain_id} not found")

        strain_row = strain_info.iloc[0]
        logger.info(
            f"Computing recommendations for strain {strain_id}: "
            f"{strain_row.get('genus', '')} {strain_row.get('species', '')}"
        )

        # Get demand scores for this strain
        demand_scores = self.genome_prior_builder.get_demand_scores(strain_id)

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
            score, components = self._compute_match_score(supply, demand_scores)

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
        demand: dict[str, float]
    ) -> tuple[float, dict[str, float]]:
        """Compute match score between peptone supply and strain demand.

        Args:
            supply: Supply scores for a peptone
            demand: Demand scores for a strain

        Returns:
            Tuple of (total_score, component_scores)
        """
        components = {}

        # Weight defaults
        w = {
            "faa_abundance": self.weights.get("faa_abundance", 1.0),
            "taa_abundance": self.weights.get("taa_abundance", 0.8),
            "mw_low": self.weights.get("mw_low", 1.2),
            "mw_medium": self.weights.get("mw_medium", 1.0),
            "mw_high": self.weights.get("mw_high", 0.6),
            "vitamin_b": self.weights.get("vitamin_b", 0.5),
            "nucleotides": self.weights.get("nucleotides", 0.4),
            "aa_match": self.weights.get("aa_biosynthesis_gap", 1.5),
            "transporter": self.weights.get("transporter_bonus", 0.5),
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

        # 4. AA-specific matching (supply essential AAs the strain can't make)
        essential_aa = ["His", "Ile", "Leu", "Lys", "Met", "Phe", "Thr", "Trp", "Val"]
        aa_match_score = 0

        for aa in essential_aa:
            supply_key = f"supply_{aa}"
            demand_key = f"demand_{aa}"

            if supply_key in supply.index and demand_key in demand:
                # High demand + high supply = good match
                aa_match = supply[supply_key] * demand[demand_key]
                aa_match_score += aa_match

        aa_match_score = (aa_match_score / len(essential_aa)) * w["aa_match"]
        components["aa_match"] = aa_match_score

        # 5. Transporter bonus for peptide-rich peptones
        # Strains with more transporters benefit more from medium MW peptides
        transporter_score = transporter_bonus * (
            supply.get("supply_medium_mw", 0) + supply.get("supply_low_mw", 0) * 0.5
        ) * w["transporter"]
        components["transporter_bonus"] = transporter_score

        # Total score
        total_score = sum(components.values())

        # Normalize to 0-100 scale
        max_possible = sum(w.values()) * 2  # Rough estimate
        normalized_score = (total_score / max_possible) * 100

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
        self.supply_scores = self.extractor.compute_supply_scores()

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

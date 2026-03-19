"""Blend Optimizer for PeptoMatch.

Optimizes peptone blend ratios (2-3 components) using scipy,
integrated with PeptoMatch's genome-driven scoring system.
"""

import logging
from dataclasses import dataclass, field
from itertools import combinations
from typing import Optional

import numpy as np
import pandas as pd
from scipy.optimize import minimize

logger = logging.getLogger("peptone_recommender")


@dataclass
class BlendResult:
    """Result of a single blend optimization."""

    peptones: list[str]
    ratios: list[float]
    score: float
    single_scores: dict[str, float]
    synergy: float
    complementarity: dict[str, float]
    optimization_success: bool
    description: str = ""

    def __post_init__(self):
        if not self.description:
            parts = [f"{p} {r*100:.0f}%" for p, r in zip(self.peptones, self.ratios)]
            self.description = " + ".join(parts)


class BlendOptimizer:
    """Optimize peptone blends for a target strain.

    Uses the same scoring logic as PeptoneRecommender but evaluates
    blended peptone profiles (weighted sum of composition features).
    """

    def __init__(
        self,
        supply_scores: pd.DataFrame,
        min_ratio: float = 0.1,
        max_ratio: float = 0.8,
    ):
        """Initialize blend optimizer.

        Args:
            supply_scores: Supply score DataFrame from CompositionFeatureExtractor
            min_ratio: Minimum ratio per component (default 10%)
            max_ratio: Maximum ratio per component (default 80%)
        """
        self.supply_scores = supply_scores
        self.min_ratio = min_ratio
        self.max_ratio = max_ratio

    def optimize_blend(
        self,
        peptone_names: list[str],
        demand_scores: dict[str, float],
        scoring_fn,
        n_restarts: int = 5,
    ) -> BlendResult:
        """Optimize blend ratios for given peptones to maximize match score.

        Args:
            peptone_names: List of peptone names (2-3)
            demand_scores: Strain demand scores dict
            scoring_fn: Function(supply_series, demand_dict) -> (score, components)
            n_restarts: Number of random restarts for optimization

        Returns:
            BlendResult with optimal ratios and score
        """
        n = len(peptone_names)
        if n < 2 or n > 3:
            raise ValueError("Blend requires 2-3 peptones")

        # Validate peptone names exist
        for name in peptone_names:
            if name not in self.supply_scores.index:
                raise ValueError(f"Peptone '{name}' not found in supply scores")

        # Get individual supply vectors
        supply_vectors = {
            name: self.supply_scores.loc[name] for name in peptone_names
        }

        # Get individual scores for reference
        single_scores = {}
        for name, supply in supply_vectors.items():
            score, _ = scoring_fn(supply, demand_scores)
            single_scores[name] = score

        # Objective: maximize score (minimize negative score)
        def objective(ratios):
            blended_supply = self._blend_supply(supply_vectors, ratios)
            score, _ = scoring_fn(blended_supply, demand_scores)
            return -score

        # Constraints and bounds
        constraints = [{"type": "eq", "fun": lambda x: np.sum(x) - 1.0}]
        bounds = [(self.min_ratio, self.max_ratio)] * n

        # Multi-start optimization
        best_result = None
        best_score = -np.inf

        # Start with equal ratios + random restarts
        starts = [np.array([1.0 / n] * n)]
        rng = np.random.RandomState(42)
        for _ in range(n_restarts - 1):
            raw = rng.dirichlet(np.ones(n))
            # Clip to bounds
            clipped = np.clip(raw, self.min_ratio, self.max_ratio)
            clipped /= clipped.sum()
            starts.append(clipped)

        for x0 in starts:
            try:
                result = minimize(
                    objective,
                    x0,
                    method="SLSQP",
                    bounds=bounds,
                    constraints=constraints,
                    options={"maxiter": 500, "ftol": 1e-8},
                )
                if result.success and -result.fun > best_score:
                    best_score = -result.fun
                    best_result = result
            except Exception as e:
                logger.debug(f"Optimization restart failed: {e}")
                continue

        if best_result is None:
            # Fallback to equal ratios
            ratios = [1.0 / n] * n
            blended = self._blend_supply(supply_vectors, ratios)
            score, _ = scoring_fn(blended, demand_scores)
            success = False
        else:
            ratios = best_result.x.tolist()
            score = best_score
            success = best_result.success

        # Normalize ratios (ensure sum=1)
        ratio_sum = sum(ratios)
        ratios = [r / ratio_sum for r in ratios]

        # Round to 1 decimal %
        ratios = [round(r, 3) for r in ratios]
        # Fix rounding error
        diff = 1.0 - sum(ratios)
        ratios[0] += diff

        # Calculate synergy (blend score vs best single score)
        best_single = max(single_scores.values())
        synergy = score - best_single

        # Calculate complementarity details
        complementarity = self._compute_complementarity(
            supply_vectors, peptone_names
        )

        return BlendResult(
            peptones=peptone_names,
            ratios=ratios,
            score=score,
            single_scores=single_scores,
            synergy=synergy,
            complementarity=complementarity,
            optimization_success=success,
        )

    def find_best_blends(
        self,
        demand_scores: dict[str, float],
        scoring_fn,
        peptone_filter: Optional[list[str]] = None,
        max_components: int = 3,
        top_k: int = 5,
    ) -> list[BlendResult]:
        """Find the best peptone blends from all possible combinations.

        Args:
            demand_scores: Strain demand scores
            scoring_fn: Scoring function
            peptone_filter: Optional list of peptone names to consider
            max_components: Maximum blend components (2 or 3)
            top_k: Number of top blends to return

        Returns:
            List of BlendResult sorted by score (descending)
        """
        available = self.supply_scores.index.tolist()
        if peptone_filter:
            available = [p for p in available if p in peptone_filter]

        if len(available) < 2:
            logger.warning("Not enough peptones for blending")
            return []

        all_results: list[BlendResult] = []

        # Try 2-component blends
        for combo in combinations(available, 2):
            try:
                result = self.optimize_blend(
                    list(combo), demand_scores, scoring_fn
                )
                all_results.append(result)
            except Exception as e:
                logger.debug(f"Failed to optimize {combo}: {e}")

        # Try 3-component blends if allowed
        if max_components >= 3 and len(available) >= 3:
            for combo in combinations(available, 3):
                try:
                    result = self.optimize_blend(
                        list(combo), demand_scores, scoring_fn
                    )
                    all_results.append(result)
                except Exception as e:
                    logger.debug(f"Failed to optimize {combo}: {e}")

        # Sort by score descending
        all_results.sort(key=lambda r: r.score, reverse=True)

        return all_results[:top_k]

    def _blend_supply(
        self,
        supply_vectors: dict[str, pd.Series],
        ratios,
    ) -> pd.Series:
        """Create blended supply score as weighted sum.

        Args:
            supply_vectors: Dict of peptone_name -> supply Series
            ratios: Mixing ratios (array-like)

        Returns:
            Blended supply Series
        """
        names = list(supply_vectors.keys())
        blended = pd.Series(0.0, index=supply_vectors[names[0]].index)

        for name, ratio in zip(names, ratios):
            blended += supply_vectors[name] * ratio

        return blended

    def _compute_complementarity(
        self,
        supply_vectors: dict[str, pd.Series],
        peptone_names: list[str],
    ) -> dict[str, float]:
        """Compute complementarity metrics between peptones.

        Args:
            supply_vectors: Supply score vectors
            peptone_names: Peptone names

        Returns:
            Complementarity metrics dict
        """
        metrics = {}

        vectors = [supply_vectors[n].values for n in peptone_names]

        # 1. Profile diversity (mean pairwise distance)
        distances = []
        for i in range(len(vectors)):
            for j in range(i + 1, len(vectors)):
                # Use only non-NaN values
                v1 = np.nan_to_num(vectors[i])
                v2 = np.nan_to_num(vectors[j])
                dist = np.linalg.norm(v1 - v2)
                distances.append(dist)
        metrics["profile_diversity"] = float(np.mean(distances)) if distances else 0

        # 2. Coverage score (do the peptones collectively cover weak areas?)
        stacked = np.array([np.nan_to_num(v) for v in vectors])
        max_per_feature = stacked.max(axis=0)
        min_per_feature = stacked.min(axis=0)

        # High coverage = each feature has at least one good supplier
        nonzero_mask = max_per_feature > 0
        if nonzero_mask.any():
            metrics["feature_coverage"] = float(
                np.mean(max_per_feature[nonzero_mask] > 0.3)
            )
        else:
            metrics["feature_coverage"] = 0

        # 3. Gap-filling score (does one peptone fill another's weakness?)
        gap_fill_scores = []
        for i in range(len(vectors)):
            weak = np.nan_to_num(vectors[i]) < 0.2
            if not weak.any():
                continue
            others = [np.nan_to_num(vectors[j]) for j in range(len(vectors)) if j != i]
            if others:
                others_max = np.max(others, axis=0)
                gap_fill = np.mean(others_max[weak])
                gap_fill_scores.append(gap_fill)
        metrics["gap_filling"] = float(np.mean(gap_fill_scores)) if gap_fill_scores else 0

        # 4. Overall complementarity (composite)
        metrics["overall"] = (
            metrics["profile_diversity"] * 0.3
            + metrics["feature_coverage"] * 0.4
            + metrics["gap_filling"] * 0.3
        )

        return metrics

    def get_blend_composition(
        self,
        blend_result: BlendResult,
        composition_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """Get the actual blended composition values.

        Args:
            blend_result: BlendResult from optimization
            composition_df: Raw composition DataFrame

        Returns:
            DataFrame with individual + blended composition
        """
        rows = []
        for name, ratio in zip(blend_result.peptones, blend_result.ratios):
            if name in composition_df.index:
                row = composition_df.loc[name].copy()
                rows.append(row)

        if not rows:
            return pd.DataFrame()

        individual_df = pd.DataFrame(rows)

        # Compute blended values (numeric columns only)
        numeric_cols = individual_df.select_dtypes(include=[np.number]).columns
        blended = pd.Series(0.0, index=numeric_cols, name="BLEND")

        for i, (name, ratio) in enumerate(
            zip(blend_result.peptones, blend_result.ratios)
        ):
            if i < len(rows):
                blended += individual_df.iloc[i][numeric_cols].fillna(0) * ratio

        # Add blend row
        result_df = individual_df.copy()
        result_df.loc["BLEND"] = blended

        return result_df

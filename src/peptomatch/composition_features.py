"""Extract and compute features from peptone composition data."""

import logging
import re
from typing import Any, Optional

import numpy as np
import pandas as pd

logger = logging.getLogger("peptone_recommender")


class CompositionFeatureExtractor:
    """Extract meaningful features from raw peptone composition data."""

    # Column prefix patterns for feature categories
    CATEGORY_PATTERNS = {
        "faa": r"^faa_",           # Free amino acids
        "taa": r"^taa_",           # Total amino acids
        "mw": r"^mw_",             # Molecular weight distribution
        "mineral": r"^mineral_",    # Minerals
        "vitamin": r"^vitb_",       # Vitamin B
        "nucleotide": r"^nucleotide_",  # Nucleotides
        "sugar": r"^sugar_",        # Sugars
        "orgacid": r"^orgacid_",    # Organic acids
        "general": r"^general_",    # General properties
    }

    # Amino acid name standardization
    AA_NAME_MAP = {
        "aspartic acid": "Asp",
        "glutamic acid": "Glu",
        "asparagine": "Asn",
        "glutamine": "Gln",
        "histidine": "His",
        "isoleucine": "Ile",
        "leucine": "Leu",
        "lysine": "Lys",
        "methionine": "Met",
        "phenylalanine": "Phe",
        "threonine": "Thr",
        "tryptophan": "Trp",
        "valine": "Val",
        "alanine": "Ala",
        "arginine": "Arg",
        "cysteine": "Cys",
        "cystine": "Cys2",
        "glycine": "Gly",
        "proline": "Pro",
        "serine": "Ser",
        "tyrosine": "Tyr",
        "hydroxyproline": "Hyp",
        "citruline": "Cit",
        "ornithine": "Orn",
        "gaba": "GABA",
    }

    def __init__(self, composition_df: pd.DataFrame, config: Optional[dict] = None):
        """Initialize feature extractor.

        Args:
            composition_df: Raw composition DataFrame (Sample_name as index)
            config: Optional configuration dictionary
        """
        self.raw_df = composition_df
        self.config = config or {}

        # Categorize columns
        self.column_categories = self._categorize_columns()

        logger.info(f"Initialized feature extractor with {len(self.raw_df)} samples")

    def _categorize_columns(self) -> dict[str, list[str]]:
        """Categorize columns by their prefix.

        Returns:
            Dictionary mapping category to list of column names
        """
        categories = {cat: [] for cat in self.CATEGORY_PATTERNS}

        for col in self.raw_df.columns:
            col_lower = str(col).lower()
            for cat, pattern in self.CATEGORY_PATTERNS.items():
                if re.match(pattern, col_lower):
                    categories[cat].append(col)
                    break

        # Log category counts
        for cat, cols in categories.items():
            if cols:
                logger.debug(f"Category '{cat}': {len(cols)} columns")

        return categories

    def get_columns_by_category(self, category: str) -> list[str]:
        """Get column names for a specific category.

        Args:
            category: Category name (faa, taa, mw, mineral, vitamin, nucleotide, sugar, orgacid, general)

        Returns:
            List of column names in that category
        """
        return self.column_categories.get(category, [])

    def compute_faa_features(self) -> pd.DataFrame:
        """Compute free amino acid abundance features.

        Returns:
            DataFrame with FAA-derived features
        """
        faa_cols = self.get_columns_by_category("faa")

        if not faa_cols:
            logger.warning("No FAA columns found")
            return pd.DataFrame(index=self.raw_df.index)

        faa_df = self.raw_df[faa_cols]

        features = pd.DataFrame(index=self.raw_df.index)

        # Total FAA
        features["faa_total"] = faa_df.sum(axis=1)

        # Mean FAA (excluding zeros)
        features["faa_mean"] = faa_df.replace(0, np.nan).mean(axis=1)

        # Essential AA sum (if identifiable)
        essential_aa = ["Ile", "Leu", "Lys", "Met", "Phe", "Thr", "Trp", "Val", "His"]
        essential_cols = [c for c in faa_cols if self._extract_aa_code(c) in essential_aa]
        if essential_cols:
            features["faa_essential"] = self.raw_df[essential_cols].sum(axis=1)

        # Aromatic AA (Phe, Tyr, Trp)
        aromatic_aa = ["Phe", "Tyr", "Trp"]
        aromatic_cols = [c for c in faa_cols if self._extract_aa_code(c) in aromatic_aa]
        if aromatic_cols:
            features["faa_aromatic"] = self.raw_df[aromatic_cols].sum(axis=1)

        # BCAA (Leu, Ile, Val)
        bcaa = ["Leu", "Ile", "Val"]
        bcaa_cols = [c for c in faa_cols if self._extract_aa_code(c) in bcaa]
        if bcaa_cols:
            features["faa_bcaa"] = self.raw_df[bcaa_cols].sum(axis=1)

        # Sulfur AA (Met, Cys)
        sulfur_aa = ["Met", "Cys", "Cys2"]
        sulfur_cols = [c for c in faa_cols if self._extract_aa_code(c) in sulfur_aa]
        if sulfur_cols:
            features["faa_sulfur"] = self.raw_df[sulfur_cols].sum(axis=1)

        # Individual AA for later matching
        for col in faa_cols:
            aa_code = self._extract_aa_code(col)
            if aa_code:
                features[f"faa_{aa_code}"] = self.raw_df[col]

        return features

    def compute_taa_features(self) -> pd.DataFrame:
        """Compute total amino acid abundance features.

        Returns:
            DataFrame with TAA-derived features
        """
        taa_cols = self.get_columns_by_category("taa")

        if not taa_cols:
            logger.warning("No TAA columns found")
            return pd.DataFrame(index=self.raw_df.index)

        taa_df = self.raw_df[taa_cols]

        features = pd.DataFrame(index=self.raw_df.index)

        # Total TAA
        features["taa_total"] = taa_df.sum(axis=1)

        # Mean TAA
        features["taa_mean"] = taa_df.replace(0, np.nan).mean(axis=1)

        # Essential AA
        essential_aa = ["Ile", "Leu", "Lys", "Met", "Phe", "Thr", "Trp", "Val", "His"]
        essential_cols = [c for c in taa_cols if self._extract_aa_code(c) in essential_aa]
        if essential_cols:
            features["taa_essential"] = self.raw_df[essential_cols].sum(axis=1)

        # Individual AA
        for col in taa_cols:
            aa_code = self._extract_aa_code(col)
            if aa_code:
                features[f"taa_{aa_code}"] = self.raw_df[col]

        return features

    def compute_mw_features(self) -> pd.DataFrame:
        """Compute molecular weight distribution features.

        Returns:
            DataFrame with MW-derived features
        """
        mw_cols = self.get_columns_by_category("mw")

        if not mw_cols:
            logger.warning("No MW columns found")
            return pd.DataFrame(index=self.raw_df.index)

        features = pd.DataFrame(index=self.raw_df.index)

        # Average MW if present
        avg_col = [c for c in mw_cols if "avg" in c.lower()]
        if avg_col:
            features["mw_avg"] = self.raw_df[avg_col[0]]

        # Percentage columns
        for col in mw_cols:
            col_lower = col.lower()
            if "lt250" in col_lower or "<250" in col_lower:
                features["mw_pct_low"] = self.raw_df[col]
            elif "gt1000" in col_lower or ">1000" in col_lower:
                features["mw_pct_high"] = self.raw_df[col]
            elif "pct" in col_lower:
                # Sum intermediate fractions as "medium"
                if "mw_pct_medium" not in features.columns:
                    features["mw_pct_medium"] = self.raw_df[col]
                else:
                    features["mw_pct_medium"] += self.raw_df[col]

        # Calculate low MW preference score (higher = more small peptides/FAA)
        if "mw_pct_low" in features.columns and "mw_pct_high" in features.columns:
            features["mw_low_high_ratio"] = (
                features["mw_pct_low"] / (features["mw_pct_high"] + 0.1)
            )

        return features

    def compute_vitamin_features(self) -> pd.DataFrame:
        """Compute vitamin B features.

        Returns:
            DataFrame with vitamin-derived features
        """
        vit_cols = self.get_columns_by_category("vitamin")

        if not vit_cols:
            logger.warning("No vitamin columns found")
            return pd.DataFrame(index=self.raw_df.index)

        features = pd.DataFrame(index=self.raw_df.index)

        # Individual vitamins
        for col in vit_cols:
            col_lower = col.lower()
            for b_num in ["b1", "b2", "b3", "b6", "b9", "b12"]:
                if b_num in col_lower:
                    features[f"vitamin_{b_num}"] = self.raw_df[col]
                    break

        # Total vitamin B score (normalized)
        vit_df = self.raw_df[vit_cols]
        features["vitamin_total"] = vit_df.sum(axis=1)

        return features

    def compute_mineral_features(self) -> pd.DataFrame:
        """Compute mineral features.

        Returns:
            DataFrame with mineral-derived features
        """
        min_cols = self.get_columns_by_category("mineral")

        if not min_cols:
            logger.warning("No mineral columns found")
            return pd.DataFrame(index=self.raw_df.index)

        features = pd.DataFrame(index=self.raw_df.index)

        # Individual minerals
        for col in min_cols:
            col_lower = col.lower()
            if "na" in col_lower:
                features["mineral_Na"] = self.raw_df[col]
            elif "k" in col_lower:
                features["mineral_K"] = self.raw_df[col]
            elif "mg" in col_lower:
                features["mineral_Mg"] = self.raw_df[col]
            elif "ca" in col_lower:
                features["mineral_Ca"] = self.raw_df[col]

        # Total minerals
        features["mineral_total"] = self.raw_df[min_cols].sum(axis=1)

        return features

    def compute_nucleotide_features(self) -> pd.DataFrame:
        """Compute nucleotide features.

        Returns:
            DataFrame with nucleotide-derived features
        """
        nuc_cols = self.get_columns_by_category("nucleotide")

        if not nuc_cols:
            logger.warning("No nucleotide columns found")
            return pd.DataFrame(index=self.raw_df.index)

        features = pd.DataFrame(index=self.raw_df.index)

        # Total nucleotides
        features["nucleotide_total"] = self.raw_df[nuc_cols].sum(axis=1)

        # Individual nucleotides
        for col in nuc_cols:
            col_lower = col.lower()
            for nuc in ["amp", "gmp", "ump", "imp", "cmp"]:
                if nuc in col_lower:
                    features[f"nucleotide_{nuc.upper()}"] = self.raw_df[col]
                    break

        return features

    def compute_general_features(self) -> pd.DataFrame:
        """Compute general composition features (TN, AN, ash, etc.).

        Returns:
            DataFrame with general features
        """
        gen_cols = self.get_columns_by_category("general")

        if not gen_cols:
            logger.warning("No general columns found")
            return pd.DataFrame(index=self.raw_df.index)

        features = pd.DataFrame(index=self.raw_df.index)

        for col in gen_cols:
            col_lower = col.lower()
            if "tn" in col_lower or "total_n" in col_lower:
                features["general_TN"] = self.raw_df[col]
            elif "an" in col_lower or "amino_n" in col_lower:
                features["general_AN"] = self.raw_df[col]
            elif "ash" in col_lower:
                features["general_ash"] = self.raw_df[col]
            elif "moisture" in col_lower:
                features["general_moisture"] = self.raw_df[col]
            elif "sugar" in col_lower:
                features["general_sugar"] = self.raw_df[col]

        # AN/TN ratio (indicator of hydrolysis degree)
        if "general_TN" in features.columns and "general_AN" in features.columns:
            features["general_AN_TN_ratio"] = (
                features["general_AN"] / (features["general_TN"] + 0.01)
            )

        return features

    def compute_all_features(self) -> pd.DataFrame:
        """Compute all features and return combined DataFrame.

        Returns:
            DataFrame with all computed features
        """
        feature_dfs = [
            self.compute_faa_features(),
            self.compute_taa_features(),
            self.compute_mw_features(),
            self.compute_vitamin_features(),
            self.compute_mineral_features(),
            self.compute_nucleotide_features(),
            self.compute_general_features(),
        ]

        # Combine all features
        all_features = pd.concat(feature_dfs, axis=1)

        # Remove duplicate columns if any
        all_features = all_features.loc[:, ~all_features.columns.duplicated()]

        logger.info(f"Computed {len(all_features.columns)} total features")

        return all_features

    def compute_supply_scores(self, weights: Optional[dict] = None) -> pd.DataFrame:
        """Compute peptone supply scores for matching.

        Args:
            weights: Optional weight dictionary for different feature categories

        Returns:
            DataFrame with supply scores per peptone
        """
        weights = weights or {}

        features = self.compute_all_features()
        scores = pd.DataFrame(index=features.index)

        # FAA abundance score (normalized)
        if "faa_total" in features.columns:
            faa_max = features["faa_total"].max()
            if faa_max > 0:
                scores["supply_faa"] = features["faa_total"] / faa_max
            else:
                scores["supply_faa"] = 0

        # TAA abundance score
        if "taa_total" in features.columns:
            taa_max = features["taa_total"].max()
            if taa_max > 0:
                scores["supply_taa"] = features["taa_total"] / taa_max
            else:
                scores["supply_taa"] = 0

        # MW distribution score (favor low MW)
        if "mw_pct_low" in features.columns:
            scores["supply_low_mw"] = features["mw_pct_low"] / 100

        if "mw_pct_medium" in features.columns:
            scores["supply_medium_mw"] = features["mw_pct_medium"] / 100

        if "mw_pct_high" in features.columns:
            scores["supply_high_mw"] = features["mw_pct_high"] / 100

        # Vitamin score
        if "vitamin_total" in features.columns:
            vit_max = features["vitamin_total"].max()
            if vit_max > 0:
                scores["supply_vitamin"] = features["vitamin_total"] / vit_max
            else:
                scores["supply_vitamin"] = 0

        # Nucleotide score
        if "nucleotide_total" in features.columns:
            nuc_max = features["nucleotide_total"].max()
            if nuc_max > 0:
                scores["supply_nucleotide"] = features["nucleotide_total"] / nuc_max
            else:
                scores["supply_nucleotide"] = 0

        # Individual AA scores (for specific deficiency matching)
        all_aa = ["His", "Ile", "Leu", "Lys", "Met", "Phe", "Thr", "Trp", "Val",
                  "Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "Pro", "Ser", "Tyr"]
        for aa in all_aa:
            faa_col = f"faa_{aa}"
            if faa_col in features.columns:
                aa_max = features[faa_col].max()
                if aa_max > 0:
                    scores[f"supply_{aa}"] = features[faa_col] / aa_max
                else:
                    scores[f"supply_{aa}"] = 0

        return scores

    def _extract_aa_code(self, column_name: str) -> Optional[str]:
        """Extract amino acid code from column name.

        Args:
            column_name: Column name like 'faa_Leucine' or 'taa_Aspartic acid'

        Returns:
            Standardized AA code or None
        """
        # Remove prefix
        name = re.sub(r'^(faa_|taa_)', '', str(column_name).lower())

        # Look up in mapping
        return self.AA_NAME_MAP.get(name.strip())

    def get_peptone_profile(self, sample_name: str) -> dict[str, Any]:
        """Get detailed composition profile for a specific peptone.

        Args:
            sample_name: Name of the peptone sample

        Returns:
            Dictionary with composition profile
        """
        if sample_name not in self.raw_df.index:
            raise ValueError(f"Sample '{sample_name}' not found")

        features = self.compute_all_features()
        row = features.loc[sample_name]

        profile = {
            "sample_name": sample_name,
            "faa": {k.replace("faa_", ""): v for k, v in row.items() if k.startswith("faa_")},
            "taa": {k.replace("taa_", ""): v for k, v in row.items() if k.startswith("taa_")},
            "mw": {k.replace("mw_", ""): v for k, v in row.items() if k.startswith("mw_")},
            "vitamins": {k.replace("vitamin_", ""): v for k, v in row.items() if k.startswith("vitamin_")},
            "minerals": {k.replace("mineral_", ""): v for k, v in row.items() if k.startswith("mineral_")},
            "nucleotides": {k.replace("nucleotide_", ""): v for k, v in row.items() if k.startswith("nucleotide_")},
            "general": {k.replace("general_", ""): v for k, v in row.items() if k.startswith("general_")},
        }

        return profile

    def filter_peptones(
        self,
        min_faa: Optional[float] = None,
        min_taa: Optional[float] = None,
        max_mw_avg: Optional[float] = None,
        min_aa: Optional[dict[str, float]] = None
    ) -> pd.DataFrame:
        """Filter peptones based on criteria.

        Args:
            min_faa: Minimum total FAA value
            min_taa: Minimum total TAA value
            max_mw_avg: Maximum average MW
            min_aa: Dict of {AA_code: min_value} for specific AA requirements

        Returns:
            Filtered DataFrame with sample names
        """
        features = self.compute_all_features()
        mask = pd.Series(True, index=features.index)

        if min_faa is not None and "faa_total" in features.columns:
            mask &= features["faa_total"] >= min_faa

        if min_taa is not None and "taa_total" in features.columns:
            mask &= features["taa_total"] >= min_taa

        if max_mw_avg is not None and "mw_avg" in features.columns:
            mask &= features["mw_avg"] <= max_mw_avg

        if min_aa:
            for aa, min_val in min_aa.items():
                col = f"faa_{aa}"
                if col in features.columns:
                    mask &= features[col] >= min_val

        return features[mask]

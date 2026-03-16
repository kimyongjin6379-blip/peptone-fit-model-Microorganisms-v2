"""Data loading functions for Excel files."""

import logging
import re
from pathlib import Path
from typing import Any, Optional, Union

import numpy as np
import pandas as pd

from .utils import clean_numeric_value, find_column_by_pattern, safe_path

logger = logging.getLogger("peptone_recommender")


def load_composition_data(
    filepath: Union[str, Path],
    sheet_name: str = "data",
    non_numeric_handling: str = "zero",
    default_loq: float = 0.1
) -> pd.DataFrame:
    """Load and clean peptone composition data from Excel.

    Args:
        filepath: Path to the composition Excel file
        sheet_name: Name of the sheet containing data
        non_numeric_handling: Method for handling non-numeric values ('zero', 'half_loq', 'nan')
        default_loq: Default LOQ value when using half_loq method

    Returns:
        Cleaned DataFrame with Sample_name as index and numeric features as columns

    Raises:
        FileNotFoundError: If the file doesn't exist
        ValueError: If required columns are missing
    """
    filepath = safe_path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"Composition file not found: {filepath}")

    logger.info(f"Loading composition data from {filepath}")

    try:
        df = pd.read_excel(filepath, sheet_name=sheet_name)
    except Exception as e:
        raise ValueError(f"Failed to read sheet '{sheet_name}' from {filepath}: {e}")

    # Find Sample_name column
    sample_col = find_column_by_pattern(
        df.columns.tolist(),
        [r'sample_?name', r'샘플명', r'펩톤명', r'peptone_?name']
    )

    if sample_col is None:
        # Fall back to third column (index 2) based on observed structure
        if len(df.columns) > 2:
            sample_col = df.columns[2]
            logger.warning(f"Sample_name column not found, using column '{sample_col}'")
        else:
            raise ValueError("Could not identify Sample_name column in composition data")

    # Set Sample_name as index
    df = df.set_index(sample_col)
    df.index.name = "Sample_name"

    # Identify numeric feature columns (exclude metadata columns)
    metadata_patterns = [
        r'^sample_?id$', r'^material_?type$', r'^raw_?material$',
        r'^manufacturer$', r'^제조사$', r'^원료$'
    ]

    feature_columns = []
    for col in df.columns:
        col_str = str(col).lower()
        is_metadata = any(re.match(p, col_str, re.IGNORECASE) for p in metadata_patterns)
        if not is_metadata:
            feature_columns.append(col)

    logger.info(f"Found {len(feature_columns)} potential feature columns")

    # Clean numeric values
    df_clean = df[feature_columns].copy()

    for col in df_clean.columns:
        df_clean[col] = df_clean[col].apply(
            lambda x: clean_numeric_value(x, non_numeric_handling, default_loq)
        )

    # Convert to numeric, coercing errors
    df_clean = df_clean.apply(pd.to_numeric, errors='coerce')

    # Fill remaining NaN with 0 if method is 'zero'
    if non_numeric_handling == "zero":
        df_clean = df_clean.fillna(0)

    # Remove rows with all NaN
    df_clean = df_clean.dropna(how='all')

    logger.info(f"Loaded {len(df_clean)} peptone samples with {len(df_clean.columns)} features")

    return df_clean


def load_strain_table(
    filepath: Union[str, Path],
    sheet_name: Optional[str] = None
) -> pd.DataFrame:
    """Load and parse strain list from Excel with complex headers.

    Args:
        filepath: Path to the strain list Excel file
        sheet_name: Name of the sheet to read. If None, auto-detect.

    Returns:
        DataFrame with columns: strain_id, genus, species, strain_name, GCF, notes

    Raises:
        FileNotFoundError: If the file doesn't exist
        ValueError: If strain data cannot be parsed
    """
    filepath = safe_path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"Strain file not found: {filepath}")

    logger.info(f"Loading strain table from {filepath}")

    # Load Excel file
    xl = pd.ExcelFile(filepath)

    # Auto-detect sheet if not specified
    if sheet_name is None:
        # Try to find a non-empty sheet
        for sheet in xl.sheet_names:
            df_test = pd.read_excel(xl, sheet_name=sheet, header=None, nrows=30)
            if not df_test.empty and df_test.shape[1] > 5:
                sheet_name = sheet
                break

        if sheet_name is None:
            sheet_name = xl.sheet_names[0]

    logger.info(f"Using sheet: {sheet_name}")

    # Load raw data without header
    df_raw = pd.read_excel(xl, sheet_name=sheet_name, header=None)

    # Find the header row containing relevant column names
    header_row = None
    gcf_col_idx = None

    # Patterns to find GCF column
    gcf_patterns = [r'ncbi', r'accession', r'gcf', r'assembly']

    for idx, row in df_raw.iterrows():
        for col_idx, cell in enumerate(row):
            cell_str = str(cell).lower() if pd.notna(cell) else ""
            if any(p in cell_str for p in gcf_patterns):
                header_row = idx
                gcf_col_idx = col_idx
                break
        if header_row is not None:
            break

    # Find genus/species columns by looking for "Microorganisms" or actual genus names
    data_start_row = None
    genus_col_idx = None
    species_col_idx = None

    # Common genus patterns for bacteria
    genus_patterns = [
        r'lactiplantibacillus', r'lactobacillus', r'levilactobacillus',
        r'lacticaseibacillus', r'lentilactobacillus', r'latilactobacillus',
        r'bifidobacterium', r'streptococcus', r'bacillus', r'escherichia'
    ]

    for idx in range(max(0, (header_row or 0) - 2), min(len(df_raw), (header_row or 0) + 10)):
        row = df_raw.iloc[idx]
        for col_idx, cell in enumerate(row):
            cell_str = str(cell).lower() if pd.notna(cell) else ""
            if any(re.search(p, cell_str) for p in genus_patterns):
                if data_start_row is None:
                    data_start_row = idx
                    genus_col_idx = col_idx
                    # Species is typically the next column
                    if col_idx + 1 < len(row):
                        species_col_idx = col_idx + 1
                break

    if data_start_row is None:
        # Fallback: look for numbered rows (column 0 contains 1, 2, 3...)
        for idx in range(len(df_raw)):
            cell = df_raw.iloc[idx, 0]
            if pd.notna(cell) and str(cell).strip() == "1":
                data_start_row = idx
                break

    if data_start_row is None:
        raise ValueError("Could not find data start row in strain table")

    logger.info(f"Data starts at row {data_start_row}")

    # Extract data rows
    strain_rows = []
    strain_id_counter = 1

    for idx in range(data_start_row, len(df_raw)):
        row = df_raw.iloc[idx]

        # Check if this is a valid data row (first column should be a number or genus name)
        first_cell = row.iloc[0]
        if pd.isna(first_cell):
            continue

        first_str = str(first_cell).strip()

        # Skip if it's clearly a header/label row
        if any(s in first_str.lower() for s in ['합계', 'total', '기타', 'other']):
            continue

        # Determine strain ID
        try:
            strain_id = int(float(first_str))
        except ValueError:
            strain_id = strain_id_counter

        strain_id_counter = strain_id + 1

        # Extract genus and species
        genus = None
        species = None

        if genus_col_idx is not None and genus_col_idx < len(row):
            genus_val = row.iloc[genus_col_idx]
            if pd.notna(genus_val) and str(genus_val).strip():
                genus = str(genus_val).strip()

        if species_col_idx is not None and species_col_idx < len(row):
            species_val = row.iloc[species_col_idx]
            if pd.notna(species_val) and str(species_val).strip():
                species = str(species_val).strip()

        # If genus is empty but species has value, this might be continuation
        if genus is None and species is None:
            continue

        # Extract strain name (usually in column with KCTC/KCCM patterns)
        strain_name = None
        for col_idx in range(len(row)):
            cell = row.iloc[col_idx]
            if pd.notna(cell):
                cell_str = str(cell).strip()
                if re.search(r'K[A-Z]{2,4}\s*\d+', cell_str, re.IGNORECASE):
                    strain_name = cell_str
                    break

        # Extract GCF
        gcf = None
        if gcf_col_idx is not None and gcf_col_idx < len(row):
            gcf_val = row.iloc[gcf_col_idx]
            if pd.notna(gcf_val):
                gcf_str = str(gcf_val).strip()
                if gcf_str.startswith("GCF_"):
                    gcf = gcf_str

        # Only add if we have at least genus or species
        if genus or species:
            strain_rows.append({
                "strain_id": strain_id,
                "genus": genus,
                "species": species,
                "strain_name": strain_name,
                "GCF": gcf,
                "notes": None
            })

    if not strain_rows:
        raise ValueError("No valid strain data found in file")

    df_strains = pd.DataFrame(strain_rows)

    # Forward fill genus for strains where it's empty (they belong to same genus)
    df_strains['genus'] = df_strains['genus'].ffill()

    # Create combined strain identifier
    df_strains['full_name'] = df_strains.apply(
        lambda r: f"{r['genus'] or ''} {r['species'] or ''} {r['strain_name'] or ''}".strip(),
        axis=1
    )

    logger.info(f"Loaded {len(df_strains)} strains")

    return df_strains


def save_processed_data(
    composition_df: pd.DataFrame,
    strain_df: pd.DataFrame,
    output_dir: Union[str, Path] = "outputs"
) -> dict[str, Path]:
    """Save processed data to CSV files.

    Args:
        composition_df: Cleaned composition DataFrame
        strain_df: Strain table DataFrame
        output_dir: Output directory path

    Returns:
        Dictionary mapping data type to output file path
    """
    output_dir = safe_path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    paths = {}

    # Save composition data
    comp_path = output_dir / "clean_peptone_composition.csv"
    composition_df.to_csv(comp_path, encoding='utf-8-sig')
    paths['composition'] = comp_path
    logger.info(f"Saved composition data to {comp_path}")

    # Save strain table
    strain_path = output_dir / "strain_table.csv"
    strain_df.to_csv(strain_path, index=False, encoding='utf-8-sig')
    paths['strains'] = strain_path
    logger.info(f"Saved strain table to {strain_path}")

    return paths


def validate_data_integrity(
    composition_df: pd.DataFrame,
    strain_df: pd.DataFrame
) -> dict[str, Any]:
    """Validate loaded data for completeness and consistency.

    Args:
        composition_df: Composition data DataFrame
        strain_df: Strain table DataFrame

    Returns:
        Validation report dictionary
    """
    report = {
        "composition": {
            "n_samples": len(composition_df),
            "n_features": len(composition_df.columns),
            "missing_ratio": composition_df.isna().sum().sum() / composition_df.size,
            "samples": composition_df.index.tolist()
        },
        "strains": {
            "n_strains": len(strain_df),
            "n_with_gcf": strain_df['GCF'].notna().sum(),
            "n_missing_gcf": strain_df['GCF'].isna().sum(),
            "genera": strain_df['genus'].dropna().unique().tolist()
        },
        "warnings": []
    }

    # Check for high missing data
    if report["composition"]["missing_ratio"] > 0.3:
        report["warnings"].append(
            f"High missing data ratio in composition: {report['composition']['missing_ratio']:.1%}"
        )

    # Check for strains without GCF
    if report["strains"]["n_missing_gcf"] > 0:
        report["warnings"].append(
            f"{report['strains']['n_missing_gcf']} strains without GCF accession"
        )

    return report

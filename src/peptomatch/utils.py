"""Utility functions for the peptone recommendation system."""

import logging
import re
from pathlib import Path
from typing import Any, Optional, Union

import yaml


def setup_logging(level: str = "INFO") -> logging.Logger:
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    return logging.getLogger("peptomatch")


def load_config(config_path: Optional[Union[str, Path]] = None) -> dict[str, Any]:
    if config_path is None:
        current = Path(__file__).parent
        while current != current.parent:
            if (current / "config" / "config.yaml").exists():
                config_path = current / "config" / "config.yaml"
                break
            current = current.parent
        else:
            config_path = Path.cwd() / "config" / "config.yaml"

    config_path = Path(config_path)
    if not config_path.exists():
        return {}

    with open(config_path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f) or {}

    # Resolve relative paths in data section relative to project root
    project_root = config_path.parent.parent
    for key in ("composition_file", "strain_file"):
        if key in config.get("data", {}):
            p = Path(config["data"][key])
            if not p.is_absolute():
                config["data"][key] = str(project_root / p)

    return config


def get_project_root() -> Path:
    current = Path(__file__).parent
    while current != current.parent:
        if (current / "config" / "config.yaml").exists():
            return current
        current = current.parent
    return Path.cwd()


def safe_path(path: Union[str, Path]) -> Path:
    return Path(path)


def clean_numeric_value(value: Any, method: str = "zero", default_loq: float = 0.1) -> float:
    import numpy as np

    if value is None:
        return 0.0 if method != "nan" else np.nan

    if isinstance(value, (int, float)):
        if np.isnan(value):
            return 0.0 if method != "nan" else np.nan
        return float(value)

    value_str = str(value).strip().upper()

    non_detect_patterns = [
        r'^N\.?D\.?$', r'^<\s*LOQ$', r'^<\s*\d+',
        r'^\-$', r'^$', r'^미량$',
    ]

    for pattern in non_detect_patterns:
        if re.match(pattern, value_str):
            if method == "zero":
                return 0.0
            elif method == "half_loq":
                num_match = re.search(r'<\s*(\d+\.?\d*)', value_str)
                if num_match:
                    return float(num_match.group(1)) / 2
                return default_loq / 2
            else:
                return np.nan

    try:
        return float(value_str.replace(",", ""))
    except ValueError:
        return 0.0 if method != "nan" else np.nan


def normalize_column_name(name: str) -> str:
    if not isinstance(name, str):
        return str(name)
    normalized = name.lower().strip()
    normalized = re.sub(r'\s+', '_', normalized)
    normalized = re.sub(r'[^\w_]', '', normalized)
    return normalized


def find_column_by_pattern(
    columns: list[str],
    patterns: list[str],
    case_sensitive: bool = False
) -> Optional[str]:
    flags = 0 if case_sensitive else re.IGNORECASE
    for col in columns:
        for pattern in patterns:
            if re.search(pattern, str(col), flags):
                return col
    return None


def ensure_output_dir(subdir: Optional[str] = None) -> Path:
    root = get_project_root()
    output_dir = root / "outputs"
    if subdir:
        output_dir = output_dir / subdir
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir

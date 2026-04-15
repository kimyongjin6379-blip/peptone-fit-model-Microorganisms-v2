"""SQLite-based growth curve database for PeptoMatch integration.

Stores growth curve data from growth-curve-app, resolves team aliases,
and computes growth metrics for ML/FBA calibration.
"""

import json
import logging
import sqlite3
from datetime import datetime
from pathlib import Path
from typing import Any, Optional

logger = logging.getLogger("peptomatch")

# ── Alias tables ────────────────────────────────────────────────

STRAIN_ALIASES = {
    "LP": {"genus": "Lactobacillus", "species": "plantarum"},
    "LR": {"genus": "Lactobacillus", "species": "rhamnosus"},
    "LA": {"genus": "Lactobacillus", "species": "acidophilus"},
    "LC": {"genus": "Lactobacillus", "species": "casei"},
    "LPC": {"genus": "Lactobacillus", "species": "paracasei"},
    "LS": {"genus": "Lactobacillus", "species": "salivarius"},
    "EF": {"genus": "Enterococcus", "species": "faecalis"},
    "STT": {"genus": "Streptococcus", "species": "thermophilus"},
    "BS": {"genus": "Bacillus", "species": "subtilis"},
    "BC": {"genus": "Bacillus", "species": "coagulans"},
    "EC": {"genus": "Escherichia", "species": "coli"},
}

PEPTONE_ALIASES = {
    "1": "SOY-1",
    "N": "SOY-N+",
    "L": "SOY-L",
    "B": "SOY-B",
    "SP": "SOY-P",
    "W": "WHEAT-1",
    "R": "RICE-1",
    "P": "PEA-1",
    "PP": "PPR Type4",
}

# ── Schema ──────────────────────────────────────────────────────

CREATE_TABLES_SQL = """
CREATE TABLE IF NOT EXISTS experiments (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    experiment_date TEXT,
    strain_name TEXT,
    media_type TEXT DEFAULT 'peptone_screening',
    goal TEXT,
    source_filename TEXT,
    processed_at TEXT,
    notes TEXT,
    base_medium_preset TEXT,
    base_medium_custom_name TEXT,
    base_medium_composition_json TEXT,
    composition_groups_json TEXT
);

CREATE TABLE IF NOT EXISTS growth_curves (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    experiment_id INTEGER NOT NULL,
    group_code TEXT,
    peptone_name TEXT,
    peptone_pct REAL,
    peptone_1 TEXT,
    ratio_1 REAL,
    peptone_2 TEXT,
    ratio_2 REAL,
    strain_name TEXT,
    time_hours_json TEXT,
    mean_od_json TEXT,
    sd_od_json TEXT,
    n_replicates INTEGER,
    condition_name TEXT,
    variation_desc TEXT,
    variation_overrides_json TEXT,
    composition_json TEXT,
    composition_group_id TEXT,
    FOREIGN KEY (experiment_id) REFERENCES experiments(id)
);

CREATE TABLE IF NOT EXISTS growth_metrics (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    growth_curve_id INTEGER UNIQUE NOT NULL,
    max_od REAL,
    final_od REAL,
    lag_time_h REAL,
    mu_max REAL,
    doubling_time_h REAL,
    auc REAL,
    t_max_od_h REAL,
    computed_at TEXT,
    FOREIGN KEY (growth_curve_id) REFERENCES growth_curves(id)
);

CREATE TABLE IF NOT EXISTS strain_aliases (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    alias TEXT UNIQUE NOT NULL,
    genus TEXT,
    species TEXT,
    full_name TEXT
);

CREATE TABLE IF NOT EXISTS peptone_aliases (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    alias TEXT UNIQUE NOT NULL,
    canonical_name TEXT NOT NULL
);
"""


def _resolve_strain(raw: str) -> str:
    """Resolve strain alias to full name. Returns original if no match."""
    key = raw.strip().upper()
    if key in STRAIN_ALIASES:
        info = STRAIN_ALIASES[key]
        return f"{info['genus']} {info['species']}"
    return raw.strip()


def _resolve_peptone(raw: str) -> str:
    """Resolve peptone alias to canonical name. Returns original if no match."""
    key = raw.strip().upper()
    if key in PEPTONE_ALIASES:
        return PEPTONE_ALIASES[key]
    # Also try original case
    if raw.strip() in PEPTONE_ALIASES:
        return PEPTONE_ALIASES[raw.strip()]
    return raw.strip()


class GrowthDB:
    """Growth curve database with alias resolution and metrics storage."""

    def __init__(self, db_path: Optional[Path] = None):
        self.db_path = db_path or Path("data/growth_data.db")
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self.conn = sqlite3.connect(str(self.db_path), check_same_thread=False)
        self.conn.row_factory = sqlite3.Row
        self._init_tables()

    def _init_tables(self):
        """Create tables and seed alias data."""
        self.conn.executescript(CREATE_TABLES_SQL)
        self.conn.commit()
        self._migrate_schema()
        self._seed_aliases()

    def _migrate_schema(self):
        """Add new columns to existing tables if missing (idempotent)."""
        migrations = [
            ("experiments", "base_medium_preset", "TEXT"),
            ("experiments", "base_medium_custom_name", "TEXT"),
            ("experiments", "base_medium_composition_json", "TEXT"),
            ("experiments", "composition_groups_json", "TEXT"),
            ("growth_curves", "condition_name", "TEXT"),
            ("growth_curves", "variation_desc", "TEXT"),
            ("growth_curves", "variation_overrides_json", "TEXT"),
            ("growth_curves", "composition_json", "TEXT"),
            ("growth_curves", "composition_group_id", "TEXT"),
        ]
        for table, column, coltype in migrations:
            cur = self.conn.execute(f"PRAGMA table_info({table})")
            existing = [row["name"] for row in cur.fetchall()]
            if column not in existing:
                try:
                    self.conn.execute(
                        f"ALTER TABLE {table} ADD COLUMN {column} {coltype}"
                    )
                    logger.info(f"Migrated: added {table}.{column}")
                except sqlite3.OperationalError as e:
                    logger.warning(f"Migration skipped for {table}.{column}: {e}")
        self.conn.commit()

    def _seed_aliases(self):
        """Insert default aliases if not already present."""
        for alias, info in STRAIN_ALIASES.items():
            full_name = f"{info['genus']} {info['species']}"
            try:
                self.conn.execute(
                    "INSERT OR IGNORE INTO strain_aliases (alias, genus, species, full_name) VALUES (?, ?, ?, ?)",
                    (alias, info["genus"], info["species"], full_name),
                )
            except sqlite3.IntegrityError:
                pass

        for alias, canonical in PEPTONE_ALIASES.items():
            try:
                self.conn.execute(
                    "INSERT OR IGNORE INTO peptone_aliases (alias, canonical_name) VALUES (?, ?)",
                    (alias, canonical),
                )
            except sqlite3.IntegrityError:
                pass

        self.conn.commit()

    def close(self):
        self.conn.close()

    # ── Ingest from growth-curve-app ────────────────────────────

    def ingest(self, payload: dict) -> dict:
        """Ingest growth data from growth-curve-app POST payload.

        Supports three flavors:
        1) Peptone screening (legacy):
           {metadata, sample_map, chart_data, source_filename}
        2) Media optimization (legacy, per-SM variation overrides):
           {metadata, experiment_type:"media_optimization",
            base_medium: {preset, composition},
            variations: [{code, strain, description, overrides}],
            chart_data, source_filename}
        3) Media optimization (v2, composition_groups):
           {metadata, experiment_type:"media_optimization",
            base_medium: {preset, custom_name, composition},
            composition_groups: [{id, name, strain, description,
                                  composition:[{name,value,unit,category}],
                                  applied_samples:["SM1","SM2",...]}],
            variations: [...],  # server-side expanded from composition_groups
            chart_data, source_filename}

        Returns dict with experiment_id and curve_count.
        """
        metadata = payload.get("metadata", {})
        sample_map = payload.get("sample_map", []) or []
        chart_data = payload.get("chart_data", {})
        source_filename = payload.get("source_filename", "unknown")
        experiment_type = payload.get("experiment_type") or metadata.get("media_type") or "peptone_screening"
        base_medium = payload.get("base_medium") or {}
        variations = payload.get("variations") or []
        composition_groups = payload.get("composition_groups") or []

        # Resolve strain alias
        raw_strain = metadata.get("strain", "")
        resolved_strain = _resolve_strain(raw_strain)

        # Insert experiment (with optional base_medium fields)
        now = datetime.now().isoformat()
        base_preset = base_medium.get("preset") if base_medium else None
        base_custom_name = (base_medium.get("custom_name") or "").strip() if base_medium else ""
        base_comp_json = json.dumps(base_medium.get("composition", [])) if base_medium else None
        comp_groups_json = json.dumps(composition_groups) if composition_groups else None
        cur = self.conn.execute(
            """INSERT INTO experiments
               (experiment_date, strain_name, media_type, goal, source_filename, processed_at,
                base_medium_preset, base_medium_custom_name, base_medium_composition_json,
                composition_groups_json)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            (
                metadata.get("experiment_date", ""),
                resolved_strain,
                experiment_type,
                metadata.get("goal", ""),
                source_filename,
                now,
                base_preset,
                base_custom_name or None,
                base_comp_json,
                comp_groups_json,
            ),
        )
        experiment_id = cur.lastrowid

        # Build lookup by group_code for BOTH sample_map (legacy) and variations (media opt)
        sample_lookup = {}
        for entry in sample_map:
            code = entry.get("code", "")
            if code:
                sample_lookup[code] = entry

        variation_lookup = {}
        for var in variations:
            code = var.get("code", "")
            if code:
                variation_lookup[code] = var

        # Insert growth curves from chart_data series
        time_hours = chart_data.get("time_hours", [])
        series_list = chart_data.get("series", [])
        curve_count = 0

        for series in series_list:
            group_name = series.get("name", "")
            mean_values = series.get("mean", [])
            sd_values = series.get("sd", [])

            # Prefer explicit group_code from growth-curve-app's extract_chart_data().
            # The `name` field is a display label like "SOY-1 (SM1)" and can't be
            # reliably reverse-parsed. Fall back to parsing only if group_code
            # is absent (older payloads).
            group_code = series.get("group_code") or ""
            if not group_code and group_name:
                # Legacy fallback: best-effort parse from name
                group_code = group_name.split(" ")[0].split("-")[0]

            # Look up sample info + variation info
            sample_info = sample_lookup.get(group_code, {})
            variation_info = variation_lookup.get(group_code, {})

            # Resolve peptone aliases (peptone screening path)
            raw_peptone1 = sample_info.get("peptone_1", "") or sample_info.get("name", "")
            raw_peptone2 = sample_info.get("peptone_2", "")
            peptone_1 = _resolve_peptone(raw_peptone1) if raw_peptone1 else ""
            peptone_2 = _resolve_peptone(raw_peptone2) if raw_peptone2 else ""

            # Strain per curve (variation strain > sample_map strain > experiment strain)
            curve_strain = variation_info.get("strain", "") or sample_info.get("strain", "")
            if curve_strain:
                curve_strain = _resolve_strain(curve_strain)
            else:
                curve_strain = resolved_strain

            # Determine display peptone name
            ratio_1 = sample_info.get("ratio_1", 100)
            ratio_2 = sample_info.get("ratio_2", 0)
            if peptone_2 and ratio_2:
                peptone_name = f"{peptone_1}_{ratio_1}%+{peptone_2}_{ratio_2}%"
            else:
                peptone_name = peptone_1

            # Variation fields (media_optimization path)
            condition_name = (variation_info.get("condition_name") or "").strip() or None
            variation_desc = variation_info.get("description", "") or ""
            variation_overrides = variation_info.get("overrides", {}) or {}
            variation_overrides_json = json.dumps(variation_overrides) if variation_overrides else None

            # Full per-SM composition (v2 composition_groups flow).
            # Prefer variation.composition; fall back to chart_data series.composition
            # (growth-curve-app now enriches each series with its composition).
            composition_list = (
                variation_info.get("composition")
                or series.get("composition")
                or []
            )
            composition_json = json.dumps(composition_list) if composition_list else None
            composition_group_id = (
                variation_info.get("group_id")
                or series.get("group_id")
                or None
            )

            self.conn.execute(
                """INSERT INTO growth_curves
                   (experiment_id, group_code, peptone_name, peptone_pct,
                    peptone_1, ratio_1, peptone_2, ratio_2, strain_name,
                    time_hours_json, mean_od_json, sd_od_json, n_replicates,
                    condition_name, variation_desc, variation_overrides_json,
                    composition_json, composition_group_id)
                   VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                (
                    experiment_id,
                    group_code,
                    peptone_name,
                    sample_info.get("peptone_pct", 0),
                    peptone_1,
                    ratio_1,
                    peptone_2,
                    ratio_2,
                    curve_strain,
                    json.dumps(time_hours),
                    json.dumps(mean_values),
                    json.dumps(sd_values),
                    sample_info.get("n_replicates", 3),
                    condition_name,
                    variation_desc,
                    variation_overrides_json,
                    composition_json,
                    composition_group_id,
                ),
            )
            curve_count += 1

        self.conn.commit()
        logger.info(
            f"Ingested experiment {experiment_id}: {curve_count} curves from {source_filename}"
        )

        # Auto-compute metrics
        self._compute_metrics_for_experiment(experiment_id)

        return {"experiment_id": experiment_id, "curve_count": curve_count}

    # ── Metrics computation ─────────────────────────────────────

    def _compute_metrics_for_experiment(self, experiment_id: int):
        """Compute growth metrics for all curves in an experiment."""
        cur = self.conn.execute(
            "SELECT id, time_hours_json, mean_od_json FROM growth_curves WHERE experiment_id = ?",
            (experiment_id,),
        )
        for row in cur.fetchall():
            curve_id = row["id"]
            try:
                time_hours = json.loads(row["time_hours_json"])
                mean_od = json.loads(row["mean_od_json"])
                if not time_hours or not mean_od:
                    continue
                metrics = self._calc_metrics(time_hours, mean_od)
                self.conn.execute(
                    """INSERT OR REPLACE INTO growth_metrics
                       (growth_curve_id, max_od, final_od, lag_time_h, mu_max,
                        doubling_time_h, auc, t_max_od_h, computed_at)
                       VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                    (
                        curve_id,
                        metrics["max_od"],
                        metrics["final_od"],
                        metrics["lag_time_h"],
                        metrics["mu_max"],
                        metrics["doubling_time_h"],
                        metrics["auc"],
                        metrics["t_max_od_h"],
                        datetime.now().isoformat(),
                    ),
                )
            except Exception as e:
                logger.warning(f"Metrics computation failed for curve {curve_id}: {e}")

        self.conn.commit()

    @staticmethod
    def _calc_metrics(time_hours: list, mean_od: list) -> dict:
        """Calculate growth metrics from time-series data."""
        import math

        max_od = max(mean_od)
        final_od = mean_od[-1]
        t_max_od_h = time_hours[mean_od.index(max_od)]

        # AUC (trapezoidal rule)
        auc = 0.0
        for i in range(1, len(time_hours)):
            dt = time_hours[i] - time_hours[i - 1]
            auc += (mean_od[i - 1] + mean_od[i]) / 2 * dt

        # mu_max: max specific growth rate from ln(OD) slope
        mu_max = 0.0
        best_lag = 0.0
        for i in range(1, len(mean_od)):
            if mean_od[i] > 0 and mean_od[i - 1] > 0:
                dt = time_hours[i] - time_hours[i - 1]
                if dt > 0:
                    mu = (math.log(mean_od[i]) - math.log(mean_od[i - 1])) / dt
                    if mu > mu_max:
                        mu_max = mu

        # Doubling time
        doubling_time_h = math.log(2) / mu_max if mu_max > 0 else None

        # Lag time: first time OD exceeds initial OD + 2*SD (simplified)
        od_threshold = mean_od[0] * 1.2 if mean_od[0] > 0 else 0.05
        lag_time_h = 0.0
        for i, od in enumerate(mean_od):
            if od > od_threshold:
                lag_time_h = time_hours[i]
                break

        return {
            "max_od": round(max_od, 4),
            "final_od": round(final_od, 4),
            "lag_time_h": round(lag_time_h, 2),
            "mu_max": round(mu_max, 4),
            "doubling_time_h": round(doubling_time_h, 2) if doubling_time_h else None,
            "auc": round(auc, 4),
            "t_max_od_h": round(t_max_od_h, 2),
        }

    # ── Query methods ───────────────────────────────────────────

    def get_experiments(self, media_type: Optional[str] = None) -> list[dict]:
        """List experiments, optionally filtered by media_type."""
        if media_type:
            cur = self.conn.execute(
                "SELECT * FROM experiments WHERE media_type = ? ORDER BY id DESC",
                (media_type,),
            )
        else:
            cur = self.conn.execute("SELECT * FROM experiments ORDER BY id DESC")
        return [dict(r) for r in cur.fetchall()]

    def get_curves_with_metrics(self, experiment_id: Optional[int] = None) -> list[dict]:
        """Get growth curves joined with metrics."""
        query = """
            SELECT gc.*, gm.max_od, gm.final_od, gm.lag_time_h, gm.mu_max,
                   gm.doubling_time_h, gm.auc, gm.t_max_od_h
            FROM growth_curves gc
            LEFT JOIN growth_metrics gm ON gm.growth_curve_id = gc.id
        """
        if experiment_id:
            query += " WHERE gc.experiment_id = ?"
            cur = self.conn.execute(query, (experiment_id,))
        else:
            cur = self.conn.execute(query + " ORDER BY gc.id DESC")
        return [dict(r) for r in cur.fetchall()]

    def get_ml_training_data(self) -> list[dict]:
        """Get peptone screening data formatted for ML training."""
        cur = self.conn.execute(
            """SELECT gc.strain_name, gc.peptone_name, gc.peptone_1, gc.ratio_1,
                      gc.peptone_2, gc.ratio_2, gc.peptone_pct,
                      gm.max_od, gm.final_od, gm.mu_max, gm.auc, gm.lag_time_h
               FROM growth_curves gc
               JOIN growth_metrics gm ON gm.growth_curve_id = gc.id
               JOIN experiments e ON e.id = gc.experiment_id
               WHERE e.media_type = 'peptone_screening'
               ORDER BY gc.id"""
        )
        return [dict(r) for r in cur.fetchall()]

    def get_fba_validation_data(self) -> list[dict]:
        """Get media optimization data for FBA validation."""
        cur = self.conn.execute(
            """SELECT gc.strain_name, gc.peptone_name, gc.group_code,
                      gm.max_od, gm.final_od, gm.mu_max, gm.auc,
                      e.experiment_date, e.notes
               FROM growth_curves gc
               JOIN growth_metrics gm ON gm.growth_curve_id = gc.id
               JOIN experiments e ON e.id = gc.experiment_id
               WHERE e.media_type = 'media_optimization'
               ORDER BY gc.id"""
        )
        return [dict(r) for r in cur.fetchall()]

    def count_experiments(self) -> int:
        cur = self.conn.execute("SELECT COUNT(*) FROM experiments")
        return cur.fetchone()[0]

    def count_curves(self) -> int:
        cur = self.conn.execute("SELECT COUNT(*) FROM growth_curves")
        return cur.fetchone()[0]

    # ── Deletion ────────────────────────────────────────────────

    def delete_experiment(self, experiment_id: int) -> dict:
        """Delete one experiment and all its curves/metrics.

        Returns counts of deleted rows: {"experiments", "curves", "metrics"}.
        """
        # Collect curve ids first so we can delete their metrics.
        cur = self.conn.execute(
            "SELECT id FROM growth_curves WHERE experiment_id = ?",
            (experiment_id,),
        )
        curve_ids = [r[0] for r in cur.fetchall()]

        n_metrics = 0
        if curve_ids:
            placeholders = ",".join("?" * len(curve_ids))
            n_metrics = self.conn.execute(
                f"DELETE FROM growth_metrics WHERE growth_curve_id IN ({placeholders})",
                curve_ids,
            ).rowcount

        n_curves = self.conn.execute(
            "DELETE FROM growth_curves WHERE experiment_id = ?",
            (experiment_id,),
        ).rowcount

        n_exp = self.conn.execute(
            "DELETE FROM experiments WHERE id = ?",
            (experiment_id,),
        ).rowcount

        self.conn.commit()
        logger.info(
            f"Deleted experiment {experiment_id}: {n_exp} exp, {n_curves} curves, {n_metrics} metrics"
        )
        return {
            "experiments": int(n_exp),
            "curves":      int(n_curves),
            "metrics":     int(n_metrics),
        }

    def reset_all(self) -> dict:
        """Truncate all growth data tables. Alias tables are preserved."""
        n_metrics = self.conn.execute("DELETE FROM growth_metrics").rowcount
        n_curves  = self.conn.execute("DELETE FROM growth_curves").rowcount
        n_exp     = self.conn.execute("DELETE FROM experiments").rowcount
        # Reset autoincrement counters
        self.conn.execute(
            "DELETE FROM sqlite_sequence WHERE name IN ('experiments','growth_curves','growth_metrics')"
        )
        self.conn.commit()
        logger.warning(
            f"RESET ALL growth data: {n_exp} exp, {n_curves} curves, {n_metrics} metrics"
        )
        return {
            "experiments": int(n_exp),
            "curves":      int(n_curves),
            "metrics":     int(n_metrics),
        }

    def get_summary(self) -> dict:
        """DB summary stats."""
        return {
            "total_experiments": self.count_experiments(),
            "total_curves": self.count_curves(),
            "peptone_screening": len(self.get_experiments("peptone_screening")),
            "media_optimization": len(self.get_experiments("media_optimization")),
        }

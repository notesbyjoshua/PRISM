import sqlite3
from contextlib import contextmanager
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATABASE_PATH = PROJECT_ROOT / "backend" / "data" / "prism.db"


def initialize_database() -> None:
    DATABASE_PATH.parent.mkdir(parents=True, exist_ok=True)
    with connect() as connection:
        connection.executescript(
            """
            CREATE TABLE IF NOT EXISTS analysis_results (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                analysis_type TEXT NOT NULL,
                comparison TEXT NOT NULL,
                disease TEXT,
                feature TEXT NOT NULL,
                n_disease INTEGER,
                n_healthy INTEGER,
                hedges_g REAL,
                hedges_g_ci_low REAL,
                hedges_g_ci_high REAL,
                cliffs_delta REAL,
                roc_auc REAL,
                q_value REAL,
                pds_score REAL,
                high_priority INTEGER NOT NULL DEFAULT 0,
                rank_stability REAL,
                robust_hedges_g REAL,
                robust_q_value REAL,
                raw_json TEXT NOT NULL,
                is_demo INTEGER NOT NULL DEFAULT 0
            );

            CREATE INDEX IF NOT EXISTS idx_results_feature
                ON analysis_results(feature);
            CREATE INDEX IF NOT EXISTS idx_results_disease
                ON analysis_results(disease);
            CREATE INDEX IF NOT EXISTS idx_results_priority
                ON analysis_results(high_priority);

            CREATE TABLE IF NOT EXISTS feature_specificity (
                feature TEXT PRIMARY KEY,
                diseases_tested INTEGER NOT NULL,
                diseases_with_effect INTEGER NOT NULL,
                specificity_score REAL,
                weighted_effect_score REAL,
                median_abs_hedges_g REAL,
                raw_json TEXT NOT NULL,
                is_demo INTEGER NOT NULL DEFAULT 0
            );

            CREATE TABLE IF NOT EXISTS app_metadata (
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL
            );
            """
        )

        # Lightweight migration for databases created before PDS was added.
        columns = {
            row[1]
            for row in connection.execute("PRAGMA table_info(analysis_results)")
        }
        if "pds_score" not in columns:
            connection.execute(
                "ALTER TABLE analysis_results ADD COLUMN pds_score REAL"
            )


@contextmanager
def connect():
    connection = sqlite3.connect(DATABASE_PATH)
    connection.row_factory = sqlite3.Row
    try:
        yield connection
        connection.commit()
    finally:
        connection.close()

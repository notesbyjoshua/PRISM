import json
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from .database import PROJECT_ROOT, connect, initialize_database

RESULTS_ROOT = PROJECT_ROOT / "stats_results"
RESULT_FILES = {
    "combined": "stats_healthy_vs_disease.csv",
    "pairwise": "stats_each_disease_vs_healthy.csv",
    "kruskal": "stats_kruskal_all_groups.csv",
}


def clean_value(value):
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return None
    if isinstance(value, np.generic):
        return value.item()
    return value


def first_value(row, *columns):
    for column in columns:
        if column in row and pd.notna(row[column]):
            return clean_value(row[column])
    return None


def normalized_result(row, analysis_type):
    disease = first_value(row, "disease")
    comparison = first_value(row, "comparison") or (
        "all_groups" if analysis_type == "kruskal" else "all_disease_vs_healthy"
    )
    q_value = first_value(row, "mannwhitney_q_fdr_global", "mannwhitney_q_fdr", "kruskal_q_fdr")
    return (
        analysis_type,
        comparison,
        disease,
        first_value(row, "feature"),
        first_value(row, "n_disease", "n_total"),
        first_value(row, "n_healthy"),
        first_value(row, "hedges_g", "epsilon_squared"),
        first_value(row, "hedges_g_ci_95_low"),
        first_value(row, "hedges_g_ci_95_high"),
        first_value(row, "cliffs_delta"),
        first_value(row, "roc_auc"),
        q_value,
        first_value(row, "phenotype_difference_score_pdf", "pds_score"),
        int(bool(first_value(row, "high_priority") or False)),
        first_value(row, "top_k_bootstrap_frequency"),
        first_value(row, "robust_hedges_g", "robust_epsilon_squared"),
        first_value(row, "robust_mannwhitney_q_fdr_global", "robust_mannwhitney_q_fdr", "robust_kruskal_q_fdr"),
        json.dumps({key: clean_value(value) for key, value in row.items()}),
        0,
    )


def find_results_directory(results_root: Path = RESULTS_ROOT) -> Path | None:
    """Locate analysis CSVs in the current or legacy output layout."""
    candidates = [results_root / "significance_testing", results_root]
    for candidate in candidates:
        if any((candidate / filename).is_file() for filename in RESULT_FILES.values()):
            return candidate
    return None


def import_csv_results(results_dir: Path | None = None) -> bool:
    initialize_database()
    results_dir = results_dir or find_results_directory()
    if results_dir is None:
        return False
    available = {
        analysis_type: results_dir / filename
        for analysis_type, filename in RESULT_FILES.items()
        if (results_dir / filename).is_file()
    }
    if not available:
        return False
    with connect() as connection:
        connection.execute("DELETE FROM analysis_results")
        connection.execute("DELETE FROM feature_specificity")
        for analysis_type, path in available.items():
            frame = pd.read_csv(path)
            rows = [normalized_result(row, analysis_type) for row in frame.to_dict("records")]
            connection.executemany(
                """
                INSERT INTO analysis_results (
                    analysis_type, comparison, disease, feature, n_disease,
                    n_healthy, hedges_g, hedges_g_ci_low, hedges_g_ci_high,
                    cliffs_delta, roc_auc, q_value, pds_score,
                    high_priority, rank_stability, robust_hedges_g, robust_q_value,
                    raw_json, is_demo
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                rows,
            )
        specificity_path = results_dir / "feature_specificity.csv"
        if specificity_path.is_file():
            frame = pd.read_csv(specificity_path)
            rows = []
            for row in frame.to_dict("records"):
                rows.append(
                    (
                        row["feature"],
                        int(row["n_diseases_tested"]),
                        int(row["n_diseases_with_significant_moderate_or_large_effect"]),
                        clean_value(row.get("feature_specificity_score")),
                        clean_value(row.get("specificity_weighted_effect_score")),
                        clean_value(row.get("median_abs_hedges_g_among_significant")),
                        json.dumps({key: clean_value(value) for key, value in row.items()}),
                        0,
                    )
                )
            connection.executemany(
                """
                INSERT INTO feature_specificity VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """,
                rows,
            )
        connection.execute("INSERT OR REPLACE INTO app_metadata(key, value) VALUES ('data_source', 'analysis')")
        connection.execute(
            "INSERT OR REPLACE INTO app_metadata(key, value) VALUES ('imported_at', ?)",
            (datetime.now(timezone.utc).isoformat(),),
        )
    return True


def seed_demo_data() -> None:
    initialize_database()
    demo = [
        (
            "pairwise",
            "Williams syndrome vs healthy",
            "Williams syndrome",
            "chin_height",
            82,
            410,
            1.18,
            0.88,
            1.49,
            0.54,
            0.84,
            2.1e-8,
            91.4,
            1,
            0.94,
            1.12,
            4.2e-8,
        ),
        (
            "pairwise",
            "Noonan syndrome vs healthy",
            "Noonan syndrome",
            "inter_canthal_distance",
            74,
            410,
            0.81,
            0.55,
            1.07,
            0.39,
            0.76,
            3.4e-5,
            74.8,
            1,
            0.78,
            0.79,
            5.1e-5,
        ),
        (
            "pairwise",
            "Kabuki syndrome vs healthy",
            "Kabuki syndrome",
            "eye_fissure_slant_mean",
            68,
            410,
            -0.96,
            -1.24,
            -0.67,
            -0.45,
            0.79,
            7.7e-7,
            82.6,
            1,
            0.87,
            -0.91,
            1.2e-6,
        ),
        (
            "pairwise",
            "Cornelia de Lange syndrome vs healthy",
            "Cornelia de Lange syndrome",
            "philtrum_length",
            91,
            410,
            -1.31,
            -1.58,
            -1.04,
            -0.59,
            0.86,
            8.5e-11,
            95.1,
            1,
            0.97,
            -1.26,
            2.4e-10,
        ),
        (
            "pairwise",
            "Noonan syndrome vs healthy",
            "Noonan syndrome",
            "face_aspect_ratio",
            74,
            410,
            0.43,
            0.18,
            0.68,
            0.22,
            0.64,
            0.031,
            46.2,
            0,
            0.31,
            0.41,
            0.044,
        ),
        (
            "pairwise",
            "Williams syndrome vs healthy",
            "Williams syndrome",
            "mouth_width",
            82,
            410,
            0.69,
            0.42,
            0.96,
            0.34,
            0.71,
            0.0008,
            67.9,
            1,
            0.62,
            0.66,
            0.0012,
        ),
    ]
    with connect() as connection:
        existing = connection.execute("SELECT COUNT(*) FROM analysis_results").fetchone()[0]
        if existing:
            return
        for values in demo:
            raw = json.dumps({"note": "Illustrative demo row"})
            connection.execute(
                """
                INSERT INTO analysis_results (
                    analysis_type, comparison, disease, feature, n_disease,
                    n_healthy, hedges_g, hedges_g_ci_low, hedges_g_ci_high,
                    cliffs_delta, roc_auc, q_value, pds_score,
                    high_priority, rank_stability, robust_hedges_g, robust_q_value, raw_json, is_demo
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 1)
                """,
                (*values, raw),
            )
        specificity = [
            ("philtrum_length", 312, 14, 0.955, 1.19, 1.25),
            ("inter_canthal_distance", 305, 11, 0.964, 0.82, 0.85),
            ("chin_height", 318, 38, 0.881, 0.97, 1.10),
            ("eye_fissure_slant_mean", 298, 27, 0.909, 0.88, 0.97),
            ("face_aspect_ratio", 320, 220, 0.313, 0.65, 0.20),
        ]
        for row in specificity:
            connection.execute(
                "INSERT INTO feature_specificity VALUES (?, ?, ?, ?, ?, ?, ?, 1)",
                (*row, json.dumps({"note": "Illustrative demo row"})),
            )
        connection.execute("INSERT OR REPLACE INTO app_metadata(key, value) VALUES ('data_source', 'demo')")


def ensure_data() -> None:
    initialize_database()
    if import_csv_results():
        return
    seed_demo_data()


if __name__ == "__main__":
    imported = import_csv_results()
    print("Imported stats_results CSVs." if imported else "No result CSVs found.")

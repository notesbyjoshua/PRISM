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
    candidates = [results_root / "significance_testing"]
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
    config=json.loads((PROJECT_ROOT/'ACTIVE_ANALYSIS.json').read_text())
    if config.get('status')!='VALIDATED_ACTIVE' or config.get('feature_count')!=120 or config.get('synthetic_rows_in_primary')!=0:
        raise ValueError('Only the validated current real-source 120-feature study can be imported')
    if 'pairwise' not in available:
        raise ValueError('A complete current primary result table is required')
    primary=pd.read_csv(available['pairwise'])
    if primary.feature.nunique()!=120 or primary.disease.nunique()!=config['eligible_syndromes'] or len(primary)!=config['primary_hypotheses']:
        raise ValueError('The import does not match the current primary schema or family')
    if int(primary.mannwhitney_q_fdr_global.lt(.05).sum())!=config['primary_significant']:
        raise ValueError('The import does not match the validated primary results')
    features=set(primary.feature)
    if any(set(pd.read_csv(path,usecols=['feature']).feature)!=features for path in available.values()):
        raise ValueError('Every imported analysis must use the same current 120 features')
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




def ensure_data() -> None:
    initialize_database()
    if import_csv_results():
        return
    raise FileNotFoundError("Validated current 120-feature results are unavailable; unvalidated fallback data will not be substituted")


if __name__ == "__main__":
    imported = import_csv_results()
    print("Imported stats_results CSVs." if imported else "No result CSVs found.")

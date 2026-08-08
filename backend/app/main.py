import json
import math
from contextlib import asynccontextmanager
from functools import lru_cache
from pathlib import Path

import pandas as pd
from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse

from .database import connect
from .import_results import ensure_data, import_csv_results


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DISEASE_PHENOTYPES = PROJECT_ROOT / "data" / "phenotypes_all_disease.csv"
HEALTHY_PHENOTYPES = PROJECT_ROOT / "data" / "phenotypes_all_healthy.csv"
DISTRIBUTION_PLOTS = PROJECT_ROOT / "stats_results" / "data_distribution"


def serialize_row(row):
    result = dict(row)
    for key, value in result.items():
        if isinstance(value, float) and not math.isfinite(value):
            result[key] = None
    if "high_priority" in result:
        result["high_priority"] = bool(result["high_priority"])
    if "is_demo" in result:
        result["is_demo"] = bool(result["is_demo"])
    return result


@lru_cache(maxsize=256)
def distribution_values(path_string: str, feature: str, disease: str | None) -> list[float]:
    """Read one measurement column and cache it for interactive distribution plots."""
    path = Path(path_string)
    header = pd.read_csv(path, nrows=0).columns
    if feature not in header:
        raise ValueError("Unknown feature")
    columns = [feature, "frontal_ok"]
    if disease is not None:
        columns.append("disease")
    frame = pd.read_csv(path, usecols=columns)
    frame = frame[frame["frontal_ok"] == True]  # noqa: E712
    if disease is not None:
        frame = frame[frame["disease"] == disease]
    values = pd.to_numeric(frame[feature], errors="coerce").dropna()
    return [float(value) for value in values if math.isfinite(float(value))]


@asynccontextmanager
async def lifespan(_app: FastAPI):
    ensure_data()
    yield


app = FastAPI(
    title="PRISM API",
    description="Rare-disease facial phenotype statistics API",
    version="0.1.0",
    lifespan=lifespan,
)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000", "http://127.0.0.1:3000"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/api/health")
def health():
    return {"status": "ok"}


@app.get("/api/summary")
def summary():
    with connect() as connection:
        totals = connection.execute(
            """
            SELECT COUNT(*) result_count,
                   COUNT(DISTINCT feature) feature_count,
                   COUNT(DISTINCT disease) disease_count,
                   SUM(high_priority) high_priority_count,
                   AVG(roc_auc) average_auc,
                   MAX(is_demo) is_demo
            FROM analysis_results
            """
        ).fetchone()
        source = connection.execute(
            "SELECT value FROM app_metadata WHERE key = 'data_source'"
        ).fetchone()
    return {**serialize_row(totals), "data_source": source[0] if source else "unknown"}


@app.get("/api/diseases")
def diseases():
    with connect() as connection:
        rows = connection.execute(
            "SELECT DISTINCT disease FROM analysis_results WHERE disease IS NOT NULL ORDER BY disease"
        ).fetchall()
    return [row[0] for row in rows]


@app.get("/api/results")
def results(
    disease: str | None = None,
    search: str | None = None,
    high_priority: bool | None = None,
    analysis_type: str = "pairwise",
    sort_by: str = Query("pds_score", pattern="^(q_value|pds_score|roc_auc|hedges_g|rank_stability)$"),
    descending: bool = True,
    limit: int = Query(50, ge=1, le=500),
    offset: int = Query(0, ge=0),
):
    clauses = ["analysis_type = ?"]
    parameters: list[object] = [analysis_type]
    if disease:
        clauses.append("disease = ?")
        parameters.append(disease)
    if search:
        clauses.append("feature LIKE ?")
        parameters.append(f"%{search}%")
    if high_priority is not None:
        clauses.append("high_priority = ?")
        parameters.append(int(high_priority))
    where = " AND ".join(clauses)
    direction = "DESC" if descending else "ASC"
    null_order = "roc_auc IS NULL" if sort_by == "roc_auc" else f"{sort_by} IS NULL"

    with connect() as connection:
        total = connection.execute(
            f"SELECT COUNT(*) FROM analysis_results WHERE {where}", parameters
        ).fetchone()[0]
        rows = connection.execute(
            f"""
            SELECT id, analysis_type, comparison, disease, feature, n_disease,
                   n_healthy, hedges_g, hedges_g_ci_low, hedges_g_ci_high,
                   cliffs_delta, roc_auc, q_value, pds_score, high_priority,
                   rank_stability, robust_hedges_g, robust_q_value, is_demo
            FROM analysis_results
            WHERE {where}
            ORDER BY {null_order}, {sort_by} {direction}
            LIMIT ? OFFSET ?
            """,
            [*parameters, limit, offset],
        ).fetchall()
    return {"items": [serialize_row(row) for row in rows], "total": total}


@app.get("/api/specificity")
def specificity(limit: int = Query(20, ge=1, le=200)):
    with connect() as connection:
        rows = connection.execute(
            """
            SELECT feature, diseases_tested, diseases_with_effect,
                   specificity_score, weighted_effect_score,
                   median_abs_hedges_g, is_demo
            FROM feature_specificity
            ORDER BY weighted_effect_score DESC
            LIMIT ?
            """,
            (limit,),
        ).fetchall()
    return [serialize_row(row) for row in rows]


@app.get("/api/features/{feature}")
def feature_detail(feature: str):
    with connect() as connection:
        rows = connection.execute(
            """
            SELECT * FROM analysis_results
            WHERE feature = ?
            ORDER BY pds_score IS NULL, pds_score DESC
            """,
            (feature,),
        ).fetchall()
        specificity_row = connection.execute(
            "SELECT * FROM feature_specificity WHERE feature = ?", (feature,)
        ).fetchone()
    if not rows:
        raise HTTPException(status_code=404, detail="Feature not found")
    return {
        "feature": feature,
        "results": [serialize_row(row) for row in rows],
        "specificity": serialize_row(specificity_row) if specificity_row else None,
    }


@app.get("/api/distributions")
def distributions(disease: str, feature: str):
    """Return the underlying healthy and disease measurements for a selected result."""
    if not DISEASE_PHENOTYPES.is_file() or not HEALTHY_PHENOTYPES.is_file():
        raise HTTPException(status_code=404, detail="Phenotype source CSVs were not found")
    try:
        healthy = distribution_values(str(HEALTHY_PHENOTYPES), feature, None)
        affected = distribution_values(str(DISEASE_PHENOTYPES), feature, disease)
    except ValueError as error:
        raise HTTPException(status_code=404, detail=str(error)) from error
    if not affected:
        raise HTTPException(status_code=404, detail="No matching disease measurements were found")
    return {"feature": feature, "disease": disease, "healthy": healthy, "disease_values": affected}


@app.get("/api/distribution-plots")
def distribution_plots():
    """List dataset-level distribution charts generated by data_distribution.py."""
    plot_details = {
        "01_rows_per_disease.png": ("Rows per disease", "Cohort size for every disease, ordered from largest to smallest."),
        "02_top_20_diseases.png": ("Top 20 diseases", "The twenty disease cohorts with the most images."),
        "03_rows_per_ethnicity.png": ("Rows per ethnicity", "Image counts across the dataset's ethnicity categories."),
    }
    plots = []
    if DISTRIBUTION_PLOTS.is_dir():
        for path in sorted(DISTRIBUTION_PLOTS.glob("*.png")):
            if path.name.startswith("04_scatter_"):
                pair = path.stem.removeprefix("04_scatter_").replace("_vs_", " vs. ").replace("_", " ").title()
                title, description = pair, "Disease-level image counts compared between two ethnicity groups."
                category = "Ethnicity comparisons"
            elif path.name in plot_details:
                title, description = plot_details[path.name]
                category = "Dataset overview"
            else:
                continue
            plots.append({"filename": path.name, "title": title, "description": description, "category": category})
    return plots


@app.get("/api/distribution-plots/{filename}")
def distribution_plot_image(filename: str):
    """Serve a generated plot without exposing arbitrary filesystem paths."""
    if Path(filename).name != filename or not filename.endswith(".png"):
        raise HTTPException(status_code=404, detail="Plot not found")
    path = DISTRIBUTION_PLOTS / filename
    if not path.is_file():
        raise HTTPException(status_code=404, detail="Plot not found")
    return FileResponse(path, media_type="image/png")


@app.post("/api/admin/reload")
def reload_results():
    if not import_csv_results():
        raise HTTPException(status_code=404, detail="No statistics CSV files found")
    return {"status": "reloaded"}

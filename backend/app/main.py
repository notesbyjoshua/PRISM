import json
import math
from contextlib import asynccontextmanager

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware

from .database import connect
from .import_results import ensure_data, import_csv_results


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
    sort_by: str = Query("q_value", pattern="^(q_value|roc_auc|hedges_g|rank_stability)$"),
    descending: bool = False,
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
                   cliffs_delta, roc_auc, q_value, high_priority,
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
            "SELECT * FROM analysis_results WHERE feature = ? ORDER BY q_value",
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


@app.post("/api/admin/reload")
def reload_results():
    if not import_csv_results():
        raise HTTPException(status_code=404, detail="No statistics CSV files found")
    return {"status": "reloaded"}


import json
import shutil
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]


def export(kind="all"):
    config=json.loads((ROOT/"ACTIVE_ANALYSIS.json").read_text())
    if config["status"]!="VALIDATED_ACTIVE" or config["feature_count"]!=120:
        raise ValueError("Validated current 120-feature analysis required")
    run=ROOT/config["private_run"]
    source=run/"ACTIVE_PUBLICATION"
    target=ROOT/"publication"
    if kind in {"figures","all"}:
        shutil.copytree(source/"figures",target/"figures",dirs_exist_ok=True)
    if kind in {"supplementary","all"}:
        shutil.copytree(source/"supplementary",target/"supplementary",dirs_exist_ok=True)
    if kind in {"schema","all"}:
        out=target/"feature_definitions";out.mkdir(exist_ok=True)
        shutil.copyfile(source/"supplementary/Supplementary_Table_S6_exact_feature_definitions.csv",out/"current_120_feature_formulas.csv")
    if kind in {"correlation","all"}:
        out=target/"feature_correlation";out.mkdir(exist_ok=True)
        for name in ["feature_correlation_matrix.csv","feature_clusters.csv","redundancy_thresholds.csv"]:
            shutil.copyfile(run/"RESULTS_120/extensions"/name,out/name)
    if kind in {"clinical","all"}:
        shutil.copyfile(source/"supplementary/clinical_highlights_120.csv",target/"current_120_clinical_examples.csv")
    if kind in {"audit","provenance","all"}:
        table=pd.read_csv(ROOT/config["results_directory"]/"stats_each_disease_vs_healthy.csv")
        assert table.feature.nunique()==120 and len(table)==config["primary_hypotheses"]
        assert int(table.mannwhitney_q_fdr_global.lt(.05).sum())==config["primary_significant"]
        audit={"status":"PASS","feature_count":120,"primary_hypotheses":len(table),"primary_significant":config["primary_significant"],"source_run":config["private_run"]}
        (target/"CURRENT_120_PUBLICATION_AUDIT.json").write_text(json.dumps(audit,indent=2)+"\n")
    print("Current 120-feature publication export verified:",kind)

import json
import shutil
import tempfile
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]


def publication_table(frame, concordance=False):
    frame=frame.drop(columns=['historical_mapping','definition_status'],errors='ignore').copy()
    if concordance:
        frame=frame[[c for c in frame if not any(term in c for term in ['bootstrap','_ci_','_finite','bias_exceeds','interval_excludes'])]]
        frame['uncertainty_status']='excluded_biased_bootstrap_method; descriptive_point_estimates_only'
    frame=frame.rename(columns=lambda c:c.replace('healthy','reference'))
    for column in frame.select_dtypes(include=['object','string']):
        frame[column]=frame[column].replace({'eligible49':'other_eligible_syndromes','historical_fixed_anchor_pairs_120_features':'prespecified_anchor_pairs_120_features'})
    if 'comparison' in frame:
        frame['comparison']=frame.comparison.str.replace('healthy','FairFace reference',regex=False)
    return frame


def copy_tables(source,target):
    target.mkdir(parents=True,exist_ok=True)
    for path in source.glob('*.csv'):
        copy_table(path,target/path.name.replace('healthy','reference'))


def copy_table(source,target):
    original=pd.read_csv(source)
    current=publication_table(original,concordance=source.name=='concordance_summary.csv')
    with tempfile.NamedTemporaryFile(dir=target.parent,delete=False) as handle:
        temporary=Path(handle.name)
    try:
        if original.equals(current):
            shutil.copyfile(source,temporary)
        else:
            current.to_csv(temporary,index=False)
        temporary.replace(target)
    finally:
        temporary.unlink(missing_ok=True)


def export(kind="all"):
    config=json.loads((ROOT/"ACTIVE_ANALYSIS.json").read_text())
    if config["status"]!="VALIDATED_ACTIVE" or config["feature_count"]!=120:
        raise ValueError("Validated current 120-feature analysis required")
    run=ROOT/config["private_run"]
    source=run/"ACTIVE_PUBLICATION"
    target=ROOT/"publication"
    if kind in {"figures","all"}:
        shutil.copytree(source/"figures",target/"figures",dirs_exist_ok=True)
        for path in (target/"figures").glob('*.svg'):
            path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
    if kind in {"supplementary","all"}:
        copy_tables(source/"supplementary",target/"supplementary")
    if kind in {"schema","all"}:
        out=target/"feature_definitions";out.mkdir(exist_ok=True)
        copy_table(source/"supplementary/Supplementary_Table_S6_exact_feature_definitions.csv",out/"current_120_feature_formulas.csv")
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

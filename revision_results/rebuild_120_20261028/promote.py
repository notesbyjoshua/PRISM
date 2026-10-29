import hashlib
import json
import os
import shutil
import sqlite3
import subprocess
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
RUN=HERE/'local_only/run_20261028'
os.umask(0o077)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    validation=json.loads((RUN/'VALIDATION_FULL.json').read_text())
    assert validation['status']=='VALIDATED' and validation['primary_tests']==5760 and validation['synthetic_primary_rows']==0
    numbers=json.loads((RUN/'REPORT/REPORT_NUMBERS.json').read_text())
    assert numbers['primary_significant']==validation['primary_global_BH_significant']
    snapshot=RUN/'EXECUTED_CODE'
    if not snapshot.exists():
        snapshot.mkdir(mode=0o700)
        shutil.copytree(ROOT/'stats_testing',snapshot/'stats_testing',ignore=shutil.ignore_patterns('__pycache__'))
        shutil.copytree(ROOT/'revision_results/reviewer_followup_2026',snapshot/'reviewer_followup_source',
            ignore=shutil.ignore_patterns('local_only','figures','__pycache__','*.csv','*.json','*.pdf','*.png','*.svg','*.log'))
    retired=RUN/'archive/retired_active_exports';retired.mkdir(parents=True,exist_ok=True,mode=0o700)
    stage=RUN/'PROMOTION_STAGE';stage.mkdir(exist_ok=True)
    stats=stage/'stats_results';stats.mkdir(exist_ok=True)
    shutil.copytree(RUN/'RESULTS_120/primary_balanced',stats/'significance_testing',dirs_exist_ok=True,ignore=shutil.ignore_patterns('*.log'))
    shutil.copytree(RUN/'RESULTS_120/demographics',stats/'ethnicity_analysis',dirs_exist_ok=True,ignore=shutil.ignore_patterns('*.log'))
    shutil.copytree(RUN/'RESULTS_120/extensions',stats/'extensions_120',dirs_exist_ok=True,ignore=shutil.ignore_patterns('*_STATUS.json'))
    for name in ['reference_all_eligible','reference_race_sex','reference_fixed_historical','alternative_photo']:
        shutil.copytree(RUN/'RESULTS_120'/name,stats/'sensitivities'/name,dirs_exist_ok=True,ignore=shutil.ignore_patterns('*.log','JOB_STATUS.json'))
    (stats/'README.md').write_text('# Active 120-feature results\n\nValidated current analysis: 48 eligible syndromes, 5,760 primary hypotheses, 955 global-BH discoveries. Reference: 3,409 FairFace paths; disease pool: 1,845 recorded patients. All synthetic records are excluded from primary inference. These are aggregate working outputs, not an approved public release. Prior active results are preserved privately in the immutable historical archive.\n')
    publication=stage/'publication';shutil.copytree(RUN/'ACTIVE_PUBLICATION',publication,dirs_exist_ok=True)
    shutil.copytree(ROOT/'publication/historical_facekit',publication/'historical_facekit',dirs_exist_ok=True)
    scripts=publication/'scripts';scripts.mkdir(exist_ok=True)
    historic_definition=(ROOT/'publication/scripts/generate_feature_definitions.py').read_text()
    historic_definition=historic_definition.replace('def main():','def historical_main():',1)
    marker='if __name__ == "__main__": main()'
    if marker in historic_definition:
        historic_definition=historic_definition.replace(marker,'if __name__ == "__main__":\n    from current_120 import export\n    export("schema")')
    else:
        position=historic_definition.rfind('if __name__')
        historic_definition=historic_definition[:position]+'if __name__ == "__main__":\n    from current_120 import export\n    export("schema")\n'
    (scripts/'generate_feature_definitions.py').write_text(historic_definition)
    helper='''import json
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
        (target/"CURRENT_120_PUBLICATION_AUDIT.json").write_text(json.dumps(audit,indent=2)+"\\n")
    print("Current 120-feature publication export verified:",kind)
'''
    (scripts/'current_120.py').write_text(helper)
    for filename,kind in [('generate_figures.py','figures'),('generate_supplementary.py','supplementary'),
        ('generate_candidate_evidence.py','clinical'),('feature_correlation.py','correlation'),
        ('generate_provenance.py','provenance'),('verify_reproducible_outputs.py','audit')]:
        (scripts/filename).write_text(f'from current_120 import export\n\nif __name__ == "__main__":\n    export("{kind}")\n')
    for name in ['manuscript_results_audit.py','publication_lint.py']:
        (publication/name).write_text('from scripts.current_120 import export\n\nif __name__ == "__main__":\n    export("audit")\n')
    shutil.copyfile(ROOT/'publication/scripts/data_privacy_audit.py',scripts/'data_privacy_audit.py')
    manuscript=stage/'manuscript';shutil.copytree(RUN/'MANUSCRIPT_120',manuscript,dirs_exist_ok=True)
    shutil.copytree(RUN/'ACTIVE_PUBLICATION/supplementary',manuscript/'supplementary',dirs_exist_ok=True)
    (manuscript/'README.md').write_text('# Active 120-feature manuscript draft\n\nmain.tex and figures use the validated current analysis. Primary: 48 syndromes, 120 features, 5,760 tests, 955 discoveries. No real patient photographs. The previous manuscript and reviewer response are archived privately. Author review of clinical wording, IRB identifier, funding and final release terms remains required before submission. See the private run MANUSCRIPT_HANDOFF.md and CHANGELOG_120_MIGRATION.md. The live Google Doc was not edited.\n')
    (manuscript/'RESPONSE_TO_REVIEWER.md').write_text('# Current analysis revision notes\n\nThe former response refers to a retired analysis and is preserved in the read-only archive. Use the current 120-feature manuscript and source-linked handoff for any new response. Primary counts: 1,845 recorded patients, 3,409 reference paths, 48 eligible syndromes, 5,760 hypotheses and 955 global-BH discoveries. No external reviewer replies were sent.\n')
    pds=stage/'pds_results_v2';pds.mkdir(exist_ok=True)
    (pds/'README.md').write_text('# Active current score source\n\nThe historical score exports are archived. Current descriptive phenotype-difference scores are in stats_results/significance_testing/stats_each_disease_vs_healthy.csv and are derived from the new 120-feature analysis. Scores are not calibrated diagnostic probabilities.\n')
    changed=[]
    for name in ['stats_results','publication','manuscript','pds_results_v2']:
        target=ROOT/name;backup=retired/name
        if backup.exists():
            raise FileExistsError('Already promoted; refuse to overwrite the retired active copy')
        if target.exists():target.rename(backup)
        try:
            (stage/name).rename(target)
        except BaseException:
            if backup.exists() and not target.exists():backup.rename(target)
            raise
        changed.append(name)
    database=ROOT/'backend/data/prism.db'
    if database.exists():
        target=retired/'prepromotion_prism.db'
        source_connection=sqlite3.connect(f'file:{database}?mode=ro',uri=True)
        backup_connection=sqlite3.connect(target)
        source_connection.backup(backup_connection);backup_connection.close();source_connection.close();target.chmod(0o400)
    relative=str(RUN.relative_to(ROOT))
    config={'status':'VALIDATED_ACTIVE','analysis_id':'PRISM_120_20261028','feature_count':120,
        'private_run':relative,'facekit_commit':'5764bdd4af9360b55a046a16579fb5a890a8182f',
        'model_sha256':'64184e229b263107bc2b804c6625db1341ff2bb731874b0bcc2fe6544e0bc9ff',
        'disease_csv':relative+'/COHORTS/primary_phenotypes.csv','reference_csv':relative+'/COHORTS/reference_balanced_phenotypes.csv',
        'disease_master_csv':relative+'/COHORTS/primary_master.csv','reference_master_csv':relative+'/COHORTS/reference_balanced_master.csv',
        'feature_schema':relative+'/COHORTS/exact_feature_names.csv','results_directory':'stats_results/significance_testing',
        'primary_patients':1845,'reference_images':3409,'eligible_syndromes':48,'primary_hypotheses':5760,
        'primary_significant':numbers['primary_significant'],'synthetic_rows_in_primary':0,
        'historical_archive':relative+'/archive/historical_125','retired_active_exports':relative+'/archive/retired_active_exports',
        'promotion_basis':'User authorized validated new 120-feature primary analysis and removal of old analysis from active outputs',
        'publication_release':'not_performed; private/raw data and image publication rights remain restricted'}
    temporary=ROOT/'ACTIVE_ANALYSIS.json.tmp';temporary.write_text(json.dumps(config,indent=2)+'\n');os.replace(temporary,ROOT/'ACTIVE_ANALYSIS.json')
    for path in retired.rglob('*'):
        path.chmod(0o500 if path.is_dir() else 0o400)
    retired.chmod(0o500)
    record={'status':'PROMOTED','active_analysis':'PRISM_120_20261028','replaced_active_directories':changed,
        'archive':str(retired),'old_data_csvs':'unchanged; historical tracked per-row inputs flagged for separate privacy-remediation review',
        'live_Google_Doc':'not_edited','public_upload_or_commit':'not_performed'}
    (RUN/'PROMOTION_RECORD.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2),flush=True)


if __name__=='__main__':
    main()

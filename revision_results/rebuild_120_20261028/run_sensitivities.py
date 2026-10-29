import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
RUN=HERE/'local_only/run_20261028'
os.umask(0o077)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--dry-run',action='store_true')
    args=parser.parse_args()
    cohorts=RUN/'COHORTS';results=RUN/'RESULTS_120'
    jobs=[('reference_all_eligible','primary','reference_all_eligible'),
          ('reference_race_sex','primary','reference_race_sex'),
          ('reference_fixed_historical','primary','reference_fixed_historical'),
          ('alternative_photo','alternate','reference_balanced'),
          ('fixed_common_new','fixed_common_new','reference_fixed_historical'),
          ('fixed_common_old','fixed_common_old','fixed_common_reference_old')]
    if args.dry_run:
        print(json.dumps({'jobs':jobs,'bootstrap':500,'rank_bootstrap':200,'seed':2026},indent=2));return
    for name,disease,reference in jobs:
        out=results/name;out.mkdir(exist_ok=True,parents=True)
        disease_path=cohorts/f'{disease}_phenotypes.csv';reference_path=cohorts/f'{reference}_phenotypes.csv'
        executed=RUN/'EXECUTED_CODE/stats_testing'
        code=executed if executed.exists() else ROOT/'stats_testing'
        signature={'disease_sha256':sha(disease_path),'reference_sha256':sha(reference_path),
            'analysis_script_sha256':sha(code/'statistical_analysis.py'),
            'bootstrap_code_sha256':sha(code/'bootstrap_utils.py'),'bootstrap':500,'rank':200,'seed':2026}
        state=out/'JOB_STATUS.json'
        if state.exists():
            previous=json.loads(state.read_text())
            if previous['signature']!=signature:
                raise ValueError('Input or code changed; choose a new run')
            if previous['status']=='COMPLETE':
                assert all(sha(out/file)==digest for file,digest in previous['outputs'].items())
                print(name+': verified completed checkpoint',flush=True);continue
        state.write_text(json.dumps({'status':'RUNNING','signature':signature},indent=2))
        command=[sys.executable,'-u',str(code/'statistical_analysis.py'),
            '--disease-csv',str(disease_path),'--healthy-csv',str(reference_path),'--output-dir',str(out),
            '--bootstrap-iterations','500','--rank-bootstrap-iterations','200','--bootstrap-seed','2026']
        with (out/'execution.log').open('w') as log:
            completed=subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
        status={'status':'COMPLETE' if completed.returncode==0 else 'FAILED','signature':signature,'exit_code':completed.returncode}
        if completed.returncode==0:
            status['outputs']={p.name:sha(p) for p in out.glob('*.csv')}
        state.write_text(json.dumps(status,indent=2)+'\n')
        print(name+': '+status['status'],flush=True)
        if completed.returncode:
            raise RuntimeError('Statistical job failed; see private execution log')


if __name__=='__main__':
    main()

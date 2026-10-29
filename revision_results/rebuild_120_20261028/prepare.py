import argparse
import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'revision_results/rebuild_120_20261028/local_only/run_20261028'


def main():
    parser=argparse.ArgumentParser(description='Verify or privately copy validated current FaceKit cohorts')
    parser.add_argument('--run',type=Path,default=SOURCE)
    parser.add_argument('--dry-run',action='store_true')
    args=parser.parse_args()
    if args.dry_run:
        print('Current 120-feature cohorts only: primary, alternate photograph, expanded balanced/all-eligible/category-by-sex references; no image processing or inferential execution')
        return
    features=pd.read_csv(SOURCE/'COHORTS/exact_feature_names.csv').feature.tolist()
    if len(features)!=120 or len(set(features))!=120:
        raise ValueError('Validated current schema required')
    for name in ['primary','alternate','reference_balanced','reference_all_eligible','reference_race_sex']:
        frame=pd.read_csv(SOURCE/f'COHORTS/{name}_master.csv',dtype={'patient_id':'string'})
        if not np.isfinite(frame[features].to_numpy(float)).all() or not frame.frontal_ok.eq(True).all():
            raise ValueError('Current accepted cohorts must be finite and quality-passing')
        if name in ['primary','alternate'] and not frame.patient_id.is_unique:
            raise ValueError('One current image per recorded patient is required')
    target=args.run.resolve()
    if target!=SOURCE:
        if 'local_only' not in target.parts and not target.is_relative_to(ROOT/'local_review'):
            raise ValueError('Person-level inputs must remain in private ignored storage')
        if subprocess.run(['git','check-ignore','-q',str(target/'provenance.json')],cwd=ROOT).returncode:
            raise ValueError('Person-level outputs must be Git-ignored')
        target.mkdir(parents=True,exist_ok=False)
        shutil.copytree(SOURCE/'COHORTS',target/'COHORTS')
        shutil.copy2(SOURCE/'ANALYSIS_PLAN.json',target/'ANALYSIS_PLAN.json')
    print(json.dumps({'status':'PASS','feature_count':120,'primary_patients':1845,'reference_paths':3409,'synthetic_primary_rows':0}))


if __name__=='__main__':
    main()

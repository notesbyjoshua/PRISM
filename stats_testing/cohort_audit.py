import argparse
import ast
import hashlib
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MAPPING_NAMES = {'norm_label', 'map_gmdb_layer1', 'map_gmdb_layer2', 'map_fairface_layer1', 'map_fairface_layer2'}


def read_mapping_functions(path):
    tree = ast.parse(path.read_text())
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in MAPPING_NAMES]
    if {n.name for n in nodes} != MAPPING_NAMES:
        raise ValueError(f'Mapping functions missing in {path}')
    scope = {'pd': pd, 'np': np, 're': re}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), 'exec'), scope)
    return scope


def clean_id(values):
    values = values.astype('string').str.strip()
    return values.mask(values.str.lower().isin(['', 'nan', 'none', 'null', 'unknown', 'na', 'n/a']))


def quality_status(values):
    def classify(value):
        if pd.isna(value) or str(value).strip() == '':
            return 'missing'
        if isinstance(value, (bool, np.bool_)):
            return 'pass' if value else 'fail'
        if isinstance(value, str) and value.strip().lower() in {'true', 'false'}:
            return 'pass' if value.strip().lower() == 'true' else 'fail'
        return 'malformed'
    return values.map(classify)


def patient_summary(frame, patient_column='patient_id'):
    ids = clean_id(frame[patient_column]) if patient_column in frame else pd.Series(pd.NA, index=frame.index)
    counts = ids.value_counts()
    return {'rows': len(frame), 'known_patients': int(ids.nunique()), 'missing_patient_id_rows': int(ids.isna().sum()),
            'patients_with_repeats': int(counts.gt(1).sum()), 'rows_from_repeated_patients': int(counts[counts.gt(1)].sum()),
            'excess_rows_among_known_patients': int((counts - 1).sum())}


def age_diagnostics(frame):
    year = pd.to_numeric(frame['age_year'], errors='coerce')
    month = pd.to_numeric(frame['age_month'], errors='coerce')
    numeric = np.isfinite(year) & np.isfinite(month)
    valid = numeric & year.ge(0) & month.between(0, 11) & year.mod(1).eq(0) & month.mod(1).eq(0)
    status = pd.Series('missing_or_nonnumeric', index=frame.index)
    status.loc[numeric & ~valid] = 'invalid_components'
    status.loc[valid] = 'nonzero_components_unverified_semantics'
    status.loc[valid & year.eq(0) & month.eq(0)] = 'zero_zero_unresolved'
    age = (year + month / 12).where(valid & ~status.eq('zero_zero_unresolved'))
    band = pd.cut(age, [-np.inf, 3, 10, 18, np.inf], right=False, labels=['under_3', '3_to_under_10', '10_to_under_18', '18_plus'])
    return pd.DataFrame({'age_status': status, 'conditional_age_years': age,
                         'conditional_age_band': band.astype('string').fillna(status)})


def feature_counts(frame, features):
    numeric = frame[features].apply(pd.to_numeric, errors='coerce')
    finite = pd.DataFrame(np.isfinite(numeric), index=frame.index, columns=features)
    return pd.DataFrame({'feature': features, 'rows': len(frame), 'valid_n': finite.sum().to_numpy(),
                         'missing_n': numeric.isna().sum().to_numpy(),
                         'infinite_n': np.isinf(numeric).sum().to_numpy(),
                         'invalid_fraction': 1 - finite.mean().to_numpy()})




def main():
    from analysis_config import active_configuration
    parser=argparse.ArgumentParser(description='Read-only current FaceKit cohort audit')
    parser.add_argument('--output',type=Path,default=ROOT/'local_review'/datetime.now().strftime('audit_%Y%m%d_%H%M%S'))
    parser.add_argument('--dry-run',action='store_true')
    args=parser.parse_args()
    config=active_configuration()
    disease=pd.read_csv(config['disease_master_csv'],dtype={'patient_id':'string'})
    reference=pd.read_csv(config['reference_master_csv'])
    features=pd.read_csv(config['feature_schema']).feature.tolist()
    assert len(features)==120 and disease.patient_id.is_unique and reference.file.is_unique
    assert disease.frontal_ok.eq(True).all() and reference.frontal_ok.eq(True).all()
    assert np.isfinite(disease[features].to_numpy(float)).all() and np.isfinite(reference[features].to_numpy(float)).all()
    result={'status':'PASS','feature_count':120,'selected_patients':len(disease),'reference_paths':len(reference),'eligible_syndromes':int(disease.disease.value_counts().ge(10).sum())}
    print(json.dumps(result,indent=2))
    if not args.dry_run:
        args.output.mkdir(parents=True,exist_ok=False)
        (args.output/'current_120_cohort_audit.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    main()

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from cohort_audit import ROOT, age_diagnostics, clean_id, quality_status


def deduplicate_reference(frame):
    frame = frame.copy()
    frame['reference_key'] = clean_id(frame['file']).str.replace('\\', '/', regex=False)
    if frame.reference_key.isna().any():
        raise ValueError('Reference image paths are missing')
    check = [c for c in frame if c not in {'reference_key'}]
    conflicts = frame.groupby('reference_key')[check].nunique(dropna=False).gt(1).any(axis=1)
    if conflicts.any():
        raise ValueError(f'{conflicts.sum()} repeated reference paths have conflicting records')
    frame['selection_reason'] = np.where(frame.reference_key.duplicated(), 'repeated_reference_path', 'selected')
    manifest = frame[['reference_key', 'selection_reason']].copy()
    selected = frame[frame.selection_reason.eq('selected')].drop(columns='selection_reason')
    return selected, manifest


def select_patients(frame, age_scope='all', selection='pose'):
    frame = frame.copy()
    frame['patient_id'] = clean_id(frame.patient_id)
    frame['image_id'] = clean_id(frame.image_id)
    if frame.image_id.isna().any() or frame.image_id.duplicated().any():
        raise ValueError('Disease image IDs must be nonmissing and unique')
    frame = pd.concat([frame, age_diagnostics(frame)], axis=1)
    pose = frame[['pose_yaw', 'pose_pitch', 'pose_roll']].apply(pd.to_numeric, errors='coerce').abs()
    frame['pose_score'] = (pose / [15, 15, 10]).max(axis=1).where(np.isfinite(pose).all(axis=1), np.inf)
    frame['selection_reason'] = 'eligible'
    frame.loc[~quality_status(frame.frontal_ok).eq('pass'), 'selection_reason'] = 'quality_not_pass'
    frame.loc[frame.patient_id.isna(), 'selection_reason'] = 'missing_patient_id'
    if age_scope == '3-9':
        eligible_age = frame.conditional_age_years.ge(3) & frame.conditional_age_years.lt(10)
        mask = frame.selection_reason.eq('eligible') & ~eligible_age
        frame.loc[mask, 'selection_reason'] = 'outside_or_unresolved_age_3_9'
    eligible = frame[frame.selection_reason.eq('eligible')]
    for column in ['internal_syndrome_name', 'image_dataset_ethnicity_category', 'gender']:
        if eligible.groupby('patient_id')[column].nunique(dropna=False).gt(1).any():
            raise ValueError(f'Conflicting patient metadata: {column}')
    order = ['patient_id', 'pose_score', 'image_id'] if selection == 'pose' else ['patient_id', 'image_id']
    chosen = eligible.sort_values(order, kind='stable').drop_duplicates('patient_id').index
    frame.loc[frame.selection_reason.eq('eligible'), 'selection_reason'] = 'additional_patient_image'
    frame.loc[chosen, 'selection_reason'] = 'selected'
    manifest = frame[['image_id', 'patient_id', 'internal_syndrome_name', 'age_status',
                      'conditional_age_years', 'pose_score', 'selection_reason']].copy()
    selected = frame.loc[chosen].sort_values('patient_id').drop(columns='selection_reason')
    return selected, manifest




def prepare(output):
    import shutil
    import subprocess
    from analysis_config import active_configuration
    config = active_configuration()
    output = output.resolve()
    if subprocess.run(['git', 'check-ignore', '-q', str(output / 'provenance.json')], cwd=ROOT).returncode:
        raise ValueError('Current person-level cohort exports must stay in a Git-ignored private directory')
    output.mkdir(parents=True, exist_ok=False)
    source = Path(config['disease_csv']).parent
    for name in ['primary_master.csv', 'primary_phenotypes.csv', 'alternate_master.csv', 'alternate_phenotypes.csv',
                 'reference_balanced_master.csv', 'reference_balanced_phenotypes.csv', 'exact_feature_names.csv']:
        shutil.copy2(source / name, output / name)
    (output / 'provenance.json').write_text(json.dumps({'active_analysis': config['analysis_id'],
        'feature_count': 120, 'source_run': config['private_run'], 'historical_values_used': False}, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    prepare(parser.parse_args().output)

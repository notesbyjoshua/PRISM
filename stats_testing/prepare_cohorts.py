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
    disease_path, reference_path = ROOT / 'data/master_dataset_disease.csv', ROOT / 'data/master_dataset_healthy.csv'
    disease = pd.read_csv(disease_path, dtype={'patient_id': 'string', 'image_id': 'string'}, low_memory=False)
    reference = pd.read_csv(reference_path)
    if not quality_status(reference.frontal_ok).eq('pass').all():
        raise ValueError('Reference quality requires review')
    reference, reference_manifest = deduplicate_reference(reference)
    features = list(reference.loc[:, 'eb_thickness_r':'cheek_area_asym'])
    if len(features) != 125:
        raise ValueError('Expected 125 FaceKit features')
    output.mkdir(parents=True, exist_ok=False)
    reference.to_csv(output / 'reference_master.csv', index=False)
    reference_manifest.to_csv(output / 'reference_selection_manifest.csv', index=False)
    summary = []
    for name, scope, selection in [('all_ages', 'all', 'pose'), ('age_3_9', '3-9', 'pose'),
                                   ('age_3_9_id_sensitivity', '3-9', 'image_id')]:
        cohort, manifest = select_patients(disease, scope, selection)
        folder = output / name
        folder.mkdir()
        cohort.to_csv(folder / 'disease_master.csv', index=False)
        manifest.to_csv(folder / 'selection_manifest.csv', index=False)
        columns = ['disease', 'image_id', 'frontal_ok', 'pose_yaw', 'pose_pitch', 'pose_roll'] + features
        cohort[columns].to_csv(folder / 'disease_phenotypes.csv', index=False)
        summary.append({'cohort': name, 'patient_rows': len(cohort), 'syndromes': cohort.disease.nunique(),
                        'reference_unique_images': len(reference), 'selection': selection, 'age_scope': scope})
    ref_features = reference.copy()
    ref_features['disease'] = 'healthy'
    ref_features['image_id'] = ref_features['file']
    ref_features[columns].to_csv(output / 'reference_phenotypes.csv', index=False)
    pd.DataFrame(summary).to_csv(output / 'cohort_summary.csv', index=False)
    provenance = {'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [disease_path, reference_path]},
                  'selection': 'one eligible image per patient; lowest max(abs(yaw)/15, abs(pitch)/15, abs(roll)/10), then lexical image ID',
                  'missing_pose': 'infinite score, retained only if no better eligible photo exists',
                  'age': '3 <= year + month/12 < 10; integer components; month 0..11; zero-zero and invalid/missing excluded from childhood cohort',
                  'reference': 'identical repeated file paths removed; no new images or participants inferred',
                  'scope': 'current-data reanalysis; age bin overlap does not establish exact-age balance',
                  'software': {'pandas': pd.__version__, 'numpy': np.__version__}}
    (output / 'cohort_provenance.json').write_text(json.dumps(provenance, indent=2))
    print(pd.DataFrame(summary).to_string(index=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    prepare(parser.parse_args().output)

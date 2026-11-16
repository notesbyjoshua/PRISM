import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, norm
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / 'local_review/prism_robustness_20261110'
SEED = 20261110
MIN_N = 10
KEY = ['canonical_syndrome_id', 'feature']
FORBIDDEN = {'patient_id', 'image_id', 'image_path', 'file', 'key', 'image_sha256', 'recorded_patient_group_key',
             'landmarks', 'embedding', 'source_image_key', 'image_path_key', 'reference_key', 'patient_index'}
CASES = [('Cornelia de Lange syndrome', 'philtrum_length'), ('Cornelia de Lange syndrome', 'philtrum_width'),
         ('Noonan syndrome', 'inter_pupillary_distance'), ('WILLIAMS-BEUREN SYNDROME; WBS', 'eb_thickness_mean'),
         ('WILLIAMS-BEUREN SYNDROME; WBS', 'outer_canthal_distance'), ('Rubinstein-Taybi syndrome', 'malar_bulge_mean')]


class ResearchContractError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise ResearchContractError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def private_output(path):
    path = Path(path).resolve()
    require(path.is_relative_to(ROOT / 'local_review'), 'New research outputs must remain in private local_review/')
    require(subprocess.run(['git', 'check-ignore', '-q', str(path / 'privacy-check.json')], cwd=ROOT).returncode == 0,
            'Private output is not Git-ignored')
    path.mkdir(parents=True, exist_ok=True, mode=0o700)
    return path


def write_json(value, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False, default=str) + '\n')
    temporary.replace(path)


def write_aggregate(frame, path):
    require(not FORBIDDEN.intersection(frame.columns), 'Individual fields prohibited in research exports')
    for column in frame.select_dtypes(include=['object', 'str']):
        require(not frame[column].map(lambda value: isinstance(value, (list, dict, np.ndarray))).any(),
                'Nested individual values prohibited in research exports')
    path = Path(path)
    require(path.resolve().is_relative_to(ROOT / 'local_review'), 'Aggregate candidates require private review first')
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    temporary = path.with_suffix(path.suffix + '.tmp')
    frame.to_csv(temporary, index=False, float_format='%.17g', lineterminator='\n')
    temporary.replace(path)


def bh(values):
    values = np.asarray(values, float)
    require(np.isfinite(values).all() and ((values >= 0) & (values <= 1)).all(), 'Invalid planned correction p-values')
    return multipletests(values, method='fdr_bh')[1]


def corrected(frame, p_column, supported, family):
    result = frame.copy()
    supported = np.asarray(supported, bool)
    result['supported'] = supported
    result['correction_p'] = np.where(supported, result[p_column], 1.)
    result['exploratory_q'] = np.where(supported, bh(result.correction_p), np.nan)
    result['correction_family'] = family
    result['planned_hypotheses'] = len(result)
    return result


def effect(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    n, m = len(x), len(y)
    if min(n, m) < 2:
        return np.full(x.shape[1:] or (), np.nan)
    pooled = ((n - 1) * x.var(axis=0, ddof=1) + (m - 1) * y.var(axis=0, ddof=1)) / (n + m - 2)
    with np.errstate(divide='ignore', invalid='ignore'):
        result = (1 - 3 / (4 * (n + m) - 9)) * (x.mean(axis=0) - y.mean(axis=0)) / np.sqrt(pooled)
    return np.where(pooled > 1e-24, result, np.nan)


def u_test(x, y):
    test = mannwhitneyu(x, y, axis=0, alternative='two-sided', method='asymptotic', use_continuity=True)
    return test.statistic, test.pvalue


def age_years(frame):
    year, month = [pd.to_numeric(frame[name], errors='coerce') for name in ['age_year', 'age_month']]
    valid = year.ge(0) & month.between(0, 11) & year.mod(1).eq(0) & month.mod(1).eq(0)
    return (year + month / 12).where(valid)


def clean_ids(series):
    return series.astype('string').str.strip().str.replace(r'\.0$', '', regex=True).replace({'': pd.NA, 'nan': pd.NA, 'None': pd.NA})


def canonical_map():
    taxonomy = pd.read_csv(ROOT / 'data/syndrome_folder_map.csv', dtype={'internal_syndrome_id': 'string'})
    require(taxonomy.facekit_cohort.is_unique and taxonomy.internal_syndrome_id.is_unique, 'Nonunique canonical syndrome map')
    return taxonomy.rename(columns={'facekit_cohort': 'syndrome', 'internal_syndrome_id': 'canonical_syndrome_id'})[
        ['syndrome', 'canonical_syndrome_id']]


def canonicalize(frame, label='syndrome'):
    require(label in frame, 'Expected exact current cohort label is absent')
    result = frame.rename(columns={label: 'syndrome'}).merge(canonical_map(), on='syndrome', how='left', validate='many_to_one')
    require(result.canonical_syndrome_id.notna().all(), 'Unmapped current syndrome labels; no fuzzy or display-label join permitted')
    return result


def source_table(name, label='syndrome'):
    frame = pd.read_csv(ROOT / 'publication/supplementary' / name)
    return canonicalize(frame, label)


def primary_table():
    return source_table('S1_primary_all_5760_hypotheses.csv', 'disease')


def keyed_merge(left, right, prefix, columns, require_complete=True):
    require(not left.duplicated(KEY).any() and not right.duplicated(KEY).any(), 'Duplicated canonical hypothesis keys')
    right = right[KEY + columns].rename(columns={name: prefix + name for name in columns})
    result = left.copy().merge(right, on=KEY, how='left', validate='one_to_one', indicator=True)
    if require_complete:
        require(result['_merge'].eq('both').all(), 'A required canonical hypothesis grid is incomplete')
    return result.drop(columns='_merge').copy()


def features():
    frame = pd.read_csv(ROOT / 'publication/feature_definitions/current_120_feature_formulas.csv')
    require(len(frame) == 120 and frame.feature.is_unique and frame.order.tolist() == list(range(120)), 'Current ordered feature schema failed')
    other = pd.read_csv(ROOT / 'publication/supplementary/Supplementary_Table_S6_exact_feature_definitions.csv')
    require(frame.equals(other), 'Current S6/formula exports differ')
    return frame


def load_private():
    config = json.loads((ROOT / 'ACTIVE_ANALYSIS.json').read_text())
    run = ROOT / config['private_run']
    paths = {'primary': run / 'COHORTS/primary.parquet', 'reference': run / 'COHORTS/reference_balanced.parquet',
             'all_gmdb': run / 'COHORTS/all_gmdb_candidates.parquet', 'alternate': run / 'COHORTS/alternate.parquet',
             'primary_master': ROOT / config['disease_master_csv'], 'reference_master': ROOT / config['reference_master_csv'],
             'measured_gmdb': ROOT / 'local_review/feature_extraction_120/20261028_current120/gmdb/gmdb_120_all_records_measurements.parquet'}
    require(all(path.is_file() for path in paths.values()), 'BLOCKED_PRIVATE_INPUT: authorized current source tables unavailable')
    schema = features().feature.tolist()
    real, reference, raw, alternate = [pd.read_parquet(paths[name]) for name in ['primary', 'reference', 'all_gmdb', 'alternate']]
    for frame in [real, raw, alternate]:
        frame['patient_id'], frame['image_id'] = clean_ids(frame.patient_id), clean_ids(frame.image_id)
    for frame in [real, alternate]:
        require(frame.patient_id.notna().all() and frame.patient_id.is_unique and frame.image_id.is_unique,
                'One nonmissing recorded patient and selected photograph required')
        require(frame.quality_passed.eq(True).all() and np.isfinite(frame[schema]).all().all(), 'Nonfinite or non-QC selected cohort')
        require(age_years(frame).between(3, 10, inclusive='left').all(), 'Selected real ages fail frozen eligibility')
    require(len(real) == 1845 and len(reference) == 3409 and reference.image_id.is_unique, 'Frozen private sample counts differ')
    require(np.isfinite(reference[schema]).all().all() and reference.quality_passed.eq(True).all(), 'Invalid selected reference')
    for name, frame in [('primary_master', real), ('reference_master', reference)]:
        master = pd.read_csv(paths[name], dtype={'patient_id': 'string', 'image_id': 'string'})
        require(len(master) == len(frame) and clean_ids(master.image_id).eq(clean_ids(frame.image_id)).all(),
                'Frozen CSV/Parquet source ordering differs')
        require(np.allclose(master[schema], frame[schema], rtol=1e-12, atol=1e-14), 'Frozen CSV/Parquet features differ')
    require(real.disease.value_counts().ge(10).sum() == 48 and real.disease.value_counts().loc[lambda v: v >= 10].sum() == 872,
            'Frozen 48-cohort / 872-patient primary subset differs')
    mapping = canonical_map()
    own = real[['disease', 'internal_syndrome_id']].drop_duplicates().rename(columns={'disease': 'syndrome'})
    own['internal_syndrome_id'] = own.internal_syndrome_id.astype('string')
    own = own.merge(mapping, on='syndrome', validate='one_to_one')
    require(own.internal_syndrome_id.eq(own.canonical_syndrome_id).all(), 'Current private canonical IDs disagree with shared taxonomy')
    measured = pd.read_parquet(paths['measured_gmdb'])
    measured['image_id'] = clean_ids(measured.image_id)
    require(measured.image_id.is_unique, 'Measured source image keys are nonunique')
    indices = pd.Index(measured.image_id).get_indexer(real.image_id)
    require((indices >= 0).all(), 'Selected measurements are missing from current measured source')
    selected = measured.iloc[indices].reset_index(drop=True)
    require(np.allclose(selected[schema], real[schema], rtol=1e-12, atol=1e-14), 'Selected/current source measurement equality failed')
    for frame, commit, model in [(selected, 'facekit_commit', 'model_sha256'), (reference, 'facekit_commit', 'model_sha256')]:
        require(frame[commit].eq(config['facekit_commit']).all() and frame[model].eq(config['model_sha256']).all(), 'Measured-source/model pins differ')
    for name in ['width', 'height', 'screen_detected_faces', 'landmark_count', 'derotated']:
        real[name] = selected[name].to_numpy()
    for frame in [real, reference, alternate]:
        for axis in ['yaw', 'pitch', 'roll']:
            frame['abs_' + axis] = frame['pose_' + axis].abs()
        frame['pose_score'] = np.maximum.reduce([frame.abs_yaw / 15, frame.abs_pitch / 15, frame.abs_roll / 10])
    require(reference.facekit_commit.eq(config['facekit_commit']).all(), 'Synthetic/reference source separation failed')
    return real, reference, raw, alternate, paths


def case_rows(matrix):
    pairs = set(CASES)
    mask = [(str(s).split('__')[0], feature) in pairs for s, feature in zip(matrix.syndrome, matrix.feature)]
    result = matrix.loc[mask].copy()
    require(len(result) == len(CASES), 'Predeclared case coverage failed')
    return result


def offline():
    os.umask(0o077)
    def hook(event, args):
        if event in {'socket.connect', 'socket.getaddrinfo', 'http.client.connect', 'urllib.Request'}:
            raise ResearchContractError('Protected-data analysis is offline')
    sys.addaudithook(hook)


def protect_baseline(out, create=False):
    path = out / '00_frozen_file_hashes.json'
    if create and not path.exists():
        tracked = subprocess.check_output(['git', 'ls-files'], cwd=ROOT, text=True).splitlines()
        tracked.append('data/dataset_creation/reloadHFCache.py')
        write_json({name: sha(ROOT / name) for name in tracked if (ROOT / name).is_file()}, path)
    require(path.is_file(), 'Run mandatory preflight before later milestones')
    expected = json.loads(path.read_text())
    require(all((ROOT / name).is_file() and sha(ROOT / name) == digest for name, digest in expected.items()),
            'A frozen baseline or user-owned file changed; stop for review')
    private_manifest = out / '00_baseline_validation.json'
    if private_manifest.exists():
        sources = json.loads(private_manifest.read_text())['protected_input_sources']
        require(all((ROOT / item['path']).is_file() and sha(ROOT / item['path']) == item['sha256'] for item in sources.values()),
                'A protected private source changed; stop for review')
    return len(expected)

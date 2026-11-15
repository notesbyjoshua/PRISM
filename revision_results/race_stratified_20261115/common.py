import hashlib
import importlib.metadata
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PRIVATE = ROOT / 'local_review/race_stratified_20261115'
EXTRACTION = ROOT / 'local_review/feature_extraction_120/20261028_current120'
RACES = ('African', 'Asian', 'European')
SEED = 20261115
os.umask(0o077)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def private_output(path):
    path = Path(path).resolve()
    if not path.is_relative_to((ROOT / 'local_review').resolve()):
        raise ValueError('All research outputs require a directory inside local_review')
    result = subprocess.run(['git', 'check-ignore', '-q', str(path / 'privacy_probe.csv')], cwd=ROOT)
    if result.returncode or subprocess.check_output(['git', 'ls-files', '--', str(path)], cwd=ROOT).strip():
        raise ValueError('Output must be Git-ignored and contain no tracked files')
    if path.exists() and any(p.is_symlink() for p in path.rglob('*')):
        raise ValueError('Private output may not contain symlinks to other locations')
    path.mkdir(parents=True, exist_ok=True, mode=0o700)
    return path


def normalize(value):
    if pd.isna(value):
        return ''
    return re.sub(r'\s+', ' ', re.sub(r'[_/\-]', ' ', str(value).strip().lower())).strip()


def map_group(value, source):
    value = normalize(value)
    if source == 'GMDB':
        mapping = {'european': 'European', 'europe': 'European', 'caucasian': 'European', 'white': 'European',
                   'african': 'African', 'africa': 'African', 'black': 'African', 'asian': 'Asian', 'asia': 'Asian'}
    elif source == 'FairFace':
        mapping = {'black': 'African', 'african': 'African', 'white': 'European', 'european': 'European'}
        for label in ['east asian', 'central asian', 'west asian', 'western asian', 'middle eastern', 'middle east',
                      'indian', 'south asian', 'southeast asian', 'south east asian', 'se asian', 'asian']:
            mapping[label] = 'Asian'
    else:
        raise ValueError('Unknown source')
    return mapping.get(value)


def validate_features(frame, features):
    if len(features) != 120 or len(set(features)) != 120:
        raise ValueError('Exactly 120 distinct features required')
    if [column for column in frame if column in set(features)] != features:
        raise ValueError('Feature order differs from the pinned schema')
    if not np.isfinite(frame[features].to_numpy(float)).all():
        raise ValueError('Nonfinite feature measurements are not eligible')
    return frame[features].astype(float)


def validate_patients(frame):
    for field in ['patient_id', 'image_id']:
        values = frame[field].astype('string')
        if values.isna().any() or values.str.strip().eq('').any() or values.duplicated().any():
            raise ValueError(f'Unique nonmissing {field} required; duplicate images cannot inflate cells')
    if frame.disease.isna().any() or frame.disease.astype(str).str.strip().eq('').any():
        raise ValueError('Missing exact disease identifier')


def screen_reference_content(reference):
    if reference.image_sha256.isna().any() or not reference.file.is_unique:
        raise ValueError('FairFace requires content hashes and unique paths')
    conflicts = reference.groupby('image_sha256').race.nunique().gt(1)
    conflict_rows = reference.image_sha256.isin(conflicts.index[conflicts])
    dedup = reference.loc[~conflict_rows].sort_values('file').drop_duplicates('image_sha256').copy()
    return dedup, {'fairface_conflicting_label_hash_sets': int(conflicts.sum()),
                   'fairface_conflicting_label_rows_excluded': int(conflict_rows.sum()),
                   'fairface_repeated_content_rows_removed': len(reference) - int(conflict_rows.sum()) - len(dedup)}


def protected_hashes():
    paths = [ROOT / 'ACTIVE_ANALYSIS.json']
    for directory in ['stats_results', 'manuscript', 'publication', 'stats_testing',
                      'revision_results/rebuild_120_20261028', 'revision_results/reviewer_followup_2026', 'ml_pilot']:
        paths.extend(p for p in (ROOT / directory).rglob('*') if p.is_file()
                     and not any(part in {'__pycache__', '.pytest_cache', 'local_only'} for part in p.parts))
    return {str(p.relative_to(ROOT)): sha(p) for p in sorted(set(paths))}


def provenance(sources, args):
    code = [p for p in HERE.rglob('*') if p.suffix in {'.py', '.txt', '.md'} and '__pycache__' not in p.parts]
    return {'study_id': 'PRISM_RACE_STRATIFIED_20261115', 'status': 'EXPERIMENTAL_PRIVATE',
            'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
            'code_sha256': {str(p.relative_to(ROOT)): sha(p) for p in sorted(code)},
            'source_sha256': {str(p.relative_to(ROOT)): sha(p) for p in sorted(set(sources))},
            'versions': {name: importlib.metadata.version(name) for name in
                         ['numpy', 'pandas', 'scipy', 'scikit-learn', 'matplotlib', 'pyarrow', 'joblib', 'threadpoolctl']},
            'python': sys.version, 'parameters': args,
            'sensitivity_status': 'NOT RUN: deferred by direct user instruction',
            'public_export_status': 'NOT RUN: private review only; no publishing capability in these entrypoints'}


def load_sources():
    config_path = ROOT / 'ACTIVE_ANALYSIS.json'
    config = json.loads(config_path.read_text())
    if (config['status'] != 'VALIDATED_ACTIVE' or config['feature_count'] != 120
            or config['synthetic_rows_in_primary'] != 0):
        raise ValueError('Validated real-source 120-feature active configuration required')
    cohorts = ROOT / config['private_run'] / 'COHORTS'
    paths = [config_path] + [ROOT / config[k] for k in
                            ['disease_csv', 'reference_csv', 'disease_master_csv', 'reference_master_csv', 'feature_schema']]
    paths += [cohorts / name for name in ['primary.parquet', 'reference_all_eligible.parquet',
                                        'reference_all_eligible_master.csv', 'private_primary_selection_manifest.csv',
                                        'syndrome_eligibility.csv']]
    paths += [EXTRACTION / 'gmdb/gmdb_120_all_records_measurements.parquet']
    paths += [EXTRACTION / source / 'provenance.json' for source in ['gmdb', 'fairface']]
    definition_path = ROOT / 'publication/feature_definitions/current_120_feature_formulas.csv'
    paths.append(definition_path)
    missing = [str(p.relative_to(ROOT)) for p in paths if not p.is_file()]
    if missing:
        raise FileNotFoundError('BLOCKED: missing authorized inputs: ' + ', '.join(missing))
    features = pd.read_csv(ROOT / config['feature_schema']).feature.tolist()
    definitions = pd.read_csv(definition_path)
    if definitions.feature.tolist() != features or not definitions.source_commit.eq(config['facekit_commit']).all():
        raise ValueError('Feature definitions must match exact extraction schema and commit')
    for source in ['gmdb', 'fairface']:
        recorded = json.loads((EXTRACTION / source / 'provenance.json').read_text())
        if (recorded['facekit_commit'] != config['facekit_commit'] or recorded['model_sha256'] != config['model_sha256']
                or recorded['features'] != features or recorded['status'] != 'FULL'):
            raise ValueError('Pinned extraction provenance mismatch')
    disease = pd.read_csv(ROOT / config['disease_master_csv'], float_precision='round_trip',
                          dtype={'patient_id': 'string', 'image_id': 'string'}).copy()
    reference = pd.read_csv(cohorts / 'reference_all_eligible_master.csv', float_precision='round_trip',
                            dtype={'image_id': 'string'}).copy()
    for field in ['disease_csv', 'reference_csv', 'reference_master_csv']:
        validate_features(pd.read_csv(ROOT / config[field], float_precision='round_trip'), features)
    for frame in [disease, reference]:
        validate_features(frame, features)
        if not frame.frontal_ok.eq(True).all():
            raise ValueError('All rows must pass frontal QC')
    validate_patients(disease)
    if not disease.conditional_age_years.between(3, 10, inclusive='left').all():
        raise ValueError('GMDB age eligibility failed')
    if (not reference.age.eq('3-9').all() or not reference.file.is_unique
            or reference.file.isna().any() or not reference.image_id.is_unique):
        raise ValueError('FairFace age/path eligibility failed')
    gmdb_qc = pd.read_parquet(cohorts / 'primary.parquet')
    fairface_qc = pd.read_parquet(cohorts / 'reference_all_eligible.parquet')
    for master, original in [(disease, gmdb_qc), (reference, fairface_qc)]:
        original.image_id = original.image_id.astype('string')
        ordered = original.set_index('image_id').loc[master.image_id]
        if not ordered.quality_passed.eq(True).all() or not ordered.reason.eq('accepted').all():
            raise ValueError('Nonaccepted source observation')
        np.testing.assert_array_equal(master[features].to_numpy(float), ordered[features].to_numpy(float))
        fields = ['patient_id', 'disease', 'conditional_age_years'] if 'patient_id' in master else ['file', 'race', 'age']
        for field in fields:
            if master[field].astype(str).tolist() != ordered[field].astype(str).tolist():
                raise ValueError('Source identity, label or age mismatch: ' + field)
    if (not fairface_qc.facekit_commit.eq(config['facekit_commit']).all()
            or not fairface_qc.model_sha256.eq(config['model_sha256']).all()):
        raise ValueError('FairFace extraction pin mismatch')
    raw = pd.read_parquet(EXTRACTION / 'gmdb/gmdb_120_all_records_measurements.parquet')
    raw.image_id = raw.image_id.astype('string')
    original = raw.set_index('image_id').loc[disease.image_id]
    if (not original.quality_passed.eq(True).all() or not original.facekit_commit.eq(config['facekit_commit']).all()
            or not original.model_sha256.eq(config['model_sha256']).all()):
        raise ValueError('GMDB extraction provenance mismatch')
    if (not original.dataset.eq('GMDB_restored').all() or not fairface_qc.dataset.eq('FairFace').all()
            or original.patient_id.astype(str).tolist() != disease.patient_id.astype(str).tolist()):
        raise ValueError('Real source/patient linkage failed')
    np.testing.assert_array_equal(disease[features].to_numpy(float), original[features].to_numpy(float))
    disease['image_sha256'] = original.image_sha256.to_numpy()
    if disease.image_sha256.isna().any() or disease.image_sha256.duplicated().any():
        raise ValueError('GMDB duplicate content or missing content hashes; resolve selection before inference')
    selection = pd.read_csv(cohorts / 'private_primary_selection_manifest.csv',
                            dtype={'patient_id': 'string', 'image_id': 'string'})
    chosen = selection.loc[selection.selection_reason.eq('selected')]
    if set(zip(chosen.patient_id, chosen.image_id)) != set(zip(disease.patient_id, disease.image_id)):
        raise ValueError('GMDB selected-image manifest mismatch')
    reference['image_sha256'] = fairface_qc.set_index('image_id').loc[reference.image_id].image_sha256.to_numpy()
    dedup, content_audit = screen_reference_content(reference)
    audit = {'gmdb_selected_patients': len(disease), 'fairface_accepted_paths_before_content_qc': len(reference),
             **content_audit,
             'fairface_after_content_qc': len(dedup), 'gmdb_duplicate_content': 0,
             'gmdb_selection_reasons': selection.selection_reason.value_counts().to_dict(),
             'gmdb_age_semantics': 'Existing conditional age interpretation retained; not independently reverified',
             'unknown_identity_linkage': 'Recorded GMDB patients only; FairFace unique people unverified'}
    disease['canonical_race'] = disease.image_dataset_ethnicity_category.map(lambda v: map_group(v, 'GMDB'))
    dedup['canonical_race'] = dedup.race.map(lambda v: map_group(v, 'FairFace'))
    audit['gmdb_source_labels'] = disease.image_dataset_ethnicity_category.value_counts(dropna=False).to_dict()
    audit['fairface_source_labels_before_content_qc'] = reference.race.value_counts(dropna=False).to_dict()
    audit['fairface_source_labels_after_content_qc'] = dedup.race.value_counts(dropna=False).to_dict()
    audit['gmdb_unmapped_excluded'] = int(disease.canonical_race.isna().sum())
    audit['fairface_unmapped_excluded'] = int(dedup.canonical_race.isna().sum())
    return disease, dedup, features, definitions, audit, paths


def cohort_counts(disease):
    validate_patients(disease)
    totals = disease.groupby('disease').size()
    mapped = disease.loc[disease.canonical_race.isin(RACES)]
    counts = mapped.groupby(['disease', 'canonical_race']).size()
    index = pd.MultiIndex.from_product([sorted(disease.disease.unique()), RACES], names=['disease', 'canonical_race'])
    table = counts.reindex(index, fill_value=0).rename('n_current_gmdb').reset_index()
    table['unique_patient_count'] = table.n_current_gmdb
    table['image_count'] = table.n_current_gmdb
    table['eligible'] = table.n_current_gmdb.ge(5)
    table['reason'] = np.where(table.eligible, 'eligible', 'fewer_than_five_patients')
    table['previous_study_ge10'] = table.disease.map(totals).ge(10)
    table['previous_total_patients'] = table.disease.map(totals)
    table['umbrella_label'] = table.disease.str.split('__').str[0].eq('Mucopolysaccharidoses')
    denominator = int(table.loc[table.eligible, 'n_current_gmdb'].sum())
    table['percent_of_all_mapped_eligible_cases'] = 100 * table.n_current_gmdb / denominator if denominator else np.nan
    labels = mapped.groupby(['disease', 'canonical_race']).image_dataset_ethnicity_category.agg(
        lambda values: json.dumps(values.value_counts().to_dict(), sort_keys=True))
    table['source_label_counts'] = [labels.get((d, r), '{}') for d, r in zip(table.disease, table.canonical_race)]
    return table


def select_ml(counts):
    eligible = counts.loc[counts.eligible].sort_values(['n_current_gmdb', 'disease', 'canonical_race'],
                                                     ascending=[False, True, True])
    global_top = eligible.head(5).copy()
    candidates = []
    for race in RACES:
        top = eligible.loc[eligible.canonical_race.eq(race)].head(5)
        if len(top) == 5 and top.n_current_gmdb.min() >= 10:
            candidates.append((race, top))
    preferred = [item for item in candidates if item[1].n_current_gmdb.min() >= 15]
    candidates = preferred or candidates
    candidates.sort(key=lambda item: (-item[1].n_current_gmdb.sum(), -item[1].n_current_gmdb.min(), item[0]))
    selected = candidates[0][1].copy() if candidates else eligible.head(0).copy()
    return global_top, selected


def apportion(counts, capacities, maximum=3409):
    counts = np.array([counts.get(race, 0) for race in RACES], dtype=np.int64)
    capacity = np.array([capacities.get(race, 0) for race in RACES], dtype=np.int64)
    if counts.sum() <= 0 or (counts < 0).any() or (capacity < 0).any() or maximum < 1:
        raise ValueError('Invalid quota inputs')
    denominator = int(counts.sum())
    for size in range(min(maximum, int(capacity.sum())), 0, -1):
        numerators = size * counts
        quota = numerators // denominator
        residual = size - int(quota.sum())
        order = sorted(range(3), key=lambda i: (-int(numerators[i] % denominator), RACES[i]))
        quota[order[:residual]] += 1
        if (quota <= capacity).all():
            assert quota.sum() == size
            assert np.max(np.abs(quota / size - counts / denominator)) <= 1 / size + 1e-12
            return size, dict(zip(RACES, map(int, quota)))
    raise ValueError('No positive feasible proportional reference size')


def sample_reference(reference, quotas, seed):
    rng = np.random.default_rng(seed)
    pieces = []
    for race in RACES:
        pool = reference.loc[reference.canonical_race.eq(race)].sort_values('file')
        if quotas[race] > len(pool):
            raise ValueError('Insufficient within-category reference capacity')
        selected = pool.iloc[rng.choice(len(pool), quotas[race], replace=False)]
        if not selected.canonical_race.eq(race).all():
            raise ValueError('Cross-category reference contamination')
        pieces.append(selected)
    sampled = pd.concat(pieces, ignore_index=True)
    if not sampled.file.is_unique or len(sampled) != sum(quotas.values()):
        raise ValueError('Reference sampling is not without replacement')
    return sampled


def matched_frames(disease, reference, syndrome, race):
    cases = disease.loc[disease.disease.eq(syndrome) & disease.canonical_race.eq(race)]
    controls = reference.loc[reference.canonical_race.eq(race)]
    assert cases.canonical_race.eq(race).all() and controls.canonical_race.eq(race).all()
    return cases, controls

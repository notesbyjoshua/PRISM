import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, wasserstein_distance

ROOT = Path(__file__).resolve().parents[1]
PUBLIC = ROOT / 'publication/synthetic_consistency'
SEED = 20261109
MIN_N = 10
EPS = 1e-12
ANCHORS = [('Cornelia de Lange syndrome__fb4a389d', 'philtrum_length'),
           ('Noonan syndrome__c2df3f82', 'inter_pupillary_distance'),
           ('WILLIAMS-BEUREN SYNDROME; WBS__93e4a782', 'eb_thickness_mean'),
           ('WILLIAMS-BEUREN SYNDROME; WBS__93e4a782', 'outer_canthal_distance'),
           ('Rubinstein-Taybi syndrome__ad6c2718', 'malar_bulge_mean')]
PROHIBITED = {'patient_id', 'image_id', 'image_sha256', 'key', 'image_path', 'file', 'landmarks', 'embedding',
              'recorded_patient_group_key', 'source_image_key', 'image_path_key', 'reference_key'}
REAL_METADATA = {'image_id', 'patient_id', 'disease', 'internal_syndrome_name', 'image_dataset_ethnicity_category',
                 'image_dataset_ethnicity_sub_category', 'gender', 'age_year', 'age_month', 'conditional_age_years',
                 'broad_group', 'detailed_group', 'sex_label', 'frontal_ok', 'pose_yaw', 'pose_pitch', 'pose_roll'}
REFERENCE_METADATA = {'image_id', 'file', 'gender', 'race', 'age', 'broad_group', 'detailed_group', 'sex_label',
                      'frontal_ok', 'pose_yaw', 'pose_pitch', 'pose_roll'}
SYNTHETIC_METADATA = {'image_id', 'disease', 'disease_abbr', 'model', 'verification', 'age_range', 'gender', 'race',
                      'key', 'source_dataset', 'synthetic', 'dataset_revision', 'recorded_patient_group_key',
                      'source_linkage_status', 'analytic_role', 'reason', 'quality_passed', 'feature_finite',
                      'width', 'height', 'exif_orientation', 'screen_detected_faces', 'pose_yaw', 'pose_pitch',
                      'pose_roll', 'frontal_ok', 'derotated', 'landmark_count', 'image_sha256', 'extractor_commit',
                      'model_hash'}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def validate_schema(frame, features, metadata):
    observed = [c for c in frame if c not in metadata]
    if len(features) != 120 or len(set(features)) != 120 or observed != features:
        raise ValueError('Exactly 120 ordered current features required; unexpected/missing/renamed columns')
    if not np.isfinite(frame[features].to_numpy(float)).all():
        raise ValueError('Nonfinite accepted measurements; stop rather than silently dropping cells')


def validate_real(frame):
    if frame.patient_id.isna().any() or not frame.patient_id.is_unique or not frame.image_id.is_unique:
        raise ValueError('Real cohort must contain one image per nonmissing recorded patient')
    if not frame.frontal_ok.eq(True).all():
        raise ValueError('Real cohort is not quality passing')
    if 'synthetic' in frame and frame.synthetic.ne(False).any():
        raise ValueError('Synthetic rows are prohibited in the real primary cohort')
    if not frame.conditional_age_years.between(3, 10, inclusive='left').all():
        raise ValueError('Real cohort ages must satisfy 3 <= age < 10')


def reviewed_mapping(mapping, real_labels, synthetic_labels):
    if mapping.pdidb_disease.duplicated().any():
        raise ValueError('Multiply mapped synthetic labels require explicit resolution')
    accepted = mapping.loc[mapping.mapping_status.eq('reviewed_exact_name')]
    if accepted.gmdb_syndrome.isna().any() or accepted.gmdb_syndrome.duplicated().any():
        raise ValueError('Accepted label mapping must be one-to-one')
    if not set(accepted.gmdb_syndrome).issubset(real_labels):
        raise ValueError('Reviewed mapping targets are absent from the current GMDB cohort')
    result = mapping.set_index('pdidb_disease').reindex(sorted(synthetic_labels)).reset_index()
    result['mapping_status'] = result.mapping_status.fillna('excluded_unmatched')
    result['review_basis'] = result.review_basis.fillna('No explicit reviewed mapping; no fuzzy match attempted')
    return result


def summary(values):
    values = np.asarray(values, float)
    finite = values[np.isfinite(values)]
    if len(finite) == 0:
        return dict(n=0, missing=len(values), mean=np.nan, sd=np.nan, median=np.nan, q1=np.nan, q3=np.nan)
    q1, median, q3 = np.quantile(finite, [.25, .5, .75])
    return dict(n=len(finite), missing=len(values) - len(finite), mean=finite.mean(),
                sd=finite.std(ddof=1) if len(finite) > 1 else np.nan, median=median, q1=q1, q3=q3)


def hedges_g(group, reference):
    x, y = summary(group), summary(reference)
    if min(x['n'], y['n']) < 2:
        return np.nan
    pooled = np.sqrt(((x['n'] - 1) * x['sd']**2 + (y['n'] - 1) * y['sd']**2) / (x['n'] + y['n'] - 2))
    return (x['mean'] - y['mean']) / pooled * (1 - 3 / (4 * (x['n'] + y['n']) - 9)) if pooled > EPS else np.nan


def feature_metrics(real, synthetic, reference):
    a, b, ref = [summary(v) for v in [real, synthetic, reference]]
    row = {f'{source}_{k}': v for source, values in [('real', a), ('synthetic', b), ('reference', ref)]
           for k, v in values.items()}
    row['mean_shift'] = b['mean'] - a['mean']
    row['median_shift'] = b['median'] - a['median']
    row['real_sd_defined'] = bool(a['sd'] > EPS)
    row['synthetic_sd_defined'] = bool(b['sd'] > EPS)
    row['reference_sd_defined'] = bool(ref['sd'] > EPS)
    row['small_n'] = min(a['n'], b['n']) < MIN_N
    row['standardized_mean_shift'] = row['mean_shift'] / a['sd'] if row['real_sd_defined'] else np.nan
    row['standardized_median_shift'] = row['median_shift'] / a['sd'] if row['real_sd_defined'] else np.nan
    row['variance_ratio'] = b['sd']**2 / a['sd']**2 if row['real_sd_defined'] else np.nan
    finite_real, finite_synthetic = [np.asarray(v)[np.isfinite(v)] for v in [real, synthetic]]
    row['wasserstein'] = wasserstein_distance(finite_real, finite_synthetic) if len(finite_real) and len(finite_synthetic) else np.nan
    row['wasserstein_real_sd'] = row['wasserstein'] / a['sd'] if row['real_sd_defined'] else np.nan
    row['real_g'], row['synthetic_g'] = hedges_g(real, reference), hedges_g(synthetic, reference)
    for source, values in [('real', a), ('synthetic', b)]:
        row[f'{source}_reference_sd_contrast'] = (values['mean'] - ref['mean']) / ref['sd'] if ref['sd'] > EPS else np.nan
    row['g_difference'] = row['synthetic_g'] - row['real_g']
    return row


def correlations(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 3 or np.ptp(x) <= EPS or np.ptp(y) <= EPS:
        return np.nan, np.nan
    return float(pearsonr(x, y).statistic), float(spearmanr(x, y).statistic)


def profile_summary(frame):
    valid = frame.loc[np.isfinite(frame.real_g) & np.isfinite(frame.synthetic_g)]
    strong = valid.loc[valid.real_g.abs().ge(.2)]
    r, rho = correlations(valid.real_g, valid.synthetic_g)
    cr, cs = correlations(frame.real_reference_sd_contrast, frame.synthetic_reference_sd_contrast)
    ratios = frame.variance_ratio.dropna()
    return {'features': len(frame), 'defined_profile_features': len(valid),
            'median_absolute_standardized_shift': frame.standardized_mean_shift.abs().median(),
            'median_variance_ratio': ratios.median(), 'variance_ratio_defined_features': len(ratios),
            'fraction_synthetic_variance_lower': ratios.lt(1).mean() if len(ratios) else np.nan,
            'median_wasserstein_real_sd': frame.wasserstein_real_sd.median(),
            'direction_agreement': np.mean(np.sign(valid.real_g) == np.sign(valid.synthetic_g)) if len(valid) else np.nan,
            'sign_reversals': int((valid.real_g * valid.synthetic_g < 0).sum()),
            'strong_real_features': len(strong),
            'strong_direction_agreement': np.mean(np.sign(strong.real_g) == np.sign(strong.synthetic_g)) if len(strong) else np.nan,
            'pearson_g': r, 'spearman_g': rho, 'pearson_reference_sd': cr, 'spearman_reference_sd': cs,
            'median_absolute_g_difference': valid.g_difference.abs().median(),
            'age_alignment': 'unavailable: all PDIDB age_range missing; real ages 3 <= age < 10',
            'interpretation': 'conditional within-source fidelity; correlated features; no independent replication'}


def save_aggregate(frame, path):
    if PROHIBITED.intersection(frame.columns):
        raise ValueError('Individual identifiers/measurements prohibited in aggregate output')
    if any(frame[c].map(lambda v: isinstance(v, (list, dict, np.ndarray))).any() for c in frame.select_dtypes(include=['object', 'str'])):
        raise ValueError('Nested per-record values prohibited in aggregate output')
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, float_format='%.17g', lineterminator='\n')


def load_inputs(config_path):
    config = json.loads(config_path.read_text())
    if (config['analysis_id'] != 'PRISM_120_20261028' or config['feature_count'] != 120
            or config['synthetic_rows_in_primary'] != 0 or config['status'] != 'VALIDATED_ACTIVE'):
        raise ValueError('Frozen current real-source primary study required')
    run = ROOT / config['private_run']
    paths = {key: ROOT / config[key] for key in ['disease_master_csv', 'reference_master_csv', 'feature_schema']}
    paths.update(synthetic_all=run / 'PDIDB/pdidb_synthetic_120_all_records.parquet',
                 synthetic_qc=run / 'PDIDB/pdidb_synthetic_120_quality_passed.parquet',
                 provenance=run / 'PDIDB/provenance.json', real_selected=run / 'COHORTS/primary.parquet',
                 real_measured=ROOT / 'local_review/feature_extraction_120/20261028_current120/gmdb/gmdb_120_all_records_measurements.parquet',
                 reference_selected=run / 'COHORTS/reference_balanced.parquet')
    for name, path in paths.items():
        if not path.is_file():
            raise FileNotFoundError(f'Authorized current private input unavailable ({name}): {path}')
    features = pd.read_csv(paths['feature_schema']).feature.tolist()
    definitions = pd.read_csv(ROOT / 'publication/feature_definitions/current_120_feature_formulas.csv')
    provenance = json.loads(paths['provenance'].read_text())
    if features != definitions.feature.tolist() or provenance['features'] != features:
        raise ValueError('Feature definitions/order disagree across sources')
    if (provenance['facekit_commit'] != config['facekit_commit'] or provenance['model_sha256'] != config['model_sha256']
            or not definitions.source_commit.eq(config['facekit_commit']).all()
            or not definitions.extractor_sha256.eq(provenance['extractor_sha256']).all()
            or sha(ROOT / 'FaceKit/src/facekit/core/geometric/extractor.py') != provenance['extractor_sha256']):
        raise ValueError('Pinned source/model/extractor provenance differs')
    real = pd.read_csv(paths['disease_master_csv'], dtype={'patient_id': 'string', 'image_id': 'string'})
    reference = pd.read_csv(paths['reference_master_csv'], dtype={'image_id': 'string'})
    all_synthetic = pd.read_parquet(paths['synthetic_all'])
    synthetic = pd.read_parquet(paths['synthetic_qc'])
    if len(all_synthetic) != 1727 or len(synthetic) != 1563 or all_synthetic.quality_passed.sum() != 1563:
        raise ValueError('PDIDB counts disagree with frozen 1,727 total / 1,563 QC-passing extraction')
    pd.testing.assert_frame_equal(synthetic.reset_index(drop=True),
                                  all_synthetic.loc[all_synthetic.quality_passed].reset_index(drop=True), check_dtype=False)
    if (not all_synthetic.synthetic.eq(True).all() or not all_synthetic.source_dataset.eq('PDIDB').all()
            or not all_synthetic.extractor_commit.eq(config['facekit_commit']).all()
            or not all_synthetic.model_hash.eq(config['model_sha256']).all()
            or not all_synthetic.dataset_revision.eq(provenance['dataset_revision']).all()
            or not all_synthetic.key.is_unique or not synthetic.quality_passed.eq(True).all()):
        raise ValueError('Synthetic source separation, uniqueness or row-level pins failed')
    for frame, metadata in [(real, REAL_METADATA), (reference, REFERENCE_METADATA), (synthetic, SYNTHETIC_METADATA)]:
        validate_schema(frame, features, metadata)
    validate_real(real)
    if len(real) != 1845 or len(reference) != 3409 or not reference.image_id.is_unique:
        raise ValueError('Frozen real/reference cohort counts or unique reference paths differ')
    if real.disease.value_counts().ge(MIN_N).sum() != 48:
        raise ValueError('Primary eligibility no longer contains exactly 48 syndrome cohorts')
    selected = pd.read_parquet(paths['real_selected'])
    selected['image_id'] = selected.image_id.astype('string')
    selected = selected.set_index('image_id').loc[real.image_id]
    np.testing.assert_allclose(selected[features], real[features], rtol=1e-12, atol=1e-14)
    if not selected.quality_passed.eq(True).all() or not selected.patient_id.astype('string').to_numpy().tolist() == real.patient_id.tolist():
        raise ValueError('Real master disagrees with frozen selected-patient source')
    measured = pd.read_parquet(paths['real_measured'])
    measured['image_id'] = measured.image_id.astype('string')
    measured = measured.set_index('image_id').loc[real.image_id]
    np.testing.assert_allclose(measured[features], real[features], rtol=1e-12, atol=1e-14)
    if not measured.facekit_commit.eq(config['facekit_commit']).all() or not measured.model_sha256.eq(config['model_sha256']).all():
        raise ValueError('Real measured-source pins disagree')
    ref_selected = pd.read_parquet(paths['reference_selected'])
    ref_selected['image_id'] = ref_selected.image_id.astype('string')
    ref_selected = ref_selected.set_index('image_id').loc[reference.image_id]
    np.testing.assert_allclose(ref_selected[features], reference[features], rtol=1e-12, atol=1e-14)
    if not ref_selected.facekit_commit.eq(config['facekit_commit']).all() or not ref_selected.model_sha256.eq(config['model_sha256']).all():
        raise ValueError('Reference measured-source pins disagree')
    return config, paths, definitions, provenance, real, reference, all_synthetic, synthetic, measured


def audit_sources(mapping, real, all_synthetic, synthetic, measured, features):
    rows, metadata = [], []
    for item in mapping.itertuples(index=False):
        a = real.loc[real.disease.eq(item.gmdb_syndrome)]
        b = synthetic.loc[synthetic.disease.eq(item.pdidb_disease)]
        total = all_synthetic.loc[all_synthetic.disease.eq(item.pdidb_disease)]
        supported = item.mapping_status == 'reviewed_exact_name' and min(len(a), len(b)) >= MIN_N
        rows.append({'pdidb_disease': item.pdidb_disease, 'gmdb_syndrome': item.gmdb_syndrome,
                     'mapping_status': item.mapping_status, 'review_basis': item.review_basis,
                     'real_n': len(a) if item.mapping_status == 'reviewed_exact_name' else None,
                     'synthetic_total_n': len(total), 'synthetic_qc_n': len(b),
                     'matched_features': len(features) if item.mapping_status == 'reviewed_exact_name' else 0,
                     'descriptive_eligible': supported, 'reason': 'supported' if supported else item.mapping_status,
                     'real_age_rule': '3 <= age < 10', 'synthetic_age_available': int(b.age_range.notna().sum()),
                     'synthetic_gender_available': int(b.gender.notna().sum()), 'synthetic_race_available': int(b.race.notna().sum()),
                     'metadata_alignment': 'age/sex/race alignment cannot be established'})
        for source, frame, fields in [('GMDB_selected', a, ['sex_label', 'broad_group']),
                                      ('PDIDB_total', total, ['age_range', 'gender', 'race', 'model', 'verification']),
                                      ('PDIDB_QC', b, ['age_range', 'gender', 'race', 'model', 'verification'])]:
            for field in fields:
                for value, n in frame[field].fillna('<missing>').value_counts().sort_index().items():
                    metadata.append({'syndrome': item.pdidb_disease, 'source': source, 'field': field,
                                     'recorded_value': str(value), 'n': int(n), 'source_total_n': len(frame)})
    represented = set(mapping.loc[mapping.mapping_status.eq('reviewed_exact_name'), 'gmdb_syndrome'])
    for label, count in real.disease.value_counts().sort_index().items():
        if count >= MIN_N and label not in represented:
            rows.append({'pdidb_disease': None, 'gmdb_syndrome': label, 'mapping_status': 'excluded_no_synthetic_label',
                         'review_basis': 'Current eligible primary syndrome not represented by an explicit synthetic label',
                         'real_n': count, 'synthetic_total_n': 0, 'synthetic_qc_n': 0, 'matched_features': 0,
                         'descriptive_eligible': False, 'reason': 'no synthetic label', 'real_age_rule': '3 <= age < 10',
                         'synthetic_age_available': 0, 'synthetic_gender_available': 0, 'synthetic_race_available': 0,
                         'metadata_alignment': 'no synthetic counterpart'})
    for source, frame, fields in [('GMDB_all_selected', real, ['sex_label', 'broad_group']),
                                  ('PDIDB_all', all_synthetic, ['age_range', 'gender', 'race', 'model', 'verification'])]:
        for field in fields:
            for value, n in frame[field].fillna('<missing>').value_counts().sort_index().items():
                metadata.append({'syndrome': 'all selected/source records', 'source': source, 'field': field,
                                 'recorded_value': str(value), 'n': int(n), 'source_total_n': len(frame)})
    hashes = measured.image_sha256.dropna()
    duplicates = {'real_selected_duplicate_byte_groups': int(hashes.value_counts().gt(1).sum()),
                  'synthetic_all_duplicate_byte_groups': int(all_synthetic.image_sha256.dropna().value_counts().gt(1).sum()),
                  'exact_synthetic_to_selected_real_byte_overlap': int(all_synthetic.image_sha256.isin(set(hashes)).sum()),
                  'interpretation': 'encoded-byte audit only; not a privacy/identity test'}
    return pd.DataFrame(rows), pd.DataFrame(metadata), duplicates


def analyze(real, reference, synthetic, definitions, mapping):
    rows, summaries = [], []
    features = definitions.feature.tolist()
    for item in mapping.loc[mapping.mapping_status.eq('reviewed_exact_name')].itertuples(index=False):
        a = real.loc[real.disease.eq(item.gmdb_syndrome)]
        b = synthetic.loc[synthetic.disease.eq(item.pdidb_disease)]
        variants = [('pooled', b)] + [(str(tag), group) for tag, group in b.groupby('model', dropna=False, sort=True)]
        for tag, group in variants:
            if min(len(a), len(group)) < MIN_N:
                continue
            for sensitivity in ['unmodified', 'winsorized_01_99'] if tag == 'pooled' else ['unmodified']:
                cells = []
                for feature in features:
                    x, y, ref = [frame[feature].to_numpy(float) for frame in [a, group, reference]]
                    if sensitivity == 'winsorized_01_99':
                        x, y, ref = [np.clip(v, *np.quantile(v, [.01, .99])) for v in [x, y, ref]]
                    row = dict(syndrome=item.gmdb_syndrome, pdidb_disease=item.pdidb_disease,
                               model_tag=tag, sensitivity=sensitivity, feature=feature)
                    row.update(feature_metrics(x, y, ref))
                    cells.append(row)
                frame = pd.DataFrame(cells)
                summaries.append(dict(syndrome=item.gmdb_syndrome, model_tag=tag, sensitivity=sensitivity,
                                      real_n=len(a), synthetic_n=len(group), reference_n=len(reference), **profile_summary(frame)))
                rows.extend(cells)
    if not rows:
        raise ValueError('No explicitly matched cohorts meet the predeclared descriptive sample support')
    detail = pd.DataFrame(rows).merge(definitions[['feature', 'units', 'formula', 'anatomical_region']],
                                      on='feature', validate='many_to_one')
    return detail, pd.DataFrame(summaries)


def run(config_path, output_dir, aggregate_dir=PUBLIC):
    output_dir = output_dir.resolve()
    if not output_dir.is_relative_to(ROOT / 'local_review') and 'local_only' not in output_dir.relative_to(ROOT).parts:
        raise ValueError('Private intermediates must stay in ignored local_review/ or local_only/')
    if subprocess.run(['git', 'check-ignore', '-q', str(output_dir / 'inputs.json')], cwd=ROOT).returncode:
        raise ValueError('Private output directory must be Git-ignored')
    if aggregate_dir.resolve() != PUBLIC.resolve():
        raise ValueError('Aggregate destination must use the reviewed publication/synthetic_consistency directory')
    output_dir.mkdir(parents=True, exist_ok=True, mode=0o700)
    config, paths, definitions, provenance, real, reference, all_synthetic, synthetic, measured = load_inputs(config_path)
    mapping_path = PUBLIC / 'reviewed_label_mapping.csv'
    mapping = reviewed_mapping(pd.read_csv(mapping_path), set(real.disease), set(all_synthetic.disease))
    audit, metadata, duplicates = audit_sources(mapping, real, all_synthetic, synthetic, measured, definitions.feature.tolist())
    if synthetic.age_range.notna().any():
        raise ValueError('New age metadata requires an explicit age parser and a declared age-matched sensitivity')
    inputs = {name: {'path': str(path.relative_to(ROOT)), 'sha256': sha(path)} for name, path in paths.items()}
    (output_dir / 'private_input_manifest.json').write_text(json.dumps(inputs, indent=2) + '\n')
    detail, fidelity = analyze(real, reference, synthetic, definitions, mapping)
    principal = detail.loc[detail.model_tag.eq('pooled') & detail.sensitivity.eq('unmodified')].copy()
    figures = principal.loc[pd.MultiIndex.from_frame(principal[['syndrome', 'feature']]).isin(ANCHORS)]
    gaps = principal.assign(absolute_standardized_shift=principal.standardized_mean_shift.abs())
    gaps = gaps.sort_values(['absolute_standardized_shift', 'syndrome', 'feature'], ascending=[False, True, True]).head(20)
    clusters = pd.read_csv(ROOT / 'stats_results/extensions_120/feature_clusters.csv')
    clustered = principal.merge(clusters[['feature', 'cluster_id']], on='feature', validate='many_to_one')
    cluster_summary = clustered.assign(absolute_standardized_shift=clustered.standardized_mean_shift.abs()).groupby(
        ['syndrome', 'cluster_id'], sort=True).agg(features=('feature', 'size'),
        median_absolute_standardized_shift=('absolute_standardized_shift', 'median'),
        median_variance_ratio=('variance_ratio', 'median')).reset_index()
    main = fidelity.loc[fidelity.model_tag.eq('pooled') & fidelity.sensitivity.eq('unmodified')]
    sensitivity_rows = []
    for label, subset in [('all_supported', main), ('real_n_at_least_20', main.loc[main.real_n.ge(20)]),
                           ('real_n_at_least_30', main.loc[main.real_n.ge(30)])]:
        status = 'predeclared support' if label == 'all_supported' else 'exploratory cutoff after initial pooled results'
        sensitivity_rows.append(dict(sensitivity=label, selection_status=status, syndromes=len(subset),
                                     median_absolute_standardized_shift=subset.median_absolute_standardized_shift.median(),
                                     median_variance_ratio=subset.median_variance_ratio.median(),
                                     median_direction_agreement=subset.direction_agreement.median(),
                                     median_pearson_g=subset.pearson_g.median()))
    equal_cluster = cluster_summary.groupby('syndrome', sort=True).agg(
        clusters=('cluster_id', 'size'), median_cluster_absolute_shift=('median_absolute_standardized_shift', 'median'),
        median_cluster_variance_ratio=('median_variance_ratio', 'median')).reset_index()
    proxy = []
    for metric in ['median_absolute_standardized_shift', 'median_variance_ratio', 'direction_agreement', 'pearson_g']:
        rho = correlations(main.real_n, main[metric])[1] if len(main) >= 5 else np.nan
        proxy.append({'size_proxy': 'age-restricted available selected GMDB n; NOT training n', 'metric': metric,
                      'syndromes': len(main), 'descriptive_spearman': rho})
    exports = {'cohort_and_mapping_audit': audit, 'metadata_audit': metadata, 'feature_fidelity_summary': detail,
               'syndrome_fidelity_summary': fidelity, 'anchor_distribution_source': figures,
               'effect_scatter_source': principal[['syndrome', 'feature', 'real_n', 'synthetic_n', 'real_g', 'synthetic_g']],
               'largest_gaps_exploratory': gaps, 'cluster_fidelity_summary': cluster_summary,
               'equal_cluster_weight_summary': equal_cluster, 'low_real_n_sensitivity': pd.DataFrame(sensitivity_rows),
               'available_size_proxy_summary': pd.DataFrame(proxy)}
    for name, frame in exports.items():
        save_aggregate(frame, aggregate_dir / (name + '.csv'))
    manifest = {'analysis': 'conditional_PDidb_GMDB_distributional_consistency_120', 'primary_analysis_id': config['analysis_id'],
                'primary_unchanged': {'syndromes': 48, 'hypotheses': 5760, 'discoveries': 955, 'synthetic_rows': 0},
                'facekit_commit': config['facekit_commit'], 'model_sha256': config['model_sha256'],
                'extractor_sha256': provenance['extractor_sha256'], 'pdidb_revision': provenance['dataset_revision'],
                'feature_count': 120, 'seed_reserved': SEED, 'resampling': 'none; descriptive point estimates only',
                'real_selected_n': len(real), 'reference_n': len(reference), 'synthetic_total_n': len(all_synthetic),
                'synthetic_qc_n': len(synthetic), 'supported_syndromes': len(main),
                'supported_real_n': int(main.real_n.sum()), 'supported_synthetic_n': int(main.synthetic_n.sum()),
                'excluded_synthetic_labels': audit.loc[~audit.descriptive_eligible & audit.pdidb_disease.notna(), ['pdidb_disease', 'reason']].to_dict('records'),
                'primary_syndromes_without_synthetic_label': int(audit.mapping_status.eq('excluded_no_synthetic_label').sum()),
                'missingness': {field: int(all_synthetic[field].isna().sum()) for field in ['age_range', 'gender', 'race']},
                'age_matched_sensitivity': 'unavailable; all synthetic ages missing; no imputation',
                'model_tags_qc': synthetic.model.value_counts().to_dict(), 'duplicate_audit': duplicates,
                'generation_origin': {'mentor_statement_date': '2026-11-09', 'statement': 'StyleGAN trained on GMDB; PDIDB sampled its learned distribution',
                                      'verification': 'mentor context; not verified for each row/tag/checkpoint',
                                      'training_membership_and_counts': 'not available for this PDIDB release',
                                      'paper_generator_release_equivalence': 'unverified'},
                'input_sha256': {name: value['sha256'] for name, value in inputs.items()},
                'config_sha256': sha(config_path), 'mapping_sha256': sha(mapping_path),
                'plan_sha256': sha(PUBLIC / 'predeclared_analysis_plan.md'),
                'code_sha256': sha(Path(__file__)),
                'code_base_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                'publication_status': 'aggregate-only candidate; existing PI review applies; no release permission inferred',
                'outputs': {name + '.csv': sha(aggregate_dir / (name + '.csv')) for name in exports}}
    (aggregate_dir / 'synthetic_consistency_manifest.json').write_text(json.dumps(manifest, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: manifest[k] for k in ['feature_count', 'supported_syndromes', 'supported_real_n', 'supported_synthetic_n']}))
    return manifest


if __name__ == '__main__':
    os.umask(0o077)
    def offline(event, args):
        if event in {'socket.connect', 'socket.getaddrinfo', 'http.client.connect', 'urllib.Request'}:
            raise RuntimeError('Patient-data analysis runs entirely offline')
    sys.addaudithook(offline)
    parser = argparse.ArgumentParser(description='Offline conditional synthetic fidelity; never alters primary results')
    parser.add_argument('--active-config', type=Path, default=ROOT / 'ACTIVE_ANALYSIS.json')
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    run(args.active_config, args.output_dir)

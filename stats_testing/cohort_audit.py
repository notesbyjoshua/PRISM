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


def audit(root, output, dry_run=False):
    paths = list((root / 'data').glob('*.csv'))
    paths += list((root / 'stats_testing').glob('*.py'))
    paths += list((root / 'data/dataset_creation').glob('*.py'))
    paths += list((root / 'stats_results').rglob('*.csv'))
    paths += [root / 'FaceKit/src/facekit/core/geometric/extractor.py']
    hashes = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
    disease = pd.read_csv(root / 'data/master_dataset_disease.csv', dtype={'patient_id': 'string'}, low_memory=False)
    healthy = pd.read_csv(root / 'data/master_dataset_healthy.csv', low_memory=False)
    features = list(healthy.loc[:, 'eb_thickness_r':'cheek_area_asym'])
    if len(features) != 125 or any(f not in disease for f in features):
        raise ValueError('Expected the verified 125 shared FaceKit columns')
    tables, stages, labels, feature_rows, ages, conflicts, layers = {}, [], [], [], [], [], []
    mappings = read_mapping_functions(root / 'stats_testing/ancestry_analysis.py')
    interaction_mappings = read_mapping_functions(root / 'stats_testing/ancestry_disease_interaction.py')
    filtered = {}
    for source, frame, image_column in [('GMDB', disease, 'image_id'), ('FairFace', healthy, 'file')]:
        frame = frame.copy()
        frame['_patient'] = clean_id(frame['patient_id']) if 'patient_id' in frame else pd.NA
        frame['_quality'] = quality_status(frame['frontal_ok']) if 'frontal_ok' in frame else 'column_absent'
        tables[f'{source}_quality_flags'] = frame.groupby('_quality', dropna=False).size().rename('rows').reset_index()
        for stage, subset in [('loaded', frame), ('frontal_pass', frame[frame['_quality'].eq('pass')])]:
            summary = patient_summary(subset)
            ids = clean_id(subset[image_column])
            stages.append({'source': source, 'stage': stage, **summary, 'unique_image_keys': ids.nunique(),
                           'missing_image_keys': ids.isna().sum(), 'duplicate_image_rows': ids.dropna().duplicated().sum()})
            for column in ['disease', 'internal_syndrome_name', 'race', 'image_dataset_ethnicity_category',
                           'image_dataset_ethnicity_sub_category', 'gender', 'age', 'age_year', 'age_month', 'age_note']:
                if column in subset:
                    for value, n in subset[column].fillna('<missing>').value_counts(dropna=False).items():
                        labels.append({'source': source, 'stage': stage, 'field': column, 'value': value, 'rows': n})
        frame = frame[frame['_quality'].eq('pass')].copy()
        filtered[source] = frame
        if source == 'GMDB':
            frame = pd.concat([frame, age_diagnostics(frame)], axis=1)
            filtered[source] = frame
            for column in ['internal_syndrome_name', 'image_dataset_ethnicity_category', 'image_dataset_ethnicity_sub_category', 'gender']:
                counts = frame.groupby('_patient')[column].nunique()
                conflicts.append({'field': column, 'patients_with_conflicting_nonmissing_values': counts.gt(1).sum()})
            grouped = frame.groupby('_patient')
            repeats = grouped.agg(rows=('image_id', 'size'), observed_ages=('conditional_age_years', 'count'),
                                  min_conditional_age=('conditional_age_years', 'min'),
                                  max_conditional_age=('conditional_age_years', 'max'))
            repeats['conditional_span_years'] = repeats.max_conditional_age - repeats.min_conditional_age
            tables['patient_age_spans_LOCAL_ONLY'] = repeats[repeats.rows.gt(1)].reset_index()
            tables['image_patient_manifest_LOCAL_ONLY'] = frame[['image_id', '_patient', 'internal_syndrome_name',
                                                                 'age_status', 'conditional_age_years']].copy()
            for prefix in ['', 'image_dataset_']:
                if prefix:
                    for column in ['patient_id', 'age_year', 'age_month', 'gender', 'internal_syndrome_name']:
                        a, b = frame[column].astype('string'), frame[prefix + column].astype('string')
                        conflicts.append({'field': f'{column}_versus_{prefix}{column}',
                                          'patients_with_conflicting_nonmissing_values': pd.NA,
                                          'discordant_rows': int((a.fillna('<missing>') != b.fillna('<missing>')).sum())})
        for layer in [1, 2]:
            if source == 'FairFace':
                frame[f'layer{layer}'] = frame.race.map(mappings[f'map_fairface_layer{layer}'])
            elif layer == 1:
                frame['layer1'] = frame.image_dataset_ethnicity_category.map(mappings['map_gmdb_layer1'])
            else:
                pairs = zip(frame.image_dataset_ethnicity_category, frame.image_dataset_ethnicity_sub_category)
                frame['layer2'] = [mappings['map_gmdb_layer2'](a, b) for a, b in pairs]
                pairs = zip(frame.image_dataset_ethnicity_category, frame.image_dataset_ethnicity_sub_category)
                other = pd.Series([interaction_mappings['map_gmdb_layer2'](a, b) for a, b in pairs], index=frame.index)
                tables['mapping_discrepancies'] = pd.DataFrame([{'rows_differing_A_B_vs_C': int((frame.layer2.fillna('<excluded>') != other.fillna('<excluded>')).sum())}])
            for group, subset in frame.groupby(f'layer{layer}', dropna=False):
                layers.append({'source': source, 'layer': layer, 'group': group, **patient_summary(subset),
                               'unique_image_keys': subset[image_column].nunique()})
        for grouping in [None, 'layer1', 'layer2', 'internal_syndrome_name']:
            if grouping and grouping not in frame:
                continue
            subsets = [('all', frame)] if grouping is None else frame.groupby(grouping, dropna=False)
            for group, subset in subsets:
                f = feature_counts(subset, features)
                f.insert(0, 'group', group)
                f.insert(0, 'grouping', grouping or 'overall')
                f.insert(0, 'source', source)
                feature_rows.append(f)
                for column in ['conditional_age_band', 'age_status', 'age', 'gender']:
                    if column in subset:
                        for value, n in subset[column].fillna('<missing>').value_counts().items():
                            ages.append({'source': source, 'grouping': grouping or 'overall', 'group': group,
                                         'field': column, 'value': value, 'rows': n})
    tables['cohort_flow'] = pd.DataFrame(stages)
    tables['label_counts'] = pd.DataFrame(labels)
    tables['layer_counts'] = pd.DataFrame(layers)
    tables['feature_counts'] = pd.concat(feature_rows, ignore_index=True)
    tables['age_sex_distributions'] = pd.DataFrame(ages)
    tables['metadata_consistency'] = pd.DataFrame(conflicts)
    config = pd.read_csv(root / 'stats_results/significance_testing/analysis_configuration.csv').iloc[0]
    min_n, max_missing = int(config.min_group_n), float(config.max_missing_fraction)
    hcounts = feature_counts(filtered['FairFace'], features).set_index('feature')
    eligibility, syndrome_rows = [], []
    for syndrome, subset in filtered['GMDB'].groupby('disease', dropna=False):
        counts = feature_counts(subset, features).set_index('feature')
        n_eligible = 0
        for feature in features:
            d, h = counts.loc[feature], hcounts.loc[feature]
            reasons = []
            if d.valid_n < min_n:
                reasons.append('disease_valid_n_below_saved_min')
            if h.valid_n < min_n:
                reasons.append('reference_valid_n_below_saved_min')
            if d.invalid_fraction > max_missing:
                reasons.append('disease_missing_fraction_above_saved_max')
            if h.invalid_fraction > max_missing:
                reasons.append('reference_missing_fraction_above_saved_max')
            n_eligible += not reasons
            eligibility.append({'syndrome': syndrome, 'feature': feature, 'image_n_disease': int(d.valid_n),
                                'image_n_reference': int(h.valid_n), 'eligible_legacy_image_rule': not reasons,
                                'reason': ';'.join(reasons) or 'eligible_image_rule_only'})
        syndrome_rows.append({'syndrome': syndrome, **patient_summary(subset), 'eligible_features': n_eligible,
                              'reason': 'eligible_image_rule_only' if n_eligible else 'no_feature_passes_saved_rules'})
    tables['feature_eligibility'] = pd.DataFrame(eligibility)
    tables['syndrome_eligibility'] = pd.DataFrame(syndrome_rows)
    saved = pd.read_csv(root / 'stats_results/significance_testing/stats_each_disease_vs_healthy.csv')
    sig = saved.mannwhitney_q_fdr_global.lt(.05)
    actual = {'disease_images': len(filtered['GMDB']), 'reference_rows': len(filtered['FairFace']),
              'syndromes': filtered['GMDB'].disease.nunique(), 'tested_syndromes': saved.disease.nunique(),
              'comparisons': len(saved), 'significant': int(sig.sum()),
              'significant_abs_g_ge_0_8': int((sig & saved.hedges_g.abs().ge(.8)).sum())}
    historical = dict(zip(actual, [9411, 1000, 579, 356, 44500, 9002, 4028]))
    reconciliation = [{'metric': k, 'handoff': historical[k], 'current': v, 'difference': v - historical[k]} for k, v in actual.items()]
    result_rows = []
    for layer in ['layer1_broad_ethnicity', 'layer2_subcategory']:
        folder = root / 'stats_results/ethnicity_analysis' / layer
        row = {'layer': layer}
        for name in ['analysis_A_healthy_kruskal_all_groups', 'analysis_A_healthy_pairwise_groups',
                     'analysis_B_matched_disease_vs_healthy', 'analysis_B_matched_significant',
                     'analysis_C_all_interactions', 'analysis_C_high_confidence_interactions']:
            row[name] = len(pd.read_csv(folder / f'{name}.csv'))
        result_rows.append(row)
        for field, target in [('analysis_B_matched_significant', 8338 if layer.startswith('layer1') else 7844),
                              ('analysis_C_high_confidence_interactions', 397 if layer.startswith('layer1') else 348)]:
            reconciliation.append({'metric': f'{layer}_{field}', 'handoff': target, 'current': row[field],
                                   'difference': row[field] - target})
    tables['saved_layer_results'] = pd.DataFrame(result_rows)
    tables['legacy_reconciliation'] = pd.DataFrame(reconciliation)
    saved_pairs = set(zip(saved.disease, saved.feature))
    audited_pairs = {(r['syndrome'], r['feature']) for r in eligibility if r['eligible_legacy_image_rule']}
    tables['eligibility_reconciliation'] = pd.DataFrame([{'saved_only': len(saved_pairs - audited_pairs),
                                                          'audit_only': len(audited_pairs - saved_pairs)}])
    provenance = {'created_utc': datetime.now(timezone.utc).isoformat(), 'python': __import__('sys').version,
                  'numpy': np.__version__, 'pandas': pd.__version__, 'input_sha256': hashes,
                  'git_head': subprocess.check_output(['git', '-C', str(root), 'rev-parse', 'HEAD'], text=True).strip(),
                  'facekit_head': subprocess.check_output(['git', '-C', str(root / 'FaceKit'), 'rev-parse', 'HEAD'], text=True).strip(),
                  'policies': {'analysis_unit': 'descriptive image audit; no patient selection',
                               'age': 'conditional year + month/12 only; zero-zero unresolved; months outside 0..11 flagged',
                               'eligibility': 'saved image minimum and missingness; finite values required',
                               'min_n': min_n, 'max_missing': max_missing,
                               'reference_identity': 'file path identifies image, not participant',
                               'mapping': 'existing A/B function definitions; no analysis script import'}}
    if dry_run:
        return tables, provenance
    output.mkdir(parents=True, exist_ok=False)
    for name, table in tables.items():
        table.to_csv(output / f'{name}.csv', index=False)
    (output / 'provenance.json').write_text(json.dumps(provenance, indent=2))
    changed = [name for name, digest in hashes.items() if hashlib.sha256((root / name).read_bytes()).hexdigest() != digest]
    if changed:
        raise RuntimeError(f'Inputs changed during audit: {changed}')
    return tables, provenance


def main():
    parser = argparse.ArgumentParser(description='Read-only cohort audit; outputs can contain local patient identifiers.')
    parser.add_argument('--output', type=Path, default=ROOT / 'local_review' / datetime.now().strftime('audit_%Y%m%d_%H%M%S'))
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    tables, _ = audit(ROOT, args.output, args.dry_run)
    print(tables['cohort_flow'].to_string(index=False))
    print(tables['legacy_reconciliation'].to_string(index=False))
    print('Dry run; no outputs written' if args.dry_run else f'Wrote audit to {args.output.resolve()}')


if __name__ == '__main__':
    main()

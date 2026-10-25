#!/usr/bin/env python3
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import pandas as pd
import scipy
from scipy.stats import mannwhitneyu, rankdata, norm, spearmanr, pearsonr
import statsmodels
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
FROZEN = ROOT / 'revision_results/final_20260914'
sys.path.insert(0, str(ROOT / 'stats_testing'))
from cohort_audit import age_diagnostics, clean_id, quality_status, read_mapping_functions
from prepare_cohorts import select_patients, deduplicate_reference

SEED = 20261024
LAYERS = {'broad': 'layer1_broad_ethnicity', 'detailed': 'layer2_subcategory'}
HIGHLIGHTS = [('Cornelia de Lange syndrome__fb4a389d', 'philtrum_length'), ('Noonan syndrome__c2df3f82', 'inter_pupillary_distance'),
              ('WILLIAMS-BEUREN SYNDROME; WBS__93e4a782', 'eb_thickness_mean'), ('Rubinstein-Taybi syndrome__ad6c2718', 'malar_bulge_mean')]
INPUTS = set()


def check(condition, message):
    if not condition:
        raise ValueError(message)


def read(path, **kwargs):
    path = Path(path)
    INPUTS.add(path)
    return pd.read_csv(path, **kwargs)


def write(frame, name):
    frame.to_csv(OUT / name, index=False)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load():
    features = read(FROZEN / 'exact_feature_names.csv').feature.tolist()
    check(len(features) == len(set(features)) == 125, 'Historical schema mismatch')
    disease = read(ROOT / 'data/master_dataset_disease.csv', dtype={'patient_id': 'string', 'image_id': 'string'}, low_memory=False)
    reference = read(ROOT / 'data/master_dataset_healthy.csv', low_memory=False)
    primary, _ = select_patients(disease, '3-9', 'pose')
    alternate, _ = select_patients(disease, '3-9', 'image_id')
    check(quality_status(reference.frontal_ok).eq('pass').all(), 'Reference stored frontal flags changed')
    reference_unique, _ = deduplicate_reference(reference)
    for frame in [primary, alternate, reference_unique]:
        check(frame.loc[:, features[0]:features[-1]].columns.tolist() == features, 'Feature order mismatch')
        frame[features] = frame[features].apply(pd.to_numeric, errors='coerce').replace([np.inf, -np.inf], np.nan)
    check(primary.patient_id.notna().all() and primary.patient_id.is_unique, 'Recorded patients not unique')
    check(primary.patient_id.tolist() == alternate.patient_id.tolist(), 'Patient universes differ')
    frozen = read(FROZEN / 'age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz')
    targets = sorted(frozen.disease.unique())
    check(len(targets) == 49 and len(frozen) == 6125, 'Frozen test family mismatch')
    check(set(targets).issubset(primary.disease), 'Frozen target labels absent')
    primary, reference_unique = primary.copy(), reference_unique.copy()
    maps = read_mapping_functions(ROOT / 'stats_testing/ancestry_analysis.py')
    for frame in [primary, reference_unique]:
        frame['sex_label'] = frame.gender.astype('string').str.strip().str.lower().map(
            {'male': 'Male', 'm': 'Male', 'female': 'Female', 'f': 'Female'})
    primary['broad_group'] = primary.image_dataset_ethnicity_category.map(maps['map_gmdb_layer1'])
    primary['detailed_group'] = [maps['map_gmdb_layer2'](a, b) for a, b in zip(
        primary.image_dataset_ethnicity_category, primary.image_dataset_ethnicity_sub_category)]
    for layer in LAYERS:
        reference_unique[f'{layer}_group'] = reference_unique.race.map(maps[f'map_fairface_layer{1 if layer == "broad" else 2}'])
    return disease, reference, primary, alternate, reference_unique, frozen, targets, features


def audit(disease, reference, primary, alternate, ref, frozen, targets, features):
    metadata = []
    observed = {'patient_id', 'image_id', 'internal_syndrome_name', 'disease', 'age_year', 'age_month', 'gender',
                'ethnicity_category', 'ethnicity_sub_category', 'image_dataset_ethnicity_category',
                'image_dataset_ethnicity_sub_category', 'gene_names', 'hgvs', 'test_type', 'hpo_terms',
                'present_features', 'absent_features', 'patient_name', 'frontal_ok', 'pose_yaw', 'pose_pitch', 'pose_roll',
                'file', 'age', 'race', 'broad_group', 'detailed_group', 'conditional_age_years', 'sex_label'}
    absent = {'image_width', 'image_height', 'resolution', 'expression', 'occlusion', 'gaze_annotation',
              'landmarks', 'genetic_confirmation'}
    for source, frame in [('GMDB_original', disease), ('GMDB_selected', primary),
                          ('FairFace_original', reference), ('FairFace_deduplicated', ref)]:
        for field in sorted(observed | absent | set(features)):
            exists = field in frame
            if exists:
                values = frame[field].astype('string').str.strip()
                available = values.notna() & ~values.str.lower().isin(['', 'nan', 'none', 'null', 'unknown', 'na', 'n/a', '[]', '{}'])
                n = int(available.sum())
            else:
                n = 0
            provenance = 'source dataset annotation; not independently verified'
            if field in features:
                provenance = 'historical FaceKit derived normalized geometry; not a clinical annotation'
            elif field.startswith('pose_') or field == 'frontal_ok':
                provenance = 'stored detector-derived estimate/quality flag; execution settings not established per image'
            elif field in ['gene_names', 'hgvs', 'test_type', 'hpo_terms', 'present_features', 'absent_features']:
                provenance = 'supplied metadata; populated field does not establish confirmation or clinical validity'
            elif field in ['broad_group', 'detailed_group', 'sex_label', 'conditional_age_years']:
                provenance = 'deterministically derived from recorded annotations using frozen rules; not independent measurements'
            elif not exists:
                provenance = 'not a column in supplied master; not measured here'
            metadata.append({'source': source, 'field': field, 'present': exists, 'n_rows': len(frame),
                             'available_n': n, 'missing_or_unknown_n': len(frame)-n, 'provenance': provenance})
    write(pd.DataFrame(metadata), 'metadata_availability.csv')
    inventory = []
    tracked = subprocess.check_output(['git', 'ls-files', 'data/*.csv'], cwd=ROOT, text=True).splitlines()
    for name in tracked:
        header = read(ROOT / name, nrows=0).columns.tolist()
        inventory.append({'tracked_file': name, 'columns': header, 'potential_identifier_columns':
                          [c for c in header if any(x in c for x in ['patient', 'image_id', 'hgvs', 'hpo', 'relative_path'])]})
    (OUT / 'local_only/metadata_inventory.json').write_text(json.dumps(inventory, indent=2))
    frontal = quality_status(disease.frontal_ok).eq('pass')
    ages = age_diagnostics(disease).conditional_age_years
    eligible = frontal & ages.ge(3) & ages.lt(10)
    flow = [('GMDB', 'original_image_records', len(disease)), ('GMDB', 'stored_frontal_pass', int(frontal.sum())),
            ('GMDB', 'age_eligible_images', int(eligible.sum())), ('GMDB', 'selected_recorded_patients', len(primary)),
            ('GMDB', 'missing_original_recorded_patient_id', int(clean_id(disease.patient_id).isna().sum())),
            ('GMDB', 'original_distinct_recorded_patients', clean_id(disease.patient_id).nunique()),
            ('GMDB', 'additional_eligible_images_removed', int(eligible.sum())-len(primary)),
            ('GMDB', 'retained_syndrome_labels', primary.disease.nunique()), ('GMDB', 'main_eligible_syndromes', len(targets)),
            ('FairFace', 'original_image_records', len(reference)), ('FairFace', 'distinct_image_paths', len(ref))]
    expected = [10942, 9411, 2101, 1850, 0, 8328, 251, 379, 49, 1000, 942]
    flow = pd.DataFrame(flow, columns=['source', 'stage', 'count'])
    flow['frozen_target'] = expected
    flow['passed'] = flow['count'].eq(flow.frozen_target)
    write(flow, 'cohort_flow.csv')
    check(flow.passed.all(), 'Cohort flow failed')
    counts = primary.groupby('disease').size().rename('n_selected_patients').reset_index().rename(columns={'disease': 'syndrome'})
    counts['main_eligible'] = counts.syndrome.isin(targets)
    counts['target_tests'] = np.where(counts.main_eligible, 125, 0)
    write(counts.loc[counts.main_eligible], 'syndrome_counts.csv')
    write(counts, 'all_syndrome_counts.csv')
    rows = []
    for source, frame in [('GMDB', primary), ('FairFace', ref)]:
        for field in ['broad_group', 'detailed_group', 'sex_label']:
            for group, count in frame[field].fillna('<unmapped>').value_counts().items():
                rows.append({'source': source, 'variable': field, 'group': group, 'n': count,
                             'denominator': len(frame), 'unit': 'recorded patient' if source == 'GMDB' else 'image path'})
    write(pd.DataFrame(rows), 'group_denominators.csv')
    selected_ids = primary.set_index('patient_id').image_id
    alt_ids = alternate.set_index('patient_id').image_id
    changed = selected_ids.ne(alt_ids)
    check(changed.sum() == 94, 'Changed photograph baseline failed')
    by_syndrome = primary.assign(photo_changed=changed.to_numpy()).groupby('disease').agg(
        selected_patients=('patient_id', 'size'), patients_photo_changed=('photo_changed', 'sum')).reset_index()
    write(by_syndrome.rename(columns={'disease': 'syndrome'}), 'photo_change_by_syndrome.csv')
    rows = []
    for group, frame in [('GMDB_original', disease), ('GMDB_primary', primary), ('FairFace_original', reference), ('FairFace_unique', ref)]:
        if group == 'GMDB_original' or 'conditional_age_years' in frame:
            vals = age_diagnostics(frame).conditional_age_years if group == 'GMDB_original' else frame.conditional_age_years
            rows.append({'source': group, 'age_type': 'recorded year+month/12', 'n': len(frame), 'available': vals.notna().sum(),
                         'minimum': vals.min(), 'median': vals.median(), 'maximum': vals.max()})
        elif group.startswith('FairFace'):
            rows.append({'source': group, 'age_type': 'category 3-9 only', 'n': len(frame), 'available': frame.age.notna().sum(),
                         'minimum': np.nan, 'median': np.nan, 'maximum': np.nan})
    write(pd.DataFrame(rows), 'age_availability.csv')
    dictionary = read(ROOT / 'publication/supplementary/Supplementary_Table_S6_exact_feature_definitions.csv')
    check(dictionary.feature.tolist() == features and dictionary.extraction_status.eq('EXACT_FROM_SOURCE').all(), 'Dictionary mismatch')
    write(dictionary, 'feature_dictionary.csv')
    formula_rows = []
    for syndrome, feature in HIGHLIGHTS:
        row = dictionary.loc[dictionary.feature.eq(feature)].iloc[0]
        formula_rows.append({'syndrome': syndrome, 'feature': feature, 'numerator': row.numerator_definition,
                             'denominator': row.denominator_definition, 'measurement_type': row.measurement_type,
                             'landmarks': row.landmark_indices_or_named_points, 'status': 'not identifiable from normalized columns alone'})
    write(pd.DataFrame(formula_rows), 'highlighted_components.csv')
    print('Phase 0: all 11 baseline cohort checks passed; local metadata inventory written', flush=True)


def effect(x, y):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    nx, ny = np.isfinite(x).sum(axis=0), np.isfinite(y).sum(axis=0)
    pooled = ((nx-1)*np.nanvar(x, axis=0, ddof=1)+(ny-1)*np.nanvar(y, axis=0, ddof=1))/(nx+ny-2)
    with np.errstate(invalid='ignore', divide='ignore'):
        return (1-3/(4*(nx+ny)-9))*(np.nanmean(x, axis=0)-np.nanmean(y, axis=0))/np.sqrt(pooled)


def bootstrap_moments(values, iterations, rng):
    n = len(values)
    indices = rng.integers(n, size=(iterations, n))
    weights = np.zeros((iterations, n))
    for i, sample in enumerate(indices):
        weights[i] = np.bincount(sample, minlength=n)
    finite = np.isfinite(values)
    clean = np.where(finite, values, 0)
    count = weights @ finite.astype(float)
    total = weights @ clean
    total2 = weights @ (clean*clean)
    with np.errstate(invalid='ignore', divide='ignore'):
        mean = total/count
        variance = (total2-total*total/count)/(count-1)
    return count, mean, np.maximum(variance, 0)


def bootstrap_effect(x, y, iterations, seed):
    rng = np.random.default_rng(seed)
    nx, mx, vx = bootstrap_moments(x, iterations, rng)
    ny, my, vy = bootstrap_moments(y, iterations, rng)
    with np.errstate(invalid='ignore', divide='ignore'):
        pooled = ((nx-1)*vx+(ny-1)*vy)/(nx+ny-2)
        values = (1-3/(4*(nx+ny)-9))*(mx-my)/np.sqrt(pooled)
    values[(nx < 2) | (ny < 2) | ~np.isfinite(values)] = np.nan
    return values


def eligible_feature_counts(nx, ny, total_x, total_y):
    return (nx >= 10) & (ny >= 10) & (10*(total_x-nx) <= 3*total_x) & (10*(total_y-ny) <= 3*total_y)


def contrast(frame, targets, features, comparator='all', iterations=0):
    rows = []
    for index, syndrome in enumerate(targets):
        target, other = comparator_frames(frame, syndrome, targets, comparator)
        x, y = target[features].to_numpy(float), other[features].to_numpy(float)
        nx, ny = np.isfinite(x).sum(axis=0), np.isfinite(y).sum(axis=0)
        missing_x, missing_y = (len(x)-nx)/len(x), (len(y)-ny)/len(y)
        estimable = eligible_feature_counts(nx, ny, len(x), len(y))
        g = effect(x, y)
        estimable &= np.isfinite(g)
        test = mannwhitneyu(x, y, axis=0, alternative='two-sided', method='asymptotic', use_continuity=True, nan_policy='omit')
        estimates = bootstrap_effect(x, y, iterations, [SEED, index]) if iterations else None
        for j, feature in enumerate(features):
            low, high = np.nanquantile(estimates[:, j], [.025, .975]) if iterations else [np.nan, np.nan]
            rows.append({'syndrome': syndrome, 'feature': feature, 'comparator': comparator,
                         'target_total': len(x), 'comparator_total': len(y), 'n_target': nx[j], 'n_comparator': ny[j],
                         'target_missing_fraction': missing_x[j], 'comparator_missing_fraction': missing_y[j],
                         'estimable': estimable[j] and np.isfinite(g[j]), 'hedges_g': g[j] if estimable[j] else np.nan,
                         'ci_95_lower': low, 'ci_95_upper': high, 'bootstrap_iterations': iterations,
                         'bootstrap_finite': int(np.isfinite(estimates[:, j]).sum()) if iterations else 0,
                         'mannwhitney_u': test.statistic[j], 'p_value': test.pvalue[j] if estimable[j] else 1.0})
        if (index+1) % 10 == 0 or index == len(targets)-1:
            print(f'GMDB {comparator}: {index+1}/{len(targets)} syndrome targets computed', flush=True)
    result = pd.DataFrame(rows)
    result['q_global_bh'] = bh(result.p_value)
    result['significant'] = result.q_global_bh.lt(.05) & result.estimable
    check(len(result) == len(targets)*len(features) == 6125, 'BH family size wrong')
    return result


def comparator_frames(frame, syndrome, targets, comparator):
    target = frame.loc[frame.disease.eq(syndrome)]
    other = frame.loc[~frame.disease.eq(syndrome)]
    if comparator == 'eligible49':
        other = other.loc[other.disease.isin(targets)]
    check(not other.disease.eq(syndrome).any(), 'Target in comparator')
    check(not set(target.patient_id) & set(other.patient_id), 'Patient overlap across groups')
    return target, other


def bh(values):
    values = np.asarray(values, dtype=float)
    check(np.isfinite(values).all() and ((values >= 0) & (values <= 1)).all(), 'Invalid p-values')
    return multipletests(values, method='fdr_bh')[1] if len(values) else np.array([])


def comparison(result, frozen):
    base = frozen[['disease', 'feature', 'n_disease', 'n_healthy', 'hedges_g', 'mannwhitney_q_fdr_global']].rename(
        columns={'disease': 'syndrome', 'hedges_g': 'fairface_g', 'mannwhitney_q_fdr_global': 'fairface_q'})
    joined = result.merge(base, on=['syndrome', 'feature'], validate='one_to_one')
    joined['fairface_significant'] = joined.fairface_q.lt(.05)
    joined['same_direction'] = np.sign(joined.hedges_g).eq(np.sign(joined.fairface_g))
    joined['discovery_category'] = np.select([joined.significant & joined.fairface_significant,
        joined.significant, joined.fairface_significant], ['both', 'within_GMDB_only', 'FairFace_only'], default='neither')
    write(joined, 'within_gmdb_vs_fairface_comparison.csv')
    highlights = pd.concat([joined.loc[joined.syndrome.eq(s) & joined.feature.eq(f)] for s, f in HIGHLIGHTS])
    check(len(highlights) == 4, 'Missing highlighted comparison')
    write(highlights, 'highlighted_examples.csv')
    rows = []
    for name, values in [('within_GMDB', joined.hedges_g), ('FairFace', joined.fairface_g)]:
        temporary = joined.assign(significant=joined.significant if name == 'within_GMDB' else joined.fairface_significant,
                                  effect=values)
        for feature, frame in temporary.groupby('feature'):
            rows.append({'comparison': name, 'feature': feature, 'syndromes_tested': len(frame),
                         'positive_significant': int((frame.significant & frame.effect.gt(0)).sum()),
                         'negative_significant': int((frame.significant & frame.effect.lt(0)).sum()),
                         'nonsignificant': int((~frame.significant).sum())})
    write(pd.DataFrame(rows), 'feature_direction_summary.csv')
    return joined


def adjusted(primary, features):
    rows, overlap = [], []
    for syndrome, feature in HIGHLIGHTS:
        frame = primary.copy()
        frame['target'] = frame.disease.eq(syndrome).astype(float)
        frame['female'] = frame.sex_label.map({'Male': 0., 'Female': 1.})
        frame['age_centered'] = frame.conditional_age_years-6
        frame = frame.dropna(subset=[feature, 'female', 'age_centered'])
        for target, subset in frame.groupby('target'):
            overlap.append({'syndrome': syndrome, 'feature': feature, 'group': 'target' if target else 'other GMDB',
                            'n': len(subset), 'age_min': subset.conditional_age_years.min(),
                            'age_q25': subset.conditional_age_years.quantile(.25), 'age_median': subset.conditional_age_years.median(),
                            'age_q75': subset.conditional_age_years.quantile(.75), 'age_max': subset.conditional_age_years.max(),
                            'female_n': int(subset.female.eq(1).sum()), 'male_n': int(subset.female.eq(0).sum())})
        design = sm.add_constant(frame[['target', 'age_centered', 'female']], has_constant='add')
        check(np.linalg.matrix_rank(design) == 4, 'Adjusted model rank deficient')
        model = sm.OLS(frame[feature], design).fit(cov_type='HC3', use_t=False)
        interval = model.conf_int().loc['target']
        leverage = model.get_influence().hat_matrix_diag
        groups = frame.groupby('target').conditional_age_years.agg(['min', 'max'])
        rows.append({'syndrome': syndrome, 'feature': feature, 'n_complete': len(frame),
                     'n_target': int(frame.target.sum()), 'n_other': int(frame.target.eq(0).sum()),
                     'excluded': len(primary)-len(frame), 'target_beta': model.params['target'],
                     'hc3_se': model.bse['target'], 'ci_95_lower': interval.iloc[0], 'ci_95_upper': interval.iloc[1],
                     'p_value': model.pvalues['target'], 'maximum_leverage': leverage.max(), 'design_rank': 4,
                     'age_overlap_low': max(groups['min']), 'age_overlap_high': min(groups['max']),
                     'condition_number': model.condition_number})
    result = pd.DataFrame(rows)
    result['q_four_models'] = bh(result.p_value)
    write(result, 'within_gmdb_age_sex_adjusted.csv')
    write(pd.DataFrame(overlap), 'adjusted_model_overlap.csv')


def photo_sensitivity(primary):
    p = read(FROZEN / 'age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz')
    a = read(FROZEN / 'age_3_9_id_sensitivity/significance_testing/stats_each_disease_vs_healthy.csv.gz')
    columns = ['disease', 'feature', 'hedges_g', 'mannwhitney_q_fdr_global']
    joined = p[columns].merge(a[columns], on=['disease', 'feature'], suffixes=('_primary', '_alternative'), validate='one_to_one')
    joined['signed_difference'] = joined.hedges_g_alternative-joined.hedges_g_primary
    joined['absolute_difference'] = joined.signed_difference.abs()
    joined['same_direction'] = np.sign(joined.hedges_g_primary).eq(np.sign(joined.hedges_g_alternative))
    joined['primary_significant'] = joined.mannwhitney_q_fdr_global_primary.lt(.05)
    joined['alternative_significant'] = joined.mannwhitney_q_fdr_global_alternative.lt(.05)
    joined['discordant_call'] = joined.primary_significant.ne(joined.alternative_significant)
    joined['near_q_005'] = joined.mannwhitney_q_fdr_global_primary.between(.025, .10) | joined.mannwhitney_q_fdr_global_alternative.between(.025, .10)
    check(joined.primary_significant.sum() == 811 and joined.alternative_significant.sum() == 834 and
          joined.discordant_call.sum() == 105 and (joined.primary_significant & joined.alternative_significant).sum() == 770,
          'Frozen photo overlap failed')
    write(joined.rename(columns={'disease': 'syndrome'}), 'sensitivity_effect_sizes.csv')
    rows = []
    for name, frame in [('all', joined), ('near_q_005', joined.loc[joined.near_q_005])]:
        rows.append({'subset': name, 'n_hypotheses': len(frame), 'median_absolute_difference': frame.absolute_difference.median(),
                     'q90_absolute_difference': frame.absolute_difference.quantile(.9), 'q95_absolute_difference': frame.absolute_difference.quantile(.95),
                     'q99_absolute_difference': frame.absolute_difference.quantile(.99), 'max_absolute_difference': frame.absolute_difference.max(),
                     'mean_signed_difference': frame.signed_difference.mean(), 'pearson_r': frame.hedges_g_primary.corr(frame.hedges_g_alternative),
                     'spearman_rho': frame.hedges_g_primary.corr(frame.hedges_g_alternative, method='spearman'),
                     'sign_agreement_fraction': frame.same_direction.mean(), 'discordant_calls': int(frame.discordant_call.sum())})
    write(pd.DataFrame(rows), 'sensitivity_effect_size_summary.csv')
    for grouping in ['disease', 'feature']:
        grouped = joined.groupby(grouping).agg(n=('feature', 'size'), median_absolute_difference=('absolute_difference', 'median'),
            maximum_absolute_difference=('absolute_difference', 'max'), discordant_calls=('discordant_call', 'sum'),
            sign_agreement=('same_direction', 'mean')).reset_index().rename(columns={'disease': 'syndrome'})
        write(grouped, f'photo_sensitivity_by_{"syndrome" if grouping == "disease" else "feature"}.csv')
    return joined


def restricted_families():
    summaries, details = [], []
    for layer, directory in LAYERS.items():
        for analysis, filename, p_col in [('B', 'analysis_B_matched_disease_vs_healthy', 'p_value'),
                                          ('C', 'analysis_C_all_interactions', 'interaction_p_value')]:
            data = read(FROZEN / f'age_3_9/ethnicity_analysis/{directory}/{filename}.csv.gz')
            for threshold in [2, 5, 10, 20]:
                if analysis == 'B':
                    mask = data.n_disease.ge(threshold) & data.n_healthy_matched.ge(20)
                else:
                    mask = data.n_disease_group_1.ge(threshold) & data.n_disease_group_2.ge(threshold)
                    mask &= data.n_healthy_group_1.ge(20) & data.n_healthy_group_2.ge(20)
                frame = data.loc[mask].copy()
                frame['q_restricted_family'] = bh(frame[p_col])
                old = frame.global_fdr_q_value.lt(.05)
                new = frame.q_restricted_family.lt(.05)
                if threshold == 2:
                    check(np.allclose(frame.q_restricted_family, frame.global_fdr_q_value), 'Original B/C BH not reproduced')
                low = 'hedges_g_ci_lower' if analysis == 'B' else 'interaction_ci_lower'
                high = 'hedges_g_ci_upper' if analysis == 'B' else 'interaction_ci_upper'
                effect_col = 'hedges_g' if analysis == 'B' else 'interaction_beta_group2_minus_group1'
                summaries.append({'analysis': analysis, 'layer': layer, 'minimum_disease_n': threshold, 'minimum_reference_n': 20,
                                  'hypotheses': len(frame), 'significant_recorrected': int(new.sum()),
                                  'significant_original_q_filtered': int(old.sum()), 'syndromes': frame.syndrome_name.nunique(),
                                  'median_interval_width': (frame[high]-frame[low]).median(),
                                  'median_absolute_effect': frame[effect_col].abs().median(),
                                  'effect_change_from_refit': 'not refit; each retained estimate/CI unchanged'})
                keys = ['syndrome_name', 'feature', 'group'] if analysis == 'B' else ['syndrome_name', 'feature', 'group_1', 'group_2']
                extract = frame[keys+[p_col, 'global_fdr_q_value', 'q_restricted_family', effect_col, low, high]].copy()
                extract.insert(0, 'threshold', threshold)
                extract.insert(0, 'layer', layer)
                extract.insert(0, 'analysis', analysis)
                details.append(extract)
    write(pd.DataFrame(summaries), 'subgroup_threshold_sensitivity.csv')
    write(pd.concat(details, ignore_index=True), 'subgroup_restricted_hypotheses.csv')
    print('Phase 2B: all four original global BH families reproduced; thresholds re-corrected', flush=True)


def random_photos(raw, ref, targets, features, replicates):
    ages = age_diagnostics(raw).conditional_age_years
    eligible = raw.loc[quality_status(raw.frontal_ok).eq('pass') & ages.ge(3) & ages.lt(10)].copy()
    eligible['patient_id'] = clean_id(eligible.patient_id)
    check(eligible.patient_id.notna().all(), 'Random selection includes missing patient IDs')
    groups = list(eligible.groupby('patient_id', sort=True).indices.values())
    check(len(groups) == 1850, 'Random selection patient count')
    values = eligible[features].apply(pd.to_numeric, errors='coerce').to_numpy(float)
    reference = ref[features].to_numpy(float)
    labels = eligible.disease.to_numpy()
    rng = np.random.default_rng(SEED+1)
    effects, qvalues, records = [], [], []
    for replicate in range(replicates):
        selected = np.array([indices[rng.integers(len(indices))] for indices in groups])
        chosen, diseases = values[selected], labels[selected]
        gs, ps = [], []
        for syndrome in targets:
            x = chosen[diseases == syndrome]
            finite_n = np.isfinite(x).sum(axis=0)
            check(eligible_feature_counts(finite_n, np.isfinite(reference).sum(axis=0), len(x), len(reference)).all(),
                  'Random replicate feature eligibility changed')
            gs.extend(effect(x, reference))
            ps.extend(mannwhitneyu(x, reference, axis=0, method='asymptotic', alternative='two-sided', nan_policy='omit').pvalue)
        q = bh(ps)
        effects.append(gs)
        qvalues.append(q)
        records.append({'replicate': replicate+1, 'hypotheses': len(q), 'significant': int((q < .05).sum())})
        if (replicate+1) % 20 == 0:
            print(f'Random photograph selection: {replicate+1}/{replicates} replicates', flush=True)
    gs, qs = np.array(effects), np.array(qvalues)
    rows = []
    for index, (syndrome, feature) in enumerate((s, f) for s in targets for f in features):
        frequency = (qs[:, index] < .05).mean()
        rows.append({'syndrome': syndrome, 'feature': feature, 'replicates': replicates,
                     'median_g': np.median(gs[:, index]), 'g_q025': np.quantile(gs[:, index], .025),
                     'g_q975': np.quantile(gs[:, index], .975), 'min_g': gs[:, index].min(), 'max_g': gs[:, index].max(),
                     'positive_g_frequency': (gs[:, index] > 0).mean(), 'significant_frequency': frequency,
                     'frequency_mc_standard_error': np.sqrt(frequency*(1-frequency)/replicates)})
    write(pd.DataFrame(rows), 'random_photo_selection_summary.csv')
    write(pd.DataFrame(records), 'random_photo_replicate_counts.csv')


def facekit_provenance():
    path = ROOT/'FaceKit'
    checkout, status = None, 'uninitialized'
    if (path/'.git').exists():
        try:
            top = subprocess.check_output(['git', 'rev-parse', '--show-toplevel'], cwd=path, text=True).strip()
            if Path(top).resolve() == path.resolve():
                checkout = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=path, text=True).strip()
                status = 'initialized; informational application checkout, not extraction provenance'
            else:
                status = 'not a distinct FaceKit repository'
        except subprocess.CalledProcessError:
            status = 'unreadable repository'
    recorded = subprocess.check_output(['git', 'rev-parse', 'HEAD:FaceKit'], cwd=ROOT, text=True).strip()
    return {'facekit_checkout': checkout, 'facekit_checkout_status': status, 'facekit_recorded_gitlink': recorded,
            'historical_source_sha256': sha(ROOT/'publication/historical_facekit/extractor.py')}


def provenance(args, started):
    INPUTS.update([ROOT/'stats_testing/prepare_cohorts.py', ROOT/'stats_testing/cohort_audit.py',
                  ROOT/'stats_testing/ancestry_analysis.py', ROOT/'stats_testing/statistical_analysis.py',
                  ROOT/'stats_testing/ancestry_disease_interaction.py', ROOT/'stats_testing/bootstrap_utils.py',
                  ROOT/'data/dataset_creation/fairface_facekit_pipeline.py',
                  ROOT/'publication/historical_facekit/extractor.py', ROOT/'publication/historical_facekit/SHA256SUMS'])
    data = {'command': sys.argv, 'python': sys.version, 'numpy': np.__version__, 'pandas': pd.__version__,
            'scipy': scipy.__version__, 'statsmodels': statsmodels.__version__, 'seed': SEED,
            'bootstrap_iterations': args.bootstrap, 'random_photo_replicates': args.photo_replicates,
            'random_photo_seed': SEED+1, 'git_head': subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
            'git_working_tree_dirty': bool(subprocess.check_output(['git', 'status', '--porcelain', '--untracked-files=no'], cwd=ROOT, text=True).strip()),
            **facekit_provenance(), 'phase': args.phase,
            'historical_facekit': 'e474628f3cb17a83edad59ac6e1d9e5438753032',
            'script_sha256': sha(__file__), 'elapsed_seconds': time.time()-started,
            'inputs_sha256': {str(p.relative_to(ROOT)): sha(p) for p in sorted(INPUTS)},
            'output_policy': 'Local only; no redistribution clearance inferred. No patient-level outputs.'}
    (OUT/'provenance.json').write_text(json.dumps(data,indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bootstrap', type=int, default=1000)
    parser.add_argument('--photo-replicates', type=int, default=100)
    parser.add_argument('--phase', choices=['all', 'audit', 'primary', 'robustness'], default='all')
    args = parser.parse_args()
    check(args.bootstrap >= 100 and args.photo_replicates >= 20, 'Minimum numerical settings violated')
    (OUT/'local_only').mkdir(exist_ok=True)
    INPUTS.clear()
    started = time.time()
    raw, reference, primary, alternate, ref, frozen, targets, features = load()
    if args.phase in ['all', 'audit']:
        audit(raw, reference, primary, alternate, ref, frozen, targets, features)
    if args.phase in ['all', 'primary']:
        result = contrast(primary, targets, features, iterations=args.bootstrap)
        write(result, 'within_gmdb_results.csv')
        comparison(result, frozen)
        restricted = contrast(primary, targets, features, comparator='eligible49')
        write(restricted, 'within_gmdb_eligible49_sensitivity.csv')
        adjusted(primary, features)
    if args.phase in ['all', 'robustness']:
        photo_sensitivity(primary)
        restricted_families()
        random_photos(raw, ref, targets, features, args.photo_replicates)
    provenance(args, started)
    print(f'Completed {args.phase} in {time.time()-started:.1f}s; only aggregate results written', flush=True)


if __name__ == '__main__':
    main()

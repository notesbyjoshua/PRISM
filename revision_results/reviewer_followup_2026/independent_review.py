#!/usr/bin/env python3
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import ndtr
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[2]
FEATURES = pd.read_csv(ROOT/'revision_results/final_20260914/exact_feature_names.csv').feature.tolist()
LAYERS = {'broad': 'layer1_broad_ethnicity', 'detailed': 'layer2_subcategory'}


def require(condition, message):
    if not bool(condition):
        raise AssertionError(message)


def close(actual, expected, message, rtol=1e-8, atol=1e-12):
    require(np.allclose(actual, expected, rtol=rtol, atol=atol, equal_nan=True), message)


def manual_bh(p):
    p = np.asarray(p, float)
    require(np.isfinite(p).all() and ((p >= 0) & (p <= 1)).all(), 'Invalid BH input')
    order = np.argsort(p, kind='stable')
    values = np.minimum.accumulate((len(p)*p[order]/np.arange(1, len(p)+1))[::-1])[::-1]
    result = np.empty(len(p))
    result[order] = np.minimum(values, 1)
    return result


def raw_cohorts():
    disease = pd.read_csv(ROOT/'data/master_dataset_disease.csv', dtype={'patient_id': 'string', 'image_id': 'string'}, low_memory=False)
    reference = pd.read_csv(ROOT/'data/master_dataset_healthy.csv', low_memory=False).copy()
    disease = disease.copy()
    for column in ['patient_id', 'image_id']:
        value = disease[column].str.strip()
        disease[column] = value.mask(value.str.lower().isin(['', 'nan', 'none', 'null', 'unknown', 'na', 'n/a']))
    require(disease.patient_id.notna().all(), 'Original recorded patient IDs missing')
    require(disease.image_id.notna().all() and disease.image_id.is_unique, 'Image keys not unique')
    year, month = pd.to_numeric(disease.age_year, errors='coerce'), pd.to_numeric(disease.age_month, errors='coerce')
    valid = np.isfinite(year) & np.isfinite(month) & year.ge(0) & month.between(0, 11)
    valid &= year.mod(1).eq(0) & month.mod(1).eq(0) & ~(year.eq(0) & month.eq(0))
    disease['_age'] = (year+month/12).where(valid)
    frontal = disease.frontal_ok.astype('string').str.strip().str.lower().eq('true')
    eligible = disease.loc[frontal & disease['_age'].ge(3) & disease['_age'].lt(10)].copy()
    pose = eligible[['pose_yaw', 'pose_pitch', 'pose_roll']].apply(pd.to_numeric, errors='coerce').to_numpy(float)
    scores = np.abs(pose)/np.array([15., 15., 10.])
    eligible['_pose_score'] = np.where(np.isfinite(pose).all(axis=1), np.nanmax(scores, axis=1), np.inf)
    for column in ['disease', 'internal_syndrome_name', 'gender']:
        require(eligible.groupby('patient_id')[column].nunique(dropna=False).le(1).all(), f'Patient label conflicts: {column}')
    primary = eligible.sort_values(['patient_id', '_pose_score', 'image_id'], kind='stable').drop_duplicates('patient_id')
    alternative = eligible.sort_values(['patient_id', 'image_id'], kind='stable').drop_duplicates('patient_id')
    primary, alternative = primary.sort_values('patient_id'), alternative.sort_values('patient_id')
    require(len(eligible) == 2101 and len(primary) == 1850, 'Independent childhood cohort flow differs')
    require(primary.patient_id.is_unique and primary.patient_id.tolist() == alternative.patient_id.tolist(), 'Selection patient universe differs')
    reference['_key'] = reference.file.astype('string').str.strip().str.replace('\\', '/', regex=False)
    require(reference.frontal_ok.astype('string').str.lower().eq('true').all(), 'Reference quality flag differs')
    reference = reference.drop_duplicates('_key')
    require(len(reference) == 942, 'Independent reference denominator differs')
    for frame in [primary, alternative, reference]:
        frame[FEATURES] = frame[FEATURES].apply(pd.to_numeric, errors='coerce').replace([np.inf, -np.inf], np.nan)
    return eligible, primary.copy(), alternative.copy(), reference.copy()


def direct_g(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    x, y = x[np.isfinite(x)], y[np.isfinite(y)]
    n, m = len(x), len(y)
    pooled = ((n-1)*np.var(x, ddof=1)+(m-1)*np.var(y, ddof=1))/(n+m-2)
    return (1-3/(4*(n+m)-9))*(x.mean()-y.mean())/math.sqrt(pooled)


def ranked_test(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    x, y = x[np.isfinite(x)], y[np.isfinite(y)]
    n, m = len(x), len(y)
    combined = np.concatenate([x, y])
    ranks = rankdata(combined, method='average')
    u = ranks[:n].sum()-n*(n+1)/2
    _, counts = np.unique(combined, return_counts=True)
    tie = sum(int(c)**3-int(c) for c in counts)
    variance = n*m/12*((n+m+1)-tie/((n+m)*(n+m-1)))
    p = min(1., 2*ndtr(-(abs(u-n*m/2)-.5)/math.sqrt(variance))) if variance > 0 else 1.
    return u, p


def replay_bootstrap(x, y, iterations, seed, columns):
    rng = np.random.default_rng(seed)
    xi = rng.integers(len(x), size=(iterations, len(x)))
    yi = rng.integers(len(y), size=(iterations, len(y)))
    rows = []
    for column in columns:
        xb, yb = x[:, column][xi], y[:, column][yi]
        nx, ny = np.isfinite(xb).sum(axis=1), np.isfinite(yb).sum(axis=1)
        with np.errstate(invalid='ignore', divide='ignore'):
            vx, vy = np.nanvar(xb, axis=1, ddof=1), np.nanvar(yb, axis=1, ddof=1)
            pooled = ((nx-1)*vx+(ny-1)*vy)/(nx+ny-2)
            g = (1-3/(4*(nx+ny)-9))*(np.nanmean(xb, axis=1)-np.nanmean(yb, axis=1))/np.sqrt(pooled)
        finite = g[np.isfinite(g) & (nx >= 2) & (ny >= 2)]
        rows.append((column, len(finite), *np.quantile(finite, [.025, .975])))
    return rows


def audit_tests(primary, result_dir, out):
    frozen = pd.read_csv(ROOT/'revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz')
    targets = sorted(frozen.disease.unique())
    targets_set = set(targets)
    main_eligible = primary.disease.isin(targets)
    composition, checked, boot, summaries = [], [], [], []
    settings = json.loads((result_dir/'provenance.json').read_text())
    for design, filename in [('all', 'within_gmdb_results.csv'), ('eligible49', 'within_gmdb_eligible49_sensitivity.csv')]:
        result = pd.read_csv(result_dir/filename)
        require(len(result) == 6125 and not result.duplicated(['syndrome', 'feature']).any(), f'{design} hypothesis family')
        require(set(map(tuple, result[['syndrome', 'feature']].to_numpy())) == {(s, f) for s in targets for f in FEATURES}, f'{design} hypothesis keys')
        for index, syndrome in enumerate(targets):
            target = primary.loc[primary.disease.eq(syndrome)]
            comparator = primary.loc[~primary.disease.eq(syndrome) & (main_eligible if design == 'eligible49' else True)]
            require(target.patient_id.is_unique and comparator.patient_id.is_unique, 'Repeated patient observation')
            require(not set(target.patient_id)&set(comparator.patient_id), 'Patient overlap across groups')
            require(not comparator.disease.eq(syndrome).any(), 'Target syndrome leaked into comparator')
            expected_total = int(main_eligible.sum()) if design == 'eligible49' else len(primary)
            require(len(target)+len(comparator) == expected_total, 'Comparator composition denominator')
            group_rows = result.loc[result.syndrome.eq(syndrome)].set_index('feature').loc[FEATURES]
            for other_syndrome, data in comparator.groupby('disease'):
                composition.append({'design': design, 'target_syndrome': syndrome, 'comparator_syndrome': other_syndrome,
                                    'n_patients': len(data), 'comparator_is_main_eligible': other_syndrome in targets_set})
            for feature in FEATURES:
                row = group_rows.loc[feature]
                x, y = target[feature].to_numpy(float), comparator[feature].to_numpy(float)
                nx, ny = np.isfinite(x).sum(), np.isfinite(y).sum()
                missing_x, missing_y = 1-nx/len(x), 1-ny/len(y)
                eligible = nx >= 10 and ny >= 10 and missing_x <= .30 and missing_y <= .30
                require(bool(row.estimable) == eligible, 'Estimability differs')
                close([row.target_total, row.comparator_total, row.n_target, row.n_comparator], [len(x), len(y), nx, ny], 'Denominator mismatch', rtol=0, atol=0)
                close([row.target_missing_fraction, row.comparator_missing_fraction], [missing_x, missing_y], 'Missingness differs')
                u, p = ranked_test(x, y)
                g = direct_g(x, y)
                close(g, row.hedges_g, 'Independent g mismatch')
                close(u, row.mannwhitney_u, 'Independent U mismatch')
                close(p, row.p_value, 'Independent tie/continuity-corrected p mismatch', rtol=1e-7, atol=1e-14)
                checked.append({'design': design, 'syndrome': syndrome, 'feature': feature, 'n_target': nx, 'n_comparator': ny,
                                'target_missing_fraction': missing_x, 'comparator_missing_fraction': missing_y,
                                'independent_g': g, 'independent_u': u, 'independent_p': p, 'passed': True})
            if design == 'all':
                require(group_rows.bootstrap_finite.eq(1000).all() and group_rows.bootstrap_iterations.eq(1000).all(), 'Bootstrap draw denominator')
                require(np.isfinite(group_rows[['ci_95_lower', 'ci_95_upper']]).all().all() and group_rows.ci_95_lower.le(group_rows.ci_95_upper).all(), 'Invalid bootstrap intervals')
                cases = {0, 61, 124, int(np.argmin(abs(group_rows.q_global_bh.to_numpy()-.05))),
                         int(np.argmax(abs(group_rows.hedges_g.to_numpy())))}
                for column, finite, lower, upper in replay_bootstrap(target[FEATURES].to_numpy(float), comparator[FEATURES].to_numpy(float),
                        1000, [int(settings['seed']), index], sorted(cases)):
                    row = group_rows.iloc[column]
                    close([lower, upper], [row.ci_95_lower, row.ci_95_upper], 'Independent patient bootstrap CI mismatch', rtol=1e-8)
                    require(finite == row.bootstrap_finite, 'Bootstrap finite-replicate count')
                    boot.append({'syndrome': syndrome, 'feature': FEATURES[column], 'draws': finite, 'independent_lower': lower,
                                 'independent_upper': upper, 'saved_lower': row.ci_95_lower, 'saved_upper': row.ci_95_upper, 'passed': True})
            else:
                require(group_rows.bootstrap_iterations.eq(0).all() and group_rows[['ci_95_lower', 'ci_95_upper']].isna().all().all(), 'Secondary intervals must be explicitly unavailable')
        independent = pd.DataFrame([r for r in checked if r['design'] == design]).set_index(['syndrome', 'feature'])
        aligned = independent.loc[pd.MultiIndex.from_frame(result[['syndrome', 'feature']])]
        q = manual_bh(aligned.independent_p.to_numpy())
        close(q, result.q_global_bh, 'Full-family independent BH differs', rtol=1e-7, atol=1e-13)
        require(np.array_equal(q < .05, result.significant), 'Discovery flags differ')
        summaries.append({'design': design, 'target_syndromes': 49, 'hypotheses': 6125, 'target_n_min': int(result.n_target.min()),
                          'target_n_max': int(result.n_target.max()), 'comparator_n_min': int(result.n_comparator.min()),
                          'comparator_n_max': int(result.n_comparator.max()), 'significant': int((q < .05).sum()),
                          'comparator_pool_n': int(main_eligible.sum()) if design == 'eligible49' else len(primary),
                          'pool_syndromes': 49 if design == 'eligible49' else primary.disease.nunique(),
                          'pool_members_outside_main_49': int((~main_eligible).sum()) if design == 'all' else 0})
    pd.DataFrame(checked).to_csv(out/'independent_all_hypotheses.csv', index=False)
    pd.DataFrame(boot).to_csv(out/'independent_bootstrap_intervals.csv', index=False)
    pd.DataFrame(composition).to_csv(out/'comparator_composition.csv', index=False)
    pd.DataFrame(summaries).to_csv(out/'comparator_designs.csv', index=False)
    return summaries, len(boot)


def audit_adjustment(primary, result_dir, out):
    saved = pd.read_csv(result_dir/'within_gmdb_age_sex_adjusted.csv')
    require(len(saved) == 4, 'Adjusted family must have exactly four models')
    frame = primary.copy()
    frame['_female'] = frame.gender.astype('string').str.strip().str.lower().map({'male': 0., 'm': 0., 'female': 1., 'f': 1.})
    rows = []
    for row in saved.itertuples():
        data = frame.dropna(subset=[row.feature, '_female', '_age'])
        target = data.disease.eq(row.syndrome).to_numpy(float)
        x = np.column_stack([np.ones(len(data)), target, data['_age'].to_numpy()-6, data['_female'].to_numpy()])
        y = data[row.feature].to_numpy(float)
        inverse = np.linalg.inv(x.T@x)
        beta = inverse@x.T@y
        leverage = np.einsum('ij,jk,ik->i', x, inverse, x)
        residual = y-x@beta
        hc3 = inverse@(x.T@(((residual/(1-leverage))**2)[:, None]*x))@inverse
        se = math.sqrt(hc3[1, 1])
        p = 2*ndtr(-abs(beta[1]/se))
        close([beta[1], se, beta[1]-1.959963984540054*se, beta[1]+1.959963984540054*se, p],
              [row.target_beta, row.hc3_se, row.ci_95_lower, row.ci_95_upper, row.p_value], 'Independent adjusted HC3 fit differs')
        close([len(data), target.sum(), len(primary)-len(data)], [row.n_complete, row.n_target, row.excluded], 'Adjusted sample-size mismatch')
        rows.append({'syndrome': row.syndrome, 'feature': row.feature, 'n_complete': len(data), 'target_n': int(target.sum()),
                     'excluded': len(primary)-len(data), 'independent_beta': beta[1], 'independent_hc3_se': se,
                     'independent_p': p, 'passed': True})
    close(manual_bh([r['independent_p'] for r in rows]), saved.q_four_models, 'Four-model BH mismatch')
    pd.DataFrame(rows).to_csv(out/'independent_adjusted_models.csv', index=False)


def audit_photos(eligible, primary, alternate, result_dir, out):
    require(primary.patient_id.tolist() == alternate.patient_id.tolist(), 'Photo-choice patient sets differ')
    changed = primary.image_id.to_numpy() != alternate.image_id.to_numpy()
    age_delta = alternate['_age'].to_numpy()-primary['_age'].to_numpy()
    expected = primary.assign(_changed=changed).groupby('disease').agg(
        n_patients=('patient_id', 'size'), photos_changed=('_changed', 'sum')).reset_index().rename(columns={'disease': 'syndrome'})
    saved = pd.read_csv(result_dir/'photo_change_by_syndrome.csv').set_index('syndrome').loc[expected.syndrome]
    close(saved.patients_photo_changed, expected.photos_changed, 'Affected syndrome photo counts')
    frozen = ROOT/'revision_results/final_20260914'
    first = pd.read_csv(frozen/'age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz')
    second = pd.read_csv(frozen/'age_3_9_id_sensitivity/significance_testing/stats_each_disease_vs_healthy.csv.gz')
    common = first.merge(second, on=['disease', 'feature'], suffixes=('_p', '_a'), validate='one_to_one')
    common['_absolute_delta'] = abs(common.hedges_g_a-common.hedges_g_p)
    common['_same_sign'] = np.sign(common.hedges_g_a) == np.sign(common.hedges_g_p)
    common['_discordant'] = common.mannwhitney_q_fdr_global_p.lt(.05) != common.mannwhitney_q_fdr_global_a.lt(.05)
    rows = []
    for syndrome, group in common.groupby('disease'):
        rows.append({'syndrome': syndrome, 'hypotheses': len(group), 'median_absolute_g_difference': group['_absolute_delta'].median(),
                     'q95_absolute_g_difference': group['_absolute_delta'].quantile(.95), 'max_absolute_g_difference': group['_absolute_delta'].max(),
                     'pearson_g': group.hedges_g_p.corr(group.hedges_g_a), 'sign_agreement': group['_same_sign'].mean(),
                     'discordant_significance_calls': int(group['_discordant'].sum()),
                     'primary_significant': int(group.mannwhitney_q_fdr_global_p.lt(.05).sum()),
                     'alternative_significant': int(group.mannwhitney_q_fdr_global_a.lt(.05).sum())})
    by_syndrome = expected.merge(pd.DataFrame(rows), on='syndrome', how='left', validate='one_to_one')
    unchanged = by_syndrome.loc[by_syndrome.photos_changed.eq(0) & by_syndrome.hypotheses.notna()]
    require(unchanged.max_absolute_g_difference.eq(0).all(), 'g changed in syndrome with unchanged photographs')
    by_syndrome.to_csv(out/'photo_changes_and_effect_agreement_by_syndrome.csv', index=False)
    result = {'photos_changed': int(changed.sum()), 'syndromes_changed': int(expected.photos_changed.gt(0).sum()),
              'main_49_syndromes_changed': int(by_syndrome.loc[by_syndrome.hypotheses.notna(), 'photos_changed'].gt(0).sum()),
              'changed_photos_with_changed_recorded_age': int((changed & (abs(age_delta) > 1e-12)).sum()),
              'median_absolute_age_change_among_changed_photos': float(np.median(abs(age_delta[changed]))),
              'primary_significant': int(common.mannwhitney_q_fdr_global_p.lt(.05).sum()),
              'alternative_significant': int(common.mannwhitney_q_fdr_global_a.lt(.05).sum()),
              'shared_significant': int((common.mannwhitney_q_fdr_global_p.lt(.05)&common.mannwhitney_q_fdr_global_a.lt(.05)).sum()),
              'significance_discordance_in_photo_unchanged_syndromes': int(unchanged.discordant_significance_calls.sum()),
              'median_absolute_g_difference': common['_absolute_delta'].median(), 'q95_absolute_g_difference': common['_absolute_delta'].quantile(.95),
              'maximum_absolute_g_difference': common['_absolute_delta'].max(), 'pearson_g': common.hedges_g_p.corr(common.hedges_g_a),
              'spearman_g': common.hedges_g_p.corr(common.hedges_g_a, method='spearman'), 'sign_agreement': common['_same_sign'].mean()}
    (out/'photo_review_summary.json').write_text(json.dumps(result, indent=2))
    return result


def audit_thresholds(result_dir, out):
    summaries = pd.read_csv(result_dir/'subgroup_threshold_sensitivity.csv')
    all_rows = pd.read_csv(result_dir/'subgroup_restricted_hypotheses.csv', low_memory=False)
    checks = []
    for summary in summaries.itertuples():
        analysis, layer, threshold = summary.analysis, summary.layer, summary.minimum_disease_n
        file = 'analysis_B_matched_disease_vs_healthy' if analysis == 'B' else 'analysis_C_all_interactions'
        original = pd.read_csv(ROOT/f'revision_results/final_20260914/age_3_9/ethnicity_analysis/{LAYERS[layer]}/{file}.csv.gz')
        if analysis == 'B':
            mask = original.n_disease.ge(threshold)&original.n_healthy_matched.ge(20)
            keys, p, effect, lower, upper = ['syndrome_name', 'feature', 'group'], 'p_value', 'hedges_g', 'hedges_g_ci_lower', 'hedges_g_ci_upper'
        else:
            mask = original.n_disease_group_1.ge(threshold)&original.n_disease_group_2.ge(threshold)
            mask &= original.n_healthy_group_1.ge(20)&original.n_healthy_group_2.ge(20)
            keys, p, effect, lower, upper = ['syndrome_name', 'feature', 'group_1', 'group_2'], 'interaction_p_value', 'interaction_beta_group2_minus_group1', 'interaction_ci_lower', 'interaction_ci_upper'
        subset = original.loc[mask]
        q = manual_bh(subset[p])
        require(len(subset) == summary.hypotheses and (q < .05).sum() == summary.significant_recorrected, 'Restricted-family denominator/discovery mismatch')
        saved = all_rows.loc[all_rows.analysis.eq(analysis)&all_rows.layer.eq(layer)&all_rows.threshold.eq(threshold)]
        aligned = saved.merge(subset[keys+[p, effect, lower, upper]], on=keys, suffixes=('_saved', '_original'), validate='one_to_one')
        require(len(aligned) == len(subset), 'Restricted hypothesis keys mismatch')
        for column in [p, effect, lower, upper]:
            close(aligned[f'{column}_saved'], aligned[f'{column}_original'], 'Filtering changed saved fit/CI')
        close(manual_bh(aligned[f'{p}_original']), aligned.q_restricted_family, 'Independent restricted BH differs')
        checks.append({'analysis': analysis, 'layer': layer, 'minimum_disease_n': threshold, 'hypotheses': len(subset),
                       'significant_recorrected': int((q < .05).sum()), 'fits_unchanged': True, 'operation': 'restricted-family BH re-correction; not refitting'})
    pd.DataFrame(checks).to_csv(out/'independent_threshold_families.csv', index=False)


def main():
    parser = argparse.ArgumentParser(description='Independent local-only statistical checks; never exports patient keys')
    parser.add_argument('--results-dir', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--output-dir', type=Path, default=Path(__file__).resolve().parent/'local_only/independent_review')
    args = parser.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    eligible, primary, alternate, reference = raw_cohorts()
    designs, bootstraps = audit_tests(primary, args.results_dir, out)
    audit_adjustment(primary, args.results_dir, out)
    photo = audit_photos(eligible, primary, alternate, args.results_dir, out)
    audit_thresholds(args.results_dir, out)
    summary = {'status': 'passed', 'point_estimate_hypotheses_verified': 12250, 'bootstrap_intervals_replayed': bootstraps,
               'bootstrap_scope': 'representative features in every one of 49 target syndromes; 1000 materialized patient resamples each',
               'adjusted_models_verified': 4, 'restricted_families_verified': 16, 'designs': designs, 'photo': photo,
               'identity_limit': 'recorded-ID independence only; cannot exclude a person assigned multiple IDs',
               'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (out/'review_checks.json').write_text(json.dumps(summary, indent=2))
    print(f'Independent review passed: 12,250 point-estimate hypotheses; {bootstraps} bootstrap intervals; 4 HC3 models; 16 threshold families')


if __name__ == '__main__':
    main()

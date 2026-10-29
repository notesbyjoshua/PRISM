#!/usr/bin/env python3
import argparse
import importlib.metadata
import json
import os
import subprocess
import sys
import textwrap
import time
import traceback
import warnings
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import pandas as pd
import scipy.stats as st
import statsmodels.api as sm
from scipy.stats import rankdata

import followup as f

SEED = 20261026
ANCHORS = [s for s, _ in f.HIGHLIGHTS]
STRATA = [('age_3_to_under_6', 3., 6.), ('age_6_to_under_10', 6., 10.)]


def label(value):
    return value.rsplit('__', 1)[0]


def plot_label(value):
    text = label(value)
    if len(text) > 58 and ';' in text:
        text = text.rsplit(';', 1)[-1].strip()
    return textwrap.shorten(text, width=58, placeholder='...')


def require_unique(frame):
    f.check(frame.patient_id.notna().all() and frame.patient_id.is_unique, 'Nonunique or missing recorded IDs')


def planned_pairs(targets, anchors=ANCHORS):
    f.check(len(targets) == len(set(targets)) and set(anchors).issubset(targets), 'Invalid target/anchor labels')
    pairs, seen = [], set()
    for anchor in anchors:
        for other in sorted(targets):
            key = frozenset([anchor, other])
            if anchor != other and key not in seen:
                pairs.append((anchor, other))
                seen.add(key)
    return pairs


def corrected(frame, family, denominator):
    f.check(len(frame) == denominator, 'Planned family rows missing')
    frame = frame.copy()
    valid = frame.status.eq('ok') & np.isfinite(frame.p_value)
    frame['correction_p'] = np.where(valid, frame.p_value, 1.)
    frame['q_global_bh'] = f.bh(frame.correction_p)
    frame.loc[~valid, 'q_global_bh'] = np.nan
    frame['significant'] = valid & frame.q_global_bh.lt(.05)
    frame['fdr_family'] = family
    frame['fdr_planned_hypotheses'] = denominator
    return frame


def adjusted_fit(frame, syndrome, feature):
    target = frame.disease.eq(syndrome)
    female = frame.sex_label.map({'Male': 0., 'Female': 1.}).to_numpy(float)
    age = pd.to_numeric(frame.conditional_age_years, errors='coerce').to_numpy(float)
    y = frame[feature].to_numpy(float)
    covariates = np.isfinite(female) & np.isfinite(age)
    complete = covariates & np.isfinite(y)
    a, b = complete & target, complete & ~target
    nt, no = int(a.sum()), int(b.sum())
    row = {'syndrome': syndrome, 'feature': feature, 'target_total': int(target.sum()),
           'other_total': int((~target).sum()), 'n_target': nt, 'n_other': no, 'n_complete': nt+no,
           'target_missing_covariates': int((target & ~covariates).sum()),
           'other_missing_covariates': int((~target & ~covariates).sum()),
           'target_missing_feature': int((target & ~np.isfinite(y)).sum()),
           'other_missing_feature': int((~target & ~np.isfinite(y)).sum()),
           'target_excluded_fraction': 1-nt/target.sum(), 'other_excluded_fraction': 1-no/(~target).sum(),
           'target_beta': np.nan, 'hc3_se': np.nan, 'ci_95_lower': np.nan, 'ci_95_upper': np.nan,
           'p_value': np.nan, 'status': 'unsupported', 'reason': '', 'warning_flags': ''}
    reasons, flags = [], []
    if nt < 10 or no < 20:
        reasons.append('sparse_complete_cases')
    if 10*(target.sum()-nt) > 3*target.sum() or 10*((~target).sum()-no) > 3*(~target).sum():
        reasons.append('exclusion_above_30_percent')
    for name, mask in [('target', a), ('other', b)]:
        row[f'{name}_female_n'] = int((mask & (female == 1)).sum())
        row[f'{name}_male_n'] = int((mask & (female == 0)).sum())
        row[f'{name}_age_min'] = float(age[mask].min()) if mask.any() else np.nan
        row[f'{name}_age_max'] = float(age[mask].max()) if mask.any() else np.nan
        if min(row[f'{name}_female_n'], row[f'{name}_male_n']) < 5:
            flags.append(f'{name}_sex_cell_below_5')
    if nt and no:
        low, high = max(age[a].min(), age[b].min()), min(age[a].max(), age[b].max())
        in_common = (age >= low) & (age <= high)
        row.update(age_overlap_low=low, age_overlap_high=high,
                   target_in_age_overlap=int((a & in_common).sum()), other_in_age_overlap=int((b & in_common).sum()))
        if high <= low or (a & in_common).sum() < 10 or (b & in_common).sum() < 10:
            reasons.append('insufficient_recorded_age_overlap')
        if max(np.quantile(age[a], .1), np.quantile(age[b], .1)) >= min(np.quantile(age[a], .9), np.quantile(age[b], .9)):
            flags.append('central_age_ranges_do_not_overlap')
        if min((a & in_common).sum()/nt, (b & in_common).sum()/no) < .5:
            flags.append('less_than_half_in_common_age_range')
    if not reasons:
        design = np.column_stack([np.ones(complete.sum()), target[complete], age[complete]-6, female[complete]])
        rank, condition = np.linalg.matrix_rank(design), np.linalg.cond(design)
        row.update(design_rank=int(rank), condition_number=float(condition))
        if rank != 4 or not np.isfinite(condition) or condition >= 1e8:
            reasons.append('rank_or_condition_failure')
        else:
            model = sm.OLS(y[complete], design).fit(cov_type='HC3', use_t=False)
            leverage = np.einsum('ij,ji->i', design, np.linalg.pinv(design))
            row['maximum_leverage'] = float(leverage.max())
            if leverage.max() >= 1-1e-8:
                reasons.append('unit_leverage')
            if np.ptp(y[complete]) == 0 or model.ssr <= np.finfo(float).eps*np.dot(y[complete], y[complete]):
                reasons.append('zero_residual_variation')
            interval = model.conf_int()[1]
            if not np.isfinite([model.params[1], model.bse[1], model.pvalues[1], *interval]).all() or model.bse[1] <= 0:
                reasons.append('nonfinite_hc3')
            if leverage.max() > .2:
                flags.append('maximum_leverage_above_0_2')
            if not reasons:
                row.update(target_beta=model.params[1], hc3_se=model.bse[1], ci_95_lower=interval[0],
                           ci_95_upper=interval[1], p_value=model.pvalues[1], status='ok')
    row['reason'], row['warning_flags'] = ';'.join(reasons), ';'.join(flags)
    return row


def full_adjusted(primary, targets, features):
    require_unique(primary)
    rows = []
    for i, syndrome in enumerate(targets):
        rows.extend(adjusted_fit(primary, syndrome, feature) for feature in features)
        if (i+1) % 10 == 0 or i+1 == len(targets):
            print(f'Adjusted screen: {i+1}/{len(targets)}', flush=True)
    return corrected(pd.DataFrame(rows), 'adjusted_all_other_GMDB', len(targets)*len(features))


def two_groups(a, b, features, iterations, seed, min_b=10):
    require_unique(a)
    require_unique(b)
    f.check(not set(a.patient_id) & set(b.patient_id), 'Recorded ID overlap in contrast')
    x, y = a[features].to_numpy(float), b[features].to_numpy(float)
    nx, ny = np.isfinite(x).sum(axis=0), np.isfinite(y).sum(axis=0)
    valid = (nx >= 10) & (ny >= min_b) & (10*(len(x)-nx) <= 3*len(x)) & (10*(len(y)-ny) <= 3*len(y))
    g, u, p = np.full(len(features), np.nan), np.full(len(features), np.nan), np.full(len(features), np.nan)
    ci, finite = np.full((2, len(features)), np.nan), np.zeros(len(features), int)
    selected = np.flatnonzero(valid)
    if len(selected):
        g[selected] = f.effect(x[:, selected], y[:, selected])
        test = st.mannwhitneyu(x[:, selected], y[:, selected], axis=0, alternative='two-sided',
                              method='asymptotic', use_continuity=True, nan_policy='omit')
        u[selected], p[selected] = test.statistic, test.pvalue
        selected = selected[np.isfinite(g[selected])]
        if len(selected):
            draws = f.bootstrap_effect(x[:, selected], y[:, selected], iterations, seed)
            finite[selected] = np.isfinite(draws).sum(axis=0)
            for j, col in enumerate(selected):
                if finite[col]:
                    ci[:, col] = np.nanquantile(draws[:, j], [.025, .975])
    rows = []
    for j, feature in enumerate(features):
        ok = bool(valid[j] and np.isfinite(g[j]) and np.isfinite(p[j]))
        reason = '' if ok else 'sparse_or_missing' if not valid[j] else 'zero_variance_or_nonfinite'
        flags = 'bootstrap_finite_below_95_percent' if ok and finite[j] < .95*iterations else ''
        rows.append({'feature': feature, 'total_a': len(a), 'total_b': len(b), 'n_a': nx[j], 'n_b': ny[j],
                     'missing_a': len(a)-nx[j], 'missing_b': len(b)-ny[j],
                     'status': 'ok' if ok else 'unsupported', 'reason': reason, 'warning_flags': flags,
                     'hedges_g': g[j] if ok else np.nan, 'mannwhitney_u': u[j], 'p_value': p[j] if ok else np.nan,
                     'ci_95_lower': ci[0, j], 'ci_95_upper': ci[1, j], 'bootstrap_finite': finite[j],
                     'bootstrap_iterations': iterations})
    return rows


def pairwise(primary, targets, features, iterations):
    require_unique(primary)
    groups = {s: primary.loc[primary.disease.eq(s)] for s in targets}
    pairs, rows = planned_pairs(targets), []
    for i, (a, b) in enumerate(pairs):
        part = two_groups(groups[a], groups[b], features, iterations, [SEED, 1, i])
        rows.extend(dict(syndrome_a=a, syndrome_b=b, pair_index=i, **row) for row in part)
        if (i+1) % 30 == 0 or i+1 == len(pairs):
            print(f'Pairwise: {i+1}/{len(pairs)}', flush=True)
    return corrected(pd.DataFrame(rows), 'targeted_pairwise_186_pairs', len(pairs)*len(features))


def age_stratified(primary, targets, features, iterations):
    require_unique(primary)
    rows = []
    for k, (name, low, high) in enumerate(STRATA):
        frame = primary.loc[primary.conditional_age_years.ge(low) & primary.conditional_age_years.lt(high)]
        for i, syndrome in enumerate(targets):
            a, b = frame.loc[frame.disease.eq(syndrome)], frame.loc[~frame.disease.eq(syndrome)]
            part = two_groups(a, b, features, iterations, [SEED, 2, k, i], min_b=20)
            rows.extend(dict(stratum=name, age_lower_inclusive=low, age_upper_exclusive=high,
                             syndrome=syndrome, **row) for row in part)
        print(f'Age sensitivity: {name} completed', flush=True)
    return corrected(pd.DataFrame(rows), 'age_strata_two_bins', 2*len(targets)*len(features))


def from_moments(n_a, sum_a, ss_a, n_b, sum_b, ss_b):
    with np.errstate(divide='ignore', invalid='ignore'):
        pooled = (ss_a-sum_a**2/n_a+ss_b-sum_b**2/n_b)/(n_a+n_b-2)
        value = (1-3/(4*(n_a+n_b)-9))*(sum_a/n_a-sum_b/n_b)/np.sqrt(np.maximum(pooled, 0))
    return np.where((n_a >= 2) & (n_b >= 2) & np.isfinite(value), value, np.nan)


def moment_groups(values, starts):
    finite = np.isfinite(values)
    clean = np.where(finite, values, 0.)
    return tuple(np.add.reduceat(v, starts, axis=0) for v in [finite.astype(float), clean, clean**2])


def shared_effects(values, starts, target_indices, reference, comparator_indices):
    n, total, ss = moment_groups(values, starts)
    nt, tt, stt = n[target_indices], total[target_indices], ss[target_indices]
    nr = np.isfinite(reference).sum(axis=0)
    tr, sr = np.nansum(reference, axis=0), np.nansum(reference**2, axis=0)
    fair = from_moments(nt, tt, stt, nr, tr, sr)
    nc, tc, sc = n[comparator_indices].sum(axis=0), total[comparator_indices].sum(axis=0), ss[comparator_indices].sum(axis=0)
    within = from_moments(nt, tt, stt, nc-nt, tc-tt, sc-stt)
    return fair, within


METRICS = ['direction_agreement', 'pearson_r', 'spearman_rho', 'median_absolute_attenuation',
           'fraction_absolute_attenuated', 'median_signed_g_difference']


def agreement(x, y):
    valid = np.isfinite(x) & np.isfinite(y)
    x, y = np.asarray(x)[valid], np.asarray(y)[valid]
    if not len(x):
        return np.full(len(METRICS), np.nan)
    pearson = np.corrcoef(x, y)[0, 1] if np.ptp(x) > 0 and np.ptp(y) > 0 else np.nan
    rx, ry = rankdata(x), rankdata(y)
    spearman = np.corrcoef(rx, ry)[0, 1] if np.ptp(rx) > 0 and np.ptp(ry) > 0 else np.nan
    attenuation = np.abs(x)-np.abs(y)
    return np.array([np.mean(np.sign(x) == np.sign(y)), pearson, spearman, np.median(attenuation),
                     np.mean(attenuation > 0), np.median(y-x)])


def grouping_indices(mapping, features, n_targets):
    rows = [('overall', 'all_features', np.arange(n_targets*len(features)))]
    for level in ['anatomical_region', 'feature_family']:
        for name, group in mapping.groupby(level, sort=True):
            columns = np.array([features.index(v) for v in group.feature])
            indices = (np.arange(n_targets)[:, None]*len(features)+columns).ravel()
            rows.append((level, name, indices))
    return rows


def concordance(primary, reference, frozen, targets, features, mapping, existing, iterations):
    require_unique(primary)
    frame = primary.sort_values(['disease', 'patient_id']).copy()
    group_labels, sizes = np.unique(frame.disease, return_counts=True)
    starts = np.r_[0, np.cumsum(sizes)[:-1]]
    target_indices = np.array([np.flatnonzero(group_labels == s)[0] for s in targets])
    values, ref = frame[features].to_numpy(float), reference[features].to_numpy(float)
    indices = grouping_indices(mapping, features, len(targets))
    keys = pd.MultiIndex.from_product([targets, features], names=['syndrome', 'feature'])
    base = frozen.rename(columns={'disease': 'syndrome'}).set_index(['syndrome', 'feature']).reindex(keys)
    f.check(base.hedges_g.notna().all(), 'Incomplete frozen concordance input')
    modes = {'all': np.arange(len(group_labels)), 'eligible49': target_indices}
    results, joins, distributions = {}, {}, {}
    for mode, selected in modes.items():
        table = existing[mode].set_index(['syndrome', 'feature']).reindex(keys)
        f.check(table.hedges_g.notna().all() and table.estimable.all(), 'Concordance requires complete existing effects')
        fair, within = shared_effects(values, starts, target_indices, ref, selected)
        f.check(np.allclose(fair.ravel(), base.hedges_g, atol=1e-10, rtol=1e-9), 'Frozen g reconstruction mismatch')
        f.check(np.allclose(within.ravel(), table.hedges_g, atol=1e-10, rtol=1e-9), 'Within-GMDB g reconstruction mismatch')
        joined = keys.to_frame(index=False)
        joined['fairface_g'], joined['within_g'] = base.hedges_g.to_numpy(), table.hedges_g.to_numpy()
        joined['fairface_significant'] = base.mannwhitney_q_fdr_global.to_numpy() < .05
        joined['within_significant'] = table.significant.to_numpy()
        joined['absolute_attenuation'] = joined.fairface_g.abs()-joined.within_g.abs()
        joined = joined.merge(mapping[['feature', 'anatomical_region', 'feature_family']], on='feature', validate='many_to_one', sort=False)
        f.check(joined[['syndrome', 'feature']].equals(keys.to_frame(index=False)), 'Concordance order changed')
        joins[mode], results[mode] = joined, []
        distributions[mode] = np.full((iterations, len(indices), len(METRICS)), np.nan)
        for level, name, cols in indices:
            group = joined.iloc[cols]
            sf, sw = group.fairface_significant, group.within_significant
            both, union = int((sf & sw).sum()), int((sf | sw).sum())
            row = dict(comparator=mode, level=level, group=name, hypotheses=len(group),
                       features=group.feature.nunique(), syndromes=group.syndrome.nunique(),
                       fairface_significant=int(sf.sum()), within_significant=int(sw.sum()), both=both,
                       fairface_only=int((sf & ~sw).sum()), within_only=int((~sf & sw).sum()),
                       neither=int((~sf & ~sw).sum()), jaccard=both/union if union else np.nan,
                       fairface_discovery_retention=both/sf.sum() if sf.sum() else np.nan,
                       median_absolute_fairface_g=group.fairface_g.abs().median(),
                       median_absolute_within_g=group.within_g.abs().median())
            row.update(zip(METRICS, agreement(group.fairface_g.to_numpy(), group.within_g.to_numpy())))
            results[mode].append(row)
    rng = np.random.default_rng([SEED, 3])
    for rep in range(iterations):
        draw = np.concatenate([rng.integers(start, start+n, size=n) for start, n in zip(starts, sizes)])
        sampled_ref = ref[rng.integers(len(ref), size=len(ref))]
        for mode, selected in modes.items():
            fair, within = shared_effects(values[draw], starts, target_indices, sampled_ref, selected)
            for j, (_, _, cols) in enumerate(indices):
                distributions[mode][rep, j] = agreement(fair.ravel()[cols], within.ravel()[cols])
        if (rep+1) % 100 == 0 or rep+1 == iterations:
            print(f'Joint concordance bootstrap: {rep+1}/{iterations}', flush=True)
    rows = []
    for mode in modes:
        for j, row in enumerate(results[mode]):
            row['bootstrap_iterations'] = iterations
            for k, name in enumerate(METRICS):
                draws = distributions[mode][:, j, k]
                finite = draws[np.isfinite(draws)]
                low, high = np.quantile(finite, [.025, .975]) if len(finite) else (np.nan, np.nan)
                mean = float(np.mean(finite)) if len(finite) else np.nan
                se = float(np.std(finite, ddof=1)) if len(finite) > 1 else np.nan
                bias = mean-row[name]
                row.update({f'{name}_ci_low': low, f'{name}_ci_high': high, f'{name}_finite': len(finite),
                            f'{name}_bootstrap_mean': mean, f'{name}_bootstrap_se': se,
                            f'{name}_bootstrap_bias': bias,
                            f'{name}_bias_exceeds_se': bool(abs(bias) > se),
                            f'{name}_interval_excludes_point': bool(row[name] < low or row[name] > high)})
            rows.append(row)
    return pd.DataFrame(rows), pd.concat([v.assign(comparator=k) for k, v in joins.items()], ignore_index=True)


def redundancy(features, mapping, existing, adjusted):
    matrix = f.read(f.ROOT/'publication/feature_correlation/feature_correlation_matrix.csv').set_index('feature')
    clusters = f.read(f.ROOT/'publication/feature_correlation/feature_clusters.csv')
    f.check(matrix.index.tolist() == features and matrix.columns.tolist() == features, 'Correlation schema mismatch')
    rho = matrix.to_numpy(float)
    f.check(np.isfinite(rho).all() and np.allclose(rho, rho.T) and np.allclose(np.diag(rho), 1), 'Invalid correlation matrix')
    f.check(len(clusters) == len(features) and set(clusters.feature) == set(features), 'Invalid cluster membership')
    merged = clusters.merge(mapping[['feature', 'anatomical_region', 'feature_family']], on='feature', suffixes=('', '_dictionary'))
    f.check(merged.anatomical_region.eq(merged.anatomical_region_dictionary).all(), 'Cluster anatomy differs')
    upper = np.abs(rho[np.triu_indices(len(features), 1)])
    thresholds = pd.DataFrame([{'absolute_rho_threshold': t, 'pairs': int((upper >= t).sum()),
                               'total_pairs': len(upper), 'source': 'existing_FairFace_942_path_Spearman'} for t in [.7, .8, .9, .95]])
    annotated = existing.merge(clusters[['feature', 'cluster_id']], on='feature', validate='many_to_one')
    annotated = annotated.merge(adjusted[['syndrome', 'feature', 'status', 'significant']].rename(
        columns={'status': 'adjusted_status', 'significant': 'adjusted_significant'}), on=['syndrome', 'feature'], validate='one_to_one')
    rows = []
    for cluster, group in annotated.groupby('cluster_id'):
        subset = clusters.loc[clusters.cluster_id.eq(cluster)]
        rows.append({'cluster_id': cluster, 'features': subset.feature.nunique(),
                     'anatomical_regions': ';'.join(sorted(subset.anatomical_region.unique())),
                     'hypotheses': len(group), 'fairface_discoveries': int(group.fairface_significant.sum()),
                     'within_discoveries': int(group.within_significant.sum()),
                     'adjusted_supported': int(group.adjusted_status.eq('ok').sum()),
                     'adjusted_discoveries': int(group.adjusted_significant.sum()),
                     'median_absolute_attenuation': group.absolute_attenuation.median(),
                     'syndromes_with_within_discovery': group.loc[group.within_significant, 'syndrome'].nunique(),
                     'median_absolute_within_cluster_rho': subset.median_absolute_within_cluster_rho.iloc[0]})
    return pd.DataFrame(rows), thresholds


def age_agreement(age):
    first = age.loc[age.stratum.eq(STRATA[0][0])]
    second = age.loc[age.stratum.eq(STRATA[1][0])]
    merged = first.merge(second, on=['syndrome', 'feature'], suffixes=('_younger', '_older'), validate='one_to_one')
    rows = []
    for syndrome, group in merged.groupby('syndrome'):
        supported = group.status_younger.eq('ok') & group.status_older.eq('ok')
        subset = group.loc[supported]
        x, y = subset.hedges_g_younger.to_numpy(), subset.hedges_g_older.to_numpy()
        row = {'syndrome': syndrome, 'features_supported_both': len(subset),
               'younger_target_n': group.total_a_younger.iloc[0], 'older_target_n': group.total_a_older.iloc[0],
               'younger_supported': int(group.status_younger.eq('ok').sum()),
               'older_supported': int(group.status_older.eq('ok').sum()),
               'younger_significant': int(group.significant_younger.sum()),
               'older_significant': int(group.significant_older.sum()),
               'both_significant': int((group.significant_younger & group.significant_older).sum()),
               'direction_agreement': np.mean(np.sign(x) == np.sign(y)) if len(subset) else np.nan,
               'median_absolute_g_difference': np.median(np.abs(y-x)) if len(subset) else np.nan,
               'pearson_r': np.corrcoef(x, y)[0, 1] if len(x) > 1 and np.ptp(x) > 0 and np.ptp(y) > 0 else np.nan}
        rows.append(row)
    return pd.DataFrame(rows)


def figures(out, adjusted, pairs, concord, redundant, age, mapping):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none'})
    directory = out/'figures'
    directory.mkdir()

    def save(fig, name, caption):
        fig.text(.02, .015, caption, fontsize=8, va='bottom')
        fig.savefig(directory/f'{name}.png', dpi=350, bbox_inches='tight')
        fig.savefig(directory/f'{name}.svg', bbox_inches='tight')
        plt.close(fig)

    summary = adjusted.merge(mapping[['feature', 'anatomical_region']], on='feature', validate='many_to_one')
    summary['display'] = summary.syndrome.map(label)
    heat = summary.pivot_table(index='display', columns='anatomical_region', values='significant', aggfunc='mean')*100
    supported = summary.assign(supported=summary.status.eq('ok')).pivot_table(
        index='display', columns='anatomical_region', values='supported', aggfunc='sum')
    heat = heat.where(supported.gt(0))
    fig, ax = plt.subplots(figsize=(11, 15))
    cmap = plt.get_cmap('Blues').copy()
    cmap.set_bad('#bdbdbd')
    image = ax.imshow(heat, aspect='auto', cmap=cmap, vmin=0, vmax=max(5., heat.max().max()))
    ax.set_yticks(np.arange(len(heat)), heat.index, fontsize=7)
    ax.set_xticks(np.arange(len(heat.columns)), heat.columns, rotation=45, ha='right', fontsize=8)
    ax.set_title('Full age/sex-adjusted screen: significant features by anatomical region', loc='left')
    fig.colorbar(image, ax=ax, label='Percent of planned region features with global BH q < 0.05', shrink=.6)
    fig.subplots_adjust(left=.37, bottom=.2, top=.96)
    save(fig, 'adjusted_anatomical_summary', '6,125 planned models; unsupported models remain in denominators. Gray = unsupported. HC3; all-other-GMDB comparators.')

    counts = pairs.groupby(['syndrome_a', 'syndrome_b']).significant.sum()
    others = sorted(set(pairs.syndrome_a) | set(pairs.syndrome_b))
    heat = np.full((len(ANCHORS), len(others)), np.nan)
    for i, anchor in enumerate(ANCHORS):
        for j, other in enumerate(others):
            if (anchor, other) in counts:
                heat[i, j] = counts[anchor, other]
            elif (other, anchor) in counts:
                heat[i, j] = counts[other, anchor]
    names = [plot_label(s) for s in others]
    f.check(len(set(names)) == len(names), 'Plot labels are ambiguous')
    pd.DataFrame({'syndrome': others, 'plot_label': names}).to_csv(out/'syndrome_plot_labels.csv', index=False)
    fig, ax = plt.subplots(figsize=(10, 13))
    image = ax.imshow(heat.T, aspect='auto', cmap='viridis', vmin=0)
    ax.set_xticks(range(4), ['Cornelia de Lange', 'Noonan', 'Williams-Beuren', 'Rubinstein-Taybi'], rotation=35, ha='right', fontsize=8)
    ax.set_yticks(range(len(others)), names, fontsize=7)
    ax.set_title('Targeted pairwise comparisons: significant feature counts', loc='left')
    fig.colorbar(image, ax=ax, label='Features with global BH q < 0.05 (of 125)', shrink=.55)
    fig.subplots_adjust(bottom=.14, left=.49, right=.94, top=.94)
    save(fig, 'pairwise_discovery_counts', '186 unique pairs x 125 features; anchor-anchor cells repeat visually only.\nLabels abbreviated for display; full mapping supplied. No diagnostic specificity is inferred.')

    region = concord.loc[concord.level.eq('anatomical_region')]
    names = sorted(region.group.unique())
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
    for offset, mode, color in [(-.12, 'all', '#2166ac'), (.12, 'eligible49', '#b35806')]:
        group = region.loc[region.comparator.eq(mode)].set_index('group').reindex(names)
        for ax, metric in zip(axes, ['direction_agreement', 'median_absolute_attenuation']):
            x = group[metric].to_numpy()
            low, high = group[f'{metric}_ci_low'].to_numpy(), group[f'{metric}_ci_high'].to_numpy()
            y = np.arange(len(names))+offset
            ax.hlines(y, low, high, color=color, linewidth=1.7)
            ax.plot(x, y, 'o', color=color, label=mode)
    axes[0].set_yticks(range(len(names)), names)
    axes[0].set_xlabel('Direction agreement fraction')
    axes[1].set_xlabel('Median |FairFace g| - |GMDB g|')
    axes[1].axvline(0, color='.6', linestyle='--')
    axes[1].legend(title='GMDB comparator', fontsize=8)
    axes[0].set_title('Effect direction')
    axes[1].set_title('Absolute effect attenuation')
    fig.subplots_adjust(left=.27, right=.98, bottom=.2, top=.9, wspace=.15)
    save(fig, 'concordance_by_anatomy', 'Marginal 95% joint bootstrap percentile intervals; positive attenuation means smaller |g|.\nIntervals may exclude point estimates: see bootstrap-bias flags. Nominal population coverage is unverified.')

    fig, ax = plt.subplots(figsize=(8, 6))
    image = ax.scatter(redundant.fairface_discoveries, redundant.within_discoveries,
                       s=20+10*redundant.features, c=redundant.adjusted_discoveries, cmap='viridis', alpha=.8)
    ax.set_xlabel('Frozen FairFace significant hypotheses in cluster')
    ax.set_ylabel('Unadjusted within-GMDB significant hypotheses in cluster')
    ax.set_title('Discovery concentration in existing reference feature clusters', loc='left')
    fig.colorbar(image, ax=ax, label='Full adjusted significant hypotheses')
    fig.subplots_adjust(bottom=.18)
    save(fig, 'redundancy_discovery_concentration', 'Marker area scales with cluster size. Reference-only clusters do not redefine BH families or count independent phenotypes.')

    fig, axes = plt.subplots(4, 1, figsize=(9, 9))
    for ax, (syndrome, feature) in zip(axes, f.HIGHLIGHTS):
        group = age.loc[age.syndrome.eq(syndrome) & age.feature.eq(feature)].set_index('stratum').reindex([s[0] for s in STRATA])
        for j, row in enumerate(group.itertuples()):
            if row.status == 'ok':
                ax.hlines(j, row.ci_95_lower, row.ci_95_upper, color='#2166ac', linewidth=2)
                ax.plot(row.hedges_g, j, 'o', color='#2166ac')
            else:
                ax.text(.02, j, f'Unsupported: n={row.n_a}/{row.n_b}', transform=ax.get_yaxis_transform(), va='center')
        ax.set_yticks([0, 1], ['[3,6)', '[6,10)'])
        ax.set_ylim(-.5, 1.5)
        ax.axvline(0, color='.7', linestyle='--')
        ax.set_title(f'{label(syndrome)}: {feature.replace("_", " ")}', fontsize=9, loc='left')
    axes[-1].set_xlabel("Hedges' g (target minus all other GMDB in age stratum)")
    fig.subplots_adjust(left=.16, right=.98, bottom=.12, top=.96, hspace=.7)
    save(fig, 'age_stratified_highlights', 'Marginal bootstrap intervals. Separate global BH family: 12,250 planned tests. Stratum differences are not interaction tests.')


def report(out, adjusted, pairs, concord, redundant, thresholds, age, age_summary, primary, iterations):
    overall = concord.loc[concord.level.eq('overall')]
    model_counts = adjusted.groupby('status').size().to_dict()
    age_counts = age.groupby(['stratum', 'status']).size().to_dict()
    pair_summary = pairs.groupby(['syndrome_a', 'syndrome_b'], sort=False).agg(
        n_a=('total_a', 'first'), n_b=('total_b', 'first'), supported=('status', lambda x: x.eq('ok').sum()),
        significant=('significant', 'sum'), median_absolute_g=('hedges_g', lambda x: x.abs().median())).reset_index()
    pair_summary.to_csv(out/'pairwise_summary.csv', index=False)
    pair_summary.loc[pair_summary.syndrome_a.isin(ANCHORS) & pair_summary.syndrome_b.isin(ANCHORS)].to_csv(
        out/'anchor_pairwise_summary.csv', index=False)
    pairs.loc[pairs.syndrome_a.isin(ANCHORS) & pairs.syndrome_b.isin(ANCHORS) &
              pairs.feature.isin([v for _, v in f.HIGHLIGHTS])].to_csv(out/'anchor_pairwise_highlights.csv', index=False)
    highlights = pd.concat([adjusted.loc[adjusted.syndrome.eq(s) & adjusted.feature.eq(v)] for s, v in f.HIGHLIGHTS])
    highlights.to_csv(out/'full_screen_adjusted_highlights.csv', index=False)
    adjusted.merge(f.read(f.ROOT/'publication/supplementary/Supplementary_Table_S6_exact_feature_definitions.csv')[
        ['feature', 'anatomical_region', 'feature_family', 'units']], on='feature', validate='many_to_one').to_csv(
        out/'adjusted_annotated.csv', index=False)
    for column in ['syndrome', 'anatomical_region']:
        annotated = adjusted.merge(f.read(f.ROOT/'publication/supplementary/Supplementary_Table_S6_exact_feature_definitions.csv')[
            ['feature', 'anatomical_region']], on='feature', validate='many_to_one')
        annotated.groupby(column).agg(planned=('feature', 'size'), supported=('status', lambda s: s.eq('ok').sum()),
                                      significant=('significant', 'sum')).reset_index().to_csv(out/f'adjusted_by_{column}.csv', index=False)
    flagged = pd.concat([adjusted.assign(analysis='adjusted'), pairs.assign(analysis='pairwise'), age.assign(analysis='age_strata')], ignore_index=True)
    flagged.loc[flagged.status.ne('ok') | flagged.warning_flags.fillna('').ne('')].to_csv(out/'unsupported_and_flagged.csv', index=False)
    unknown = int(primary.sex_label.isna().sum())
    bias_flags = []
    for row in concord.itertuples():
        for metric in METRICS:
            if getattr(row, f'{metric}_bias_exceeds_se') or getattr(row, f'{metric}_interval_excludes_point'):
                bias_flags.append({'comparator': row.comparator, 'level': row.level, 'group': row.group,
                                   'metric': metric, 'estimate': getattr(row, metric),
                                   'bootstrap_bias': getattr(row, f'{metric}_bootstrap_bias'),
                                   'bootstrap_se': getattr(row, f'{metric}_bootstrap_se'),
                                   'ci_low': getattr(row, f'{metric}_ci_low'), 'ci_high': getattr(row, f'{metric}_ci_high'),
                                   'interval_excludes_point': getattr(row, f'{metric}_interval_excludes_point')})
    pd.DataFrame(bias_flags).to_csv(out/'concordance_bootstrap_warnings.csv', index=False)
    lines = ['# Historical feature analysis extension', '',
             'Aggregate derivatives for author/statistician review only. No manuscript or frozen result was edited. '
             'These results do not establish data-use authorization, clinical validity or diagnostic specificity.', '',
             '## Full age/sex-adjusted screen', '',
             f'{len(adjusted):,} planned models; {model_counts.get("ok", 0):,} supported; '
             f'{model_counts.get("unsupported", 0):,} unsupported; {int(adjusted.significant.sum()):,} global BH q<0.05.', '',
             f'The cohort contains {len(primary):,} unique recorded IDs across {primary.disease.nunique()} labels. '
             f'{unknown} selected records lack a usable recorded sex label. Complete cases are selected per feature. '
             'The minimal age-range support criterion is not proof of balance or a correctly specified common linear slope. '
             'Coefficients are feature-unit effects with HC3 normal-reference marginal intervals; they are not Hedges\' g. '
             'The existing four adjusted models retain their four-model q-values separately. Regression tests a conditional '
             'linear mean coefficient, while the unadjusted screen uses rank tests. Different discovery counts cannot isolate '
             'the effect of adjustment, and normal-reference HC3 inference is approximate for small target groups.', '',
             '## Targeted pairwise family', '',
             f'{len(pair_summary)} unique pairs, {len(pairs):,} planned tests; {int(pairs.status.eq("ok").sum()):,} supported; '
             f'{int(pairs.significant.sum()):,} global BH q<0.05. Each pair contains two disjoint recorded-ID groups. '
             f'Marginal g intervals use {iterations:,} patient-row bootstrap draws. Positive g denotes syndrome_a minus syndrome_b. '
             'These contrasts are unadjusted. All four anchors were compared to all other historical main-eligible syndromes before examining new significance. '
             'The full family includes 6 anchor-anchor pairs and 180 anchor-other pairs; results are not specificity estimates.', '',
             '## Extended concordance', '',
             '| Comparator | Direction agreement | Pearson r | Spearman rho | Both significant | FairFace only | GMDB only | Median absolute attenuation |',
             '|---|---:|---:|---:|---:|---:|---:|---:|']
    for row in overall.itertuples():
        lines.append(f'| {row.comparator} | {row.direction_agreement:.4f} | {row.pearson_r:.4f} | '
                     f'{row.spearman_rho:.4f} | {row.both} | {row.fairface_only} | {row.within_only} | {row.median_absolute_attenuation:.4f} |')
    lines += ['', 'Anatomical region and all 89 feature-family summaries with joint-bootstrap marginal confidence intervals '
              'are in concordance_summary.csv. Intervals retain shared target/control observations by resampling GMDB rows '
              'within syndrome and using the same target/reference draws for all contrasts. Syndrome sizes and photo selection '
              'are fixed. Significance overlap uses original global families, with no new significance test for concordance. '
              'Attenuation is |FairFace g| - |GMDB g|; positive values indicate a smaller standardized effect. It is not '
              'a measure of removal of technical bias. Correlated features and reused controls preclude treating test entries as independent.', '',
              f'**Bootstrap qualification:** {len(bias_flags)} group/metric summaries have absolute estimated bootstrap bias '
              'larger than their bootstrap SE or have a percentile interval excluding the observed estimate. '
              'concordance_bootstrap_warnings.csv lists each case; all summaries export bootstrap mean, SE, bias and exclusion flags. '
              'The 95% percentile intervals characterize conditional resampling uncertainty; their nominal population coverage '
              'is unverified and large-bias cases should not be presented as calibrated confidence statements. '
              'No independent-hypothesis/binomial confidence interval or silently recentered interval replaces them.', '',
              '## Redundancy', '',
              f'The existing reference matrix and {len(redundant)} clusters are reused, not re-estimated. '
              f'{int(thresholds.loc[thresholds.absolute_rho_threshold.eq(.8), "pairs"].iloc[0])} of 7,750 feature pairs have |rho|>=0.8. '
              f'{int(redundant.within_discoveries.gt(0).sum())} clusters contain unadjusted within-GMDB discoveries; '
              f'{int(redundant.adjusted_discoveries.gt(0).sum())} contain full-screen adjusted discoveries. '
              'Clusters may span anatomy, depend on the reference population and shared geometry, and do not redefine multiplicity.', '',
              '## Age-stratified sensitivity', '',
              '| Stratum | Supported tests | Unsupported tests | Significant tests |', '|---|---:|---:|---:|']
    for name, _, _ in STRATA:
        group = age.loc[age.stratum.eq(name)]
        lines.append(f'| {name} | {age_counts.get((name, "ok"), 0)} | {age_counts.get((name, "unsupported"), 0)} | {int(group.significant.sum())} |')
    lines += ['', f'{int(age_summary.features_supported_both.gt(0).sum())} of 49 target syndromes have estimable features in both strata. '
              'All sparse tests remain in a single separate 12,250-test BH denominator. These newly estimated age-specific '
              'contrasts are not the earlier threshold re-correction. Adequate counts do not imply power; sex is unadjusted '
              'within strata. Cross-stratum changes are not interaction tests, longitudinal trajectories, or developmental validation.', '',
              '## Failures, warnings and limitations', '',
              'unsupported_and_flagged.csv records every failed support check and warning. numerical_warnings.json records '
              'runtime warnings, including empty when none occurred. execution_status.json records failure tracebacks and exit '
              'status. Bootstrap finite counts accompany each interval; low-finite intervals are flagged. Missing covariate '
              'counts overlap missing feature counts and must not be added as mutually exclusive exclusions.', '',
              'Recorded-ID separation cannot establish biological independence or remove familial clustering. FairFace uses '
              'distinct paths, without verified participant linkage. No raw landmarks, clinical confirmation, calibrated image '
              'scale, validated expression/occlusion labels or clinical measurement annotations were created. No clinical '
              'validation or redistribution clearance is inferred. Covariate models omit ancestry, imaging conditions, nonlinear '
              'age and interactions. Pairwise and strata families are exploratory; family boundaries do not permit choosing '
              'the best q-value across analyses. Marginal intervals are not simultaneous or multiplicity-adjusted. '
              'BH dependence assumptions are not verified; q-values do not guarantee the realized false-discovery fraction.', '',
              '## Evidence and reuse', '',
              'analysis_plan.json precedes fits and records the ordered target/pair list, methods hash, family sizes and seeds. '
              'provenance.json retains source/input/output hashes, versions and unavailable-submodule state. Existing outputs '
              'are read only. Full adjusted, pairwise, age-stratum, concordance and redundancy tables are aggregate results. '
              'figures/ contains editable SVGs and 350-dpi PNGs. All output remains ignored and local-only.', '']
    (out/'RESULTS_FOR_AUTHOR_REVIEW.md').write_text('\n'.join(lines))
    return {'adjusted_supported': model_counts.get('ok', 0), 'adjusted_significant': int(adjusted.significant.sum()),
            'pairwise_pairs': len(pair_summary), 'pairwise_significant': int(pairs.significant.sum()),
            'age_supported_tests': int(age.status.eq('ok').sum()), 'age_significant': int(age.significant.sum()),
            'age_syndromes_supported_both': int(age_summary.features_supported_both.gt(0).sum())}


def version_record():
    return {name: importlib.metadata.version(name) for name in ['numpy', 'pandas', 'scipy', 'statsmodels', 'matplotlib', 'pytest']}


def main():
    parser = argparse.ArgumentParser(description='New exploratory analyses of the frozen historical 125-feature cohort')
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--bootstrap', type=int, default=1000)
    parser.add_argument('--concordance-bootstrap', type=int, default=1000)
    args = parser.parse_args()
    f.check(args.bootstrap >= 100 and args.concordance_bootstrap >= 100, 'Use at least 100 resamples')
    out = args.output_dir.resolve()
    local = (f.OUT/'local_only').resolve()
    f.check(out.is_relative_to(local), 'Output must stay under reviewer_followup_2026/local_only')
    f.check(not out.exists(), 'Refusing to overwrite existing outputs')
    out.mkdir(parents=True)
    started = time.time()
    status = {'started_utc': datetime.now(timezone.utc).isoformat(), 'command': sys.argv, 'status': 'running'}
    (out/'execution_status.json').write_text(json.dumps(status, indent=2))
    captured = []
    try:
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter('always')
            f.INPUTS.clear()
            _, _, primary, _, reference, frozen, targets, features = f.load()
            require_unique(primary)
            mapping = f.read(f.ROOT/'publication/supplementary/Supplementary_Table_S6_exact_feature_definitions.csv')
            f.check(mapping.feature.tolist() == features, 'Dictionary schema mismatch')
            existing = {'all': f.read(f.OUT/'within_gmdb_results.csv'),
                        'eligible49': f.read(f.OUT/'within_gmdb_eligible49_sensitivity.csv')}
            pairs = planned_pairs(targets)
            plan = {'contract_sha256': f.sha(Path(__file__).with_name('EXTENSION_METHODS.md')), 'seed': SEED,
                    'targets': targets, 'features': features, 'ordered_pairs': pairs,
                    'planned_families': {'adjusted': len(targets)*len(features), 'pairwise': len(pairs)*len(features),
                                         'age_strata': 2*len(targets)*len(features)},
                    'bootstrap_iterations': args.bootstrap, 'concordance_bootstrap_iterations': args.concordance_bootstrap,
                    'age_strata': STRATA, 'fixed_utc': datetime.now(timezone.utc).isoformat()}
            (out/'analysis_plan.json').write_text(json.dumps(plan, indent=2))
            adjusted = full_adjusted(primary, targets, features)
            adjusted.to_csv(out/'full_adjusted_screen.csv', index=False)
            pair_results = pairwise(primary, targets, features, args.bootstrap)
            pair_results.to_csv(out/'targeted_pairwise.csv', index=False)
            age = age_stratified(primary, targets, features, args.bootstrap)
            age.to_csv(out/'age_stratified.csv', index=False)
            age_summary = age_agreement(age)
            age_summary.to_csv(out/'age_stratum_agreement.csv', index=False)
            concord, entries = concordance(primary, reference, frozen, targets, features, mapping, existing, args.concordance_bootstrap)
            concord.to_csv(out/'concordance_summary.csv', index=False)
            entries.to_csv(out/'concordance_entries.csv', index=False)
            redundant, thresholds = redundancy(features, mapping, entries.loc[entries.comparator.eq('all')], adjusted)
            redundant.to_csv(out/'redundancy_cluster_discoveries.csv', index=False)
            thresholds.to_csv(out/'redundancy_thresholds.csv', index=False)
            figures(out, adjusted, pair_results, concord, redundant, age, mapping)
            summary = report(out, adjusted, pair_results, concord, redundant, thresholds, age, age_summary, primary, args.bootstrap)
            (out/'SUMMARY.json').write_text(json.dumps(summary, indent=2))
            print(json.dumps(summary, indent=2), flush=True)
        status['status'], status['exit_code'] = 'completed', 0
    except Exception:
        status['status'], status['exit_code'], status['traceback'] = 'failed', 1, traceback.format_exc()
        raise
    finally:
        warning_rows = [{'category': w.category.__name__, 'message': str(w.message),
                         'source_file': Path(w.filename).name, 'line': w.lineno} for w in captured]
        (out/'numerical_warnings.json').write_text(json.dumps(warning_rows, indent=2))
        provenance = {'git_head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=f.ROOT, text=True).strip(),
                      'git_dirty': bool(subprocess.check_output(['git', 'status', '--porcelain'], cwd=f.ROOT, text=True).strip()),
                      'python': sys.version, 'versions': version_record(), 'seed': SEED, **f.facekit_provenance(),
                      'scripts_sha256': {str(p.relative_to(f.ROOT)): f.sha(p) for p in [
                          Path(__file__), Path(f.__file__), Path(__file__).with_name('EXTENSION_METHODS.md'),
                          Path(__file__).with_name('test_extended_analysis.py'), Path(__file__).with_name('requirements.txt'),
                          f.ROOT/'stats_testing/prepare_cohorts.py', f.ROOT/'stats_testing/cohort_audit.py',
                          f.ROOT/'stats_testing/ancestry_analysis.py']},
                      'inputs_sha256': {str(p.relative_to(f.ROOT)): f.sha(p) for p in sorted(f.INPUTS) if p.is_file()},
                      'missing_inputs': [str(p.relative_to(f.ROOT)) for p in sorted(f.INPUTS) if not p.is_file()],
                      'output_policy': 'Local aggregate derivatives only; no authorization or clinical validity inferred'}
        (out/'provenance.json').write_text(json.dumps(provenance, indent=2))
        status.update(finished_utc=datetime.now(timezone.utc).isoformat(), elapsed_seconds=time.time()-started,
                      warning_count=len(warning_rows))
        (out/'execution_status.json').write_text(json.dumps(status, indent=2))
        outputs = {str(p.relative_to(out)): f.sha(p) for p in sorted(out.rglob('*')) if p.is_file()}
        (out/'output_hashes.json').write_text(json.dumps(outputs, indent=2))


if __name__ == '__main__':
    raise SystemExit('Use revision_results/rebuild_120_20261028/extensions.py for the validated current study')

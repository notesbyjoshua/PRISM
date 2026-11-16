import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import rankdata, spearmanr

from .common import (ROOT, DEFAULT_OUT, KEY, SEED, ResearchContractError, bh, canonicalize, case_rows, corrected,
                     effect, features, load_private, offline, private_output, protect_baseline, require, sha,
                     u_test, write_aggregate, write_json)

PROXIES = ['abs_yaw', 'abs_pitch', 'abs_roll', 'pose_score', 'mouth_opening', 'eye_fissure_height_mean']


def quota(n, weights, capacity):
    weights, capacity = np.asarray(weights, float), np.asarray(capacity, int)
    require(capacity.sum() >= n, 'Insufficient disjoint reference-path pool for the declared control design')
    target = n * weights / weights.sum()
    chosen = np.minimum(np.floor(target).astype(int), capacity)
    while chosen.sum() < n:
        eligible = np.flatnonzero(chosen < capacity)
        index = eligible[np.argmax((target - chosen)[eligible])]
        chosen[index] += 1
    return chosen


def control_indices(strata, pitch, sizes, rng, design):
    levels = sorted(set(strata))
    pools, references = [], []
    for level in levels:
        indices = np.flatnonzero(strata == level)
        indices = rng.permutation(indices)
        if design == 'pitch_imbalanced':
            indices = indices[np.argsort(pitch[indices], kind='stable')]
        cut = len(indices) // 2
        references.append(indices[:cut])
        pools.append(rng.permutation(indices[cut:]))
    reference = np.concatenate(references)
    weights = np.array([len(group) for group in references]) / len(reference)
    offsets = np.zeros(len(pools), int)
    groups, imbalance = [], []
    for n in sizes:
        counts = quota(n, weights, [len(group) - offset for group, offset in zip(pools, offsets)])
        group = np.concatenate([pool[offset:offset + count] for pool, offset, count in zip(pools, offsets, counts)])
        offsets += counts
        groups.append(group)
        imbalance.append(np.max(np.abs(counts / n - weights)))
    all_cases = np.concatenate(groups)
    require(len(np.unique(all_cases)) == len(all_cases) and not np.isin(all_cases, reference).any(),
            'Within-source controls leaked paths across case/reference pools')
    return groups, reference, float(max(imbalance))


def negative_controls(real, reference, names, out, repetitions):
    counts = real.disease.value_counts()
    counts = counts.loc[counts.ge(10)].sort_index()
    values = reference[names].to_numpy(float)
    strata = (reference.race.astype(str) + '|' + reference.gender.astype(str)).to_numpy()
    pitch = reference.abs_pitch.to_numpy()
    state_path = out / 'checkpoints/source_negative_controls_v2.json'
    signature = {'code': sha(Path(__file__)), 'reference_shape': list(values.shape), 'seed': SEED,
                 'primary_input_sha256': sha(ROOT / 'ACTIVE_ANALYSIS.json'), 'repetitions_per_design': repetitions,
                 'reference_bytes': sha(ROOT / json.loads((ROOT / 'ACTIVE_ANALYSIS.json').read_text())['reference_master_csv']),
                 'family': 5760, 'design_version': 'v2_randomized_within_source_pools'}
    state = json.loads(state_path.read_text()) if state_path.exists() else {'signature': signature, 'rows': []}
    require(state['signature'] == signature, 'Changed negative-control signature; select a fresh private run')
    rows = state['rows']
    completed = {(row['design'], row['repetition']) for row in rows}
    for design_index, design in enumerate(['stratified_random', 'pitch_imbalanced']):
        for rep in range(repetitions):
            if (design, rep) in completed:
                continue
            rng = np.random.default_rng([SEED, 10, design_index, rep])
            groups, reference_indices, imbalance = control_indices(strata, pitch, counts.to_numpy(), rng, design)
            comparator = values[reference_indices]
            pvalues, gs = [], []
            for indices in groups:
                _, p = u_test(values[indices], comparator)
                pvalues.extend(p)
                gs.extend(effect(values[indices], comparator))
            qs, gs = bh(pvalues), np.array(gs)
            rows.append({'design': design, 'repetition': rep, 'hypotheses': len(qs), 'reference_paths': len(reference_indices),
                         'pseudo_case_paths': int(counts.sum()), 'distinct_pseudo_groups': len(groups), 'discoveries': int((qs < .05).sum()),
                         'positive_discoveries': int(((qs < .05) & (gs > 0)).sum()),
                         'negative_discoveries': int(((qs < .05) & (gs < 0)).sum()), 'median_absolute_g': float(np.median(np.abs(gs))),
                         'max_source_label_fraction_imbalance': imbalance,
                         'case_mean_abs_pitch': float(pitch[np.concatenate(groups)].mean()),
                         'reference_mean_abs_pitch': float(pitch[reference_indices].mean()),
                         'interpretation': 'within-FairFace workflow sensitivity; not clinical false-positive rate'})
            write_json({'signature': signature, 'rows': rows}, state_path)
            if (rep + 1) % 25 == 0:
                print(f'Negative control {design}: {rep + 1}/{repetitions}', flush=True)
    frame = pd.DataFrame(rows).sort_values(['design', 'repetition'])
    write_aggregate(frame, out / 'aggregate_candidates/P0_negative_control_repetitions.csv')
    summary = frame.groupby('design', sort=True).agg(repetitions=('repetition', 'size'), mean_discoveries=('discoveries', 'mean'),
        median_discoveries=('discoveries', 'median'), discovery_q025=('discoveries', lambda x: x.quantile(.025)),
        discovery_q975=('discoveries', lambda x: x.quantile(.975)), minimum=('discoveries', 'min'), maximum=('discoveries', 'max'),
        repetitions_with_any=('discoveries', lambda x: x.gt(0).sum())).reset_index()
    write_aggregate(summary, out / 'aggregate_candidates/P0_negative_control_summary.csv')
    return summary


def quality_correlations(real, reference, names, out):
    rows = []
    for source, frame in [('GMDB', real), ('FairFace', reference)]:
        groups = [('all', frame)] + [(str(group), part) for group, part in frame.groupby('broad_group', dropna=False) if len(part) >= 100]
        for group, part in groups:
            for proxy in PROXIES:
                for feature in names:
                    supported = proxy != feature and part[proxy].nunique() > 1 and part[feature].nunique() > 1
                    rho, p = spearmanr(part[proxy], part[feature]) if supported else (np.nan, np.nan)
                    rows.append(dict(source=source, recorded_group=group, proxy=proxy, feature=feature, n=len(part), rho=rho,
                                     p_value=p, supported=supported,
                                     reason='observed association' if supported else 'constant or deterministic self-relation'))
    frame = pd.DataFrame(rows)
    frame['exploratory_q'] = np.nan
    for _, indices in frame.groupby(['source', 'recorded_group']).groups.items():
        part = corrected(frame.loc[indices], 'p_value', frame.loc[indices, 'supported'], 'source/group quality correlations; 6 x 120 planned')
        frame.loc[indices, 'exploratory_q'] = part.exploratory_q
    frame['correction_family'] = 'separate source/group family; 720 planned correlations; unsupported self/constant cells retain missing q'
    write_aggregate(frame, out / 'aggregate_candidates/P0_feature_quality_correlations.csv')
    focal = ['eb_lateral_thickness_mean', 'eb_lateral_thickness_l', 'eb_lateral_thickness_r', 'forehead_height', 'chin_height']
    quantiles = []
    for source_index, (source, part) in enumerate([('GMDB', real), ('FairFace', reference)]):
        rng = np.random.default_rng([SEED, 11, source_index])
        for proxy in ['abs_yaw', 'abs_pitch', 'abs_roll']:
            for feature in focal:
                draws = []
                values = part[[proxy, feature]].to_numpy(float)
                for _ in range(200):
                    chosen = values[rng.integers(len(values), size=len(values))]
                    draws.append(spearmanr(chosen[:, 0], chosen[:, 1]).statistic)
                quantiles.append(dict(source=source, proxy=proxy, feature=feature, n=len(part),
                    rho=spearmanr(values[:, 0], values[:, 1]).statistic, resampling_q025=np.quantile(draws, .025),
                    resampling_q975=np.quantile(draws, .975), positive_frequency=np.mean(np.array(draws) > 0),
                    interval_role='conditional descriptive resampling range; not a calibrated clinical CI'))
    write_aggregate(pd.DataFrame(quantiles), out / 'aggregate_candidates/P0_recurrent_quality_resampling.csv')


def restrictions(real, reference, names, matrix, out):
    rows = []
    for threshold in [.75, .5]:
        disease = real.loc[real.pose_score.le(threshold)]
        ref = reference.loc[reference.pose_score.le(threshold)]
        for syndrome in sorted(matrix.syndrome.unique()):
            x = disease.loc[disease.disease.eq(syndrome), names].to_numpy(float)
            y = ref[names].to_numpy(float)
            supported = min(len(x), len(y)) >= 10
            if supported:
                _, p = u_test(x, y)
                g = effect(x, y)
            else:
                p = g = np.full(120, np.nan)
            rows.extend(dict(threshold=threshold, syndrome=syndrome, feature=feature, n_disease=len(x), n_reference=len(y),
                             supported=supported, hedges_g=g[j], p_value=p[j], reason='supported' if supported else 'fewer than 10 after restriction')
                        for j, feature in enumerate(names))
    frame = canonicalize(pd.DataFrame(rows))
    frame['exploratory_q'] = np.nan
    for threshold, indices in frame.groupby('threshold').groups.items():
        part = corrected(frame.loc[indices], 'p_value', frame.loc[indices, 'supported'], f'quality restriction {threshold}; 5760 planned')
        frame.loc[indices, 'exploratory_q'] = part.exploratory_q
    frame['correction_family'] = 'separate 5760-cell exploratory family per pose restriction; unsupported display q missing'
    write_aggregate(frame, out / 'aggregate_candidates/P0_quality_restrictions.csv')
    return frame


def quality_adjusted(real, reference, matrix, out):
    rows = []
    for row in case_rows(matrix).itertuples():
        for comparator in ['FairFace', 'other_GMDB']:
            target = real.loc[real.disease.eq(row.syndrome)].copy()
            other = reference.copy() if comparator == 'FairFace' else real.loc[real.disease.ne(row.syndrome)].copy()
            frame = pd.concat([target.assign(target=1), other.assign(target=0)], ignore_index=True).copy()
            frame['female'] = frame.sex_label.map({'Male': 0., 'Female': 1.})
            fields = ['female', 'broad_group', 'pose_yaw', 'pose_pitch', 'pose_roll', row.feature]
            if comparator == 'other_GMDB':
                fields.append('conditional_age_years')
            frame = frame.dropna(subset=fields)
            n, m = int(frame.target.sum()), int(frame.target.eq(0).sum())
            item = dict(syndrome=row.syndrome, feature=row.feature, comparator=comparator, n_target=n, n_comparator=m,
                        beta=np.nan, ci_lower=np.nan, ci_upper=np.nan, p_value=np.nan, supported=False,
                        adjustment='sex + recorded broad group + signed pose' + (' + exact recorded age' if comparator == 'other_GMDB' else ''),
                        inference='source/case effect not separately identifiable; adjustment does not remove source bias')
            if min(n, m) >= 10 and m >= 20:
                design = frame[['target', 'female', 'pose_yaw', 'pose_pitch', 'pose_roll']].copy()
                if comparator == 'other_GMDB':
                    design['age_centered'] = frame.conditional_age_years - 6
                design = pd.concat([design, pd.get_dummies(frame.broad_group, prefix='group', drop_first=True, dtype=float)], axis=1)
                design = sm.add_constant(design.astype(float), has_constant='add')
                if np.linalg.matrix_rank(design) == design.shape[1]:
                    model = sm.OLS(frame[row.feature].to_numpy(float), design.to_numpy()).fit(cov_type='HC3', use_t=False)
                    ci = model.conf_int()[1]
                    item.update(beta=model.params[1], ci_lower=ci[0], ci_upper=ci[1], p_value=model.pvalues[1],
                                supported=bool(np.isfinite(model.pvalues[1])), maximum_leverage=model.get_influence().hat_matrix_diag.max())
            rows.append(item)
    result = canonicalize(pd.DataFrame(rows))
    result = corrected(result, 'p_value', result.supported, 'predeclared six cases x two comparators; 12 exploratory HC3 models')
    write_aggregate(result, out / 'aggregate_candidates/P0_quality_adjusted_cases.csv')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUT)
    parser.add_argument('--repetitions', type=int, default=100)
    args = parser.parse_args()
    offline()
    out = private_output(args.output_dir)
    protect_baseline(out)
    real, reference, raw, alternate, paths = load_private()
    names = features().feature.tolist()
    matrix = pd.read_csv(out / 'evidence/P0_exhaustive_evidence_matrix.csv', dtype={'canonical_syndrome_id': 'string'})
    summary = negative_controls(real, reference, names, out, args.repetitions)
    quality_correlations(real, reference, names, out)
    restricted = restrictions(real, reference, names, matrix, out)
    quality_adjusted(real, reference, matrix, out)
    report = out / 'P0_source_bias_audit.md'
    text = report.read_text().split('\n## Executed negative controls')[0] + '\n## Executed negative controls and sensitivities\n\n'
    text += ('Two control designs use 48 disjoint pseudo-case groups with the observed primary syndrome sizes (872 FairFace paths), '
             'a shared disjoint reference half, exact recorded race-by-sex strata at splitting, and finite-n approximate quota matching at group allocation. '
             'Paths are not verified people. The deliberate imbalance separates absolute pitch within each label stratum. '
             'Every repetition uses all 120 features and the complete 5,760-hypothesis global BH family. These outcomes estimate workflow sensitivity, '
             'not clinical false positives. Repetitions checkpoint by code, inputs, parameters and seed.\n\n')
    for row in summary.itertuples():
        text += (f'{row.design}: {row.repetitions} repetitions, median discoveries {row.median_discoveries:.6g}, '
                 f'2.5th–97.5th empirical percentiles [{row.discovery_q025:.6g}, {row.discovery_q975:.6g}], '
                 f'{row.repetitions_with_any} repetitions with any discoveries.\n\n')
    text += ('The correlation screen retains a separate 720-cell source/group family. Self-relations/constant proxies have missing estimates and q. '
             'Recurrent-feature resampling intervals are descriptive ranges, not calibrated clinical confidence intervals. '
             'Tighter pose gates have separate planned 5,760-cell exploratory families; unsupported cells retain missing display q, '
             'not evidence for a null. Twelve case-specific HC3 models have their own family; coefficients are in normalized feature units, not g. '
             'They adjust observable covariates only and cannot isolate a biological disease effect from dataset source.\n')
    report.write_text(text)
    protect_baseline(out)
    print('Source-bias controls complete; no primary outputs changed')


if __name__ == '__main__':
    main()

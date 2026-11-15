import argparse
import json
import sys

import numpy as np
import pandas as pd
from scipy.special import gammaln
from scipy.stats import kruskal, mannwhitneyu

from common import (HERE, PRIVATE, RACES, SEED, apportion, cohort_counts, load_sources, matched_frames,
                    private_output, protected_hashes, provenance, sample_reference, select_ml, sha, write_json)


def bh(p_values):
    p = np.asarray(p_values, float)
    if p.ndim != 1 or not np.isfinite(p).all() or ((p < 0) | (p > 1)).any():
        raise ValueError('BH requires every planned hypothesis, with unsupported comparisons set to p=1')
    order = np.argsort(p, kind='stable')
    q = np.empty(len(p))
    q[order] = np.minimum(1, np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1])
    return q


def two_sample(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if min(len(x), len(y)) < 2 or not np.isfinite(np.r_[x, y]).all():
        raise ValueError('Finite two-sample observations with n>=2 required')
    df = len(x) + len(y) - 2
    pooled = ((len(x) - 1) * x.var(ddof=1) + (len(y) - 1) * y.var(ddof=1)) / df
    correction = np.exp(gammaln(df / 2) - .5 * np.log(df / 2) - gammaln((df - 1) / 2))
    g = correction * (x.mean() - y.mean()) / np.sqrt(pooled) if pooled > 0 else np.nan
    ties = len(np.unique(np.r_[x, y])) != len(x) + len(y)
    method = 'exact' if min(len(x), len(y)) <= 8 and not ties else 'asymptotic'
    test = mannwhitneyu(x, y, alternative='two-sided', method=method, use_continuity=True)
    return {'mannwhitney_u': float(test.statistic), 'p': float(test.pvalue), 'hedges_g': float(g),
            'effect_status': 'defined' if pooled > 0 else 'undefined_zero_pooled_variance',
            'mw_method': method, 'has_ties': ties, 'hedges_correction': float(correction)}


def descriptive(values, prefix):
    values = np.asarray(values, float)
    if not len(values):
        return {f'{prefix}_{key}': (0 if key == 'n' else np.nan)
                for key in ['n', 'mean', 'sd', 'median', 'q25', 'q75', 'iqr', 'missing']}
    q25, median, q75 = np.quantile(values, [.25, .5, .75])
    return {f'{prefix}_n': len(values), f'{prefix}_mean': values.mean(), f'{prefix}_sd': values.std(ddof=1),
            f'{prefix}_median': median, f'{prefix}_q25': q25, f'{prefix}_q75': q75,
            f'{prefix}_iqr': q75 - q25, f'{prefix}_missing': int((~np.isfinite(values)).sum())}


def analysis_a(reference, features, definitions):
    groups = {race: reference.loc[reference.canonical_race.eq(race)] for race in RACES}
    rows = []
    for feature in features:
        values = [groups[race][feature].to_numpy(float) for race in RACES]
        supported = min(map(len, values)) >= 20
        constant = supported and np.ptp(np.concatenate(values)) == 0
        if supported and not constant:
            test = kruskal(*values, nan_policy='raise')
            h, p = float(test.statistic), float(test.pvalue)
        else:
            h, p = (0., 1.) if constant else (np.nan, 1.)
        row = {'feature': feature, 'H': h, 'p': p, 'supported': supported,
               'status': 'constant_all_groups' if constant else ('tested' if supported else 'unsupported_reference')}
        for race, value in zip(RACES, values):
            row.update(descriptive(value, race))
        rows.append(row)
    result = pd.DataFrame(rows)
    result['q'] = bh(result.p)
    result['significant'] = result.supported & result.q.lt(.05)
    return result.merge(definitions[['feature', 'formula', 'units', 'anatomical_region']], on='feature', validate='one_to_one')


def primary_tests(disease, reference, counts, features, definitions):
    rows = []
    for cell in counts.loc[counts.eligible].itertuples():
        cases, controls = matched_frames(disease, reference, cell.disease, cell.canonical_race)
        assert cases.patient_id.nunique() == cell.n_current_gmdb == len(cases)
        supported = len(controls) >= 20
        for feature in features:
            x, y = cases[feature].to_numpy(float), controls[feature].to_numpy(float)
            row = {'disease': cell.disease, 'canonical_race': cell.canonical_race, 'feature': feature,
                   'supported': supported, 'support_status': 'tested' if supported else 'reference_n_below_20',
                   **descriptive(x, 'case'), **descriptive(y, 'reference')}
            if supported:
                row.update(two_sample(x, y))
            else:
                row.update({'mannwhitney_u': np.nan, 'p': 1., 'hedges_g': np.nan,
                            'effect_status': 'unsupported_reference', 'mw_method': 'not_tested',
                            'has_ties': False, 'hedges_correction': np.nan})
            rows.append(row)
    result = pd.DataFrame(rows)
    assert len(result) == int(counts.eligible.sum()) * 120
    result['q'] = bh(result.p)
    result['significant'] = result.supported & result.q.lt(.05)
    result['large_effect'] = result.hedges_g.abs().ge(.8)
    result['effect_direction'] = np.select([result.hedges_g.gt(0), result.hedges_g.lt(0), result.hedges_g.eq(0)],
                                           ['case_higher', 'case_lower', 'equal_mean'], default='undefined')
    return result.merge(definitions[['feature', 'formula', 'units', 'anatomical_region']], on='feature', validate='many_to_one')


def covariate_audit(disease, reference):
    rows = []
    for source, frame in [('GMDB', disease), ('FairFace', reference)]:
        for race in RACES:
            group = frame.loc[frame.canonical_race.eq(race)]
            for field in ['conditional_age_years', 'sex_label', 'pose_pitch', 'pose_yaw', 'pose_roll']:
                row = {'source': source, 'canonical_race': race, 'field': field, 'total': len(group)}
                if field not in group:
                    row.update({'available': 0, 'unavailable_reason': 'Only source age category 3-9; exact age unavailable'})
                elif field == 'sex_label':
                    row.update({'available': int(group[field].isin(['Female', 'Male']).sum()),
                                'label_counts': json.dumps(group[field].value_counts(dropna=False).to_dict())})
                else:
                    values = pd.to_numeric(group[field], errors='coerce').dropna()
                    row.update({'available': len(values), **descriptive(values, 'value')})
                rows.append(row)
    return pd.DataFrame(rows)


def summary_tables(primary):
    rows = []
    for (disease, race), frame in primary.groupby(['disease', 'canonical_race'], sort=True):
        significant = frame.loc[frame.significant]
        top = significant.sort_values(['q', 'feature']).head(3).feature.tolist()
        rows.append({'disease': disease, 'canonical_race': race, 'n_patients': int(frame.case_n.iloc[0]),
                     'n_reference': int(frame.reference_n.iloc[0]), 'planned_tests': len(frame),
                     'supported_tests': int(frame.supported.sum()), 'significant_features': len(significant),
                     'positive_effects': int(significant.hedges_g.gt(0).sum()),
                     'negative_effects': int(significant.hedges_g.lt(0).sum()),
                     'significant_large_effects': int(significant.large_effect.sum()),
                     'top_features_by_q_unassessed_stability': '|'.join(top)})
    regions = primary.groupby(['canonical_race', 'anatomical_region']).agg(
        planned=('feature', 'size'), tested=('supported', 'sum'), significant=('significant', 'sum')).reset_index()
    regions['significant_rate'] = regions.significant / regions.tested.replace(0, np.nan)
    return pd.DataFrame(rows), regions


def figures(output, a, primary, cells, regions):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 4.5), layout='constrained')
    ranked = a.sort_values(['q', 'feature']).head(30).iloc[::-1]
    ax.barh(ranked.feature, -np.log10(ranked.q.clip(lower=1e-300)), color='#28748c')
    ax.axvline(-np.log10(.05), color='#b54143', linestyle='--')
    ax.set(xlabel='−log10(BH q), 120-feature family', title='FairFace recorded-category associations (top 30)')
    ax.tick_params(axis='y', labelsize=7)
    fig.savefig(output / 'analysis_A.png', dpi=180)
    plt.close(fig)
    pivot = regions.pivot(index='anatomical_region', columns='canonical_race', values='significant_rate')
    ax = pivot.plot.bar(figsize=(10, 5), color=['#28748c', '#de9f49', '#7565a1'], rot=35)
    ax.set(ylabel='Significant / supported planned tests', title='Within-category discovery rates by region', xlabel='')
    ax.figure.tight_layout()
    ax.figure.savefig(output / 'anatomical_region_rates.png', dpi=180)
    plt.close(ax.figure)
    fig, ax = plt.subplots(figsize=(8, 4), layout='constrained')
    for race in RACES:
        values = cells.loc[cells.canonical_race.eq(race), 'significant_features']
        ax.hist(values, bins=np.arange(0, 125, 5), alpha=.55, label=f'{race} ({len(values)} cells)')
    ax.set(xlabel='Significant features per eligible syndrome-category cell', ylabel='Cells')
    ax.legend()
    fig.savefig(output / 'discovery_count_distribution.png', dpi=180)
    plt.close(fig)
    selected = primary.loc[primary.significant].sort_values(['q', 'disease', 'canonical_race', 'feature']).head(24)
    if len(selected):
        features = selected.feature.drop_duplicates().tolist()[:12]
        keys = selected[['disease', 'canonical_race']].drop_duplicates().head(14)
        frame = primary.merge(keys, on=['disease', 'canonical_race'])
        frame = frame.loc[frame.feature.isin(features)].copy()
        frame['cell'] = frame.disease + ' / ' + frame.canonical_race
        matrix = frame.pivot(index='cell', columns='feature', values='hedges_g').reindex(columns=features)
        significance = frame.pivot(index='cell', columns='feature', values='significant').reindex_like(matrix)
        fig, ax = plt.subplots(figsize=(12, max(4, .4 * len(matrix))), layout='constrained')
        values = ax.imshow(matrix, aspect='auto', cmap='RdBu_r', vmin=-3, vmax=3)
        ax.set_xticks(range(len(matrix.columns)), matrix.columns, rotation=50, ha='right', fontsize=8)
        labels = []
        for value in matrix.index:
            disease_id, race = value.rsplit(' / ', 1)
            name = disease_id.split('__')[0]
            labels.append((name if len(name) <= 48 else name[:45] + '…') + ' / ' + race)
        ax.set_yticks(range(len(matrix)), labels, fontsize=7)
        for i, j in np.argwhere(significance.to_numpy()):
            ax.text(j, i, '*', ha='center', va='center', color='black')
        fig.colorbar(values, ax=ax, label="Hedges' g (case − reference)")
        ax.set_title('Selected associations ranked by q; * global BH q < .05; stability untested\nDisplay clipped at ±3')
        fig.savefig(output / 'selected_effect_heatmap.png', dpi=180)
        plt.close(fig)


def selection_audit(output, counts):
    global_top, selected = select_ml(counts)
    columns = ['disease', 'canonical_race', 'n_current_gmdb', 'umbrella_label']
    global_top[columns].to_csv(output / 'global_top_five_cells.csv', index=False)
    selected[columns].to_csv(output / 'ml_selected_classes.csv', index=False)
    lines = ['# ML selection audit', '', 'Selection uses current patient counts only, before feature tests or fitting.',
             'Global ties: descending n, exact disease ID, then category name. Within a category: descending n, disease ID.',
             'Prefer five classes each n>=15; otherwise require each n>=10. Choose maximum total, then maximum minimum n, then category name.',
             '', '## Global five largest cells', '']
    lines += [f'- {r.disease} / {r.canonical_race}: {r.n_current_gmdb}' for r in global_top.itertuples()]
    lines += ['', '## Selected five different syndromes within one category', '']
    lines += [f'- {r.disease} / {r.canonical_race}: {r.n_current_gmdb}' for r in selected.itertuples()]
    if selected.empty:
        lines += ['', 'BLOCKED: No common category has five distinct syndromes with at least ten patients each.',
                  'Fallback: retain descriptive cell ranking; do not train a mixed-category five-class model.']
    else:
        lines += ['', 'COMPLETED: feasible selection. The global ranking mixes categories; the classifier uses the shared-category ranking.',
                  'No diagnoses were merged. Mucopolysaccharidoses is explicitly flagged as an umbrella wherever present.',
                  'Other source labels may hide subtypes; molecular homogeneity is unverified.']
    (output / 'ml_selection_audit.md').write_text('\n'.join(lines) + '\n')
    return selected


def main():
    parser = argparse.ArgumentParser(description='Private experimental race-stratified primary and demographic analyses')
    parser.add_argument('--seed', type=int, default=SEED)
    parser.add_argument('--output-dir', type=str, default=str(PRIVATE))
    parser.add_argument('--reference-size', type=int, default=3409)
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    output = private_output(args.output_dir)
    before = protected_hashes()
    try:
        disease, reference, features, definitions, audit, sources = load_sources()
    except (FileNotFoundError, ValueError, AssertionError) as error:
        write_json(output / 'statistics_status.json', {'status': 'BLOCKED', 'reason': str(error)})
        print(str(error), file=sys.stderr)
        return 2
    counts = cohort_counts(disease)
    counts.to_csv(output / 'cohort_and_race_audit.csv', index=False)
    selected = selection_audit(output, counts)
    eligible = counts.loc[counts.eligible]
    composition = eligible.groupby('canonical_race').n_current_gmdb.sum().reindex(RACES, fill_value=0).to_dict()
    pool = reference.canonical_race.value_counts().reindex(RACES, fill_value=0).to_dict()
    size, quotas = apportion(composition, pool, args.reference_size)
    summary = {'status': 'AUDIT_COMPLETE', 'available_sources': True, 'eligible_cells': len(eligible),
               'eligible_syndromes': int(eligible.disease.nunique()), 'eligible_patients': int(eligible.n_current_gmdb.sum()),
               'planned_primary_family': len(eligible) * 120, 'analysis_A_family': 120,
               'composition': composition, 'pool': pool, 'reference_size': size, 'quotas': quotas,
               'ml_feasible': not selected.empty, 'sensitivity': 'NOT RUN: user deferred'}
    write_json(output / 'source_cohort_audit.json', {**audit, **summary})
    print(json.dumps(summary, indent=2), flush=True)
    manifest = provenance(sources, {'seed': args.seed, 'maximum_reference_size': args.reference_size,
                                    'case_minimum': 5, 'reference_minimum': 20, 'alpha': .05,
                                    'output_dir': str(output), 'dry_run': args.dry_run})
    manifest.update({'features': features, 'eligibility': summary, 'exclusions': audit,
                     'crosswalk': {'GMDB': audit['gmdb_source_labels'], 'FairFace': audit['fairface_source_labels_after_content_qc'],
                                  'canonical_groups': list(RACES),
                                  'broad_Asian_includes': ['East Asian', 'Southeast Asian', 'Indian', 'Middle Eastern'],
                                  'decision': 'Retain documented inclusive Asian crosswalk; narrower definition deferred by user'},
                     'hypothesis_rules': 'One global primary BH family, all eligible cells x 120; unsupported references p=1; no missing accepted measurements',
                     'mw_rule': 'Per-feature SciPy auto-equivalent: exact if minimum n<=8 and no ties, otherwise asymptotic with tie and continuity corrections',
                     'hedges_rule': 'Exact gamma-function J(df) times pooled-SD standardized mean difference, case minus reference',
                     'composition_denominator': 'Unique GMDB patients in eligible syndrome-category cells only',
                     'content_qc_rule': 'Exclude all conflicting-original-race hashes; keep lexical-first file for other byte-identical FairFace images',
                     'protected_before': before})
    write_json(output / ('dry_run_manifest.json' if args.dry_run else 'statistics_manifest.json'), manifest)
    if args.dry_run:
        assert protected_hashes() == before
        return 0
    definitions.to_csv(output / 'exact_120_feature_definitions.csv', index=False)
    sampled = sample_reference(reference, quotas, args.seed)
    sampled[['image_id', 'file', 'image_sha256', 'race', 'canonical_race']].to_csv(
        output / 'PRIVATE_reference_selection.csv', index=False)
    denominator = sum(composition.values())
    sampling = {'X': size, 'seed': args.seed, 'maximum_X': args.reference_size,
                'rule': 'Largest exactly feasible integer size; largest remainder; alphabetical category ties; no replacement',
                'denominator': 'Unique patients belonging to eligible syndrome-category cells',
                'quotas': quotas, 'pool_sizes': pool,
                'target_proportions': {r: composition[r] / denominator for r in RACES},
                'achieved_proportions': {r: quotas[r] / size for r in RACES},
                'proportion_discrepancy': {r: quotas[r] / size - composition[r] / denominator for r in RACES},
                'all_1845_context_composition': disease.canonical_race.fillna('Excluded').value_counts().to_dict(),
                'selection_sha256': sha(output / 'PRIVATE_reference_selection.csv'),
                'unique_paths': bool(sampled.file.is_unique), 'unique_content': bool(sampled.image_sha256.is_unique)}
    write_json(output / 'fairface_sampling_summary.json', sampling)
    a = analysis_a(reference, features, definitions)
    a.to_csv(output / 'analysis_A_all_120.csv', index=False)
    a.loc[a.significant].sort_values(['q', 'feature']).to_csv(output / 'analysis_A_significant.csv', index=False)
    primary = primary_tests(disease, sampled, counts, features, definitions)
    primary.to_csv(output / 'primary_all_tests.csv', index=False)
    primary.loc[primary.significant].sort_values(['q', 'disease', 'canonical_race', 'feature']).to_csv(
        output / 'primary_significant_aggregate.csv', index=False)
    cells, regions = summary_tables(primary)
    cells.to_csv(output / 'syndrome_race_summary.csv', index=False)
    regions.to_csv(output / 'anatomical_region_summary.csv', index=False)
    covariate_audit(disease, reference).to_csv(output / 'covariate_availability_descriptive.csv', index=False)
    figures(output, a, primary, cells, regions)
    assert protected_hashes() == before, 'Protected analysis or manuscript changed during execution'
    status = {**summary, 'status': 'COMPLETED', 'analysis_A_significant': int(a.significant.sum()),
              'primary_significant': int(primary.significant.sum()),
              'primary_significant_large_effect': int((primary.significant & primary.large_effect).sum()),
              'protected_files_unchanged': True}
    write_json(output / 'statistics_status.json', status)
    output_names = ['cohort_and_race_audit.csv', 'global_top_five_cells.csv', 'ml_selected_classes.csv',
                    'ml_selection_audit.md', 'source_cohort_audit.json', 'exact_120_feature_definitions.csv',
                    'PRIVATE_reference_selection.csv', 'fairface_sampling_summary.json', 'analysis_A_all_120.csv',
                    'analysis_A_significant.csv', 'primary_all_tests.csv', 'primary_significant_aggregate.csv',
                    'syndrome_race_summary.csv', 'anatomical_region_summary.csv', 'covariate_availability_descriptive.csv',
                    'analysis_A.png', 'anatomical_region_rates.png', 'discovery_count_distribution.png',
                    'selected_effect_heatmap.png', 'statistics_status.json']
    manifest['output_sha256'] = {name: sha(output / name) for name in output_names if (output / name).is_file()}
    manifest['completed'] = status
    write_json(output / 'statistics_manifest.json', manifest)
    print(json.dumps(status, indent=2), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

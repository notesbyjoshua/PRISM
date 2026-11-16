import argparse
import hashlib
import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from revision_results.race_stratified_20261115 import common

OUTPUT = ROOT / 'local_review/race_stratified_literature_20261115'
BASELINE = ROOT / 'local_review/race_stratified_20261115'


def load_primary():
    manifest = json.loads((BASELINE / 'statistics_manifest.json').read_text())
    path = BASELINE / 'primary_all_tests.csv'
    if manifest['study_id'] != 'PRISM_RACE_STRATIFIED_20261115' or common.sha(path) != manifest['output_sha256'][path.name]:
        raise ValueError('Verified new main-analysis output required')
    frame = pd.read_csv(path, float_precision='round_trip')
    features = manifest['features']
    if len(features) != 120 or len(frame) != manifest['eligibility']['planned_primary_family']:
        raise ValueError('Exact planned 120-feature primary family required')
    if frame[['disease', 'canonical_race', 'feature']].duplicated().any():
        raise ValueError('Duplicate association keys')
    return frame, manifest


def comparisons(primary, tolerance=1e-12):
    rows = []
    for (disease, feature), group in primary.groupby(['disease', 'feature'], sort=True):
        group = group.loc[group.supported & np.isfinite(group.hedges_g)].sort_values('canonical_race')
        for a, b in itertools.combinations(group.itertuples(), 2):
            sign_a = int(np.sign(a.hedges_g)) if abs(a.hedges_g) > tolerance else 0
            sign_b = int(np.sign(b.hedges_g)) if abs(b.hedges_g) > tolerance else 0
            opposite = sign_a * sign_b < 0
            significance = int(a.q < .05) + int(b.q < .05)
            rows.append({'disease': disease, 'feature': feature, 'race_1': a.canonical_race, 'race_2': b.canonical_race,
                         'g_1': a.hedges_g, 'g_2': b.hedges_g, 'p_1': a.p, 'p_2': b.p, 'q_1': a.q, 'q_2': b.q,
                         'case_n_1': a.case_n, 'case_n_2': b.case_n,
                         'reference_n_1': a.reference_n, 'reference_n_2': b.reference_n,
                         'case_mean_1': a.case_mean, 'case_mean_2': b.case_mean,
                         'reference_mean_1': a.reference_mean, 'reference_mean_2': b.reference_mean,
                         'case_minus_reference_mean_1': a.case_mean - a.reference_mean,
                         'case_minus_reference_mean_2': b.case_mean - b.reference_mean,
                         'g_difference_1_minus_2_descriptive': a.hedges_g - b.hedges_g,
                         'case_mean_difference_1_minus_2': a.case_mean - b.case_mean,
                         'reference_mean_difference_1_minus_2': a.reference_mean - b.reference_mean,
                         'opposite_g_sign': opposite,
                         'significance_tier': ['neither_significant', 'one_significant', 'both_significant'][significance],
                         'both_large_effect': abs(a.hedges_g) >= .8 and abs(b.hedges_g) >= .8,
                         'anatomical_region': a.anatomical_region, 'formula': a.formula, 'units': a.units,
                         'interpretation': 'descriptive sign comparison; no syndrome-by-category interaction test'})
    return pd.DataFrame(rows)


def markdown_table(frame):
    lines = ['| ' + ' | '.join(frame.columns) + ' |', '| ' + ' | '.join(['---'] * len(frame.columns)) + ' |']
    for row in frame.itertuples(index=False, name=None):
        values = [f'{value:.5g}' if isinstance(value, float) else str(value) for value in row]
        lines.append('| ' + ' | '.join(value.replace('|', '\\|') for value in values) + ' |')
    return '\n'.join(lines)


def build(output):
    primary, manifest = load_primary()
    result = comparisons(primary)
    flips = result.loc[result.opposite_g_sign].sort_values(['significance_tier', 'disease', 'feature', 'race_1', 'race_2'])
    both = flips.loc[flips.significance_tier.eq('both_significant')]
    any_significant = flips.loc[flips.significance_tier.ne('neither_significant')]
    feature_summary = flips.groupby(['disease', 'feature'], as_index=False).agg(
        category_pair_reversals=('race_1', 'size'),
        at_least_one_significant=('significance_tier', lambda values: sum(values != 'neither_significant')),
        both_significant=('significance_tier', lambda values: sum(values == 'both_significant')))
    result.to_csv(output / 'all_cross_race_feature_comparisons.csv', index=False)
    flips.to_csv(output / 'all_direction_reversals.csv', index=False)
    any_significant.to_csv(output / 'direction_reversals_at_least_one_significant.csv', index=False)
    both.to_csv(output / 'direction_reversals_both_significant.csv', index=False)
    feature_summary.to_csv(output / 'unique_syndrome_features_changing_direction.csv', index=False)
    columns = ['disease', 'feature', 'race_1', 'g_1', 'q_1', 'case_n_1', 'race_2', 'g_2', 'q_2', 'case_n_2', 'significance_tier']
    summary = {'all_pairwise_comparisons': len(result), 'opposite_sign_comparisons': len(flips),
               'both_significant_reversals': len(both), 'at_least_one_significant_reversals': len(any_significant),
               'unique_syndrome_features_with_reversal': len(flips[['disease', 'feature']].drop_duplicates()),
               'syndromes_with_reversal': int(flips.disease.nunique()),
               'sign_tolerance': 1e-12, 'global_q_threshold': .05,
               'primary_family': len(primary), 'primary_sha256': common.sha(BASELINE / 'primary_all_tests.csv')}
    common.write_json(output / 'direction_reversal_summary.json', summary)
    lines = ['# Features changing direction across recorded categories', '',
             'This report exhaustively compares every available pair of eligible recorded-category cells for the same exact '
             'syndrome and feature in the new main analysis. It uses all planned main rows, including nonsignificant rows, '
             'and retains the original global primary q-values. No new inferential family or interaction test is introduced.', '',
             f'There are {len(result):,} available cross-category comparisons and {len(flips):,} opposite-sign comparisons '
             f'(absolute g must exceed 1e-12 on each side). Of these, {len(both)} have q<.05 on both sides and '
             f'{len(any_significant)} have q<.05 on at least one side. All cells contain at least five distinct recorded GMDB '
             'patients and sufficiently supported same-category references.', '',
             'A positive g means the syndrome group has a higher feature mean than its same-category FairFace reference; '
             'a negative g means a lower mean. Opposite signs can arise through different case means, different reference '
             'means, or both. The complete CSV includes these means and their differences, sample counts, formulas, '
             'units and both original q-values. Standardized g magnitudes also depend on group-specific SDs.', '',
             'Significance in one category and nonsignificance in another does not establish a category interaction. '
             'Even two significant opposite directions are descriptive evidence requiring a formal heterogeneity test '
             'before attributing the difference to category. The labels are recorded dataset annotations, not verified '
             'genetic ancestry. Sample size, source acquisition, broad Asian harmonization and reference composition remain limitations.', '',
             '## Opposite signs with both original q-values below .05', '']
    lines += [markdown_table(both[columns]) if len(both) else 'None. No supported two-significant-direction reversal is claimed.']
    lines += ['', '## Complete list with at least one significant result', '', markdown_table(any_significant[columns]), '',
              '## Opposite signs with neither result significant', '',
              'These are supplied in all_direction_reversals.csv and are not evidence of two established effects. '
              'The all_cross_race_feature_comparisons.csv file also includes every same-sign comparison, so the denominator is auditable.', '',
              '## Exhaustive appendix: every opposite-sign comparison', '',
              'This includes nonsignificant exploratory sign changes. They must not be described as established category-specific effects. '
              'The unique_syndrome_features_changing_direction.csv file collapses pairwise comparisons to one row per syndrome/feature.', '',
              markdown_table(flips[columns]), '',
              '## Provenance', '', f'Primary baseline: `{manifest["study_id"]}`. Family size: {len(primary):,}. '
              f'Primary file SHA-256: `{summary["primary_sha256"]}`. '
              'No patient rows, photographs or individual measurements are included in these aggregate reports. '
              'All outputs remain in ignored private storage.']
    (output / 'DIRECTION_REVERSALS_REPORT.md').write_text('\n'.join(lines) + '\n')
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, default=OUTPUT)
    args = parser.parse_args()
    print(json.dumps(build(common.private_output(args.output_dir)), indent=2))

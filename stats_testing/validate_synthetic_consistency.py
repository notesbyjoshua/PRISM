import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from synthetic_consistency import ROOT, load_inputs, sha


def wasserstein_cdf(x, y):
    x, y = np.sort(x), np.sort(y)
    support = np.sort(np.concatenate([x, y]))
    cdf_difference = np.searchsorted(x, support[:-1], side='right') / len(x)
    cdf_difference -= np.searchsorted(y, support[:-1], side='right') / len(y)
    return np.sum(np.diff(support) * np.abs(cdf_difference))


def validate():
    config, paths, definitions, provenance, real, reference, all_synthetic, synthetic, measured = load_inputs(ROOT / 'ACTIVE_ANALYSIS.json')
    features = definitions.feature.tolist()
    exported = pd.read_csv(ROOT / 'publication/synthetic_consistency/feature_fidelity_summary.csv')
    checked, maximum = 0, 0.
    for (syndrome, synthetic_label, tag, sensitivity), output in exported.groupby(['syndrome', 'pdidb_disease', 'model_tag', 'sensitivity']):
        x = real.loc[real.disease.eq(syndrome), features].to_numpy(float)
        subset = synthetic.loc[synthetic.disease.eq(synthetic_label)]
        if tag != 'pooled':
            subset = subset.loc[subset.model.eq(tag)]
        y = subset[features].to_numpy(float)
        ref = reference[features].to_numpy(float)
        if sensitivity == 'winsorized_01_99':
            x, y, ref = [np.clip(v, *np.quantile(v, [.01, .99], axis=0)) for v in [x, y, ref]]
        expected = {}
        for name, values in [('real', x), ('synthetic', y), ('reference', ref)]:
            expected[name + '_mean'] = values.mean(axis=0)
            expected[name + '_sd'] = values.std(axis=0, ddof=1)
            expected[name + '_n'] = np.repeat(len(values), 120)
            for suffix, q in [('q1', .25), ('median', .5), ('q3', .75)]:
                expected[name + '_' + suffix] = np.quantile(values, q, axis=0)
        expected['mean_shift'] = y.mean(axis=0) - x.mean(axis=0)
        expected['median_shift'] = np.median(y, axis=0) - np.median(x, axis=0)
        expected['standardized_mean_shift'] = expected['mean_shift'] / x.std(axis=0, ddof=1)
        expected['variance_ratio'] = y.var(axis=0, ddof=1) / x.var(axis=0, ddof=1)
        expected['wasserstein'] = [wasserstein_cdf(x[:, j], y[:, j]) for j in range(120)]
        for name, values in [('real', x), ('synthetic', y)]:
            n, m = len(values), len(ref)
            pooled = np.sqrt(((n - 1) * values.var(axis=0, ddof=1) + (m - 1) * ref.var(axis=0, ddof=1)) / (n + m - 2))
            expected[name + '_g'] = (values.mean(axis=0) - ref.mean(axis=0)) / pooled * (1 - 3 / (4 * (n + m) - 9))
        output = output.set_index('feature').loc[features]
        for column, values in expected.items():
            np.testing.assert_allclose(output[column], values, rtol=1e-10, atol=1e-12)
            maximum = max(maximum, np.max(np.abs(output[column].to_numpy() - values)))
            checked += len(features)
    return {'status': 'PASS', 'syndrome_feature_rows': len(exported), 'numeric_cells_independently_checked': checked,
            'feature_summary_sha256': sha(ROOT / 'publication/synthetic_consistency/feature_fidelity_summary.csv'),
            'validation_code_sha256': sha(Path(__file__)),
            'maximum_absolute_difference': float(maximum), 'source_rows': {'real': len(real), 'synthetic_qc': len(synthetic), 'reference': len(reference)},
            'method': 'Independent vector means/SD/quantiles/effects and empirical-CDF integral Wasserstein, all model and winsorization strata',
            'primary_results_changed': False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    if not args.report.resolve().is_relative_to(ROOT / 'local_review'):
        raise ValueError('Private validation report must stay in local_review/')
    result = validate()
    args.report.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))

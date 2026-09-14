import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def run(disease_path, reference_path, output, repeats=100, seed=2026):
    disease, reference = pd.read_csv(disease_path), pd.read_csv(reference_path)
    features = list(reference.loc[:, 'eb_thickness_r':'cheek_area_asym'])
    y = reference[features].to_numpy(float)
    rng = np.random.default_rng(seed)
    sampled_n = len(y) // 2
    sampled = np.stack([y[rng.choice(len(y), sampled_n, replace=False)] for _ in range(repeats)])
    means, variances = sampled.mean(axis=1), sampled.var(axis=1, ddof=1)
    rows = []
    for syndrome, group in disease.groupby('disease'):
        if len(group) < 10:
            continue
        x = group[features].to_numpy(float)
        mean, variance = x.mean(axis=0), x.var(axis=0, ddof=1)
        full_df, sample_df = len(x) + len(y) - 2, len(x) + sampled_n - 2
        full_g = (1 - 3 / (4 * full_df - 1)) * (mean - y.mean(axis=0)) / np.sqrt(
            ((len(x) - 1) * variance + (len(y) - 1) * y.var(axis=0, ddof=1)) / full_df)
        sample_g = (1 - 3 / (4 * sample_df - 1)) * (mean - means) / np.sqrt(
            ((len(x) - 1) * variance + (sampled_n - 1) * variances) / sample_df)
        lower, upper = np.quantile(sample_g, [.025, .975], axis=0)
        for i, feature in enumerate(features):
            rows.append({'disease': syndrome, 'feature': feature, 'n_disease': len(x), 'full_reference_n': len(y),
                         'subsample_reference_n': sampled_n, 'full_g': full_g[i],
                         'subsample_g_p025': lower[i], 'subsample_g_p975': upper[i],
                         'same_direction_fraction': np.mean(np.sign(sample_g[:, i]) == np.sign(full_g[i])),
                         'repetitions': repeats, 'seed': seed})
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError(output)
    pd.DataFrame(rows).to_csv(output, index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--disease-csv', type=Path, required=True)
    parser.add_argument('--reference-csv', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(args.disease_csv, args.reference_csv, args.output)

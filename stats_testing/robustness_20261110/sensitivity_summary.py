import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from .common import DEFAULT_OUT, KEY, SEED, case_rows, features, load_private, offline, private_output, protect_baseline, write_aggregate
from .stability import patient_samples, fixed_reference_effect


def seed_rankings(real, reference, matrix, names, out):
    original = pd.read_csv(out / 'P1_ranked_candidate_evidence.csv', dtype={'canonical_syndrome_id': 'string'}, low_memory=False)
    original = original.loc[original.ranking_rule.eq('high_stability')]
    original_keys = set(map(tuple, original[KEY].to_numpy()))
    summaries, chosen_rows = [], []
    for seed_variant in [0, 1, 2]:
        parts = []
        for identifier, group in matrix.groupby('canonical_syndrome_id', sort=True):
            point = group.set_index('feature').loc[names]
            if point.ff_n_disease.iloc[0] < 20:
                continue
            syndrome = point.syndrome.iloc[0]
            target = real.loc[real.disease.eq(syndrome)]
            samples = patient_samples(target[names].to_numpy(float), target.patient_id, 200,
                                      np.random.default_rng([SEED, 50, seed_variant, int(identifier)]))
            frequencies = []
            for source, sign_column in [(reference, 'ff_hedges_g'), (real.loc[real.disease.ne(syndrome)], 'gmdb_hedges_g')]:
                draws = fixed_reference_effect(samples, source[names].to_numpy(float))
                frequencies.append(np.mean(np.sign(draws) == np.sign(point[sign_column].to_numpy()), axis=0))
            frame = point.reset_index().copy()
            frame['seed_minimum_direction_frequency'] = np.minimum.reduce([*frequencies, frame.same_direction_photo_frequency.to_numpy()])
            mask = frame.ff_mannwhitney_q_fdr_global.lt(.05) & frame.reference_sign_consistent & frame.seed_minimum_direction_frequency.ge(.95)
            parts.append(frame.loc[mask])
        selected = pd.concat(parts, ignore_index=True).sort_values(['seed_minimum_direction_frequency', 'minimum_absolute_effect'] + KEY,
                                                                 ascending=[False, False, True, True]).head(20)
        keys = set(map(tuple, selected[KEY].to_numpy()))
        summaries.append(dict(seed_variant=seed_variant, base_seed=SEED, bootstrap_draws=200, selected=len(selected),
                               shared_with_original_high_stability=len(keys & original_keys),
                               status='post hoc Monte Carlo sensitivity; same patients; no held-out performance claim'))
        selected['seed_variant'] = seed_variant
        chosen_rows.append(selected[KEY + ['syndrome', 'seed_variant', 'seed_minimum_direction_frequency']])
    write_aggregate(pd.DataFrame(summaries), out / 'aggregate_candidates/P1_ranking_seed_sensitivity.csv')
    write_aggregate(pd.concat(chosen_rows, ignore_index=True), out / 'aggregate_candidates/P1_seed_candidate_lists.csv')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    offline()
    out = private_output(args.output_dir)
    protect_baseline(out)
    real, reference, raw, alternate, paths = load_private()
    matrix = pd.read_csv(out / 'aggregate_candidates/P1_all_candidate_evidence_private.csv',
                         dtype={'canonical_syndrome_id': 'string'}, low_memory=False)
    names = features().feature.tolist()
    seed_rankings(real, reference, matrix, names, out)
    age = pd.read_csv(out / 'evidence/P0_age_stratum_long.csv', dtype={'canonical_syndrome_id': 'string'})
    focal = case_rows(matrix)[KEY + ['syndrome', 'ff_hedges_g', 'gmdb_hedges_g', 'screen_target_beta', 'screen_q_global_bh']]
    age = age.merge(focal.drop(columns='syndrome'), on=KEY, how='inner', validate='many_to_one')
    write_aggregate(age, out / 'aggregate_candidates/P1_age_case_summary.csv')
    family = matrix.groupby(['anatomical_region', 'feature_family'], sort=True).agg(features=('feature', 'nunique'),
        tested=('feature', 'size'), significant=('ff_mannwhitney_q_fdr_global', lambda p: p.lt(.05).sum())).reset_index()
    family['significant_cell_weight_per_family'] = family.significant / family.features
    family['maximum_cell_weight_per_family'] = family.tested / family.features
    family['interpretation'] = 'equal descriptive family weights; not independent tests; original BH unchanged'
    write_aggregate(family, out / 'aggregate_candidates/P1_family_weighted_discovery_summary.csv')
    region = matrix.groupby('anatomical_region', sort=True).agg(raw_features=('feature', 'nunique'),
        raw_significant=('ff_mannwhitney_q_fdr_global', lambda p: p.lt(.05).sum()), tested=('feature', 'size')).reset_index()
    representatives = matrix.loc[matrix.representative].groupby('anatomical_region', sort=True).agg(
        representative_features=('feature', 'nunique'), representative_significant=('ff_mannwhitney_q_fdr_global', lambda p: p.lt(.05).sum()))
    region = region.merge(representatives, on='anatomical_region', how='left').fillna(0)
    region['interpretation'] = 'Representative assigned to its own anatomical region; q from original complete family'
    write_aggregate(region, out / 'aggregate_candidates/P1_regional_representative_summary.csv')
    protect_baseline(out)
    print('Three-seed candidate sensitivity, existing age-family cases and family weighting summaries complete')


if __name__ == '__main__':
    main()

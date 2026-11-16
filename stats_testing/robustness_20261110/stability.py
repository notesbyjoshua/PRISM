import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, trim_mean

from stats_testing.cohort_audit import age_diagnostics, clean_id
from .common import (ROOT, DEFAULT_OUT, KEY, SEED, canonicalize, effect, features, load_private, offline,
                     private_output, protect_baseline, require, sha, u_test, bh, write_aggregate, write_json)


def patient_samples(values, identifiers, repetitions, rng):
    identifiers = np.asarray(identifiers)
    require(len(identifiers) == len(values) and len(np.unique(identifiers)) == len(values), 'Repeated real patients in bootstrap input')
    return values[rng.integers(len(values), size=(repetitions, len(values)))]


def fixed_reference_effect(samples, reference):
    n, m = samples.shape[1], len(reference)
    pooled = ((n - 1) * samples.var(axis=1, ddof=1) + (m - 1) * reference.var(axis=0, ddof=1)) / (n + m - 2)
    with np.errstate(divide='ignore', invalid='ignore'):
        return (1 - 3 / (4 * (n + m) - 9)) * (samples.mean(axis=1) - reference.mean(axis=0)) / np.sqrt(pooled)


def photo_recovery(real, reference, alternate, names, out):
    config = json.loads((ROOT / 'ACTIVE_ANALYSIS.json').read_text())
    run = ROOT / config['private_run']
    raw_path = run / 'COHORTS/all_gmdb_candidates_master.csv'
    raw = pd.read_csv(raw_path, dtype={'patient_id': 'string', 'image_id': 'string'})
    age = age_diagnostics(raw).conditional_age_years
    eligible = raw.loc[raw.frontal_ok.eq(True) & age.ge(3) & age.lt(10)].copy()
    eligible['patient_id'] = clean_id(eligible.patient_id)
    require(eligible.patient_id.notna().all() and np.isfinite(eligible[names]).all().all(), 'Photo pool has invalid recorded IDs or features')
    groups = list(eligible.groupby('patient_id', sort=True).indices.values())
    require(len(groups) == 1845, 'Original photo-pool patient count differs')
    targets = json.loads((run / 'ANALYSIS_PLAN.json').read_text())['targets']
    counts = eligible.groupby('patient_id').size()
    primary_ids = real.patient_id.astype('string')
    alternate = alternate.set_index('patient_id')
    require(set(primary_ids) == set(alternate.index), 'Alternate photograph patient set differs')
    changes = real.image_id.astype('string').to_numpy() != alternate.loc[primary_ids, 'image_id'].astype('string').to_numpy()
    photo_counts = []
    for syndrome, part in real.assign(changed=changes).groupby('disease', sort=True):
        if syndrome in targets:
            photo_counts.append(dict(syndrome=syndrome, primary_patients=len(part), changed_alternate=int(part.changed.sum()),
                                     patients_with_multiple_eligible_images=int(counts.reindex(part.patient_id).gt(1).sum()),
                                     eligible_images=int(counts.reindex(part.patient_id).sum())))
    write_aggregate(canonicalize(pd.DataFrame(photo_counts)), out / 'aggregate_candidates/P1_primary_photo_opportunities.csv')
    checkpoint = out / 'checkpoints/exact_photo_recovery.npz'
    signature = {'raw_sha256': sha(raw_path), 'reference_sha256': sha(ROOT / config['reference_master_csv']),
                 'code_sha256': sha(Path(__file__)), 'original_seed': 20261025, 'repetitions': 100, 'hypotheses': 5760}
    metadata = checkpoint.with_suffix('.json')
    if checkpoint.exists():
        require(json.loads(metadata.read_text()) == signature, 'Changed exact-photo replay signature; choose a new run')
        data = np.load(checkpoint)
        gs, ps, qs, ranks = [data[key] for key in ['g', 'p', 'q', 'rank']]
    else:
        values = eligible[names].to_numpy(float)
        labels = eligible.disease.to_numpy()
        ref = reference[names].to_numpy(float)
        rng = np.random.default_rng(20261025)
        gs, ps, qs, ranks = [], [], [], []
        for rep in range(100):
            chosen = np.array([indices[rng.integers(len(indices))] for indices in groups])
            sample, selected_labels = values[chosen], labels[chosen]
            effects, pvalues = [], []
            for syndrome in targets:
                x = sample[selected_labels == syndrome]
                require(len(x) >= 10, 'Photo replay changed primary support')
                effects.extend(effect(x, ref))
                pvalues.extend(u_test(x, ref)[1])
            effects, pvalues = np.array(effects), np.array(pvalues)
            gs.append(effects)
            ps.append(pvalues)
            qs.append(bh(pvalues))
            ranks.append(rankdata(-np.abs(effects).reshape(48, 120), axis=1, method='average').reshape(-1))
            if (rep + 1) % 25 == 0:
                print(f'Exact original photo recovery: {rep + 1}/100', flush=True)
        gs, ps, qs, ranks = [np.array(x) for x in [gs, ps, qs, ranks]]
        checkpoint.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
        np.savez_compressed(checkpoint, g=gs, p=ps, q=qs, rank=ranks)
        write_json(signature, metadata)
    recovered = pd.DataFrame([(syndrome, feature) for syndrome in targets for feature in names], columns=['syndrome', 'feature'])
    for name, values in [('median_g', np.median(gs, axis=0)), ('g_q025', np.quantile(gs, .025, axis=0)),
                         ('g_q975', np.quantile(gs, .975, axis=0)), ('min_g', gs.min(axis=0)), ('max_g', gs.max(axis=0)),
                         ('positive_g_frequency', (gs > 0).mean(axis=0)), ('significant_frequency', (qs < .05).mean(axis=0))]:
        recovered[name] = values
    frozen = pd.read_csv(ROOT / 'publication/supplementary/random_photo_sensitivity.csv')
    compared = recovered.merge(frozen, on=['syndrome', 'feature'], validate='one_to_one', suffixes=('_replay', '_frozen'))
    maximum = 0.
    for column in ['median_g', 'g_q025', 'g_q975', 'min_g', 'max_g', 'positive_g_frequency', 'significant_frequency']:
        a, b = compared[column + '_replay'].to_numpy(), compared[column + '_frozen'].to_numpy()
        require(np.allclose(a, b, rtol=1e-10, atol=1e-12), 'Exact recovery disagrees with frozen random-photo summaries')
        maximum = max(maximum, float(np.max(np.abs(a - b))))
    for name, values in [('raw_p_q025', np.quantile(ps, .025, axis=0)), ('raw_p_median', np.median(ps, axis=0)),
                         ('raw_p_q975', np.quantile(ps, .975, axis=0)), ('global_q_q025', np.quantile(qs, .025, axis=0)),
                         ('global_q_median', np.median(qs, axis=0)), ('global_q_q975', np.quantile(qs, .975, axis=0)),
                         ('rank_q025', np.quantile(ranks, .025, axis=0)), ('rank_median', np.median(ranks, axis=0)),
                         ('rank_q975', np.quantile(ranks, .975, axis=0)), ('top10_frequency', (ranks <= 10).mean(axis=0))]:
        recovered[name] = values
    write_aggregate(canonicalize(recovered), out / 'aggregate_candidates/P1_exact_photo_stability.csv')
    write_json({'status': 'PASS', 'repetitions': 100, 'summary_rows': len(recovered), 'max_frozen_summary_difference': maximum,
                'all_selected_changed_photos': int(changes.sum()), 'primary_cohort_changed_photos': int(sum(r['changed_alternate'] for r in photo_counts)),
                'primary_patients_with_multiple_images': int(sum(r['patients_with_multiple_eligible_images'] for r in photo_counts)),
                'role': 'Exact recovery of existing 100 draws; no new photo experiment'}, out / 'P1_PHOTO_REPLAY_VALIDATION.json')


def patient_stability(real, reference, names, out):
    counts = real.disease.value_counts()
    targets = sorted(counts.loc[counts.ge(10)].index)
    rows, influence, estimators = [], [], []
    all_draws = []
    for index, syndrome in enumerate(targets):
        target = real.loc[real.disease.eq(syndrome)]
        x = target[names].to_numpy(float)
        samples = patient_samples(x, target.patient_id, 200, np.random.default_rng([SEED, 30, index]))
        for source, comparator in [('FairFace', reference), ('other_GMDB', real.loc[real.disease.ne(syndrome)])]:
            y = comparator[names].to_numpy(float)
            point = effect(x, y)
            draws = fixed_reference_effect(samples, y)
            all_draws.append(draws)
            subsample_n = int(np.floor(.7 * len(x)))
            subsamples = []
            if subsample_n >= 10:
                rng = np.random.default_rng([SEED, 31, index])
                subsamples = np.array([effect(x[rng.choice(len(x), subsample_n, replace=False)], y) for _ in range(100)])
            u, _ = u_test(x, y)
            delta = 2 * u / (len(x) * len(y)) - 1
            for j, feature in enumerate(names):
                rows.append(dict(syndrome=syndrome, feature=feature, comparator=source, n_target=len(x), n_comparator=len(y),
                                 point_g=point[j], bootstrap_q025=np.quantile(draws[:, j], .025), bootstrap_q975=np.quantile(draws[:, j], .975),
                                 same_direction_frequency=np.mean(np.sign(draws[:, j]) == np.sign(point[j])), bootstrap_draws=200,
                                 subsample_n=subsample_n, subsample_supported=subsample_n >= 10,
                                 subsample_q025=np.quantile(subsamples[:, j], .025) if len(subsamples) else np.nan,
                                 subsample_q975=np.quantile(subsamples[:, j], .975) if len(subsamples) else np.nan,
                                 subsample_direction_frequency=np.mean(np.sign(subsamples[:, j]) == np.sign(point[j])) if len(subsamples) else np.nan,
                                 uncertainty='case-patient sampling conditional on fixed comparator; descriptive quantile range; not concordance CI'))
                estimators.append(dict(syndrome=syndrome, feature=feature, comparator=source, hedges_g=point[j], cliffs_delta=delta[j],
                                       median_difference=np.median(x[:, j]) - np.median(y[:, j]),
                                       trimmed_mean_difference=trim_mean(x[:, j], .1) - trim_mean(y[:, j], .1),
                                       feature_units='see exact S6; medians/trimmed means in original normalized units',
                                       direction_agreement_g_delta=bool(np.sign(point[j]) == np.sign(delta[j]))))
            if len(x) >= 20:
                leave = np.array([effect(np.delete(x, k, axis=0), y) for k in range(len(x))])
                for j, feature in enumerate(names):
                    influence.append(dict(syndrome=syndrome, feature=feature, comparator=source, n=len(x),
                                          maximum_absolute_g_change=np.max(np.abs(leave[:, j] - point[j])),
                                          direction_change_fraction=np.mean(np.sign(leave[:, j]) != np.sign(point[j])),
                                          privacy='aggregate only; no patient identity/index or individual influence exported'))
    write_aggregate(canonicalize(pd.DataFrame(rows)), out / 'aggregate_candidates/P1_patient_resampling.csv')
    write_aggregate(canonicalize(pd.DataFrame(estimators)), out / 'aggregate_candidates/P1_robust_estimators.csv')
    write_aggregate(canonicalize(pd.DataFrame(influence)), out / 'aggregate_candidates/P1_aggregate_influence.csv')
    np.savez_compressed(out / 'checkpoints/patient_effect_draws.npz', effects=np.array(all_draws))


def coverage_experiment(out, repetitions=300, bootstrap=200):
    rows = []
    for distribution_index, distribution in enumerate(['normal', 'scaled_t5']):
        for n in [10, 35]:
            rng = np.random.default_rng([SEED, 32, distribution_index, n])
            ref = rng.normal(size=(500, 1))
            correction = 1 - 3 / (4 * (n + len(ref)) - 9)
            target = correction * (.5 - ref.mean()) / np.sqrt(((n - 1) + (len(ref) - 1) * ref.var(ddof=1)) / (n + len(ref) - 2))
            covered = 0
            for rep in range(repetitions):
                values = rng.normal(size=(n, 1)) if distribution == 'normal' else rng.standard_t(5, size=(n, 1)) / np.sqrt(5 / 3)
                values += .5
                samples = patient_samples(values, np.arange(n), bootstrap, rng)
                low, high = np.quantile(fixed_reference_effect(samples, ref)[:, 0], [.025, .975])
                covered += low <= target <= high
            p = covered / repetitions
            z, denominator = 1.959963984540054, 1 + 1.959963984540054**2 / repetitions
            center = (p + z**2 / (2 * repetitions)) / denominator
            half = z * np.sqrt(p * (1 - p) / repetitions + z**2 / (4 * repetitions**2)) / denominator
            rows.append(dict(distribution=distribution, case_n=n, fixed_reference_n=500, independent_simulations=repetitions,
                             bootstrap_draws=bootstrap, target_nominal_coverage=.95, covered=covered, empirical_coverage=p,
                             coverage_mc_wilson_low=center - half, coverage_mc_wilson_high=center + half,
                             role='simulation-only conditional functional; not clinical coverage validation'))
    frame = pd.DataFrame(rows)
    write_aggregate(frame, out / 'aggregate_candidates/P1_bootstrap_coverage_simulation.csv')
    return frame


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    offline()
    out = private_output(args.output_dir)
    protect_baseline(out)
    real, reference, raw, alternate, paths = load_private()
    names = features().feature.tolist()
    photo_recovery(real, reference, alternate, names, out)
    patient_stability(real, reference, names, out)
    coverage = coverage_experiment(out)
    opportunities = pd.read_csv(out / 'aggregate_candidates/P1_primary_photo_opportunities.csv')
    replay = json.loads((out / 'P1_PHOTO_REPLAY_VALIDATION.json').read_text())
    text = ('# P1 photograph, patient and estimator stability\n\n'
            f"The deterministic alternate choice changes {replay['all_selected_changed_photos']}/1,845 selected photographs, "
            f"including {replay['primary_cohort_changed_photos']}/872 in the primary 48 cohorts. "
            f"{replay['primary_patients_with_multiple_images']}/872 primary patients have more than one eligible image. "
            'Per-syndrome opportunity counts are supplied rather than treating every patient as perturbed.\n\n'
            'Exactly the existing 100 random-photo draws were replayed with seed 20261025 and original patient-group/order rules, solely to recover per-repetition raw p, global q and feature ranks. '
            f"All seven frozen per-hypothesis summaries match (maximum absolute difference {replay['max_frozen_summary_difference']:.6g}). "
            'No new photo experiment replaces the original. Rank ranges and significance frequency are distinguished from effect-direction frequency.\n\n'
            'Case-patient bootstraps (200) preserve correlated feature rows and hold the selected empirical comparator fixed. '
            'The same case draw feeds both contrasts, but no joint concordance-bootstrap interval is computed. '
            'Seventy-percent case subsampling (100 draws) is unsupported when retained n<10; its missing range is not a tested null. '
            'Leave-one-patient influence is limited to n>=20 and exports only aggregate maxima/fractions, never individual identities or influence rows. '
            'Median and 10%-trimmed-mean contrasts use original feature units; Cliff’s delta uses rank dominance. They are not interchangeable with g.\n\n'
            'Percentile ranges are conditional sampling descriptions, not calibrated clinical confidence intervals. In normal/scaled-t5 numerical simulations with case n=10/35, '
            f"observed nominal-95% coverage ranges from {coverage.empirical_coverage.min():.1%} to {coverage.empirical_coverage.max():.1%}. "
            'Monte Carlo Wilson intervals use independent simulated datasets. These checks cannot validate clinical coverage and may reveal undercoverage, especially at small n. '
            'Existing age-stratum and adjusted estimates retain their original declared families; FairFace exact-age matching remains impossible.\n')
    (out / 'P1_stability_report.md').write_text(text)
    protect_baseline(out)
    print('Photo recovery, patient-level resampling, influence and coverage checks complete')


if __name__ == '__main__':
    main()

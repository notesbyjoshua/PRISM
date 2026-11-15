import numpy as np
import pandas as pd
from scipy.stats import beta, rankdata

from core import (KEYS, RACES, Checkpoints, bh, effect_moments, empirical_p, p_from_untied_u,
                  placebo_allocation, permuted_order, rng_for, signature, tests_direct)

REPETITION_COLUMNS = ['repetition', 'discoveries', 'African_discoveries', 'Asian_discoveries',
                      'European_discoveries', 'minimum_valid_q', 'maximum_absolute_g',
                      'overlap_with_main', 'effect_direction_changes']


def repetition_summary(study, rep, q, g, supported):
    discoveries = q < .05
    original = study.primary.q.to_numpy() < .05
    by_race = [int(discoveries[study.races[race]['features_index']].sum()) for race in RACES]
    valid = np.isfinite(g)
    return [rep, int(discoveries.sum()), *by_race, float(q[supported].min()) if supported.any() else np.nan,
            float(np.abs(g[valid]).max()) if valid.any() else np.nan, int((discoveries & original).sum()),
            int((valid & (np.sign(g) != np.sign(study.primary.hedges_g.to_numpy()))).sum())]


def support_table(study):
    rows = []
    for race, data in study.races.items():
        nref, pool = study.sampling['quotas'][race], len(data['y'])
        shared = pool >= int(data['ns'].sum()) + nref
        for cell_index, n in zip(data['rows'], data['ns']):
            rows.append({'disease': study.cells.disease.iloc[cell_index], 'canonical_race': race,
                         'n_case': int(n), 'n_reference': nref, 'pool': pool,
                         'supported': bool(pool >= n + nref and nref >= 20),
                         'design': 'shared_disjoint' if shared else 'per_pair_disjoint_fallback',
                         'shortage': max(0, int(n) + nref - pool)})
    return pd.DataFrame(rows)


def run_placebo(study, plan, output, repetitions, resume):
    checkpoint = Checkpoints(output, 'placebo', signature(plan, 'placebo'), resume)
    support = support_table(study)
    support.to_csv(output / 'fairface_placebo_support.csv', index=False)
    summaries, hits = [], np.zeros(study.total, dtype=np.int64)
    supported = np.zeros(study.total, dtype=bool)
    for rep in range(checkpoint.state['completed'], repetitions):
        p, g = np.ones(study.total), np.full(study.total, np.nan)
        for race_index, (race, data) in enumerate(study.races.items()):
            reference_n = study.sampling['quotas'][race]
            design, allocations = placebo_allocation(len(data['y']), data['ns'], reference_n,
                                                       rng_for(plan['seed'], 1, rep, race_index))
            if design == 'shared_disjoint':
                case_indices = np.concatenate([pair[0] for pair in allocations])
                reference_indices = allocations[0][1]
                case_values, reference_values = data['y'][case_indices], data['y'][reference_indices]
                scores = (rankdata(np.r_[case_values, reference_values], axis=0)[:len(case_values)]
                          - rankdata(case_values, axis=0))
                us = np.add.reduceat(scores, np.r_[0, np.cumsum(data['ns'])[:-1]], axis=0)
                ymean, yvar = reference_values.mean(0), reference_values.var(0, ddof=1)
            else:
                us = None
            for j, (cell_index, pair) in enumerate(zip(data['rows'], allocations)):
                if pair is None or reference_n < 20:
                    continue
                case, reference = pair
                assert len(case) == data['ns'][j] and len(reference) == reference_n
                assert not set(case) & set(reference) and len(set(reference)) == reference_n
                indices = np.arange(cell_index * len(study.features), (cell_index + 1) * len(study.features))
                supported[indices] = True
                x, y = data['y'][case], data['y'][reference]
                if us is not None and data['ff_untied']:
                    p[indices] = p_from_untied_u(us[j], len(x), len(y))
                    g[indices] = effect_moments(len(x), len(y), x.mean(0), x.var(0, ddof=1), ymean, yvar)
                else:
                    _, p[indices], g[indices] = tests_direct(x, y)
        q = bh(p)
        hits += q < .05
        summaries.append(repetition_summary(study, rep, q, g, supported))
        if (rep + 1) % 25 == 0 or rep + 1 == repetitions:
            checkpoint.append(rep + 1, summary=np.asarray(summaries), hits=hits)
            summaries, hits = [], np.zeros(study.total, dtype=np.int64)
            print(f'FairFace placebo: {rep + 1}/{repetitions}', flush=True)
    frame = pd.DataFrame(np.vstack(checkpoint.arrays('summary')), columns=REPETITION_COLUMNS)
    frame['planned_family'] = study.total
    frame['supported_hypotheses'] = int(support.supported.sum()) * len(study.features)
    frame['interpretation'] = 'same-source placebo; observed-primary exceedance is descriptive, not an inferential p'
    frame.to_csv(output / 'fairface_placebo_summary.csv', index=False)
    return frame


class ReferenceBank:
    def __init__(self, study):
        self.study, self.races = study, {}
        for race, data in study.races.items():
            y = data['y']
            contributions = np.empty((len(y), len(data['rows']), len(study.features)), dtype=np.float32)
            cases, tied = [], []
            for j, cell_index in enumerate(data['rows']):
                target = study.cells.disease.iloc[cell_index]
                x = data['x'][data['patients'].disease.eq(target).to_numpy()]
                cases.append(x)
                combined = np.r_[x, y]
                tied.append((np.diff(np.sort(combined, axis=0), axis=0) == 0).any(0))
                for feature in range(len(study.features)):
                    ordered = np.sort(x[:, feature])
                    left = np.searchsorted(ordered, y[:, feature], side='left')
                    right = np.searchsorted(ordered, y[:, feature], side='right')
                    contributions[:, j, feature] = len(x) - (left + right) / 2
            self.races[race] = {'contributions': contributions, 'cases': cases, 'tied': np.asarray(tied)}
        self.full_cache = {}

    def evaluate(self, selections):
        p, g = np.ones(self.study.total), np.full(self.study.total, np.nan)
        for race, data in self.study.races.items():
            indices = np.asarray(selections[race], dtype=int)
            if len(indices) != self.study.sampling['quotas'][race] or len(np.unique(indices)) != len(indices):
                raise ValueError('Reference sample count/replacement constraint failed')
            if len(indices) < 20:
                continue
            bank = self.races[race]
            if len(indices) == len(data['y']) and race in self.full_cache:
                p[data['features_index']], g[data['features_index']] = self.full_cache[race]
                continue
            y = data['y'][indices]
            u = bank['contributions'][indices].sum(0, dtype=np.float64)
            ymean, yvar = y.mean(0), y.var(0, ddof=1)
            for j, cell_index in enumerate(data['rows']):
                where = np.arange(cell_index * len(self.study.features), (cell_index + 1) * len(self.study.features))
                x = bank['cases'][j]
                if bank['tied'][j].any():
                    _, p[where], g[where] = tests_direct(x, y)
                else:
                    p[where] = p_from_untied_u(u[j], len(x), len(y))
                    g[where] = effect_moments(len(x), len(y), x.mean(0), x.var(0, ddof=1), ymean, yvar)
            if len(indices) == len(data['y']):
                self.full_cache[race] = (p[data['features_index']].copy(), g[data['features_index']].copy())
        return p, g


def prepare_fast_paths(study):
    for data in study.races.values():
        data['ff_untied'] = not (np.diff(np.sort(data['y'], axis=0), axis=0) == 0).any()


def run_reference(study, plan, output, repetitions, resume):
    checkpoint = Checkpoints(output, 'reference-resampling', signature(plan, 'reference-resampling'), resume)
    bank = ReferenceBank(study)
    original = {}
    for race, data in study.races.items():
        selected = study.reference.loc[study.reference.canonical_race.eq(race)].sort_values('file')
        primary_sample = study.baseline_sample
        original[race] = np.flatnonzero(selected.file.isin(primary_sample.file).to_numpy())
    p, g = bank.evaluate(original)
    np.testing.assert_allclose(p, study.primary.p, rtol=1e-10, atol=1e-15)
    np.testing.assert_allclose(g, study.primary.hedges_g, rtol=1e-10, atol=1e-12)
    summaries, gs, qs = [], [], []
    for rep in range(checkpoint.state['completed'], repetitions):
        selections = {race: rng_for(plan['seed'], 3, rep, i).choice(len(data['y']), study.sampling['quotas'][race], replace=False)
                      for i, (race, data) in enumerate(study.races.items())}
        p, g = bank.evaluate(selections)
        q = bh(p)
        summaries.append(repetition_summary(study, rep, q, g, np.isfinite(g)))
        gs.append(g)
        qs.append(q)
        if (rep + 1) % 25 == 0 or rep + 1 == repetitions:
            checkpoint.append(rep + 1, summary=np.asarray(summaries), g=np.asarray(gs), q=np.asarray(qs))
            summaries, gs, qs = [], [], []
            print(f'FairFace reference resampling: {rep + 1}/{repetitions}', flush=True)
    g = np.vstack(checkpoint.arrays('g'))
    q = np.vstack(checkpoint.arrays('q'))
    result = study.primary[KEYS + ['case_n', 'reference_n', 'hedges_g', 'q', 'supported']].rename(
        columns={'hedges_g': 'observed_main_g', 'q': 'observed_main_q', 'case_n': 'n_case', 'reference_n': 'n_reference'})
    result['median_g'], result['g_resampling_q025'], result['g_resampling_q975'] = np.quantile(g, [.5, .025, .975], axis=0)
    result['sign_agreement_frequency'] = (np.sign(g) == np.sign(result.observed_main_g.to_numpy())).mean(0)
    result['significant_frequency'] = (q < .05).mean(0)
    result['median_q'], result['minimum_q'], result['maximum_q'] = np.median(q, axis=0), q.min(0), q.max(0)
    result['repetitions'] = len(g)
    result['descriptive_stable'] = result.significant_frequency.ge(.8) & result.sign_agreement_frequency.ge(.95)
    result['reference_varies'] = result.canonical_race.map(
        {race: study.sampling['quotas'][race] < len(data['y']) for race, data in study.races.items()})
    result.to_csv(output / 'reference_resampling_aggregate.csv', index=False)
    frame = pd.DataFrame(np.vstack(checkpoint.arrays('summary')), columns=REPETITION_COLUMNS)
    frame['planned_family'] = study.total
    frame.to_csv(output / 'reference_resampling_repetitions.csv', index=False)
    return result, frame


def permutation_counts(study, seed, rep, blocked, counts):
    for race_index, (race, data) in enumerate(study.races.items()):
        blocks = data['blocks'] if blocked else None
        order = permuted_order(len(data['x']), rng_for(seed, 2, rep, race_index, int(blocked)), blocks)
        n = data['ns'][:, None]
        u = np.add.reduceat(data['ranks'][order], data['starts'], axis=0)[data['positions']] - n * (n + 1) / 2
        center = data['block_center'] if blocked else n * (len(data['x']) - n) / 2
        observed = np.abs(data['observed_u'] - center)
        counts[data['features_index']] += (np.abs(u - center) >= observed - 1e-10).ravel()


def run_permutation(study, plan, output, repetitions, resume):
    checkpoint = Checkpoints(output, 'permutation', signature(plan, 'permutation'), resume)
    variants = ['within_race'] + (['within_race_sex_age'] if study.blocked else [])
    counts = checkpoint.arrays('counts')[-1].copy() if checkpoint.state['chunks'] else np.zeros((len(variants), study.total), dtype=np.int64)
    progress = []
    for rep in range(checkpoint.state['completed'], repetitions):
        for variant in range(len(variants)):
            permutation_counts(study, plan['seed'], rep, bool(variant), counts[variant])
        if (rep + 1) % 1000 == 0 or rep + 1 == repetitions:
            for variant in range(len(variants)):
                p = empirical_p(counts[variant], rep + 1)
                progress.append([rep + 1, variant, int((bh(p) < .05).sum()), 1 / (rep + 2), int((counts[variant] == 0).sum())])
            checkpoint.append(rep + 1, counts=counts, progress=np.asarray(progress))
            progress = []
            print(f'GMDB label permutations: {rep + 1}/{repetitions}; {len(variants)} separate nulls', flush=True)
    completed = checkpoint.state['completed']
    tables = []
    for variant, name in enumerate(variants):
        frame = study.primary[KEYS].copy()
        ns, others, us, gs, ps, centers = [np.empty(study.total) for _ in range(6)]
        for race, data in study.races.items():
            where = data['features_index']
            ns[where] = np.repeat(data['ns'], len(study.features))
            others[where] = len(data['x']) - ns[where]
            us[where], gs[where], ps[where] = data['observed_u'].ravel(), data['observed_g'].ravel(), data['observed_p'].ravel()
            centers[where] = data['block_center'].ravel() if variant else ns[where] * others[where] / 2
        frame['n_case'], frame['n_other'], frame['observed_U'], frame['observed_hedges_g'] = ns.astype(int), others.astype(int), us, gs
        frame['p_asymptotic_or_exact_auto'] = ps
        frame['q_global_BH_asymptotic_or_exact_auto'] = bh(ps)
        frame['null_specification'] = name
        frame['null_expected_U'] = centers
        frame['supported'] = (ns >= 5) & (others >= 2)
        frame['support_status'] = np.where(frame.supported, 'tested', 'insufficient_other_GMDB_patients')
        if variant:
            for race, data in study.races.items():
                immobile = np.repeat(data['movable_case_counts'] == 0, len(study.features))
                where = data['features_index'][immobile]
                frame.loc[where, 'supported'] = False
                frame.loc[where, 'support_status'] = 'no_movable_case_labels_within_blocks'
        p = empirical_p(counts[variant], completed)
        p[~frame.supported] = 1
        frame['exceedances'] = counts[variant]
        frame['B'] = completed
        frame['p_empirical_permutation'] = p
        frame['q_global_BH_empirical'] = bh(p)
        frame['p_value_floor'] = 1 / (completed + 1)
        alpha = .05 / (2 * study.total)
        k = counts[variant]
        lower = np.zeros(study.total)
        upper = np.ones(study.total)
        lower[k > 0] = beta.ppf(alpha, k[k > 0], completed - k[k > 0] + 1)
        upper[k < completed] = beta.ppf(1 - alpha, k[k < completed] + 1, completed - k[k < completed])
        lower, upper = (1 + completed * lower) / (1 + completed), (1 + completed * upper) / (1 + completed)
        lower[~frame.supported] = upper[~frame.supported] = 1
        frame['q_MC_simultaneous_lower'] = bh(lower)
        frame['q_MC_simultaneous_upper'] = bh(upper)
        frame['MC_FDR_decision_ambiguous'] = frame.q_MC_simultaneous_lower.lt(.05) & frame.q_MC_simultaneous_upper.ge(.05)
        frame['main_g'], frame['main_q'] = study.primary.hedges_g.to_numpy(), study.primary.q.to_numpy()
        frame['same_effect_direction_as_main'] = np.sign(frame.observed_hedges_g) == np.sign(frame.main_g)
        tables.append(frame)
    result = pd.concat(tables, ignore_index=True)
    result.to_csv(output / 'gmdb_permutation_aggregate.csv', index=False)
    progress = pd.DataFrame(np.vstack(checkpoint.arrays('progress')), columns=['B', 'null_index', 'discoveries', 'p_floor', 'zero_exceedances'])
    progress['null_specification'] = progress.null_index.astype(int).map(dict(enumerate(variants)))
    progress.to_csv(output / 'gmdb_permutation_progress.csv', index=False)
    return result

import hashlib
import importlib.metadata
import json
import os
import subprocess
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import gammaln
from scipy.stats import mannwhitneyu, norm, rankdata

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from revision_results.race_stratified_20261115 import common as baseline

HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / 'local_review/race_stratified_sensitivity'
RACES = baseline.RACES
KEYS = ['disease', 'canonical_race', 'feature']
os.umask(0o077)


def write_json(path, value):
    path = Path(path)
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + '\n')
    os.replace(temporary, path)


def bh(values):
    p = np.asarray(values, float)
    if p.ndim != 1 or not len(p) or not np.isfinite(p).all() or ((p < 0) | (p > 1)).any():
        raise ValueError('Full planned BH family required; unsupported tests must have p=1')
    order = np.argsort(p, kind='stable')
    q = np.empty_like(p)
    q[order] = np.minimum(1, np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1])
    return q


def rng_for(seed, stream, repetition, race, variant=0):
    return np.random.default_rng(np.random.SeedSequence([seed * 100 + stream, variant, repetition, race]))


def effect(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if min(len(x), len(y)) < 2 or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Effects require finite measurements and n>=2')
    return effect_moments(len(x), len(y), x.mean(0), x.var(0, ddof=1), y.mean(0), y.var(0, ddof=1))


def effect_moments(n, m, xmean, xvar, ymean, yvar):
    df = n + m - 2
    pooled = ((n - 1) * xvar + (m - 1) * yvar) / df
    correction = np.exp(gammaln(df / 2) - .5 * np.log(df / 2) - gammaln((df - 1) / 2))
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(pooled > 0, correction * (xmean - ymean) / np.sqrt(pooled), np.nan)


@lru_cache(maxsize=64)
def exact_p_lookup(n, m):
    u = np.arange(n * m // 2 + 1)
    b = u[:, None] // n + (np.arange(n)[None, :] >= n - u[:, None] % n)
    x = b + np.arange(n)[None, :]
    y = np.empty((len(u), m), int)
    for i, indices in enumerate(x):
        y[i] = np.setdiff1d(np.arange(n + m), indices, assume_unique=True)
    tests = mannwhitneyu(x, y, axis=1, alternative='two-sided', method='exact')
    np.testing.assert_array_equal(tests.statistic, u)
    return tests.pvalue


def p_from_untied_u(u, n, m):
    u = np.asarray(u, float)
    if min(n, m) <= 8:
        if not np.array_equal(u, np.rint(u)):
            raise ValueError('Untied exact U must be integer')
        return exact_p_lookup(min(n, m), max(n, m))[np.minimum(u, n * m - u).astype(int)]
    z = (np.maximum(u, n * m - u) - n * m / 2 - .5) / np.sqrt(n * m * (n + m + 1) / 12)
    return np.minimum(1, 2 * norm.sf(z))


def tests_direct(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if x.ndim != 2 or y.ndim != 2 or x.shape[1] != y.shape[1]:
        raise ValueError('Expected two compatible numeric feature matrices')
    g = effect(x, y)
    joined = np.r_[x, y]
    ranks = rankdata(joined, axis=0)
    u = ranks[:len(x)].sum(0) - len(x) * (len(x) + 1) / 2
    tied = (np.diff(np.sort(joined, axis=0), axis=0) == 0).any(0)
    p = np.empty(x.shape[1])
    if (~tied).any():
        p[~tied] = p_from_untied_u(u[~tied], len(x), len(y))
    for j in np.flatnonzero(tied):
        p[j] = mannwhitneyu(x[:, j], y[:, j], method='asymptotic', alternative='two-sided').pvalue
    return u, p, g


def placebo_allocation(capacity, sizes, reference_n, rng):
    if capacity >= sum(sizes) + reference_n:
        draw = rng.permutation(capacity)[:sum(sizes) + reference_n]
        reference = draw[:reference_n]
        cases = np.split(draw[reference_n:], np.cumsum(sizes)[:-1])
        assert len(np.unique(np.r_[reference, *cases])) == reference_n + sum(sizes)
        return 'shared_disjoint', [(case, reference) for case in cases]
    pairs = []
    for n in sizes:
        if capacity < n + reference_n:
            pairs.append(None)
            continue
        draw = rng.choice(capacity, n + reference_n, replace=False)
        pairs.append((draw[:n], draw[n:]))
    return 'per_pair_disjoint_fallback', pairs


def permuted_order(size, rng, blocks=None):
    if blocks is None:
        return rng.permutation(size)
    result = np.arange(size)
    for indices in blocks:
        result[indices] = rng.permutation(indices)
    return result


def empirical_p(exceedances, repetitions):
    counts = np.asarray(exceedances)
    if repetitions < 1 or (counts < 0).any() or (counts > repetitions).any():
        raise ValueError('Invalid permutation counts')
    return (1 + counts) / (1 + repetitions)


def source_signature(paths):
    return {str(path.relative_to(ROOT)): baseline.sha(path) for path in sorted(set(paths))}


def load_baseline(directory=None):
    config_path = ROOT / 'SENSITIVITY_ANALYSIS.json'
    config = json.loads(config_path.read_text())
    directory = Path(directory) if directory else ROOT / config['baseline_private_output']
    required = ['statistics_manifest.json', 'statistics_status.json', 'primary_all_tests.csv',
                'fairface_sampling_summary.json', 'PRIVATE_reference_selection.csv', 'cohort_and_race_audit.csv']
    missing = [str(directory / name) for name in required if not (directory / name).is_file()]
    if missing:
        raise FileNotFoundError('BLOCKED: new race-stratified baseline unavailable: ' + ', '.join(missing))
    manifest = json.loads((directory / 'statistics_manifest.json').read_text())
    if manifest.get('study_id') != 'PRISM_RACE_STRATIFIED_20261115':
        raise ValueError('BLOCKED: old pooled results cannot substitute for the new baseline')
    for field in ['source_sha256', 'code_sha256']:
        for name, expected in manifest[field].items():
            if baseline.sha(ROOT / name) != expected:
                raise ValueError('BLOCKED: baseline provenance mismatch: ' + name)
    for name in required[1:]:
        if baseline.sha(directory / name) != manifest['output_sha256'][name]:
            raise ValueError('BLOCKED: baseline artifact changed: ' + name)
    disease, reference, features, definitions, audit, sources = baseline.load_sources()
    counts = baseline.cohort_counts(disease)
    eligible = counts.loc[counts.eligible].sort_values(['disease', 'canonical_race'])
    primary = pd.read_csv(directory / 'primary_all_tests.csv', float_precision='round_trip')
    ordered_keys = pd.DataFrame([dict(disease=row.disease, canonical_race=row.canonical_race, feature=feature)
                                for row in eligible.itertuples() for feature in features])
    if len(primary) != len(ordered_keys) or primary[KEYS].duplicated().any():
        raise ValueError('Baseline full planned family is inconsistent')
    primary = ordered_keys.merge(primary, on=KEYS, validate='one_to_one')
    np.testing.assert_allclose(bh(primary.p), primary.q, rtol=1e-12, atol=1e-15)
    sampling = json.loads((directory / 'fairface_sampling_summary.json').read_text())
    sampled = baseline.sample_reference(reference, sampling['quotas'], sampling['seed'])
    selection = pd.read_csv(directory / 'PRIVATE_reference_selection.csv')
    if sampled.file.tolist() != selection.file.tolist():
        raise ValueError('Locked primary sample cannot be replayed')
    for cell in eligible.itertuples():
        frame = primary.loc[primary.disease.eq(cell.disease) & primary.canonical_race.eq(cell.canonical_race)]
        if not frame.case_n.eq(cell.n_current_gmdb).all() or not frame.reference_n.eq(sampling['quotas'][cell.canonical_race]).all():
            raise ValueError('Baseline case/reference sizes do not match the locked design')
    sources += [config_path] + [directory / name for name in required]
    return disease, reference, features, definitions, audit, counts, primary, sampling, sources


class Study:
    def __init__(self, loaded):
        self.disease, self.reference, self.features, self.definitions, self.audit, self.counts, self.primary, self.sampling, self.sources = loaded
        self.cells = self.counts.loc[self.counts.eligible].sort_values(['disease', 'canonical_race']).reset_index(drop=True)
        self.total = len(self.primary)
        self.races = {}
        self.covariates = []
        for race in RACES:
            patients = self.disease.loc[self.disease.canonical_race.eq(race)].sort_values(['disease', 'patient_id'])
            reference = self.reference.loc[self.reference.canonical_race.eq(race)].sort_values('file')
            x, y = patients[self.features].to_numpy(float), reference[self.features].to_numpy(float)
            labels, starts, sizes = np.unique(patients.disease, return_index=True, return_counts=True)
            rows = self.cells.index[self.cells.canonical_race.eq(race)].to_numpy()
            positions = np.array([np.flatnonzero(labels == self.cells.disease.iloc[i])[0] for i in rows])
            feature_indices = (rows[:, None] * len(self.features) + np.arange(len(self.features))).ravel()
            ranks = rankdata(x, axis=0)
            bands = pd.cut(patients.conditional_age_years, [3, 5, 7, 10], right=False).astype(str)
            block_labels = patients.sex_label.fillna('Unknown').astype(str) + '|' + bands
            blocks = [np.flatnonzero(block_labels.to_numpy() == label) for label in sorted(block_labels.unique())]
            center = np.zeros((len(rows), len(self.features)))
            for block in blocks:
                for j, pos in enumerate(positions):
                    center[j] += np.sum(patients.disease.to_numpy()[block] == labels[pos]) * ranks[block].mean(0)
            ns = sizes[positions]
            center -= (ns * (ns + 1) / 2)[:, None]
            observed_u = np.add.reduceat(ranks, starts, axis=0)[positions] - (ns * (ns + 1) / 2)[:, None]
            observed_g, observed_p = [], []
            for cell_index, pos in zip(rows, positions):
                take = patients.disease.eq(labels[pos]).to_numpy()
                if np.sum(~take) >= 2:
                    u, p, g = tests_direct(x[take], x[~take])
                else:
                    u = observed_u[len(observed_g)]
                    p, g = np.ones(len(self.features)), np.full(len(self.features), np.nan)
                np.testing.assert_allclose(u, observed_u[len(observed_g)], atol=0)
                observed_g.append(g)
                observed_p.append(p)
                case, other = patients.loc[take], patients.loc[~take]
                if len(other) >= 2:
                    pooled_sd = np.sqrt((case.conditional_age_years.var() + other.conditional_age_years.var()) / 2)
                    difference = abs(case.conditional_age_years.mean() - other.conditional_age_years.mean())
                    age_smd = difference / pooled_sd if pooled_sd > 0 else (np.inf if difference else 0.)
                    sex_gap = max(abs(case.sex_label.fillna('Unknown').eq(sex).mean()
                                      - other.sex_label.fillna('Unknown').eq(sex).mean()) for sex in ['Female', 'Male', 'Unknown'])
                else:
                    age_smd, sex_gap = np.nan, np.nan
                self.covariates.append({'disease': labels[pos], 'canonical_race': race, 'age_abs_smd': age_smd,
                                        'sex_max_absolute_fraction_gap': sex_gap,
                                        'substantial_imbalance': bool(age_smd >= .5 or sex_gap >= .2)})
            movable = [block for block in blocks if len(np.unique(patients.disease.to_numpy()[block])) > 1]
            movable_indices = np.concatenate(movable) if movable else np.array([], dtype=int)
            case_mobility = np.array([np.sum(patients.disease.to_numpy()[movable_indices] == labels[pos]) for pos in positions])
            self.races[race] = {'x': x, 'y': y, 'patients': patients, 'rows': rows, 'positions': positions,
                                'features_index': feature_indices, 'labels': labels, 'starts': starts, 'sizes': sizes,
                                'ns': ns, 'ranks': ranks, 'blocks': movable, 'block_center': center,
                                'observed_u': observed_u, 'observed_g': np.asarray(observed_g),
                                'observed_p': np.asarray(observed_p),
                                'movable_patients': sum(map(len, movable)), 'movable_case_counts': case_mobility,
                                'all_block_count': len(blocks),
                                'immobile_patients': len(x) - sum(map(len, movable))}
        self.blocked = any(row['substantial_imbalance'] for row in self.covariates)

    def plan(self, seed, repetitions, permutation_repetitions):
        races = {}
        for race, data in self.races.items():
            n = self.sampling['quotas'][race]
            sizes = data['ns']
            races[race] = {'gmdb_all_mapped': len(data['x']), 'eligible_patients': int(sizes.sum()),
                           'eligible_cells': len(sizes), 'fairface_pool': len(data['y']), 'locked_reference_n': n,
                           'selection_fraction': n / len(data['y']),
                           'placebo_full_disjoint_shortage': max(0, int(sizes.sum()) + n - len(data['y'])),
                           'placebo_pair_unsupported': int(np.sum(sizes + n > len(data['y']))),
                           'movable_blocked_patients': data['movable_patients'],
                           'immobile_blocked_patients': data['immobile_patients'], 'all_blocks': data['all_block_count']}
        paths = list(HERE.rglob('*.py')) + list(HERE.glob('*.md')) + [ROOT / 'stats_testing/reference_sensitivity.py']
        return {'study_id': 'PRISM_THREE_RACE_SENSITIVITY_20261115', 'baseline_id': 'PRISM_RACE_STRATIFIED_20261115',
                'feature_count': len(self.features), 'features': self.features,
                'feature_order_sha256': hashlib.sha256('\n'.join(self.features).encode()).hexdigest(),
                'planned_hypotheses': self.total, 'eligible_cells': len(self.cells), 'races': races,
                'reference': self.sampling, 'seed': seed,
                'analysis_seeds': {name: seed * 100 + i for i, name in enumerate(['placebo', 'permutation', 'reference-resampling'], 1)},
                'repetition_seed': 'SeedSequence([analysis_seed, variant, zero_based_repetition, race_index])',
                'requested_repetitions': repetitions, 'requested_permutation_repetitions': permutation_repetitions,
                'fixed_permutation_milestones': [1000, 10000, 100000],
                'permutation_cap': 100000, 'blocked_enabled': self.blocked,
                'blocked_trigger': 'Any cell age absolute SMD>=0.5 or sex category fraction gap>=0.20, evaluated before permutations',
                'blocked_strata': 'recorded sex including Unknown x conditional age [3,5), [5,7), [7,10)',
                'blocked_center': 'Exact conditional expected U from fixed syndrome counts and mean feature ranks within each block',
                'nulls': {'placebo': 'same-source FairFace disjoint pseudo-cases/references; descriptive calibration',
                          'permutation': 'multiclass syndrome labels exchangeable within recorded race',
                          'permutation_blocked': 'multiclass syndrome labels exchangeable within race/sex/age-band blocks',
                          'reference-resampling': 'fixed GMDB cases; locked reference allocation; selection stability, not a null'},
                'correction': 'Full planned family separately per repetition/empirical null; unsupported p=1',
                'statistics': 'Primary auto-equivalent two-sided U; exact gamma J pooled-SD Hedges g; finite features required',
                'empirical_p': '(1 + count(abs(U_perm-center)>=abs(U_observed-center)))/(B+1); ties included with 1e-10 tolerance',
                'stability_rule': 'Descriptive q<.05 frequency>=.80 AND sign agreement>=.95; robust shortlist also main q<.05 and empirical q<.05 in both enabled GMDB nulls',
                'source_sha256': source_signature(self.sources), 'code_sha256': source_signature(paths),
                'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                'versions': {name: importlib.metadata.version(name) for name in ['numpy', 'pandas', 'scipy', 'matplotlib', 'pyarrow']},
                'privacy': 'All results, checkpoints and small cells remain Git-ignored; no publication action',
                'crosswalk': 'Reuse fixed broad primary mapping, including Middle Eastern/West Asian within Asian; unsupported groups excluded'}


class Checkpoints:
    def __init__(self, output, stage, signature, resume):
        self.directory = output / 'checkpoints' / stage
        self.directory.mkdir(parents=True, exist_ok=True)
        self.path = self.directory / 'state.json'
        self.state = {'signature': signature, 'completed': 0, 'chunks': []}
        if self.path.exists():
            if not resume:
                raise ValueError('Existing checkpoint requires --resume or a fresh output directory')
            self.state = json.loads(self.path.read_text())
            if self.state['signature'] != signature:
                raise ValueError('Checkpoint design/source/code signature changed; use a fresh run')
            expected = 0
            for chunk in self.state['chunks']:
                if chunk['start'] != expected or baseline.sha(self.directory / chunk['file']) != chunk['sha256']:
                    raise ValueError('Checkpoint corruption, gap or duplicate iterations')
                expected = chunk['stop']
            if expected != self.state['completed']:
                raise ValueError('Checkpoint completed count inconsistent')

    def append(self, stop, **arrays):
        start = self.state['completed']
        if stop <= start:
            raise ValueError('Checkpoint iterations must increase')
        name = f'{start:06d}_{stop:06d}.npz'
        temporary = self.directory / (name + '.tmp')
        with temporary.open('wb') as handle:
            np.savez(handle, **arrays)
        os.replace(temporary, self.directory / name)
        self.state['chunks'].append({'start': start, 'stop': stop, 'file': name,
                                      'sha256': baseline.sha(self.directory / name)})
        self.state['completed'] = stop
        write_json(self.path, self.state)

    def arrays(self, field):
        return [np.load(self.directory / item['file'])[field] for item in self.state['chunks']]


def signature(plan, stage):
    design = {key: value for key, value in plan.items()
              if key not in ['requested_repetitions', 'requested_permutation_repetitions', 'revision']}
    return hashlib.sha256(json.dumps({'design': design, 'stage': stage}, sort_keys=True).encode()).hexdigest()

import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from research.race_stratified_sensitivity import core as locked

ROOT = locked.ROOT
HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT/'local_review/race_stratified_pitch_sensitivity'
RACES = locked.RACES
ARMS = ('random', 'pitch_imbalanced')
SEED = 2026111704
POSE_FIELDS = ('pose_pitch', 'pose_yaw', 'pose_roll')
METRICS = [f'{group}_abs_pitch_{stat}' for group in ['case', 'reference']
           for stat in ['mean', 'median', 'q25', 'q75', 'min', 'max']] + [
    'mean_abs_pitch_gap', 'median_abs_pitch_gap', 'pitch_range_overlap', 'minimum_abs_pitch_gap',
    'case_fraction_above_reference_max', 'signed_pitch_mean_gap', 'signed_pitch_median_gap',
    'signed_yaw_mean_gap', 'abs_yaw_mean_gap', 'signed_roll_mean_gap', 'abs_roll_mean_gap',
    'female_fraction_gap', 'male_fraction_gap', 'unknown_sex_fraction_gap']


class Blocked(ValueError):
    def __init__(self, message, audit=None):
        super().__init__('BLOCKED: '+message)
        self.audit = audit


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def rng_for(seed, repetition, category, arm):
    return np.random.default_rng(np.random.SeedSequence([seed, repetition, category, arm]))


def pose_audit(frame, thresholds, verified_units):
    rows = []
    for race in RACES:
        group = frame.loc[frame.canonical_race.eq(race)]
        for field in POSE_FIELDS:
            values = pd.to_numeric(group[field], errors='coerce').to_numpy(float) if field in group else np.full(len(group), np.nan)
            finite = values[np.isfinite(values)]
            row = {'canonical_race': race, 'field': field, 'field_present': field in group,
                   'n_total': len(group), 'n_finite': len(finite), 'n_nonfinite': len(values)-len(finite),
                   'missing_fraction': 1-len(finite)/len(values) if len(values) else np.nan,
                   'units': verified_units, 'accepted_absolute_limit': thresholds.get(field.removeprefix('pose_')),
                   'distinct_values': len(np.unique(finite)), 'tied_observations': len(finite)-len(np.unique(finite))}
            for prefix, data in [('', finite), ('absolute_', np.abs(finite))]:
                for stat, function in [('min', np.min), ('max', np.max), ('mean', np.mean), ('median', np.median)]:
                    row[prefix+stat] = float(function(data)) if len(data) else np.nan
                row[prefix+'q25'], row[prefix+'q75'] = np.quantile(data, [.25, .75]) if len(data) else (np.nan, np.nan)
            rows.append(row)
    return pd.DataFrame(rows)


def pitch_halves(pitch, keys, verified_units):
    pitch, keys = np.asarray(pitch, float), np.asarray(keys, str)
    if verified_units != 'degrees':
        raise Blocked('Head-pitch units/convention have not been verified')
    if len(pitch) != len(keys) or not len(pitch) or not np.isfinite(pitch).all():
        raise Blocked('Missing/nonfinite head pitch or missing category pool')
    if len(np.unique(keys)) != len(keys) or np.any(np.char.str_len(keys) == 0):
        raise Blocked('Stable private image keys must be unique and nonmissing')
    order = np.lexsort((keys, np.abs(pitch)))
    cut = len(order)//2
    return order[:cut], order[cut:]


def require_pose(frame, thresholds, units):
    audit = pose_audit(frame, thresholds, units)
    if units != 'degrees' or not set(POSE_FIELDS).issubset(frame.columns):
        raise Blocked('Missing pose_pitch/yaw/roll or unverified pose units', audit)
    for field in POSE_FIELDS:
        values = pd.to_numeric(frame[field], errors='coerce').to_numpy(float)
        if not np.isfinite(values).all():
            raise Blocked(f'Nonfinite {field}; accepted-pool capacity must be resolved before outcomes', audit)
        if np.any(np.abs(values) > thresholds[field.removeprefix('pose_')]):
            raise Blocked(f'{field} violates the existing extraction/QC threshold; no silent exclusions', audit)
    return audit


def support_table(cells, pools, quotas):
    rows = []
    for race in RACES:
        group = cells.loc[cells.canonical_race.eq(race)]
        total, low, high = pools[race]
        nref = quotas[race]
        sizes = group.n_current_gmdb.to_numpy(int)
        nominal = (sizes+nref <= total) & (nref >= 20)
        split = (sizes <= high) & (nref <= low) & (nref >= 20)
        common = nominal & split
        shared_a = bool(common.any() and sizes[common].sum()+nref <= total)
        shared_b = bool(common.any() and sizes[common].sum() <= high and nref <= low)
        for i, cell in enumerate(group.itertuples()):
            reasons = []
            if nref < 20:
                reasons.append('reference_quota_below_20')
            if not nominal[i]:
                reasons.append('random_insufficient_disjoint_pool')
            if nref > low:
                reasons.append('low_pitch_insufficient_reference_capacity')
            if sizes[i] > high:
                reasons.append('high_pitch_insufficient_case_capacity')
            rows.append({'cell_index': cell.Index, 'disease': cell.disease, 'canonical_race': race,
                         'n_case': int(sizes[i]), 'n_reference': nref, 'pool_total': total,
                         'low_pitch_pool': low, 'high_pitch_pool': high, 'odd_element_assignment': 'high_half',
                         'nominal_random_supported': bool(nominal[i]), 'pitch_split_supported': bool(split[i]),
                         'supported': bool(common[i]), 'newly_unsupported_after_split': bool(nominal[i] and not split[i]),
                         'support_reason': 'supported' if common[i] else ';'.join(reasons),
                         'random_design': ('shared_disjoint' if shared_a else 'per_pair_disjoint_fallback') if common[i] else 'unsupported',
                         'pitch_design': ('shared_disjoint' if shared_b else 'per_pair_disjoint_fallback') if common[i] else 'unsupported',
                         'shared_random_all_cells_feasible': bool(sizes.sum()+nref <= total),
                         'shared_pitch_all_cells_feasible': bool(sizes.sum() <= high and nref <= low),
                         'cross_pair_reuse_allowed_random': bool(common[i] and not shared_a),
                         'cross_pair_reuse_allowed_pitch': bool(common[i] and not shared_b)})
    return pd.DataFrame(rows).sort_values('cell_index').reset_index(drop=True)


def allocations(capacity, low, high, sizes, reference_n, rng, arm):
    sizes = np.asarray(sizes, int)
    if arm not in ARMS or not len(sizes) or np.any(sizes < 2) or reference_n < 2:
        raise ValueError('A supported arm, positive group sizes and reference are required')
    if arm == 'random':
        design, pairs = locked.placebo_allocation(capacity, sizes, reference_n, rng)
    else:
        if len(low) < reference_n or np.any(sizes > len(high)):
            raise ValueError('Unsupported pitch contrast cannot be sampled')
        if sum(sizes) <= len(high):
            reference = rng.choice(low, reference_n, replace=False)
            cases = np.split(rng.permutation(high)[:sum(sizes)], np.cumsum(sizes)[:-1])
            design, pairs = 'shared_disjoint', [(case, reference) for case in cases]
        else:
            design = 'per_pair_disjoint_fallback'
            pairs = [(rng.choice(high, n, replace=False), rng.choice(low, reference_n, replace=False)) for n in sizes]
    for n, pair in zip(sizes, pairs):
        if pair is None:
            raise ValueError('Frozen common support was not respected')
        case, reference = pair
        if (len(case) != n or len(reference) != reference_n or len(np.unique(case)) != n
                or len(np.unique(reference)) != reference_n or np.intersect1d(case, reference).size):
            raise ValueError('Invalid case/reference sizes, replacement or overlap')
        if arm == 'pitch_imbalanced' and (not np.isin(case, high).all() or not np.isin(reference, low).all()):
            raise ValueError('Pitch half assignment violated')
    if design == 'shared_disjoint':
        cases = np.concatenate([pair[0] for pair in pairs])
        if len(np.unique(cases)) != sum(sizes) or np.intersect1d(cases, pairs[0][1]).size:
            raise ValueError('Shared reference or pseudo-cases overlap across pairs')
    return design, pairs


def pose_metrics(pose, sex, case, reference):
    x, y = pose[case], pose[reference]
    a, b = np.abs(x[:, 0]), np.abs(y[:, 0])
    result = []
    for values in [a, b]:
        result.extend([values.mean(), np.median(values), *np.quantile(values, [.25, .75]), values.min(), values.max()])
    result.extend([a.mean()-b.mean(), np.median(a)-np.median(b),
                   max(0, min(a.max(), b.max())-max(a.min(), b.min())), a.min()-b.max(), np.mean(a > b.max()),
                   x[:, 0].mean()-y[:, 0].mean(), np.median(x[:, 0])-np.median(y[:, 0]),
                   x[:, 1].mean()-y[:, 1].mean(), np.abs(x[:, 1]).mean()-np.abs(y[:, 1]).mean(),
                   x[:, 2].mean()-y[:, 2].mean(), np.abs(x[:, 2]).mean()-np.abs(y[:, 2]).mean()])
    result.extend([np.mean(sex[case] == label)-np.mean(sex[reference] == label) for label in ['Female', 'Male', 'Unknown']])
    return np.asarray(result)


def full_family(p_supported, g_supported, supported, features):
    if features != list(dict.fromkeys(features)) or len(features) != 120 or len(supported) != 9600:
        raise ValueError('Exact 120-feature, 9,600-test family required')
    if len(p_supported) != np.sum(supported) or len(g_supported) != np.sum(supported):
        raise ValueError('Supported vector does not match the frozen mask')
    p, g = np.ones(9600), np.full(9600, np.nan)
    p[supported], g[supported] = p_supported, g_supported
    q = locked.bh(p)
    discoveries = supported & (q < .05)
    return p, g, q, discoveries


def protected_hashes(sources):
    paths = set(sources)
    paths.update(ROOT/name for name in subprocess.check_output(['git', 'ls-files'], cwd=ROOT, text=True).splitlines()
                 if (ROOT/name).is_file() and not (ROOT/name).is_relative_to(HERE))
    directories = ['research', 'stats_results', 'manuscript', 'publication', 'stats_testing',
                   'revision_results/rebuild_120_20261028', 'revision_results/reviewer_followup_2026',
                   'local_review/race_stratified_20261115', 'local_review/race_stratified_sensitivity',
                   'local_review/reference_stability_diagnostic', 'local_review/sensitivity_summary_figure']
    for name in directories:
        paths.update(p for p in (ROOT/name).rglob('*') if p.is_file() and not p.is_relative_to(HERE)
                     and not any(part in {'__pycache__', '.pytest_cache', 'local_only'} for part in p.parts)
                     and p.name != '.DS_Store')
    return locked.source_signature(paths)


def load_study():
    loaded = locked.load_baseline()
    disease, reference, features, definitions, audit, counts, primary, sampling, sources = loaded
    cells = counts.loc[counts.eligible].sort_values(['disease', 'canonical_race']).reset_index(drop=True)
    status = json.loads((locked.baseline.PRIVATE/'statistics_status.json').read_text())
    expected = (67, 80, 736, 120, 9600, 982)
    observed = (cells.disease.nunique(), len(cells), int(cells.n_current_gmdb.sum()), len(features), len(primary), int(primary.q.lt(.05).sum()))
    if observed != expected or status['status'] != 'COMPLETED':
        raise Blocked(f'Locked baseline dimensions/status mismatch: {observed}')
    if sampling['quotas'] != {'African': 52, 'Asian': 190, 'European': 647} or sampling['X'] != 889:
        raise Blocked('Locked reference allocation differs')
    manifest = json.loads((locked.baseline.PRIVATE/'statistics_manifest.json').read_text())
    if manifest['features'] != features or definitions.feature.tolist() != features:
        raise Blocked('Baseline feature order mismatch')
    expected_keys = [(cell.disease, cell.canonical_race, feature) for cell in cells.itertuples() for feature in features]
    if list(primary[locked.KEYS].itertuples(index=False, name=None)) != expected_keys:
        raise Blocked('Baseline hypothesis key order mismatch')
    config = json.loads((ROOT/'ACTIVE_ANALYSIS.json').read_text())
    provenance_path = locked.baseline.EXTRACTION/'fairface/provenance.json'
    provenance = json.loads(provenance_path.read_text())
    extractor_path = ROOT/'FaceKit/src/facekit/core/geometric/extractor.py'
    if locked.baseline.sha(extractor_path) != provenance['extractor_sha256']:
        raise Blocked('Pose extractor hash differs from pinned extraction')
    pinned = subprocess.check_output(['git', '-C', str(ROOT/'FaceKit'), 'show',
                                      config['facekit_commit']+':src/facekit/core/geometric/extractor.py'])
    if hashlib.sha256(pinned).hexdigest() != provenance['extractor_sha256']:
        raise Blocked('Pose source differs from pinned FaceKit commit')
    if b'np.degrees(pitch)' not in pinned or b'pitch = np.arctan2(R[2, 1], R[2, 2])' not in pinned:
        raise Blocked('Unverified head-pitch coordinate conversion')
    thresholds = provenance['thresholds_degrees']
    availability = require_pose(reference.loc[reference.canonical_race.isin(RACES)], thresholds, 'degrees')
    raw_path = locked.baseline.EXTRACTION/'fairface/fairface_120_all_eligible_measurements.parquet'
    raw = pd.read_parquet(raw_path)
    required = {'image_path_key', 'quality_passed', 'reason', 'facekit_commit', 'model_sha256',
                'dataset', 'image_sha256', *POSE_FIELDS, *features}
    if not required.issubset(raw.columns):
        raise Blocked('Original FairFace schema missing columns: '+str(sorted(required-set(raw.columns))))
    if not raw.image_path_key.is_unique or raw.image_path_key.isna().any():
        raise Blocked('Original FairFace extraction has duplicate keys')
    ordered = raw.set_index('image_path_key').loc[reference.file]
    if (not ordered.quality_passed.eq(True).all() or not ordered.reason.eq('accepted').all()
            or not ordered.facekit_commit.eq(config['facekit_commit']).all()
            or not ordered.model_sha256.eq(config['model_sha256']).all() or not ordered.dataset.eq('FairFace').all()):
        raise Blocked('Pose extraction/QC provenance mismatch')
    np.testing.assert_array_equal(reference.image_sha256.to_numpy(), ordered.image_sha256.to_numpy())
    np.testing.assert_array_equal(reference[list(POSE_FIELDS)+features].to_numpy(float),
                                  ordered[list(POSE_FIELDS)+features].to_numpy(float))
    benchmark_path = ROOT/'research/race_stratified_sensitivity_results_20261115/summary.json'
    benchmark = json.loads(benchmark_path.read_text())
    if (benchmark['status'] != 'COMPLETED' or benchmark['baseline_id'] != 'PRISM_RACE_STRATIFIED_20261115'
            or benchmark['planned_hypotheses'] != 9600 or benchmark['placebo']['B'] != 1000
            or benchmark['observed_main_discoveries'] != 982):
        raise Blocked('Historical placebo benchmark differs from locked study')
    sources += [raw_path, extractor_path, benchmark_path, ROOT/'research/race_stratified_sensitivity/core.py',
                ROOT/'research/race_stratified_sensitivity/experiments.py']
    groups, capacities = {}, {}
    for race in RACES:
        group = reference.loc[reference.canonical_race.eq(race)].sort_values('file').reset_index(drop=True)
        low, high = pitch_halves(group.pose_pitch, group.file, 'degrees')
        groups[race] = {'frame': group, 'values': group[features].to_numpy(float),
                        'pose': group[list(POSE_FIELDS)].to_numpy(float),
                        'sex': group.sex_label.fillna('Unknown').where(group.sex_label.isin(['Female', 'Male']), 'Unknown').to_numpy(),
                        'low': low, 'high': high}
        capacities[race] = (len(group), len(low), len(high))
    support = support_table(cells, capacities, sampling['quotas'])
    supported = np.repeat(support.supported.to_numpy(bool), 120)
    if int(supported.sum()) != benchmark['placebo']['supported_hypotheses'] or not supported.any():
        raise Blocked('Historical/common support differs; review the support audit before outcomes')
    pose_provenance = {'units': 'degrees', 'signed': True,
                       'coordinate_convention': 'Pinned MediaPipe facial transformation matrix; small-angle XYZ Tait-Bryan; pitch atan2(R[2,1],R[2,2]); gimbal branch atan2(-R[1,2],R[1,1])',
                       'magnitude': 'Absolute deviation from extractor zero; sign is a coordinate convention, not independently verified physical up/down',
                       'extraction': 'Stored estimates from the pinned FaceKit/MediaPipe extraction; no re-estimation',
                       'uncertainty': 'No image-specific pose ground truth or uncertainty estimates available',
                       'thresholds_degrees': thresholds, 'additional_exclusions': 0,
                       'extractor_sha256': provenance['extractor_sha256'], 'facekit_commit': config['facekit_commit']}
    return {'cells': cells, 'features': features, 'definitions': definitions, 'primary': primary,
            'groups': groups, 'support': support, 'supported': supported, 'quotas': sampling['quotas'],
            'availability': availability, 'pose_provenance': pose_provenance, 'sources': list(set(sources)),
            'benchmark': benchmark['placebo'], 'baseline_audit': audit}


def analysis_plan(study, seed, repetitions):
    numerical_code = [HERE/'design.py', HERE/'run_pitch_imbalance.py']
    return {'study_id': 'PRISM_RACE_STRATIFIED_HEAD_PITCH_SENSITIVITY_20261117',
            'baseline_id': 'PRISM_RACE_STRATIFIED_20261115', 'seed': seed, 'repetitions_per_arm': repetitions,
            'arms': list(ARMS), 'features': study['features'], 'planned_family': 9600,
            'key_order_sha256': digest(study['primary'][locked.KEYS].values.tolist()),
            'supported_mask_sha256': digest(study['supported'].tolist()),
            'support_sha256': digest(study['support'].to_dict('records')),
            'supported_cells': int(study['support'].supported.sum()), 'supported_tests': int(study['supported'].sum()),
            'quotas': study['quotas'], 'pose_provenance': study['pose_provenance'],
            'split_rule': 'Sort absolute pose_pitch ascending; ties by lexical file; lower floor(N/2), odd element in high half',
            'allocation_rule': 'Shared reference and globally disjoint cases when feasible; otherwise per-pair disjoint without-replacement draws, reuse across pairs allowed',
            'random_arm': 'Accepted same-category full pool; no post-hoc pose matching',
            'imbalanced_arm': 'High absolute-pitch cases, low absolute-pitch reference; quotas and common mask fixed',
            'rng_rule': 'NumPy SeedSequence([seed, repetition_zero_based, category_index, arm_index]); independent substreams',
            'case_n_rule': 'Exact locked eligible syndrome-category case count; no GMDB feature vector used in either arm',
            'statistics': 'Original per-feature two-sided MW: exact if min n<=8 and untied, otherwise tie/continuity-corrected asymptotic; exact-gamma-J pooled-SD Hedges g',
            'correction': 'One BH family of 9,600 p per arm/repetition; unsupported p=1; discovery=supported & q<0.05',
            'gap_metrics': METRICS, 'optional_analyses': 'None; no reversed or pitch-matched arm',
            'source_sha256': locked.source_signature(study['sources']),
            'numerical_code_sha256': locked.source_signature(numerical_code),
            'versions': {name: locked.importlib.metadata.version(name) for name in ['numpy', 'pandas', 'scipy']},
            'privacy': 'All results and checkpoints remain ignored; no image extraction, upload, manuscript or selector edits'}


def evaluate(study, seed, repetition, arm_index):
    p, g = np.ones(9600), np.full(9600, np.nan)
    support = study['support'].loc[study['support'].supported]
    metrics = np.empty((len(support), len(METRICS)))
    position = {cell: i for i, cell in enumerate(support.cell_index)}
    selection = []
    for race_index, race in enumerate(RACES):
        rows = support.loc[support.canonical_race.eq(race)]
        if rows.empty:
            continue
        group = study['groups'][race]
        design, pairs = allocations(len(group['values']), group['low'], group['high'], rows.n_case.to_numpy(int),
                                     study['quotas'][race], rng_for(seed, repetition, race_index, arm_index), ARMS[arm_index])
        expected = rows.random_design.iloc[0] if arm_index == 0 else rows.pitch_design.iloc[0]
        if design != expected:
            raise ValueError('Allocation differs from frozen plan')
        for row, (case, reference) in zip(rows.itertuples(), pairs):
            pair_metrics = pose_metrics(group['pose'], group['sex'], case, reference)
            if arm_index == 1 and not (pair_metrics[12] > 0 and pair_metrics[13] > 0 and pair_metrics[15] >= 0):
                raise Blocked('Achieved pitch separation was not in the intended direction; no redrawing')
            metrics[position[row.cell_index]] = pair_metrics
            where = slice(row.cell_index*120, (row.cell_index+1)*120)
            _, p[where], g[where] = locked.tests_direct(group['values'][case], group['values'][reference])
            selection.append([row.cell_index, case.tolist(), reference.tolist()])
    return p[study['supported']], g[study['supported']], metrics, digest(selection)

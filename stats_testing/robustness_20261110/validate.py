import argparse
import importlib.util
import json
import sys
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from .common import (ROOT, DEFAULT_OUT, KEY, bh, features, load_private, private_output, protect_baseline,
                     require, sha, write_json)


def independent_bh(pvalues):
    pvalues = np.asarray(pvalues)
    order = np.argsort(pvalues, kind='stable')
    q = np.empty(len(pvalues))
    q[order] = np.minimum(1, np.minimum.accumulate((pvalues[order] * len(pvalues) / np.arange(1, len(pvalues) + 1))[::-1])[::-1])
    return q


def full_baseline_readonly(out):
    source = ROOT / 'revision_results/rebuild_120_20261028/validate.py'
    spec = importlib.util.spec_from_file_location('frozen_full_validation', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    expected_json, expected_csv = module.RUN / 'VALIDATION_FULL.json', module.RUN / 'validation_full_summary.csv'
    before = {str(path): sha(path) for path in [expected_json, expected_csv]}
    original_write, original_csv = Path.write_text, pd.DataFrame.to_csv
    outputs = out / 'validation'
    outputs.mkdir(exist_ok=True, mode=0o700)
    routed = []
    def route_write(path, data, *args, **kwargs):
        require(path == expected_json, 'Unexpected write attempted by frozen full validator')
        routed.append('json')
        return original_write(outputs / 'frozen_primary_full_validation.json', data, *args, **kwargs)
    def route_csv(frame, path=None, *args, **kwargs):
        require(Path(path) == expected_csv, 'Unexpected CSV write attempted by frozen full validator')
        routed.append('csv')
        return original_csv(frame, outputs / 'frozen_primary_full_validation_summary.csv', *args, **kwargs)
    old_argv = sys.argv[:]
    sys.argv = [str(source), '--scope', 'full']
    try:
        with patch.object(Path, 'write_text', route_write), patch.object(pd.DataFrame, 'to_csv', route_csv):
            module.main()
    finally:
        sys.argv = old_argv
    require(routed == ['json', 'csv'], 'Expected both full-validation reports to route to new private output')
    require(all(sha(Path(path)) == digest for path, digest in before.items()), 'Existing validation record changed')
    result = json.loads((outputs / 'frozen_primary_full_validation.json').read_text())
    require(result['status'] == 'VALIDATED' and len(result['components']) == 15, 'Frozen full private validation did not pass')
    write_json({'status': 'PASS', 'original_validator_sha256': sha(source), 'numerical_functions_modified': False,
                'write_routing_only': routed, 'existing_validation_report_hashes_unchanged': before,
                'components': len(result['components'])}, outputs / 'READONLY_REPLAY_ATTESTATION.json')


def new_outputs(out):
    real, reference, raw, alternate, paths = load_private()
    names = features().feature.tolist()
    folder = out / 'aggregate_candidates'
    points = pd.read_csv(folder / 'P1_patient_resampling.csv')
    original = pd.read_csv(out / 'evidence/P0_exhaustive_evidence_matrix.csv', low_memory=False)
    for comparator, column in [('FairFace', 'ff_hedges_g'), ('other_GMDB', 'gmdb_hedges_g')]:
        joined = points.loc[points.comparator.eq(comparator)].merge(original[['syndrome', 'feature', column]], on=['syndrome', 'feature'], validate='one_to_one')
        require(np.allclose(joined.point_g, joined[column], rtol=1e-10, atol=1e-12), 'Patient-resampling point estimate differs from frozen output')
    photo = np.load(out / 'checkpoints/exact_photo_recovery.npz')
    require(photo['p'].shape == photo['q'].shape == (100, 5760), 'Photo replay grid shape differs')
    for p, q in zip(photo['p'], photo['q']):
        require(np.allclose(independent_bh(p), q, rtol=1e-10, atol=1e-12), 'Recovered photo global BH differs')
    cases = pd.read_csv(folder / 'P1_normalization_cases.csv')
    require(len(cases) == 12 and cases.planned_hypotheses.eq(12).all(), 'Geometry family denominator differs')
    require(np.allclose(independent_bh(cases.p_value), cases.exploratory_q, rtol=1e-10, atol=1e-12), 'Geometry BH differs')
    restricted = pd.read_csv(folder / 'P0_quality_restrictions.csv')
    for threshold, part in restricted.groupby('threshold'):
        require(len(part) == 5760, 'Quality restriction lost planned cells')
        expected = independent_bh(np.where(part.supported, part.p_value, 1))
        require(part.loc[~part.supported, 'exploratory_q'].isna().all(), 'Unsupported restrictions have a displayed q')
        require(np.allclose(part.loc[part.supported, 'exploratory_q'], expected[part.supported], rtol=1e-10, atol=1e-12), 'Restriction BH differs')
    adjusted = pd.read_csv(folder / 'P0_quality_adjusted_cases.csv')
    require(len(adjusted) == 12, 'Quality-adjusted case family differs')
    require(np.allclose(independent_bh(np.where(adjusted.supported, adjusted.p_value, 1))[adjusted.supported],
                        adjusted.loc[adjusted.supported, 'exploratory_q'], rtol=1e-10, atol=1e-12), 'Adjusted-case BH differs')
    for source in ['FairFace', 'GMDB_selected']:
        coassign = pd.read_csv(folder / f'P1_coclustering_{source}.csv')[names].to_numpy()
        require(np.allclose(coassign, coassign.T) and np.allclose(np.diag(coassign), 1), 'Invalid co-clustering stability matrix')
        require(((coassign >= 0) & (coassign <= 1)).all(), 'Invalid co-clustering probabilities')
    negative = pd.read_csv(folder / 'P0_negative_control_repetitions.csv')
    require(len(negative) == 200 and negative.hypotheses.eq(5760).all(), 'Negative control family/repetitions differ')
    require(negative.groupby('design').size().eq(100).all(), 'Negative-control designs are incomplete')
    from .source_bias import control_indices
    strata = (reference.race.astype(str) + '|' + reference.gender.astype(str)).to_numpy()
    counts = real.disease.value_counts().loc[lambda x: x >= 10].sort_index().to_numpy()
    for index, design in enumerate(['stratified_random', 'pitch_imbalanced']):
        groups, ri, _ = control_indices(strata, reference.abs_pitch.to_numpy(), counts,
                                       np.random.default_rng([20261110, 10, index, 0]), design)
        pvalues = []
        values = reference[names].to_numpy(float)
        for group in groups:
            pvalues.extend(mannwhitneyu(values[group], values[ri], axis=0, method='asymptotic', alternative='two-sided').pvalue)
        expected = int((independent_bh(pvalues) < .05).sum())
        observed = negative.loc[negative.design.eq(design) & negative.repetition.eq(0), 'discoveries'].iloc[0]
        require(expected == observed, 'Negative-control first repetition did not reproduce')
    for path in out.rglob('*.csv'):
        columns = set(pd.read_csv(path, nrows=0).columns)
        require(not columns.intersection({'patient_id', 'image_id', 'image_path', 'landmarks', 'image_sha256', 'embedding', 'patient_index'}),
                'Protected individual field appeared in an aggregate draft')
    figures = list((out / 'figures').glob('*.json'))
    require(len(figures) == 15, 'Private figure set incomplete')
    for path in figures:
        item = json.loads(path.read_text())
        require(item['private_review_draft'] is True, 'Figure release classification missing')
        require(item['code_sha256'] == sha(ROOT / 'stats_testing/robustness_20261110/plots.py'), 'Figure code hash differs')
        require(all(sha(out / source) == digest for source, digest in item['source_sha256'].items()), 'Figure source changed')
        require(all(sha(out / 'figures' / output) == digest for output, digest in item['outputs'].items()), 'Figure bytes changed')
    claims = pd.read_csv(out / 'P3_claim_to_evidence_ledger.csv')
    require(len(claims) == 14 and claims.claim_id.is_unique, 'Claim ledger incomplete')
    require(all(sha(out / row.source) == row.source_sha256 for row in claims.itertuples()), 'Claim source hash differs')
    opportunities = pd.read_csv(folder / 'P1_primary_photo_opportunities.csv')
    require(opportunities.changed_alternate.sum() == 54 and opportunities.patients_with_multiple_eligible_images.sum() == 94,
            'Photograph opportunity claims differ')
    changes = pd.read_csv(folder / 'P1_random_photo_change_counts.csv').changed_primary_photos
    require((changes.min(), changes.median(), changes.max()) == (38, 53, 65), 'Random-photo count summary differs')
    require(negative.loc[negative.design.eq('stratified_random'), 'discoveries'].eq(0).all(), 'Random control summary differs')
    imbalanced = negative.loc[negative.design.eq('pitch_imbalanced'), 'discoveries']
    require((imbalanced.min(), imbalanced.median(), imbalanced.max()) == (0, 25, 73), 'Pitch-control summary differs')
    clinical = pd.read_csv(out / 'P2_phenotype_source_matrix.csv')
    require(len(clinical) == 17 and not clinical.exact_comparable_measure.any(), 'Clinical audit denominator/scope differs')
    seed = pd.read_csv(folder / 'P1_ranking_seed_sensitivity.csv')
    require(len(seed) == 3 and seed.selected.eq(20).all() and seed.shared_with_original_high_stability.eq(20).all(),
            'Seed sensitivity report differs')
    age = pd.read_csv(folder / 'P1_age_case_summary.csv')
    require(len(age) == 12 and age.status.eq('ok').all() and age.fdr_planned_hypotheses.eq(11760).all(), 'Age case support/family differs')
    age_counts = age.groupby('stratum').q_global_bh.apply(lambda values: values.lt(.05).sum()).to_dict()
    require(age_counts == {'age_3_to_under_6': 0, 'age_6_to_under_10': 3}, 'Age case interpretation differs')
    wbs = adjusted.loc[adjusted.syndrome.str.startswith('WILLIAMS-BEUREN') & adjusted.comparator.eq('other_GMDB')]
    require(len(wbs) == 2 and wbs.supported.all() and wbs.exploratory_q.gt(.05).all(), 'WBS quality-model interpretation differs')
    result = {'status': 'PASS', 'photo_BH_replays': 100, 'photo_hypotheses_each': 5760, 'patient_point_estimates': len(points),
              'geometry_family': 12, 'quality_restriction_families': 2, 'quality_adjusted_family': 12,
              'negative_control_first_replays': 2, 'aggregate_privacy_schema_checks': True, 'figure_hash_checks': len(figures),
              'claim_source_checks': len(claims), 'source_clinical_questions': len(clinical), 'seed_sensitivity_checks': len(seed),
              'frozen_files_unchanged': protect_baseline(out)}
    write_json(result, out / 'validation/NEW_RESEARCH_VALIDATION.json')
    print(json.dumps(result, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUT)
    parser.add_argument('--full-private-baseline', action='store_true')
    args = parser.parse_args()
    sys.dont_write_bytecode = True
    out = private_output(args.output_dir)
    protect_baseline(out)
    new_outputs(out)
    if args.full_private_baseline:
        full_baseline_readonly(out)
    protect_baseline(out)


if __name__ == '__main__':
    main()

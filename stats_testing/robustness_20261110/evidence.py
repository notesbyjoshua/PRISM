import argparse
import json
import platform
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import fisher_exact, pearsonr, spearmanr

from .common import (ROOT, DEFAULT_OUT, KEY, CASES, SEED, ResearchContractError, bh, canonicalize, case_rows,
                     features, keyed_merge, load_private, offline, primary_table, private_output, protect_baseline,
                     require, sha, source_table, write_aggregate, write_json)


def preflight(out):
    out = private_output(out)
    count = protect_baseline(out, create=True)
    config = json.loads((ROOT / 'ACTIVE_ANALYSIS.json').read_text())
    require(config['analysis_id'] == 'PRISM_120_20261028' and config['feature_count'] == 120
            and config['synthetic_rows_in_primary'] == 0, 'Active study selection differs from the handout')
    definitions = features()
    commit = subprocess.check_output(['git', '-C', str(ROOT / 'FaceKit'), 'rev-parse', 'HEAD'], text=True).strip()
    extractor = ROOT / 'FaceKit/src/facekit/core/geometric/extractor.py'
    require(commit == config['facekit_commit'] and definitions.source_commit.eq(commit).all(), 'FaceKit source pin differs')
    require(definitions.extractor_sha256.eq(sha(extractor)).all(), 'Pinned extractor checksum differs')
    model = Path.home() / '.cache/facekit/face_landmarker.task'
    if model.exists():
        require(sha(model) == config['model_sha256'], 'Local model checksum differs')
    primary, within = primary_table(), source_table('within_gmdb_results.csv')
    require(len(primary) == 5760 and primary.syndrome.nunique() == 48 and primary.feature.nunique() == 120,
            'Primary hypothesis counts differ')
    require(primary.groupby('canonical_syndrome_id').feature.nunique().eq(120).all(), 'Primary ordered schema coverage differs')
    require(primary.n_reference.eq(3409).all() and primary.n_disease.ge(10).all(), 'Primary sample counts differ')
    require(np.allclose(primary.mannwhitney_q_fdr_global, bh(primary.mannwhitney_p), rtol=1e-10, atol=1e-12), 'Primary BH differs')
    require(len(within) == 5760 and within.estimable.all(), 'Current within-GMDB grid differs')
    require(set(map(tuple, primary[KEY].to_numpy())) == set(map(tuple, within[KEY].to_numpy())), 'Canonical grids differ')
    a = set(map(tuple, primary.loc[primary.mannwhitney_q_fdr_global.lt(.05), KEY].to_numpy()))
    b = set(map(tuple, within.loc[within.q_global_bh.lt(.05), KEY].to_numpy()))
    counts = dict(features=120, syndromes=48, hypotheses=len(primary), primary_significant=len(a),
                  within_significant=len(b), shared=len(a & b), fairface_only=len(a - b), within_only=len(b - a))
    require(list(counts.values()) == [120, 48, 5760, 955, 222, 173, 782, 49], 'Frozen comparator baseline disagrees; stop')
    input_status, input_hashes = 'AVAILABLE_AUTHORIZED_LOCAL_MEASUREMENTS', {}
    try:
        real, reference, raw, alternate, paths = load_private()
        input_hashes = {name: {'path': str(path.relative_to(ROOT)), 'sha256': sha(path)} for name, path in paths.items()}
        private_counts = {'selected_real': len(real), 'primary_patients': 872, 'reference_paths': len(reference),
                          'one_record_per_recorded_patient': bool(real.patient_id.is_unique), 'synthetic_primary_rows': 0,
                          'source_dimensions_constant': bool(real.width.eq(224).all() and reference.width.eq(224).all())}
    except ResearchContractError as error:
        if 'BLOCKED_PRIVATE_INPUT' not in str(error):
            raise
        input_status, private_counts = 'BLOCKED_PRIVATE_INPUT', {}
    validation = {'status': 'PASS', 'git_head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  'analysis_id': config['analysis_id'], 'frozen_counts': counts, 'private_status': input_status,
                  'private_counts': private_counts, 'facekit_commit': commit, 'extractor_sha256': sha(extractor),
                  'model_sha256': config['model_sha256'], 'model_asset_present_and_verified': model.exists(),
                  'protected_input_sources': input_hashes, 'active_config_sha256': sha(ROOT / 'ACTIVE_ANALYSIS.json'),
                  'baseline_files_protected': count, 'publication_permission': 'NOT_GRANTED_FOR_NEW_OUTPUTS; private drafts only',
                  'versions': {'python': platform.python_version(), 'numpy': np.__version__, 'pandas': pd.__version__, 'scipy': scipy.__version__}}
    write_json(validation, out / '00_baseline_validation.json')
    text = f"# Preflight, November 10, 2026\n\nHEAD: `{validation['git_head']}`. Active study: `{config['analysis_id']}`. "
    text += f"Baseline checks PASS: {counts}. The active config and {count} pre-existing tracked/user-owned files are protected by checksums.\n\n"
    text += 'Exact config:\n\n```json\n' + json.dumps(config, indent=2) + '\n```\n\n'
    text += ('Current ordered S6 and feature formulas match exactly. FaceKit commit/extractor and the local model (if present) are verified. '
             'Canonical joins resolve full current internal syndrome IDs through the unique shared taxonomy, checked against private metadata; display labels are never join keys between results.\n\n')
    text += f'Private input status: {input_status}. Exact safe path classifications and digests are in 00_baseline_validation.json, without individual values.\n\n'
    text += ('Available commands: local_review/ci_120_fix/venv/bin/python (public requirements and PyArrow), .venv/bin/python (existing full validator), '
             'offline Python/SciPy/statsmodels/Matplotlib computations, existing aggregate verification and tests. No new GPU job or extractor run is required.\n\n'
             'Limitations: source and case/reference status are not separately identifiable; FairFace has categorical age and unverified identities/health. '
             'Both processed image sources are 224×224, which says nothing about original resolution/crop conventions. Detection confidence, crop margins, brightness/exposure '
             'and compatible saved landmarks are not established. Real-image transformations remain blocked pending explicit image-reprocessing permission; no pixels will be read. '
             'PDIDB remains excluded. New data/figures remain private and unpushed; manuscript and public manifest stay unchanged.\n')
    (out / '00_preflight.md').write_text(text)
    plan = (ROOT / 'stats_testing/robustness_20261110/ANALYSIS_PLAN.md').read_text()
    (out / '00_work_plan.md').write_text('# Work plan and resource estimate\n\n'
        'Existing computations reused: primary and all-other-GMDB 5,760 grids, three current alternate/reference specifications, 100 random-photo summaries, '
        'full/prespecified adjusted and age families, reference feature clusters, and the separate synthetic fidelity audit. '
        'Missing per-photo p/q/rank traces may require exact recovery of the same 100 draws, not a new photo design.\n\n'
        'New CPU-only tasks: evidence/quality inventory (seconds); 200 negative-control repetitions, correlation bootstraps and exact photo recovery (minutes); '
        'case-patient/subsample stability and geometry alternatives (minutes). Expected RAM below 2 GB and private outputs below 200 MB. '
        'Expensive repetitions checkpoint by code/input/parameter signature. Blockers: raw components/real-image perturbations require compatible saved landmarks or explicit permission. '
        'P2 new-method development is conditional; do not create a composite score without a defensible target and held-out benchmark.\n\n' + plan)
    print(json.dumps({'preflight': 'PASS', 'private_status': input_status, 'counts': counts}))
    return validation


def evidence_matrix(out):
    protect_baseline(out)
    matrix = primary_table()
    rename = {name: 'ff_' + name for name in ['n_disease', 'n_reference', 'hedges_g', 'hedges_g_ci_95_low',
              'hedges_g_ci_95_high', 'mannwhitney_p', 'mannwhitney_q_fdr_global']}
    matrix = matrix.rename(columns=rename)
    matrix['ff_supported'] = np.isfinite(matrix.ff_hedges_g) & matrix.ff_n_disease.ge(10)
    within = source_table('within_gmdb_results.csv')
    columns = ['n_target', 'n_comparator', 'estimable', 'hedges_g', 'ci_95_lower', 'ci_95_upper', 'p_value', 'q_global_bh']
    matrix = keyed_merge(matrix, within, 'gmdb_', columns)
    for name, folder in [('alternate', 'alternative_photo'), ('all_reference', 'reference_all_eligible'), ('sex_reference', 'reference_race_sex')]:
        source = canonicalize(pd.read_csv(ROOT / f'stats_results/sensitivities/{folder}/stats_each_disease_vs_healthy.csv'), 'disease')
        matrix = keyed_merge(matrix, source, name + '_', ['n_disease', 'n_healthy', 'hedges_g', 'hedges_g_ci_95_low',
                                                        'hedges_g_ci_95_high', 'mannwhitney_p', 'mannwhitney_q_fdr_global'])
    matrix = keyed_merge(matrix, source_table('random_photo_sensitivity.csv'), 'photo_',
        ['replicates', 'median_g', 'g_q025', 'g_q975', 'min_g', 'max_g', 'positive_g_frequency', 'significant_frequency'])
    adjusted = source_table('full_adjusted_screen.csv')
    matrix = keyed_merge(matrix, adjusted, 'screen_', ['n_target', 'n_other', 'n_complete', 'target_beta', 'ci_95_lower',
        'ci_95_upper', 'p_value', 'q_global_bh', 'status', 'reason', 'fdr_family', 'fdr_planned_hypotheses'])
    four = source_table('within_gmdb_age_sex_adjusted.csv')
    matrix = keyed_merge(matrix, four, 'four_', ['target_beta', 'n_complete', 'ci_95_lower', 'ci_95_upper', 'p_value', 'q_four_models'], False)
    three = pd.read_csv(ROOT / 'publication/supplementary/three_predefined_reference_adjusted_models.csv')
    exact = matrix[['syndrome', 'canonical_syndrome_id']].drop_duplicates().copy()
    exact['display'] = exact.syndrome.str.split('__').str[0]
    require(exact.display.is_unique, 'Nonunique short names in prespecified reference source')
    three = three.merge(exact.rename(columns={'display': 'syndrome', 'syndrome': 'full_syndrome'}), on='syndrome', validate='many_to_one')
    require(three.canonical_syndrome_id.notna().all(), 'Prespecified reference model lacks a canonical match')
    matrix = keyed_merge(matrix, three, 'three_', ['adjusted_disease_coefficient', 'n_total', 'ci_95_lower', 'ci_95_upper',
                                                 'raw_p', 'bh_q_across_three_predefined_models'], False)
    definitions = features().drop(columns=['order'])
    matrix = matrix.merge(definitions, on='feature', validate='many_to_one').copy()
    matrix['same_direction'] = matrix.gmdb_estimable & np.sign(matrix.ff_hedges_g).eq(np.sign(matrix.gmdb_hedges_g))
    matrix['intervals_overlap'] = ((matrix.ff_hedges_g_ci_95_low <= matrix.gmdb_ci_95_upper)
                                   & (matrix.gmdb_ci_95_lower <= matrix.ff_hedges_g_ci_95_high))
    matrix['interval_comparison'] = 'marginal comparator-specific interval overlap; dependent estimates; not a difference test'
    matrix['category'] = np.select([matrix.ff_mannwhitney_q_fdr_global.lt(.05) & matrix.gmdb_q_global_bh.lt(.05),
                                   matrix.ff_mannwhitney_q_fdr_global.lt(.05), matrix.gmdb_q_global_bh.lt(.05)],
                                  ['shared', 'FairFace-only', 'within-GMDB-only'], default='neither')
    matrix['apparently_strong'] = matrix.ff_mannwhitney_q_fdr_global.lt(.05) & matrix.ff_hedges_g.abs().ge(.8)
    matrix['magnitude_halved'] = matrix.gmdb_hedges_g.abs().lt(.5 * matrix.ff_hedges_g.abs())
    matrix['sign_reversal'] = (matrix.ff_hedges_g * matrix.gmdb_hedges_g).lt(0)
    matrix['magnitude_change'] = matrix.gmdb_hedges_g.abs() - matrix.ff_hedges_g.abs()
    matrix['reference_min_g'] = matrix[['ff_hedges_g', 'all_reference_hedges_g', 'sex_reference_hedges_g']].min(axis=1)
    matrix['reference_max_g'] = matrix[['ff_hedges_g', 'all_reference_hedges_g', 'sex_reference_hedges_g']].max(axis=1)
    matrix['reference_sign_consistent'] = matrix.reference_min_g.mul(matrix.reference_max_g).gt(0)
    matrix['same_direction_photo_frequency'] = np.where(matrix.ff_hedges_g.gt(0), matrix.photo_positive_g_frequency,
                                                       1 - matrix.photo_positive_g_frequency)
    matrix['clinical_correspondence'] = 'unmapped; not evidence of novelty'
    clinical = canonicalize(pd.read_csv(ROOT / 'publication/final_research_120/clinical_literature_audit_17.csv'))
    matrix = keyed_merge(matrix, clinical, 'clinical_', ['clinical_support', 'exact_ratio_validated', 'retain_supported_clinical_example'], False)
    matrix.loc[matrix.clinical_clinical_support.notna(), 'clinical_correspondence'] = matrix.clinical_clinical_support
    matrix = matrix.sort_values(KEY).reset_index(drop=True)
    write_aggregate(matrix, out / 'evidence/P0_exhaustive_evidence_matrix.csv')
    pairs = canonicalize(pd.read_csv(ROOT / 'publication/supplementary/targeted_pairwise.csv'), 'syndrome_a')
    from .common import canonical_map
    mapping_b = canonical_map().rename(columns={'syndrome': 'syndrome_b', 'canonical_syndrome_id': 'comparator_syndrome_id'})
    pairs = pairs.merge(mapping_b, on='syndrome_b', how='left', validate='many_to_one')
    require(pairs.comparator_syndrome_id.notna().all(), 'Targeted comparator lacks a canonical ID')
    write_aggregate(pairs, out / 'evidence/P0_targeted_specificity_long.csv')
    age = source_table('age_stratified.csv')
    write_aggregate(age, out / 'evidence/P0_age_stratum_long.csv')
    categories = matrix.groupby('category', sort=True).size().rename('hypotheses').reset_index()
    require(categories.set_index('category').hypotheses.to_dict() == {'FairFace-only': 782, 'neither': 4756,
            'shared': 173, 'within-GMDB-only': 49}, 'Comparator categories differ')
    write_aggregate(categories, out / 'aggregate_candidates/P0_overlap_categories.csv')
    summaries = []
    for factor in ['all', 'anatomical_region', 'syndrome', 'size_band']:
        if factor == 'size_band':
            matrix['size_band'] = pd.cut(matrix.ff_n_disease, [9, 19, 39, np.inf], labels=['10–19', '20–39', '40+']).astype(str)
        groups = [('all', matrix)] if factor == 'all' else matrix.groupby(factor, observed=True, sort=True)
        for label, part in groups:
            strong = part.loc[part.apparently_strong]
            summaries.append({'factor': factor, 'group': label, 'hypotheses': len(part), 'same_direction': int(part.same_direction.sum()),
                'same_direction_fraction': part.same_direction.mean(), 'pearson': pearsonr(part.ff_hedges_g, part.gmdb_hedges_g).statistic,
                'spearman': spearmanr(part.ff_hedges_g, part.gmdb_hedges_g).statistic,
                'median_absolute_g_difference': (part.ff_hedges_g - part.gmdb_hedges_g).abs().median(),
                'marginal_interval_overlap_fraction': part.intervals_overlap.mean(), 'strong_pairs': len(strong),
                'strong_halved': int(strong.magnitude_halved.sum()), 'strong_sign_reversals': int(strong.sign_reversal.sum()),
                'strong_not_significant_within': int(strong.gmdb_q_global_bh.ge(.05).sum()),
                'estimand_note': 'dependent comparator contrasts; no clinical false-positive interpretation'})
    write_aggregate(pd.DataFrame(summaries), out / 'aggregate_candidates/P0_concordance_strata.csv')
    enrichments = []
    for factor in ['feature_family', 'anatomical_region', 'syndrome']:
        for label in sorted(matrix[factor].unique()):
            inside, shared = matrix[factor].eq(label), matrix.category.eq('shared')
            cells = [int((inside & shared).sum()), int((inside & ~shared).sum()),
                     int((~inside & shared).sum()), int((~inside & ~shared).sum())]
            odds, p = fisher_exact(np.array(cells).reshape(2, 2))
            enrichments.append(dict(factor=factor, group=label, shared=cells[0], tested=int(inside.sum()),
                                    odds_ratio=odds, p_value=p, ecological=True))
    enrich = pd.DataFrame(enrichments)
    enrich['exploratory_q'] = enrich.groupby('factor').p_value.transform(lambda p: bh(p))
    enrich['warning'] = 'Exploratory; correlated feature cells and shared patients violate independent-table interpretation'
    write_aggregate(enrich, out / 'aggregate_candidates/P0_exploratory_overlap_enrichment.csv')
    discordant = matrix.loc[matrix.apparently_strong & (matrix.sign_reversal | matrix.magnitude_halved)]
    discordant = discordant.assign(discrepancy=(discordant.ff_hedges_g - discordant.gmdb_hedges_g).abs())
    discordant = discordant.sort_values(['discrepancy'] + KEY, ascending=[False, True, True]).head(20)
    write_aggregate(discordant, out / 'aggregate_candidates/P0_discordant_candidates.csv')
    write_aggregate(case_rows(matrix), out / 'aggregate_candidates/P0_predeclared_cases.csv')
    overall = summaries[0]
    report = '# P0 comparator evidence\n\n'
    report += (f"All 5,760 canonical syndrome–feature hypotheses are represented: 955 FairFace discoveries, 222 within-GMDB, 173 shared, "
               f"782 FairFace-only, 49 within-only and 4,756 neither. {overall['same_direction']}/5,760 "
               f"({overall['same_direction_fraction']:.1%}) share an effect direction across all supported pairs. "
               f"Pearson r={overall['pearson']:.6g}; Spearman rho={overall['spearman']:.6g}.\n\n"
               f"Among {overall['strong_pairs']} apparently strong primary associations (global q<.05, |g|>=.8), "
               f"{overall['strong_halved']} have within-GMDB magnitude below half, {overall['strong_sign_reversals']} reverse sign, "
               f"and {overall['strong_not_significant_within']} do not pass the separate within-GMDB q threshold. "
               'Neither loss of significance nor magnitude attenuation estimates a false-positive rate. Different reference means/variances and the syndrome-versus-other-conditions estimand can explain nonagreement; '
               'source and case status cannot be separated in the FairFace contrast.\n\n')
    report += ('Evidence matrix: evidence/P0_exhaustive_evidence_matrix.csv. The full adjusted 5,880 family and four-/three-model families are explicitly distinguished. '
               'Targeted 22,320-cell and age 11,760-cell families are retained as separate canonical long tables, not averaged into one q. '
               'Absent targeted/prespecified models remain missing, while unsupported full-screen estimates/q remain missing with their original support reason. '
               'Intervals are original marginal source intervals; overlapping intervals are descriptive and cannot establish equality or independent replication. '
               'No excluded concordance-bootstrap uncertainty is restored.\n\n'
               'The overlap enrichment tables are exploratory summaries with separate BH correction per factor; correlated hypotheses and shared cases limit Fisher-test interpretation. '
               'A small q is not evidence for independent anatomical enrichment. Candidate and discordance tables retain counterexamples. '
               'All outputs remain private review drafts; no new aggregate release is assumed.\n')
    (out / 'P0_comparator_evidence.md').write_text(report)
    return matrix


def quality_inventory(out):
    rows, balance = [], []
    try:
        real, reference, raw, alternate, paths = load_private()
    except ResearchContractError as error:
        if 'BLOCKED_PRIVATE_INPUT' not in str(error):
            raise
        (out / 'P0_source_bias_audit.md').write_text('# Source-bias audit\n\nBLOCKED_PRIVATE_INPUT. Source/case confounding is not identifiable from aggregate effects. No quality values fabricated.\n')
        return
    available = ['pose_yaw', 'pose_pitch', 'pose_roll', 'abs_yaw', 'abs_pitch', 'abs_roll', 'pose_score',
                 'width', 'height', 'screen_detected_faces', 'landmark_count', 'derotated', 'mouth_opening', 'eye_fissure_height_mean']
    missing = ['face_box', 'crop_margin', 'detection_confidence', 'brightness', 'contrast', 'exposure', 'native_resolution']
    for source, frame in [('GMDB_selected', real), ('FairFace_reference', reference)]:
        for field in available + missing:
            present = field in frame
            values = pd.to_numeric(frame[field], errors='coerce') if present else pd.Series(dtype=float)
            rows.append(dict(source=source, proxy=field, present=present, total_n=len(frame), observed_n=int(values.notna().sum()),
                             missing_n=int(len(frame) - values.notna().sum()), distinct_values=int(values.nunique()),
                             status='observed' if present else 'not measured; no pixel processing authorized'))
        for group, part in [('all', frame)] + list(frame.groupby('broad_group', dropna=False, sort=True)):
            for field in available:
                if field not in part:
                    continue
                values = pd.to_numeric(part[field], errors='coerce').astype(float)
                balance.append(dict(source=source, group=str(group), proxy=field, n=len(part), mean=values.mean(),
                                    sd=values.std(ddof=1), q1=values.quantile(.25), median=values.median(), q3=values.quantile(.75)))
    write_aggregate(pd.DataFrame(rows), out / 'aggregate_candidates/P0_proxy_inventory.csv')
    write_aggregate(pd.DataFrame(balance), out / 'aggregate_candidates/P0_acquisition_proxy_balance.csv')
    matrix = pd.read_csv(out / 'evidence/P0_exhaustive_evidence_matrix.csv', dtype={'canonical_syndrome_id': 'string'})
    ecological = []
    for syndrome, part in real.groupby('disease', sort=True):
        if len(part) >= 10:
            for field in ['abs_yaw', 'abs_pitch', 'abs_roll', 'pose_score', 'mouth_opening', 'eye_fissure_height_mean']:
                ecological.append(dict(syndrome=syndrome, proxy=field, n=len(part), cohort_mean=part[field].mean(),
                                       cohort_median=part[field].median()))
    ecological = canonicalize(pd.DataFrame(ecological))
    cohort = matrix.groupby('canonical_syndrome_id', sort=True).agg(shared=('category', lambda s: (s == 'shared').sum()),
                                                                  significant=('ff_mannwhitney_q_fdr_global', lambda s: (s < .05).sum())).reset_index()
    ecological = ecological.merge(cohort, on='canonical_syndrome_id', validate='many_to_one')
    write_aggregate(ecological, out / 'aggregate_candidates/P0_ecological_quality_overlap.csv')
    text = ('# P0 source-bias inventory and identification limits\n\n'
            'GMDB selected photographs → FaceKit → case measurements; FairFace reference paths → the SAME FaceKit → reference measurements. '
            'Dataset source and case status align in the primary contrast; no source coefficient can identify a separate disease effect. '
            'Within-GMDB compares syndromes to other conditions, not normal controls.\n\n'
            'Both sources have 224×224 processed image dimensions, signed/absolute pose and normalized pose score. The detector count and landmark count are constant after acceptance. '
            'These constants do not establish equal original resolution, margins, illumination or image restoration. Crop/face-box bounds, confidence, native resolution and exposure markers are absent. '
            'Mouth opening and eye-fissure height are expression-sensitive morphology, so controlling them may remove biology as well as acquisition variation. '
            'No new pixel-based proxies or transformations were performed: permission and compatible saved landmarks are not established.\n\n'
            'P0_proxy_inventory.csv documents presence/missingness; P0_acquisition_proxy_balance.csv reports all-source and recorded broad-group summaries. '
            'P0_ecological_quality_overlap.csv relates cohort-average proxies to overlap counts; its unit is the syndrome, not independent feature cells. '
            'These descriptive distributions cannot prove or rule out source bias. Negative controls and restricted/model sensitivities follow only after the Milestone A checkpoint.\n')
    (out / 'P0_source_bias_audit.md').write_text(text)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUT)
    parser.add_argument('--preflight-only', action='store_true')
    args = parser.parse_args()
    offline()
    try:
        out = private_output(args.output_dir)
        preflight(out)
        if not args.preflight_only:
            evidence_matrix(out)
            quality_inventory(out)
        protect_baseline(out)
        print('Milestone A complete; protected inputs and frozen baseline unchanged')
    except (ResearchContractError, FileNotFoundError, KeyError, AssertionError, TypeError, ValueError) as error:
        print('Research stage stopped: ' + (str(error) if isinstance(error, ResearchContractError) else type(error).__name__))
        raise SystemExit(1)


if __name__ == '__main__':
    main()

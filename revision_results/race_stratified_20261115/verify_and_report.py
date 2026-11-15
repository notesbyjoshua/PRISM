import argparse
import json
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.metrics import confusion_matrix
from statsmodels.stats.multitest import multipletests

from common import (PRIVATE, RACES, ROOT, apportion, cohort_counts, load_sources, matched_frames, private_output,
                    protected_hashes, sample_reference, select_ml, sha, write_json)
from run_five_syndrome_ml import metrics


def read_csv(output, name):
    return pd.read_csv(output / name, float_precision='round_trip')


def md_table(frame):
    columns = frame.columns.tolist()
    lines = ['| ' + ' | '.join(columns) + ' |', '| ' + ' | '.join(['---'] * len(columns)) + ' |']
    for row in frame.itertuples(index=False, name=None):
        cells = []
        for value in row:
            if isinstance(value, float):
                value = f'{value:.4g}'
            cells.append(str(value).replace('|', '\\|').replace('\n', ' '))
        lines.append('| ' + ' | '.join(cells) + ' |')
    return '\n'.join(lines)


def verify(output):
    disease, reference, features, _, audit, _ = load_sources()
    manifests = {name: json.loads((output / f'{name}_manifest.json').read_text()) for name in ['statistics', 'ml']}
    for manifest in manifests.values():
        assert manifest['features'] == features
        for name, digest in manifest['source_sha256'].items():
            assert sha(ROOT / name) == digest, f'Source changed: {name}'
        for name, digest in manifest['code_sha256'].items():
            assert sha(ROOT / name) == digest, f'Code changed after execution: {name}'
        for name, digest in manifest['output_sha256'].items():
            assert sha(output / name) == digest, f'Output changed: {name}'
        assert protected_hashes() == manifest['protected_before'], 'Existing analysis or manuscript changed'
    counts = cohort_counts(disease)
    stored_counts = read_csv(output, 'cohort_and_race_audit.csv')
    pd.testing.assert_frame_equal(counts, stored_counts, check_dtype=False)
    sampling = json.loads((output / 'fairface_sampling_summary.json').read_text())
    composition = counts.loc[counts.eligible].groupby('canonical_race').n_current_gmdb.sum().to_dict()
    pool = reference.canonical_race.value_counts().to_dict()
    size, quotas = apportion(composition, pool, sampling['maximum_X'])
    assert size == sampling['X'] and quotas == sampling['quotas']
    sampled = sample_reference(reference, quotas, sampling['seed'])
    selected_paths = read_csv(output, 'PRIVATE_reference_selection.csv')
    assert selected_paths.file.tolist() == sampled.file.tolist() and selected_paths.file.is_unique
    assert selected_paths.image_sha256.is_unique
    for race in RACES:
        assert int(sampled.canonical_race.eq(race).sum()) == quotas[race]
    a = read_csv(output, 'analysis_A_all_120.csv')
    primary = read_csv(output, 'primary_all_tests.csv')
    assert len(a) == 120 and len(primary) == counts.eligible.sum() * 120
    assert not primary[['disease', 'canonical_race', 'feature']].duplicated().any()
    for frame in [a, primary]:
        np.testing.assert_allclose(frame.q, multipletests(frame.p, method='fdr_bh')[1], rtol=1e-12, atol=1e-15)
        assert np.array_equal(frame.significant, frame.supported & frame.q.lt(.05))
    verified_u = 0
    for (syndrome, race), frame in primary.groupby(['disease', 'canonical_race']):
        cases, controls = matched_frames(disease, sampled, syndrome, race)
        assert frame.feature.tolist() == features
        assert frame.case_n.eq(len(cases)).all() and frame.reference_n.eq(len(controls)).all()
        for row in frame.itertuples():
            if row.supported:
                x, y = cases[row.feature].to_numpy(float), controls[row.feature].to_numpy(float)
                u = rankdata(np.r_[x, y])[:len(x)].sum() - len(x) * (len(x) + 1) / 2
                assert u == row.mannwhitney_u
                verified_u += 1
            else:
                assert row.p == 1 and not row.significant
    global_top, selected = select_ml(counts)
    classes = selected.disease.tolist()
    stored_ml = read_csv(output, 'ml_selected_classes.csv')
    assert stored_ml.disease.tolist() == classes
    ml_summary = read_csv(output, 'ml_cv_summary.csv')
    long_confusion = read_csv(output, 'ml_confusion_matrix.csv')
    for row in ml_summary.itertuples():
        oof = read_csv(output, f'PRIVATE_{row.model}_oof_predictions.csv')
        assert oof.patient_id.is_unique
        folds = read_csv(output, f'PRIVATE_{row.model}_outer_folds.csv')
        assert folds.patient_id.is_unique and set(folds.patient_id) == set(oof.patient_id)
        scores = None if row.model == 'dummy' else oof[[f'ranking_score_{i}' for i in range(5)]].to_numpy()
        recomputed = metrics(oof.true_class, oof.predicted_class, scores)
        for key, value in recomputed.items():
            np.testing.assert_allclose(value, getattr(row, key), rtol=1e-12)
        matrix = confusion_matrix(oof.true_class, oof.predicted_class, labels=range(5))
        stored = long_confusion.loc[long_confusion.model.eq(row.model)].pivot(
            index='true_disease', columns='predicted_disease', values='count').reindex(index=classes, columns=classes)
        np.testing.assert_array_equal(matrix, stored)
        inner = read_csv(output, f'PRIVATE_{row.model}_inner_folds.csv')
        for outer_fold, frame in inner.groupby('outer_fold'):
            holdout = set(folds.loc[folds.outer_fold.eq(outer_fold), 'patient_id'])
            assert not holdout & set(frame.patient_id)
            for _, inner_fold in frame.groupby('inner_fold'):
                assert not (set(inner_fold.loc[inner_fold.role.eq('train'), 'patient_id'])
                            & set(inner_fold.loc[inner_fold.role.eq('validation'), 'patient_id']))
    suites = ET.parse(output / 'tests.xml').getroot()
    test_summary = {field: sum(int(suite.get(field, 0)) for suite in suites.iter('testsuite'))
                    for field in ['tests', 'failures', 'errors', 'skipped']}
    assert test_summary['tests'] > 0 and test_summary['failures'] == test_summary['errors'] == 0
    regression = {}
    if (output / 'regression_tests.xml').exists():
        suites = ET.parse(output / 'regression_tests.xml').getroot()
        regression = {field: sum(int(suite.get(field, 0)) for suite in suites.iter('testsuite'))
                      for field in ['tests', 'failures', 'errors', 'skipped']}
        assert regression['failures'] == regression['errors'] == 0
    validation = {'status': 'PASS', 'primary_tests': len(primary), 'independent_U_checks': verified_u,
                  'primary_BH_reconstructed': True, 'analysis_A_BH_reconstructed': True,
                  'sampling_replayed': True, 'all_ML_metrics_and_confusion_reconstructed': True,
                  'protected_files_unchanged': True, 'tests': test_summary, 'existing_regression_tests': regression,
                  'source_and_execution_code_hashes_match': True}
    write_json(output / 'validation_summary.json', validation)
    return a, primary, counts, audit, sampling, ml_summary, global_top, selected, validation, manifests


def report(output, result):
    a, primary, counts, audit, sampling, ml_summary, global_top, selected, validation, manifests = result
    significant = primary.loc[primary.significant]
    cells = read_csv(output, 'syndrome_race_summary.csv')
    regions = read_csv(output, 'anatomical_region_summary.csv')
    regional = regions.groupby('anatomical_region')[['planned', 'tested', 'significant']].sum().reset_index()
    regional['rate'] = regional.significant / regional.tested
    regional = regional.sort_values(['rate', 'anatomical_region'], ascending=[False, True])
    selection_columns = ['disease', 'canonical_race', 'n_current_gmdb']
    top = significant.sort_values(['q', 'disease', 'canonical_race', 'feature']).head(12)
    statuses = [
        ('Source/QC, content duplication, counts and label audits', 'COMPLETED'),
        ('Three-category FairFace Analysis A (120-test BH)', 'COMPLETED'),
        ('Proportional reference and within-category primary family', 'COMPLETED'),
        ('Global top-five and shared-category five-class selection', 'COMPLETED'),
        ('Four-model nested-CV ML and aggregate figures', 'COMPLETED'),
        ('Reproducibility, privacy guards and synthetic fixture tests', 'COMPLETED'),
        ('100-seed reference replay and all-reference inferential sensitivity', 'NOT RUN: deferred by user'),
        ('Narrower Asian definition and tiny-cell permutation sensitivity', 'NOT RUN: deferred by user'),
        ('Bootstrap intervals, repeated CV and label permutations', 'NOT RUN: deferred by user'),
        ('Covariate-adjusted inferential alternatives', 'NOT RUN: deferred by user'),
        ('Clinical/HPO, photographs, literature novelty and manuscript editing', 'NOT RUN: outside current scope'),
        ('Aggregate publication or active-study promotion', 'NOT RUN: private experimental review'),
        ('Missing source/package/permission blockers', 'NOT APPLICABLE: none for completed components'),
    ]
    lines = ['# PRISM race-stratified redesign: private results', '',
             '## Scope and completion', '',
             'This experimental analysis asks which of the exact 120 pinned FaceKit geometry features differ between '
             'GMDB syndromes and same-recorded-category FairFace references. Cases are the current quality-passing, '
             'one-image-per-recorded-patient cohort with conditional age 3 to <10 and complete finite measurements. '
             'Each syndrome/category cell needs five distinct patients. Statistical eligibility is rebuilt across all '
             'syndromes and does not inherit the previous 48-syndrome restriction. The classifier separately predicts '
             'five syndromes within one recorded category using patient-held-out nested CV.', '',
             md_table(pd.DataFrame(statuses, columns=['Component', 'Status'])), '',
             'The direct user instruction to defer sensitivity analysis supersedes those tasks in the supplied document. '
             'Existing main, demographic, sensitivity, publication and manuscript files were hash-checked and preserved. '
             'The redesign has not been promoted into ACTIVE_ANALYSIS.json. No raw data or aggregates were published.', '',
             '## Source and cohort audit', '',
             f"The validated GMDB source has {audit['gmdb_selected_patients']:,} patients. Recorded source labels are "
             f"{json.dumps(audit['gmdb_source_labels'], sort_keys=True)}. The three-category analysis excludes "
             f"{audit['gmdb_unmapped_excluded']} Others/unknown records. All selected GMDB patient IDs, image IDs and "
             'image-content hashes are unique; source selection, exact measurements, age range, model and extractor pins pass.', '',
             f"The original accepted FairFace pool has {audit['fairface_accepted_paths_before_content_qc']:,} age-3–9 paths. "
             f"Two conflicting-label content-hash sets remove {audit['fairface_conflicting_label_rows_excluded']} rows; "
             f"{audit['fairface_repeated_content_rows_removed']} additional byte-identical copies are removed by lexical-first path. "
             f"This leaves {audit['fairface_after_content_qc']:,} observations, of which "
             f"{audit['fairface_unmapped_excluded']} unsupported-category observations are excluded. This content-QC rule "
             'was set before p-values. Source-label counts after content QC: '
             f"{json.dumps(audit['fairface_source_labels_after_content_qc'], sort_keys=True)}.", '',
             'The fixed crosswalk maps GMDB African/Asian/European to the same canonical categories, and FairFace '
             'Black→African, White→European, East Asian/Southeast Asian/Indian/Middle Eastern→Asian. '
             'The documented West/Central Asian spelling variants remain supported by code. Including Middle Eastern/West '
             'Asian within Asian is a substantive and heterogeneous harmonization decision. No equivalent granular GMDB '
             'ancestry is manufactured. The narrower-map sensitivity is deferred.', '',
             f"There are {counts.disease.nunique()} candidate syndrome identifiers, {counts.eligible.sum()} eligible cells, "
             f"{counts.loc[counts.eligible].disease.nunique()} eligible syndromes and "
             f"{counts.loc[counts.eligible].n_current_gmdb.sum()} eligible patients. Mapped patients outside eligible cells: "
             f"{counts.loc[~counts.eligible].n_current_gmdb.sum()}. The count audit contains all three category cells for "
             'every candidate syndrome, including zeroes and <5 cells, pass/fail reasons, counts, source labels and '
             'old >=10-total eligibility. Mucopolysaccharidoses is flagged as an umbrella diagnosis.', '',
             'The previous study used 1,845 patients as its source cohort, 48 >=10-total syndromes, 5,760 hypotheses '
             'and 955 associations. The redesigned eligibility, comparisons and correction family differ; discovery '
             'counts do not establish replication, improved accuracy or survival of the old q-values.', '',
             '## FairFace demographic Analysis A', '',
             f"Three-category Kruskal–Wallis tests on all content-screened mapped FairFace images yield "
             f"{int(a.significant.sum())}/120 features with BH q<.05; {int((~a.significant).sum())} are not significant. "
             'The family contains exactly 120 tests. No pairwise post-hoc tests were run. Counts, means, SDs, medians, '
             'IQRs, missingness, H, p, q, units and formulas are in analysis_A_all_120.csv.', '',
             md_table(a.sort_values(['q', 'feature']).head(12)[['feature', 'anatomical_region', 'H', 'p', 'q', 'units']]), '',
             'Nonsignificant features: ' + (', '.join(a.loc[~a.significant, 'feature']) or 'none') + '.', '',
             'These are recorded-category associations and motivate examining matched references; they do not prove '
             'a causal race effect or removal of source/acquisition confounding.', '',
             '## Reference construction and new main results', '',
             f"Composition uses only the {sum(manifests['statistics']['eligibility']['composition'].values())} patients in "
             'eligible cells. The separate full-1,845 context counts are '
             f"{json.dumps(sampling['all_1845_context_composition'], sort_keys=True)}. Largest-remainder apportionment "
             f"with alphabetical ties and exact capacity checking gives maximum feasible X={sampling['X']} "
             f"under the prespecified cap {sampling['maximum_X']}. The European reference capacity limits X; "
             f"seed={sampling['seed']}. Paths are sampled without replacement in sorted order and the one global "
             'sample is reused through same-category slices.', '']
    sampling_rows = [{'category': race, 'pool': sampling['pool_sizes'][race], 'selected': sampling['quotas'][race],
                      'target_proportion': sampling['target_proportions'][race],
                      'achieved_proportion': sampling['achieved_proportions'][race],
                      'difference': sampling['proportion_discrepancy'][race]} for race in RACES]
    lines += [md_table(pd.DataFrame(sampling_rows)), '',
              f"The global primary BH family has {len(primary):,} planned comparisons; "
              f"{int(primary.supported.sum()):,} are supported. There are {len(significant):,} q<.05 associations, "
              f"including {int(significant.large_effect.sum()):,} with |g|>=.8. Significant effects include "
              f"{int(significant.hedges_g.gt(0).sum()):,} positive and {int(significant.hedges_g.lt(0).sum()):,} negative g; "
              f"{int((~np.isfinite(primary.hedges_g)).sum())} planned effects have undefined g. "
              f"{int(cells.significant_features.gt(0).sum())}/{len(cells)} cells have at least one discovery. "
              'Every p-value was calculated from the current within-category observations. Cases-minus-reference g '
              'uses pooled sample SD and exact gamma J(df). Two-sided Mann–Whitney compares distributions; it '
              'does not generally test median equality. Exact/asymptotic method and ties are recorded per feature.', '',
              md_table(top[['disease', 'canonical_race', 'feature', 'case_n', 'reference_n', 'hedges_g', 'p', 'q']]), '',
              'Feature formulas and units for these associations:', '',
              md_table(top[['feature', 'units', 'formula']].drop_duplicates()), '',
              'All significant results are retained in primary_significant_aggregate.csv and every planned result '
              'in primary_all_tests.csv. Both files are private, including their small-cell aggregates.', '',
              'Region patterns are descriptive rates with their testing denominators:', '', md_table(regional), '',
              '## Sparse-cell support and deferred sensitivity work', '',
              'No new inferential sensitivity family was run. Reference resampling, all-reference comparisons, '
              'narrower Asian mapping and tiny-cell permutations remain NOT RUN at the user’s request. '
              'Effect stability is unassessed. The following sample-size counts describe the primary cohort only.', '']
    support = cells.assign(patient_band=pd.cut(cells.n_patients, [4, 5, 9, 100000], labels=['n=5', 'n=6–9', 'n>=10']))
    support = support.groupby('patient_band', observed=False).agg(cells=('n_patients', 'size'),
        patients=('n_patients', 'sum'), discoveries=('significant_features', 'sum')).reset_index()
    lines += [md_table(support), '',
              'Covariate availability and descriptive age/sex/pose summaries are in covariate_availability_descriptive.csv. '
              'GMDB exact ages retain the existing conditional age interpretation. FairFace has category 3–9 only, '
              'so exact-age differences cannot be assessed. Both sources have recorded binary sex and pose fields; '
              'availability does not establish comparability. Source remains perfectly associated with case/reference.', '',
              '## Five-syndrome ML', '', 'Global five largest eligible cells:', '',
              md_table(global_top[selection_columns]), '',
              'Selected five distinct syndromes in one shared category, frozen from patient counts before fitting:', '',
              md_table(selected[selection_columns]), '',
              f'The global top five mix categories. The primary classifier uses the {selected.canonical_race.iloc[0]}-annotated group and '
              'has five adequately supported classes. No healthy FairFace examples enter training. Inputs are '
              'exactly 120 numeric features; metadata, patient IDs, paths, demographics, pose, source and '
              'statistical significance are excluded. Five outer and three inner stratified folds separate tuning '
              'from held-out prediction; each patient gets exactly one OOF prediction. All imputation/scaling '
              'is training-fold-local. Full grids, selected hyperparameters and fold sizes are recorded.', '',
              md_table(ml_summary[['model', 'n', 'balanced_accuracy', 'macro_f1', 'accuracy', 'top_2_accuracy',
                                   'top_3_accuracy', 'fold_balanced_accuracy_sd_descriptive']]), '',
              'Balanced chance is 20%; the observed majority dummy is shown. Missing dummy top-k values mean '
              'no ranking metric was reported. Fold SD/range describe variation, not independent uncertainty. '
              'These are internal CV estimates; choosing a winner from these scores would require further evaluation.', '',
              'Largest directed off-diagonal confusions for each fitted model:', '']
    confusion = read_csv(output, 'ml_confusion_matrix.csv')
    wrong = confusion.loc[confusion.true_disease.ne(confusion.predicted_disease) & confusion.model.ne('dummy')]
    wrong = wrong.sort_values(['model', 'count', 'true_disease', 'predicted_disease'], ascending=[True, False, True, True])
    lines += [md_table(wrong.groupby('model', sort=False).head(3)), '',
              'The full 5×5 count matrix for every model is in ml_confusion_matrix.csv, with figures in '
              'ml_confusion_matrix.png and ml_model_comparison.png. Confusion counts do not prove shared phenotype. '
              'Per-class precision, recall, F1 and support are in ml_per_class_summary.csv. Standardized logistic '
              'coefficients are descriptive across five fitted folds; correlated normalized features prevent '
              'causal anatomical interpretation. No calibrated clinical probabilities or independent validation are claimed.', '']
    per_class = read_csv(output, 'ml_per_class_summary.csv')
    zero_recall = per_class.loc[per_class.recall.eq(0) & per_class.model.ne('dummy'), ['model', 'disease', 'support']]
    lines += ['Logistic per-class performance:', '',
              md_table(per_class.loc[per_class.model.eq('logistic'), ['disease', 'precision', 'recall', 'support']]), '',
              'Classes with zero held-out recall in fitted models:', '', md_table(zero_recall), '',
              '## Interpretation and limits', '',
              'Within-category comparisons control the measured mismatch in recorded category between each syndrome '
              'and its reference. They do not establish that results are more accurate: source acquisition remains '
              'confounded, group correspondence is broad and uncertain, ages/sex may differ, and many cells are small. '
              'No source supplies verified genetic ancestry. FairFace is not clinically certified healthy. '
              'GMDB patient independence is based on recorded identifiers and content hashes, without a comprehensive '
              'identity-linkage review. FairFace distinct paths/content do not verify distinct people. '
              'Selection, referral, subtype heterogeneity and feature normalization can affect associations and ML. '
              'Clinical/HPO comparisons and novel phenotype claims are deferred.', '',
              '## Reproducibility and verification', '',
              f"Validation: {json.dumps(validation, sort_keys=True)}.", '',
              f"Git revision: `{manifests['statistics']['revision']}`. Exact uncommitted source hashes, 120-feature order, "
              'input/output SHA-256 values, model/extractor pins and parameters are in statistics_manifest.json and '
              'ml_manifest.json. Protected-file hashes establish that existing outputs and sensitivity files were unchanged. '
              'validation_summary.json records hash matching, reconstructed BH, every supported U statistic, '
              'sampling replay and reconstructed ML scores/confusion. Tests use fabricated fixtures only.', '',
              md_table(pd.DataFrame(manifests['statistics']['versions'].items(), columns=['Dependency', 'Version'])), '',
              'Successful commands from repository root:', '', '```sh',
              '.venv/bin/python revision_results/race_stratified_20261115/run_race_stratified.py --dry-run',
              '.venv/bin/python revision_results/race_stratified_20261115/run_five_syndrome_ml.py --dry-run',
              '.venv/bin/python revision_results/race_stratified_20261115/run_race_stratified.py',
              '.venv/bin/python revision_results/race_stratified_20261115/run_five_syndrome_ml.py',
              '.venv/bin/python -m pytest revision_results/race_stratified_20261115/tests -q --junitxml=local_review/race_stratified_20261115/tests.xml',
              '.venv/bin/python revision_results/race_stratified_20261115/verify_and_report.py', '```', '',
              f'Private output directory: `{output}`. Reviewable implementation: '
              '`revision_results/race_stratified_20261115/`. No aggregate release has been approved or performed.', '',
              'Method references: [SciPy Mann–Whitney](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html) '
              'and [scikit-learn nested CV](https://scikit-learn.org/stable/auto_examples/model_selection/plot_nested_cross_validation_iris.html).']
    (output / 'RESULTS_REPORT.md').write_text('\n'.join(lines) + '\n')
    files = {p.name: sha(p) for p in sorted(output.iterdir()) if p.is_file() and p.name != 'final_artifact_hashes.json'}
    write_json(output / 'final_artifact_hashes.json', files)


def main():
    parser = argparse.ArgumentParser(description='Verify executed private results and write the research report')
    parser.add_argument('--output-dir', default=str(PRIVATE))
    args = parser.parse_args()
    output = private_output(args.output_dir)
    result = verify(output)
    report(output, result)
    print(json.dumps(result[8], indent=2))


if __name__ == '__main__':
    main()

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[2]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def close(observed, expected):
    np.testing.assert_allclose(observed, expected, rtol=1e-10, atol=1e-12, equal_nan=True)


def verify(root=ROOT):
    folder, synthetic = root / 'publication/final_research_120', root / 'publication/synthetic_consistency'
    source = root / 'publication/supplementary'
    schema = pd.read_csv(source / 'Supplementary_Table_S6_exact_feature_definitions.csv')
    primary = pd.read_csv(source / 'S1_primary_all_5760_hypotheses.csv')
    manifest = json.loads((synthetic / 'synthetic_consistency_manifest.json').read_text())
    assert manifest['feature_count'] == 120 and manifest['reference_n'] == 3409
    assert manifest['primary_unchanged'] == {'syndromes': 48, 'hypotheses': 5760, 'discoveries': 955, 'synthetic_rows': 0}
    assert manifest['synthetic_total_n'] == 1727 and manifest['synthetic_qc_n'] == 1563
    assert manifest['supported_syndromes'] == 9 and manifest['supported_real_n'] == 314 and manifest['supported_synthetic_n'] == 1424
    assert manifest['code_sha256'] == sha(root / 'stats_testing/synthetic_consistency.py')
    assert manifest['mapping_sha256'] == sha(synthetic / 'reviewed_label_mapping.csv')
    assert manifest['plan_sha256'] == sha(synthetic / 'predeclared_analysis_plan.md')
    assert manifest['report_generation_code_sha256'] == sha(root / 'publication/scripts/final_research_text.py')
    assert manifest['figure_source_manifest_sha256'] == sha(folder / 'FIGURE_SOURCE_MANIFEST.json')
    config = json.loads((root / 'ACTIVE_ANALYSIS.json').read_text())
    assert manifest['config_sha256'] == sha(root / 'ACTIVE_ANALYSIS.json')
    assert manifest['facekit_commit'] == config['facekit_commit'] and manifest['model_sha256'] == config['model_sha256']
    for name, digest in manifest['outputs'].items():
        assert sha(synthetic / name) == digest, name
    frame = pd.read_csv(synthetic / 'feature_fidelity_summary.csv')
    assert not frame.duplicated(['syndrome', 'feature', 'model_tag', 'sensitivity']).any()
    assert set(frame.feature) == set(schema.feature) and frame.groupby(['syndrome', 'model_tag', 'sensitivity']).size().eq(120).all()
    assert frame.real_n.ge(10).all() and frame.synthetic_n.ge(10).all() and frame.reference_n.eq(3409).all()
    assert frame[['real_missing', 'synthetic_missing', 'reference_missing']].eq(0).all().all()
    close(frame.mean_shift, frame.synthetic_mean - frame.real_mean)
    close(frame.median_shift, frame.synthetic_median - frame.real_median)
    close(frame.standardized_mean_shift, frame.mean_shift / frame.real_sd)
    close(frame.variance_ratio, (frame.synthetic_sd / frame.real_sd)**2)
    close(frame.wasserstein_real_sd, frame.wasserstein / frame.real_sd)
    for side in ['real', 'synthetic']:
        n, m = frame[f'{side}_n'], frame.reference_n
        sd = np.sqrt(((n - 1) * frame[f'{side}_sd']**2 + (m - 1) * frame.reference_sd**2) / (n + m - 2))
        close(frame[f'{side}_g'], (frame[f'{side}_mean'] - frame.reference_mean) / sd * (1 - 3 / (4 * (n + m) - 9)))
        close(frame[f'{side}_reference_sd_contrast'], (frame[f'{side}_mean'] - frame.reference_mean) / frame.reference_sd)
    main = frame.loc[frame.model_tag.eq('pooled') & frame.sensitivity.eq('unmodified')]
    joined = main.merge(primary, left_on=['syndrome', 'feature'], right_on=['disease', 'feature'], validate='one_to_one')
    assert len(joined) == 1080
    close(joined.real_g, joined.hedges_g)
    close(joined.real_mean, joined.mean_disease)
    close(joined.reference_mean, joined.mean_reference)
    fidelity = pd.read_csv(synthetic / 'syndrome_fidelity_summary.csv')
    for key, part in frame.groupby(['syndrome', 'model_tag', 'sensitivity']):
        row = fidelity.set_index(['syndrome', 'model_tag', 'sensitivity']).loc[key]
        close(row.median_absolute_standardized_shift, np.median(np.abs(part.standardized_mean_shift)))
        close(row.median_variance_ratio, part.variance_ratio.median())
        close(row.direction_agreement, np.mean(np.sign(part.real_g) == np.sign(part.synthetic_g)))
        close(row.pearson_g, pearsonr(part.real_g, part.synthetic_g).statistic)
        close(row.spearman_g, spearmanr(part.real_g, part.synthetic_g).statistic)
        assert row.sign_reversals == (part.real_g * part.synthetic_g < 0).sum()
    for label, part in frame.loc[frame.sensitivity.eq('unmodified')].groupby('syndrome'):
        for field in ['real_mean', 'real_sd', 'reference_mean', 'reference_sd']:
            assert part.groupby('feature')[field].nunique().eq(1).all(), (label, field)
    audit = pd.read_csv(synthetic / 'cohort_and_mapping_audit.csv')
    assert audit.descriptive_eligible.sum() == 9 and audit.mapping_status.eq('excluded_no_synthetic_label').sum() == 39
    assert audit.loc[audit.pdidb_disease.notna(), 'synthetic_qc_n'].sum() == 1563
    assert audit.loc[audit.descriptive_eligible, 'real_n'].sum() == 314
    anchors = pd.read_csv(synthetic / 'anchor_distribution_source.csv')
    assert len(anchors) == 5
    match = anchors.merge(main, on=['syndrome', 'feature'], suffixes=('_chart', '_source'), validate='one_to_one')
    for field in ['real_q1', 'real_median', 'real_q3', 'synthetic_q1', 'synthetic_median', 'synthetic_q3']:
        close(match[field + '_chart'], match[field + '_source'])
    regions = pd.read_csv(folder / 'regional_discovery_source.csv')
    merged = primary.merge(schema[['feature', 'anatomical_region']], on='feature')
    expected = merged.assign(significant=merged.mannwhitney_q_fdr_global.lt(.05)).groupby('anatomical_region').agg(
        tested=('feature', 'size'), significant=('significant', 'sum'), features=('feature', 'nunique')).reset_index()
    pd.testing.assert_frame_equal(regions[expected.columns], expected, check_dtype=False)
    top = pd.read_csv(folder / 'top_five_recurrent_features.csv')
    counts = primary.loc[primary.mannwhitney_q_fdr_global.lt(.05)].groupby('feature').size().rename('count').reset_index()
    counts = counts.sort_values(['count', 'feature'], ascending=[False, True]).head(5)
    assert counts.feature.tolist() == top.feature.tolist()
    close(counts['count'], top.significant_syndromes)
    close(top.positive + top.negative, top.significant_syndromes)
    table_manifest = json.loads((folder / 'TABLE_SOURCE_MANIFEST.json').read_text())
    assert table_manifest['code_sha256'] == sha(root / 'publication/scripts/final_research_tables.py')
    for name, digest in table_manifest['source_sha256'].items():
        assert sha(source / name) == digest
    demos = pd.read_csv(folder / 'demographic_examples.csv')
    assert len(demos) == 18 and demos.groupby(['layer', 'analysis']).size().eq(3).all()
    for row in demos.itertuples():
        original = pd.read_csv(source / row.source)
        original = original.loc[original.feature.eq(row.feature)]
        if row.analysis == 'A':
            assert len(original) == 1
            close(row.q, original.iloc[0].fdr_q_value)
            close(row.omnibus_H, original.iloc[0].kruskal_H)
        elif row.analysis == 'B':
            original = original.loc[original.syndrome_name.eq(row.syndrome) & original.group.eq(row.group_1)]
            assert len(original) == 1
            expected = original.iloc[0]
            close([row.effect, row.q, row.ci_lower, row.ci_upper, row.n_disease_1, row.n_reference_1],
                  [expected.hedges_g, expected.global_fdr_q_value, expected.hedges_g_ci_lower, expected.hedges_g_ci_upper,
                   expected.n_disease, expected.n_reference_matched])
        else:
            original = original.loc[original.syndrome_name.eq(row.syndrome) & original.group_1.eq(row.group_1) & original.group_2.eq(row.group_2)]
            assert len(original) == 1
            expected = original.iloc[0]
            close([row.effect, row.q, row.ci_lower, row.ci_upper], [expected.interaction_beta_group2_minus_group1,
                  expected.global_fdr_q_value, expected.interaction_ci_lower, expected.interaction_ci_upper])
        assert row.family_hypotheses == len(pd.read_csv(source / row.source))
    adjusted = pd.read_csv(folder / 'adjusted_examples.csv')
    assert len(adjusted) == 12 and adjusted.family_hypotheses.value_counts().to_dict() == {5880: 5, 4: 4, 3: 3}
    for row in adjusted.itertuples():
        original = pd.read_csv(source / row.source)
        selected = original.loc[original.syndrome.eq(row.syndrome) & original.feature.eq(row.feature)]
        assert len(selected) == 1
        expected = selected.iloc[0]
        if row.family_hypotheses == 3:
            close([row.beta, row.q, row.n_complete], [expected.adjusted_disease_coefficient,
                  expected.bh_q_across_three_predefined_models, expected.n_total])
        else:
            close([row.beta, row.q, row.n_complete], [expected.target_beta,
                  expected.q_four_models if row.family_hypotheses == 4 else expected.q_global_bh, expected.n_complete])
        close([row.ci_lower, row.ci_upper], [expected.ci_95_lower, expected.ci_95_upper])
    rts = adjusted.loc[adjusted.family_hypotheses.eq(4) & adjusted.feature.eq('malar_bulge_mean')].iloc[0]
    assert rts.q > .05 and rts.ci_lower < 0 < rts.ci_upper
    literature = pd.read_csv(folder / 'clinical_literature_audit_17.csv')
    assert len(literature) == 17 and not literature.exact_ratio_validated.any()
    assert literature.retain_supported_clinical_example.sum() == 4
    figures = json.loads((folder / 'FIGURE_SOURCE_MANIFEST.json').read_text())
    assert figures['script_sha256'] == sha(root / 'publication/scripts/final_research_figures.py')
    assert len(figures['figures']) == 6
    for item in figures['figures']:
        assert item['caption'] and item['alt_text']
        for path, digest in item['sources'].items():
            assert sha(root / path) == digest, path
        for ext, digest in item['outputs'].items():
            assert sha(root / f'manuscript/figures/{item["name"]}.{ext}') == digest
    for directory in [folder, synthetic]:
        for path in directory.glob('*.csv'):
            columns = set(pd.read_csv(path, nrows=0).columns)
            assert not columns & {'patient_id', 'image_id', 'image_sha256', 'image_path', 'landmarks', 'embedding', 'file'}
    text = (root / 'manuscript/FINAL_120_RESULTS_INSERTIONS.md').read_text()
    assert '0.2066604' in text and '48' in text and '5,760' in text and '955' in text
    return {'status': 'PASS', 'pooled_synthetic_feature_cells': 1080, 'demographic_examples': 18,
            'adjusted_examples': 12, 'literature_entries': 17, 'new_figures': 6,
            'scope': 'Aggregate arithmetic, fixed-reference agreement, source correspondence and figure hashes; no image reprocessing'}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))

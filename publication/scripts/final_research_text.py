import json
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / 'publication/final_research_120'
SYNTHETIC = ROOT / 'publication/synthetic_consistency'


def fmt(value):
    if isinstance(value, (float, np.floating)):
        return f'{value:.6g}' if np.isfinite(value) else '—'
    return str(value).replace('|', '/')


def table(frame):
    return '| ' + ' | '.join(frame.columns) + ' |\n| ' + ' | '.join(['---'] * len(frame.columns)) + ' |\n' + ''.join(
        '| ' + ' | '.join(fmt(value) for value in row) + ' |\n' for row in frame.itertuples(index=False, name=None))


def short(value):
    return value.split('__')[0]


def synthetic_report():
    fidelity = pd.read_csv(SYNTHETIC / 'syndrome_fidelity_summary.csv')
    principal = fidelity.loc[fidelity.model_tag.eq('pooled') & fidelity.sensitivity.eq('unmodified')]
    detail = pd.read_csv(SYNTHETIC / 'feature_fidelity_summary.csv')
    detail = detail.loc[detail.model_tag.eq('pooled') & detail.sensitivity.eq('unmodified')]
    manifest = json.loads((SYNTHETIC / 'synthetic_consistency_manifest.json').read_text())
    anchors = pd.read_csv(SYNTHETIC / 'anchor_distribution_source.csv')
    medians = principal.select_dtypes('number').median()
    text = '# Conditional PDIDB–GMDB consistency results\n\n'
    text += ('This secondary analysis uses frozen current measurements without reprocessing images or refitting the primary study. '
             'The mentor reports that StyleGAN was trained on GMDB and PDIDB sampled that learned distribution. This statement is recorded separately from provider metadata: '
             'the identities, training membership and relationship of v1/v2 checkpoints to the published generator remain unverified. '
             'Synthetic observations do not increase the real-patient n or independently confirm the primary 955 associations.\n\n')
    text += (f"## Coverage and methods\n\nThe actual files contain 1,727 synthetic records and 1,563 QC-passing rows. Nine exact-name matches support comparisons of "
             f"{manifest['supported_real_n']} selected GMDB patients and {manifest['supported_synthetic_n']} synthetic images on 120 features each (1,080 pooled syndrome-feature cells). "
             'The 139 QC-passing synthetic 22q11.2 rows remain excluded because the generic label cannot establish equivalence to DiGeorge or distal deletion labels. '
             'All accepted comparisons have at least ten observations per source; available selected real n ranges from 16 to 67. '
             'Age_range, gender and race are missing for all 1,727 synthetic rows. No age-matched or demographic-adjusted synthetic sensitivity is possible, '
             'while real GMDB ages are 3 <= age < 10. The fixed reference is the same 3,409 FairFace paths in every unmodified profile.\n\n'
             'The analysis plan was written after metadata inspection and before fidelity computation; it was not registered externally. '
             'It defines sample SD (ddof=1), synthetic-minus-real shifts, real-SD scaling, empirical Wasserstein distance, sample variance ratios and Hedges’ g against the fixed reference. '
             'The second profile scaling divides both mean contrasts by the identical reference SD. Winsorized sensitivity clips each source separately at its 1st/99th percentiles, including a once-fixed clipped reference. '
             'No p-values, inferential intervals, equivalence margins or generator-level uncertainty are computed.\n\n')
    text += ('## Main descriptive findings\n\n'
             f"Across the nine equally weighted syndrome summaries, median absolute standardized mean shift is {medians.median_absolute_standardized_shift:.3f} real SD, "
             f"median variance ratio is {medians.median_variance_ratio:.3f}, median effect-direction agreement is {medians.direction_agreement:.1%}, "
             f"and median Pearson/Spearman profile correlation is {medians.pearson_g:.3f}/{medians.spearman_g:.3f}. "
             f"Across all 1,080 cells (equal 120-feature contribution per syndrome), {(detail.variance_ratio < 1).sum()} have lower synthetic sample variance. "
             'These medians summarize correlated measurements, not independent biological replications.\n\n')
    shown = principal[['syndrome', 'real_n', 'synthetic_n', 'median_absolute_standardized_shift', 'median_variance_ratio',
                       'fraction_synthetic_variance_lower', 'direction_agreement', 'pearson_g', 'spearman_g']].copy()
    shown['syndrome'] = shown.syndrome.map(short)
    text += table(shown) + '\n'
    text += ('Cornelia de Lange has a median variance ratio of 0.839 and lower synthetic variance in 78.3% of features; KBG has a median ratio of 0.785. '
             'Angelman, Nicolaides–Baraitser and Rubinstein–Taybi have syndrome-level median ratios above one. Thus compression is heterogeneous, '
             'and a pooled near-one dispersion summary does not establish fidelity for each syndrome or measurement.\n\n')
    text += '## Five anchors, retained regardless of agreement\n\n'
    anchor_table = anchors[['syndrome', 'feature', 'real_g', 'synthetic_g', 'standardized_mean_shift', 'variance_ratio']].copy()
    anchor_table['syndrome'] = anchor_table.syndrome.map(short)
    text += table(anchor_table) + '\n'
    text += ('Philtrum and pupil-spacing effects are close in sign and magnitude. The Williams–Beuren eyebrow-area and outer-canthal effects and the Rubinstein–Taybi malar effect '
             'are attenuated in these synthetic samples, although all five retain the primary effect sign. The malar index is a projected cheek/outline construct and the eyebrow index '
             'is polygon area; neither is a direct diagnosis of a clinical trait. No significance criterion was used to choose these anchors.\n\n')
    gaps = pd.read_csv(SYNTHETIC / 'largest_gaps_exploratory.csv').head(5)
    gap_table = gaps[['syndrome', 'feature', 'standardized_mean_shift', 'variance_ratio', 'real_n', 'synthetic_n']].copy()
    gap_table['syndrome'] = gap_table.syndrome.map(short)
    text += '## Largest discrepancies (exploratory ranking)\n\n' + table(gap_table) + '\n'
    text += '## Sensitivities\n\n'
    text += table(pd.read_csv(SYNTHETIC / 'low_real_n_sensitivity.csv')) + '\n'
    sensitivity = fidelity.groupby(['model_tag', 'sensitivity'], sort=True).agg(
        syndromes=('syndrome', 'nunique'), median_shift=('median_absolute_standardized_shift', 'median'),
        median_variance_ratio=('median_variance_ratio', 'median'), median_direction_agreement=('direction_agreement', 'median'),
        median_pearson_g=('pearson_g', 'median')).reset_index()
    text += table(sensitivity) + '\n'
    text += ('Model-tag strata describe version-labeled release subsets, not independently trained generators. Every supported tag result is supplied in the feature and syndrome tables. '
             'Fixed-reference-SD profile correlations and contrasts are supplied alongside pooled-SD Hedges’ g. The frozen feature-correlation clusters have separate feature-within-cluster '
             'and equal-cluster-weighted summaries; they do not modify inferential denominators. The n >= 20 and n >= 30 real-cohort cutoffs were added after initial pooled results and are explicitly exploratory. '
             'Low-n exclusions and winsorization change descriptive summaries but do not resolve age/source mismatch.\n\n')
    text += table(pd.read_csv(SYNTHETIC / 'available_size_proxy_summary.csv')) + '\n'
    text += ('The size variable here is the available age-restricted GMDB selection, not the actual number of StyleGAN training patients/images. '
             'These nine descriptive correlations cannot establish a training-size effect. No exact duplicate encoded-byte groups or exact synthetic/selected-real byte overlaps were found; '
             'that check cannot establish patient independence, absence of memorization or privacy.\n\n'
             '## Relationship to FaceKit and limitations\n\n'
             '[FaceKit §§2.6/3.6 and Fig. 9](https://arxiv.org/html/2610.09288v1#S3.SS6) already assess synthetic phenotype preservation. '
             'Its 4,559 real/synthetic images and HPO-mapped contrasts differ from this age-restricted selection and fixed FairFace baseline. Its reported median Pearson correlation (0.866) '
             'and variance ratio (0.94) are not directly comparable benchmarks. Exact generator-release equivalence is unverified. '
             'This audit extends a question already studied in FaceKit; it does not establish a new extraction method, independent replication or general generator validity.\n\n'
             'All images, per-image measurements, hashes and sensitive linkage records remain local. The output directory contains only syndrome/feature statistics, source-level category counts, '
             'code/documentation and original aggregate plots. Existing PI publication review still applies; this pass does not perform privacy/IRB/legal release review.\n\n'
             '## Figures and reproducibility\n\n'
             'Figures 8–10 in manuscript/figures show only predeclared-anchor quartiles, aggregate effect profiles and syndrome-level fidelity. '
             'Captions and alt text are in publication/final_research_120/FIGURE_CAPTIONS.md. CSVs reproduce the plotted values exactly. '
             'Run stats_testing/synthetic_consistency.py with --active-config ACTIVE_ANALYSIS.json and a private --output-dir, followed by publication/scripts/final_research_figures.py '
             'and publication/scripts/final_research_text.py. Source hashes and the mentor/provider provenance distinction are in synthetic_consistency_manifest.json.\n')
    (SYNTHETIC / 'synthetic_consistency_results.md').write_text(text)


def insertions():
    top = pd.read_csv(TABLES / 'top_five_recurrent_features.csv')
    regional = pd.read_csv(TABLES / 'regional_discovery_source.csv')
    demos = pd.read_csv(TABLES / 'demographic_examples.csv')
    groups = pd.read_csv(TABLES / 'demographic_A_group_summaries.csv')
    adjusted = pd.read_csv(TABLES / 'adjusted_examples.csv')
    literature = pd.read_csv(TABLES / 'clinical_literature_audit_17.csv')
    text = ('# Final 120-feature results insertions\n\n'
            'Prepared November 9, 2026 for the Formal Paper tab of FaceKit Expansion Results. The Google Doc was inspected read-only; no tab or comment was edited. '
            'These replacements use the frozen current aggregate package. Exact machine-readable values and formulas are in publication/final_research_120; '
            'displayed numbers below are rounded to six significant digits. No primary models were rerun or replaced.\n\n'
            '## Primary recurrence: replace the top-five-feature placeholder\n\n'
            'Across the 48 eligible syndromes, the most recurrent globally significant measurements were bilateral mean lateral eyebrow thickness (32 syndromes), '
            'forehead height (29), chin height (27), left lateral eyebrow thickness (27) and right lateral eyebrow thickness (26). '
            'Eyebrow and forehead effects were positive; all 27 significant chin-height effects were negative. Left/right eyebrow indices share anatomy, '
            'so their recurrence is not evidence for independent phenotypes. Significance means q < 0.05 in the single 5,760-test global BH family, and all five measurements were tested in all 48 cohorts.\n\n')
    text += table(top[['feature', 'feature_name', 'anatomical_region', 'significant_syndromes', 'tested_syndromes', 'positive', 'negative']]) + '\n'
    text += ('Table caption: Five measurements ranked by number of primary significant syndrome cohorts, with lexical feature-ID tie breaking. '
             'Counts use global BH q < .05, not within-syndrome q. Forehead/chin height and lateral eyebrow thickness use the precise S6 normalized geometry; '
             'neither pixel distances nor clinical hair density are implied. Exact formulas are retained in the source CSV.\n\n'
             '## Regional results: replace region-table/figure placeholders\n\n'
             'Discovery rates varied across anatomical regions. Chin/jaw contributed 87/240 significant hypotheses (36.25%) and forehead 60/192 (31.25%), '
             'compared with 197/960 for eyebrows (20.52%) and 35/384 for cheeks/malar geometry (9.11%). '
             'Eyebrows generated more absolute discoveries while chin/jaw had a higher rate. Different regional feature counts produce unequal denominators, '
             'and shared landmarks make these correlated summaries rather than independent anatomical enrichment tests.\n\n')
    text += table(regional) + '\n'
    text += ('Table caption: Significant and tested primary syndrome-feature hypotheses in all nine regions. Regional denominators are 48 times the region’s feature count; '
             'all cells remain in the same global correction family. Figure 7 uses this exact source table.\n\n'
             '## Recorded-demographic examples: replace Analysis A/B/C placeholders\n\n'
             'Selection was declared in SELECTION_RULES.md before extracting examples: A uses three fixed readable ratios; B/C prioritize supported cell counts and use the same feature sequence, '
             'without filtering by q or effect direction. Nonsignificant examples are retained. The selected effects describe source-labeled groups and cannot be interpreted as genetic ancestry, '
             'intrinsic group characteristics or causal demographic effects.\n\n')
    text += 'The examples use these exact S6 definitions (all unitless):\n\n'
    text += table(demos[['feature', 'formula', 'units']].drop_duplicates()) + '\n'
    for layer, part in demos.groupby('layer', sort=False):
        label = 'Broad recorded groups' if layer.startswith('layer1') else 'Detailed recorded groups'
        text += f'### {label}\n\n'
        a = part.loc[part.analysis.eq('A')]
        text += ('**A: reference-only differences.** ' + '; '.join(
            f"{row.feature}: H={row.omnibus_H:.6g}, q={row.q:.6g}" for row in a.itertuples()) + '. '
            'Each omnibus q uses the layer’s 120-feature family. An omnibus association does not identify a particular group pair or a disease effect; '
            'no pairwise interval is inferred from this output. The following means/medians retain the six-significant-digit precision encoded in the source A table.\n\n')
        text += table(groups.loc[groups.layer.eq(layer), ['feature', 'group', 'n', 'mean', 'median']]) + '\n'
        text += '**B: syndrome versus matched reference.**\n\n'
        for row in part.loc[part.analysis.eq('B')].itertuples():
            text += (f"Within the recorded {row.group_1} group, {row.syndrome} ({row.n_disease_1:.0f} disease, {row.n_reference_1:.0f} reference) "
                     f"has {row.feature} Hedges’ g={row.effect:+.6g}, marginal 95% bootstrap CI [{row.ci_lower:.6g}, {row.ci_upper:.6g}], "
                     f"global layer-B q={row.q:.6g} across {row.family_hypotheses:,} tests. "
                     f"The normalized means are {row.disease_mean_1:.6g} versus {row.reference_mean_1:.6g}, a difference of {row.mean_difference:+.6g} feature units.\n\n")
        text += '**C: between-group differences in disease contrasts.**\n\n'
        for row in part.loc[part.analysis.eq('C')].itertuples():
            text += (f"For {row.syndrome}, {row.feature}, group 1={row.group_1} (disease/reference {row.n_disease_1:.0f}/{row.n_reference_1:.0f}) "
                     f"and group 2={row.group_2} ({row.n_disease_2:.0f}/{row.n_reference_2:.0f}), the group-2-minus-group-1 disease/reference interaction is "
                     f"{row.effect:+.6g} normalized feature units (95% HC3 CI [{row.ci_lower:.6g}, {row.ci_upper:.6g}]; "
                     f"global layer-C q={row.q:.6g}, {row.family_hypotheses:,} tests). "
                     'This is a difference of contrasts, not a Hedges’ g or reference-only group difference.\n\n')
    text += ('These support-selected examples show what global counts conceal. Broad-group forehead height has no omnibus BH association at q=0.0546, '
             'although the broad Asian Williams–Beuren syndrome-reference contrast is large (g=1.453). The detailed Kabuki forehead interaction is significant '
             '(q=0.0176), while the selected detailed Cornelia de Lange philtrum and Williams–Beuren pupil-spacing interactions are not. '
             'A shared reference difference therefore cannot be substituted for evidence that a disease contrast differs across groups. All C models remain exploratory and source-confounded.\n\n'
             'Table caption: Preselected A/B/C examples with source-label names, complete cell counts, exact feature formulas/units and the appropriate correction family. '
             'A reports omnibus H; B reports syndrome-minus-reference Hedges’ g and marginal bootstrap intervals; C reports a group-2-minus-group-1 interaction in normalized units and HC3 intervals. '
             'Global q is separate for every layer/analysis. Full estimates and formulas are in demographic_examples.csv; detailed A summaries are in demographic_A_group_summaries.csv.\n\n'
             '## Adjusted results: replace Table 5 placeholder\n\n')
    displayed = adjusted[['syndrome', 'feature', 'comparator', 'beta', 'ci_lower', 'ci_upper', 'n_complete', 'n_disease', 'n_comparator', 'q', 'family_hypotheses']].copy()
    displayed['syndrome'] = displayed.syndrome.map(short)
    text += table(displayed) + '\n'
    text += ('Table caption: Seven prespecified adjusted models and five anchor cells from the full adjusted screen. Coefficients and HC3 95% intervals are in the feature’s normalized units, '
             'not pooled-SD units. The first four rows use age/sex adjustment and a four-model BH family; the next three use recorded sex, broad group and yaw/pitch/roll and a three-model BH family. '
             'The remaining five retain the full 5,880-model family rather than being recorrected as a selected subset. Complete-case n and target/comparator counts are shown; all unsupported screen cells '
             'retain correction p=1 in the planned denominator. All exact family labels, formulas and source filenames are in adjusted_examples.csv.\n\n')
    for row in adjusted.loc[adjusted.family_hypotheses.eq(3)].itertuples():
        text += (f"Against FairFace, the covariate-adjusted {row.feature} coefficient for {row.syndrome} is {row.beta:+.8g} "
                 f"(95% HC3 CI [{row.ci_lower:.8g}, {row.ci_upper:.8g}]; three-model q={row.q:.8g}; "
                 f"complete-case n={row.n_complete}, disease/reference={row.n_disease}/{row.n_comparator}).\n\n")
    rts = adjusted.loc[adjusted.family_hypotheses.eq(4) & adjusted.feature.eq('malar_bulge_mean')].iloc[0]
    text += (f"The age/sex-adjusted within-GMDB Rubinstein–Taybi malar coefficient is {rts.beta:+.8g} "
             f"(95% HC3 CI [{rts.ci_lower:.8g}, {rts.ci_upper:.8g}]; four-model q={rts.q:.8g}). "
             'It is not significant at q < .05. The FairFace association must not be described as persisting under every comparator or adjustment.\n\n'
             '### Why 48 primary cohorts but 49 adjusted targets\n\n'
             'Primary eligibility is computed from the current selected-patient counts with n >= 10. The separately prespecified adjusted target list retains 49 labels, '
             'including FBXW7__cbba0000, which has nine current selected patients. All other 48 targets match the primary set. The preparation source fixed the main list from current counts '
             'and retained the earlier prespecified target list for extension bookkeeping; no additional current measurement or obsolete estimate is used to enlarge the primary family. '
             'extensions.py passes plan["targets"] to 48-cohort comparisons and plan["adjustment_targets"] to the full adjusted screen. Its support rules reject sparse target cells; '
             'the FBXW7 rows and their correction denominator remain rather than forcing the lists to coincide. The active ANALYSIS_PLAN.json, syndrome_eligibility.csv, '
             'extensions.py and full_adjusted_screen.csv establish this distinction. Archived preparation code was consulted only to audit how the target list was declared.\n\n'
             '### Comparator nonagreement\n\n'
             'The primary FairFace family has 955 associations and the all-other-GMDB family 222; 173 are shared, 782 FairFace-only and 49 GMDB-only. '
             'These comparisons estimate different contrasts. They do not imply that 782 primary results are false positives or estimate a clinical false-positive rate. '
             'Source effects, broad disease-shared morphology and comparator variance can all change which associations pass their separate global BH thresholds.\n\n'
             '## Clinical literature reassessment\n\n'
             'All 17 appendix questions were linked to a unique active syndrome-feature row and its exact S6 definition; their old estimates were not reused. '
             'All 17 active associations pass the current global threshold, but that does not validate their clinical meaning. Four have directly related trait support: '
             'Cornelia de Lange/philtrum length, Noonan/pupil spacing, Wiedemann–Steiner/inner-canthal spacing and Sotos/mid-forehead width. '
             'Even these do not validate the exact normalized ratios as clinical measurements. The other 13 remain related geometric observations, scope mismatches or unverified interpretations, '
             'and should not be retained as established clinical confirmations. The five original anchors may remain as geometric examples with their existing qualifications.\n\n')
    for row in literature.itertuples():
        text += (f"**{row.appendix_entry}. {short(row.syndrome)} / {row.feature}.** n={row.n}; g={row.hedges_g:+.6g}, "
                 f"95% CI [{row.ci_lower:.6g}, {row.ci_upper:.6g}], global q={row.global_q:.6g}. "
                 f"{row.interpretation}. Status: {row.clinical_support.replace('_', ' ')}. "
                 f"[Clinical source: {row.clinical_source_section}]({row.clinical_source}).\n\n")
    text += ('## Separate synthetic consistency: ready-to-paste framing\n\n'
             'The primary PRISM study analyzes real GMDB patient photographs against selected FairFace reference photographs. The mentor reports that PDIDB contains '
             'syndrome-conditioned images sampled from a StyleGAN model trained on GMDB. Because this learned-source relationship does not supply independent patients, '
             'PDIDB measurements were evaluated separately for conditional distributional fidelity. Checkpoint-specific generation lineage remains unverified. '
             'Nine matched syndromes cover 314 selected GMDB patients and 1,424 QC-passing synthetic images. Across syndrome summaries, median absolute mean shift was 0.166 real SD, '
             'median variance ratio 0.986, median effect-direction agreement 85.0% and median Pearson effect-profile correlation 0.896 against the shared frozen FairFace baseline. '
             'Age/sex/race are missing for all synthetic records, preventing equivalent age filtering and limiting interpretation of these similarities and discrepancies. '
             'No synthetic-derived hypothesis tests contribute support to the 955 primary associations. Full findings and sensitivities are in synthetic_consistency_results.md.\n\n'
             '## Figure captions and placement\n\n')
    text += (TABLES / 'FIGURE_CAPTIONS.md').read_text().split('\n', 1)[1] + '\n'
    text += ('Place Figures 5/6 in Methods to replace workflow/anatomy placeholders, Figure 7 beside regional Results, and Figures 8–10 in the separate synthetic subsection or supplement. '
             'These repository numbers continue after existing Figures 1–4; renumber in the Google Doc only when inserted.\n\n'
             '## Results-source map and unresolved items\n\n'
             '| Formal Paper placeholder / claim | Exact current source and insertion |\n| --- | --- |\n'
             '| Top-five recurrent features | S1 → top_five_recurrent_features.csv; primary recurrence paragraph |\n'
             '| Nine regions and regional figure | regional_summary_120.csv + independent S1 regrouping → regional_discovery_source.csv; Figure 7 |\n'
             '| Demographic A/B/C examples | Complete layer1/layer2 A/B/C CSVs → demographic_examples.csv and demographic_A_group_summaries.csv; six separate families |\n'
             '| Table 5 adjusted models | within_gmdb_age_sex_adjusted.csv; three_predefined_reference_adjusted_models.csv; full_adjusted_screen.csv → adjusted_examples.csv |\n'
             '| Comparator overlap | S1 + within_gmdb_results.csv → comparator_overlap.csv |\n'
             '| Appendix A1–A17 | S1 + S6 + cited clinical sections → clinical_literature_audit_17.csv |\n'
             '| Extraction/workflow and geometry art | ACTIVE_ANALYSIS.json + S6 + exact pinned source → Figures 5/6 |\n'
             '| New PDIDB fidelity | synthetic_consistency_manifest.json and its aggregate CSVs → Figures 8–10 and separate results report |\n\n'
             'Pending scientific judgments: confirm the nine broad syndrome mappings clinically; establish checkpoint/generator release identity and exact training membership; '
             'obtain supported age metadata before claiming age-matched synthetic fidelity; validate the exact eyebrow/malar and other related geometric constructs clinically. '
             'The Coffin–Siris long-philtrum interpretation remains unverified in the retrieved source. No new clinical novelty is asserted from absence in a general description. '
             'The separate FaceKit overlap audit defines the limited statistical contribution. No external clinical validation, diagnostic accuracy, independent synthetic replication '
             'or new extractor claim is supported. Google Docs placeholders/comments still require insertion/resolution; BIBM and PRISM tabs were not edited. '
             'The handoff excludes privacy/IRB/legal release review, history remediation and full frontend testing; existing PI items remain pending.\n')
    (ROOT / 'manuscript/FINAL_120_RESULTS_INSERTIONS.md').write_text(text)


if __name__ == '__main__':
    synthetic_report()
    insertions()
    path = SYNTHETIC / 'synthetic_consistency_manifest.json'
    manifest = json.loads(path.read_text())
    manifest['outputs']['synthetic_consistency_results.md'] = hashlib.sha256((SYNTHETIC / 'synthetic_consistency_results.md').read_bytes()).hexdigest()
    manifest['report_generation_code_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    manifest['figure_source_manifest_sha256'] = hashlib.sha256((TABLES / 'FIGURE_SOURCE_MANIFEST.json').read_bytes()).hexdigest()
    path.write_text(json.dumps(manifest, indent=2) + '\n')
    print('Generated source-tied synthetic report and manuscript insertions')

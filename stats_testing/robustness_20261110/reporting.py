import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .common import DEFAULT_OUT, private_output, protect_baseline, require, sha, write_aggregate


def build(out):
    folder = out / 'aggregate_candidates'
    matrix = pd.read_csv(out / 'evidence/P0_exhaustive_evidence_matrix.csv', low_memory=False)
    geometry = pd.read_csv(folder / 'P1_normalization_cases.csv')
    coverage = pd.read_csv(folder / 'P1_bootstrap_coverage_simulation.csv')
    age = pd.read_csv(folder / 'P1_age_case_summary.csv')
    age_effects = age.pivot(index=['canonical_syndrome_id', 'feature'], columns='stratum', values='hedges_g')
    require(len(age_effects) == 6 and age.status.eq('ok').all(), 'Predeclared age cases have incomplete support')
    age_agreement = int((np.sign(age_effects.iloc[:, 0]) == np.sign(age_effects.iloc[:, 1])).sum())
    quality = pd.read_csv(folder / 'P0_quality_adjusted_cases.csv')
    candidates = pd.read_csv(folder / 'P2_ranked_candidate_clinical_status.csv')
    same = int(matrix.same_direction.sum())
    strong = matrix.loc[matrix.ff_mannwhitney_q_fdr_global.lt(.05) & matrix.ff_hedges_g.abs().ge(.8)]
    weakened = int((strong.gmdb_hedges_g.abs() < strong.ff_hedges_g.abs() / 2).sum())
    reverse = int((np.sign(strong.gmdb_hedges_g) != np.sign(strong.ff_hedges_g)).sum())
    require((len(matrix), same, len(strong), weakened, reverse) == (5760, 4107, 482, 211, 9), 'P3 comparator claim check failed')
    claims = [
        ('R00', 'The frozen study includes 1,845 selected recorded GMDB patients, including 872 across 48 primary cohorts, and 3,409 FairFace reference paths; zero synthetic observations enter primary inference.',
         '00_baseline_validation.json', 'private_counts; frozen_counts; active config', 'Frozen baseline; locally revalidated'),
        ('R01', f'Across all 5,760 pairs, {same:,} ({same / 5760:.1%}) had the same effect direction in the two comparator designs.',
         'evidence/P0_exhaustive_evidence_matrix.csv', 'same_direction; complete supported grid', 'Descriptive; overlapping cases, different estimands; no concordance CI'),
        ('R02', f'Among {len(strong)} primary associations with global q<0.05 and |g|≥0.8, {weakened} ({weakened / len(strong):.1%}) had less than half the effect magnitude against other GMDB conditions and {reverse} reversed sign.',
         'evidence/P0_exhaustive_evidence_matrix.csv', 'ff_mannwhitney_q_fdr_global<.05 & abs(ff_hedges_g)>=.8; abs(gmdb_hedges_g)<abs(ff_hedges_g)/2; sign comparison', 'Exploratory threshold; not a clinical false-positive rate'),
        ('R03', 'Average-linkage clustering of 1−|Spearman rho| at distance 0.3 yielded 55 descriptive groups in each global source; source-partition adjusted Rand agreement was 0.857.',
         'aggregate_candidates/P1_cluster_source_comparison.csv', 'source=FairFace/GMDB_selected; clusters and adjusted_rand_to_FairFace', 'Descriptive; not independent biological dimensions or an effective testing denominator'),
        ('R04', 'Only 54 of 872 primary patients changed photographs under the deterministic alternate rule, and 94 had more than one eligible image.',
         'aggregate_candidates/P1_primary_photo_opportunities.csv', 'sum changed_alternate; sum patients_with_multiple_eligible_images', 'Limited observational photograph perturbation; not extractor-transform robustness'),
        ('R05', 'The existing 100 random-photo repetitions changed a median of 53 primary photographs per repetition (range 38–65).',
         'aggregate_candidates/P1_random_photo_change_counts.csv', 'changed_primary_photos; original seed20261025', 'Same recorded patients; exact replay of frozen runs'),
        ('R06', 'Across 100 repetitions per design, stratified random within-FairFace partitions produced zero discoveries, whereas deliberate within-stratum pitch imbalance produced a median of 25 discoveries (range 0–73).',
         'aggregate_candidates/P0_negative_control_repetitions.csv', 'design, discoveries; corrected v2 allocation;5760 global-BH tests per repetition', 'Workflow susceptibility; not clinical false-positive rate'),
        ('R09', f'The new fixed-reference case-bootstrap procedure had empirical nominal-95% coverage from {coverage.empirical_coverage.min():.1%} to {coverage.empirical_coverage.max():.1%} across four numerical regimes.',
         'aggregate_candidates/P1_bootstrap_coverage_simulation.csv', 'empirical_coverage; normal/scaled_t5 x n10/35;300datasets/regime;200draws;fixedreference500', 'Simulation diagnostic; not validation/invalidation of frozen two-sample primary bootstrap'),
        ('R10', 'The high-magnitude top-20 list shared 16–19 pairs with photo-perturbed lists (median 18).',
         'aggregate_candidates/P1_top20_photo_variation.csv', 'top20_shared_with_primary; |g|,case n>=20,100same-patient repeats', 'Exploratory ranking sensitivity; not held-out performance'),
        ('R11', f'The source review assessed 17 selected prior clinical questions and marked {len(candidates)} unique ranked candidate pairs; no exact normalized clinical-measure equivalence was established.',
         'P2_phenotype_source_matrix.csv', '17rows; exact_comparable_measure=False; candidate denominator separately in P2_ranked_candidate_clinical_status.csv', 'Selected literature audit; not all955; unmapped does not mean novel'),
        ('R12', f'All {age_agreement} predeclared case effects retained direction in both within-GMDB age strata, but none met the age-family q<0.05 threshold in ages 3–<6 and three did in ages 6–<10.',
         'aggregate_candidates/P1_age_case_summary.csv', '6cases x2strata;status=ok;hedges_g,q_global_bh;original11760-cell age family',
         'Different support/comparator distributions; not a formal age-interaction test or exact FairFace matching'),
        ('R13', 'In the new age/sex/broad-group/pose-adjusted other-GMDB models, Williams–Beuren eyebrow area and outer-canthal span did not meet the separate 12-model q<0.05 threshold.',
         'aggregate_candidates/P0_quality_adjusted_cases.csv', 'WBS eb_thickness_mean and outer_canthal_distance;comparator=other_GMDB;exploratory_q',
         'New HC3 complete-case model/family; not interchangeable with frozen5880-screen q; no source-confounding identification'),
    ]
    for identifier, feature, label in [('R07', 'outer_canthal_distance', 'Williams–Beuren outer-canthal span'),
                                       ('R08', 'philtrum_width', 'Cornelia de Lange upper-lip-peak span')]:
        part = geometry.loc[geometry.feature.eq(feature) & geometry.comparator.eq('FairFace')].set_index('definition')
        sentence = (f'{label} changed from g={part.loc["original_bizygomatic", "hedges_g"]:+.3f} under bizygomatic normalization '
                    f'to g={part.loc["height_normalized", "hedges_g"]:+.3f} under height normalization against FairFace.')
        claims.append((identifier, sentence, 'aggregate_candidates/P1_normalization_cases.csv',
                       f'feature={feature}; FairFace;original_bizygomatic vs height_normalized;new12-cell BH family',
                       'Exploratory different estimands; no independent raw component attribution; original q unchanged'))
    ledger = pd.DataFrame(claims, columns=['claim_id', 'exact_draft_sentence', 'source', 'columns_or_filter', 'inference_scope']).sort_values('claim_id')
    ledger['source_sha256'] = [sha(out / source) for source in ledger.source]
    ledger['method_version'] = 'stats_testing.robustness_20261110'
    ledger['clinical_source'] = np.where(ledger.claim_id.eq('R11'), 'PUBLIC_SOURCE_NOTES.json; URLs/sections/access states in source matrix', '')
    ledger['release'] = 'PRIVATE draft; PI/disclosure review pending'
    write_aggregate(ledger, out / 'P3_claim_to_evidence_ledger.csv')
    (out / 'P3_claim_to_evidence_ledger.md').write_text(
        '# Claim-to-evidence ledger\n\nExact draft sentences, source paths, source SHA-256, filters, method versions and inference scopes are in '
        'P3_claim_to_evidence_ledger.csv. Bracketed R identifiers resolve to those rows. All new empirical robustness/ranking claims are exploratory. '
        'Clinical claims require source scope and access-state review; the 17 selected questions do not represent all 955 discoveries.\n')
    methods = (
        'PRISM is an exploratory statistical investigation built on the pinned FaceKit 120-feature extractor. We retained the original '
        '48×120 primary family and analyzed syndrome-minus-FairFace and syndrome-minus-other-GMDB estimands separately. '
        'Complete hypothesis grids were joined by canonical syndrome identifiers and exact ordered features. Two-sided asymptotic Mann–Whitney tests '
        'used continuity correction and BH correction within each declared family; signed Hedges g described standardized mean differences. '
        'Original marginal intervals retain their original methods and do not estimate uncertainty for dependent comparator differences. '
        'We did not estimate a concordance-bootstrap confidence interval.\n\n'
        'New FairFace controls used 100 repetitions per design, 48 disjoint pseudo-case groups totaling 872 paths and a disjoint shared reference half. '
        'Splitting was stratified by recorded race/sex; finite group sizes allowed approximate quota matching at allocation. The deliberate design '
        'separated absolute pitch within strata before random group allocation. Each repetition corrected all 5,760 tests together. Paths were not '
        'assumed to identify unique people. Processed dimensions do not establish native resolution, crop or exposure.\n\n'
        'Case-patient resampling preserved each full correlated feature row and held each empirical comparator fixed: 200 bootstrap draws and '
        '100 subsamples retaining 70% of cases. Subsampling was unsupported when retained n<10. Leave-one-patient summaries were limited to n≥20 '
        'and exported only aggregate maxima/fractions. The 100 frozen photograph-selection runs were replayed with their original ordering and seed '
        '20261025, without image processing. New percentile ranges are descriptive because coverage was imperfect in numerical checks. Cliff delta, '
        'median and 10%-trimmed-mean contrasts supplement g with distinct interpretations and scales.\n\n'
        'Spearman and Pearson matrices described dependence by source and supported syndrome; average linkage at 1−|rho|=0.3 and 100 global-source '
        'resamples described cluster stability. Reference-only representatives favored non-asymmetry, fewer derived references, shorter formulas '
        'and lexical identifiers. Choosing largest displayed syndrome effects was outcome-dependent. The original BH denominator was preserved. '
        'Candidate filters were specified after P0/P1 review and are post hoc. Three additional Monte Carlo seeds evaluated the same patients, '
        'without providing independent validation.\n\n'
        'Two geometry cases×three definitions×two comparators formed a new 12-test BH family. Normalized outputs identified alternative ratios, '
        'not raw numerator/denominator values. New 12-case HC3 quality models, two 5,760-cell pose-restricted families, source/group 720-cell '
        'correlation families and six ecological proxy tests were corrected separately. Existing adjusted-screen (5,880), targeted (22,320), '
        'age-stratum (11,760) and four-/three-model results retain distinct scopes. Unsupported estimates/display q remain missing; conservative '
        'internal placeholders preserve planned multiplicity only. FairFace exact-age matching and a separate source/disease effect were '
        'not identifiable. PDIDB contributed no new inferential observations.\n')
    results = '\n\n'.join(f'{r.exact_draft_sentence} [{r.claim_id}]' for r in ledger.itertuples())
    case_text = ''
    for syndrome, feature in [('Cornelia de Lange', 'philtrum_length'), ('Noonan', 'inter_pupillary_distance'),
                              ('WILLIAMS-BEUREN', 'eb_thickness_mean'), ('Rubinstein-Taybi', 'malar_bulge_mean')]:
        r = matrix.loc[matrix.syndrome.str.startswith(syndrome) & matrix.feature.eq(feature)].iloc[0]
        case_text += (f'{r.syndrome.split("__")[0]} / {feature}: primary g={r.ff_hedges_g:+.3f}; other-GMDB g={r.gmdb_hedges_g:+.3f}, '
                      f'other-GMDB global q={r.gmdb_q_global_bh:.4g}.\n\n')
    discussion = (
        'Comparator and normalization sensitivity provide the clearest contribution. Shared associations indicate agreement between dependent '
        'questions, not independent replication; loss of significance is not a false-positive estimate. Pitch imbalance demonstrates observable '
        'within-source susceptibility without isolating measurement error from correlated appearance. A different denominator changes the '
        'anatomical estimand; the Williams–Beuren height-normalized point estimate near zero does not establish a robust reversed trait or prove '
        'the original ratio erroneous. Cornelia de Lange length, Noonan pupil spacing and Williams–Beuren eyebrow area remain directionally '
        'positive across comparators, while Rubinstein–Taybi malar geometry is markedly attenuated.\n\n'
        'Clinical descriptions are related constructs, not calibration of the exact ratios. Long-philtrum appearance is supported by the '
        '[Cornelia de Lange consensus](https://pmc.ncbi.nlm.nih.gov/articles/PMC7136165/); '
        '[Noonan GeneReviews](https://www.ncbi.nlm.nih.gov/sites/books/NBK1124/) describes eye spacing; '
        '[Williams GeneReviews](https://ncbi.nlm.nih.gov/sites/books/NBK1249/?report=classic) describes facial appearance. Upper-lip peaks need '
        'not be philtral columns and eyebrow polygon area is not hair density. Dates, sections and access limits are recorded in PUBLIC_SOURCE_NOTES.json. '
        'No reviewed ratio is an established diagnostic marker or new clinical sign.\n')
    limitations = (
        'Source and case/reference status align in the primary contrast, so acquisition and disease cannot be separated. Other GMDB conditions '
        'are not normal controls; FairFace labels do not verify health, exact ages, unique people or biological ancestry. Selection, small n, '
        'demographic missingness, correlated features and shared controls constrain interpretation. Only 94 primary patients offered multiple '
        'eligible photographs. Fixed empirical comparators omit comparator-sampling uncertainty and small-sample percentile ranges undercovered '
        'in simulation. Post hoc selection and repeated use of the same patients preclude unbiased ranking-performance claims. Compatible saved '
        'landmarks/components were unavailable and image reprocessing was not approved. The clinical review covers selected questions. No external '
        'clinical replication, novel-method superiority or independent synthetic validation was established.\n')
    text = '# Private manuscript insertion drafts\n\nMain narrative: comparator and normalization sensitivity. Existing manuscripts were not edited.\n\n'
    for title, body in [('Methods', methods), ('Results', results), ('Highlighted case checks', case_text),
                         ('Discussion', discussion), ('Limitations', limitations)]:
        text += f'## {title}\n\n{body}\n\n'
    text += ('## Provenance\n\nR00–R13 resolve to the claims ledger. Highlighted checks use evidence/P0_exhaustive_evidence_matrix.csv: '
             'ff_hedges_g, gmdb_hedges_g and gmdb_q_global_bh (the latter only the all-other-GMDB 5,760 family). '
             'Method parameters are in ANALYSIS_PLAN.md and script modules; input hashes/environment in 00_baseline_validation.json and '
             'RESEARCH_PROVENANCE.json. Clinical/PI and disclosure review are required before incorporation/export.\n')
    (out / 'P3_manuscript_insertions.md').write_text(text)
    sensitivity_text = (
        '# Age and observable-covariate case review\n\n'
        'These are six predeclared examples, not an outcome-selected age-interaction screen. '
        'Existing age-stratum g and q values retain the full 11,760-cell age-family correction. '
        'All six signs agree across ages 3–<6 and 6–<10. No younger-stratum example has q<.05; three older examples do. '
        'Age-specific other-GMDB case/comparator composition, variances and support differ, so these contrasts do not establish an age interaction. '
        'Exact-age matching against FairFace remains impossible. Full estimates and marginal intervals are in P1_age_case_summary.csv.\n\n')
    for row in age.itertuples():
        sensitivity_text += (f'- {row.syndrome.split("__")[0]} / {row.feature}, {row.stratum}: '
                             f'n={row.n_a}/{row.n_b}, g={row.hedges_g:+.4f}, age-family q={row.q_global_bh:.5g}.\n')
    sensitivity_text += (
        '\nNew quality models use sex, recorded broad group and signed pose; other-GMDB also adjusts recorded exact age. '
        'All 12 models are supported under the declared complete-case thresholds, but estimates concern selected complete cases. '
        'Beta and intervals are in feature units, not g, and cannot be compared in magnitude across different features. '
        'WBS eyebrow area and outer span weaken under this specification despite their positive point directions; their q here must not '
        'be replaced with the much smaller frozen adjusted-screen q. Source and case status are still not separately identifiable.\n\n')
    for row in quality.itertuples():
        sensitivity_text += (f'- {row.syndrome.split("__")[0]} / {row.feature}, {row.comparator}: '
                             f'n={row.n_target}/{row.n_comparator}, beta={row.beta:+.5g}, '
                             f'HC3 marginal interval [{row.ci_lower:+.5g}, {row.ci_upper:+.5g}], new 12-model q={row.exploratory_q:.5g}.\n')
    restrictions = pd.read_csv(folder / 'P0_quality_restrictions.csv')
    sensitivity_text += '\nPose restriction summaries (each planned 5,760-cell family; unsupported q remains missing):\n\n'
    for threshold, part in restrictions.groupby('threshold'):
        sensitivity_text += f'- pose_score≤{threshold}: {part.supported.sum()} supported cells, {part.exploratory_q.lt(.05).sum()} discoveries.\n'
    sensitivity_text += ('\nRestriction changes support and composition, so discovery totals are not directly comparable power-adjusted performance. '
                         'No restriction or model identifies independent biology. Retain both favorable and weakening cases.\n')
    (out / 'P1_covariate_age_sensitivity.md').write_text(sensitivity_text)
    (out / 'P3_editorial_decision_memo.md').write_text(
        '# Editorial decision memo\n\nRecommended main narrative: **reference/comparator sensitivity and observable acquisition susceptibility**. '
        'Present an exploratory statistical study built on FaceKit. Retain compact phenotype examples as illustrations; consider demographics '
        'and the separate synthetic catalog for supplements after review.\n\n'
        'Decisive evidence: 4,107/5,760 directions agree, but 211/482 strong primary effects more than halve against other GMDB conditions. '
        'The 173 shared discoveries are dependent. Pitch-imbalanced controls produce discoveries without clinical syndrome contrasts. '
        'Williams–Beuren outer-span normalization, its weakened quality-adjusted evidence, and Rubinstein–Taybi malar attenuation are important counterexamples. '
        'Limited photograph opportunities and imperfect bootstrap coverage restrict broad robustness claims. Keep negative findings central.\n\n'
        'A standalone robust phenotype-map claim is premature without comparable-measure clinical assessment and source-matched controls. '
        'A benchmark-method paper is NOT_JUSTIFIED: no distinct score or honest performance advantage was established. FaceKit already supplies '
        'the extractor and related phenotype/reference questions (see the unchanged publication/final_research_120/FACEKIT_OVERLAP_AUDIT.md). '
        'The incremental contribution is traceable comparator/geometry evidence and explicit identification limits.\n\n'
        'Potential journal types are computational clinical-genetics measurement studies or reproducibility-focused biomedical informatics papers. '
        'This is a scope suggestion, not a claim about current journal policy, eligibility or acceptance probability.\n\n'
        'GO for private author/mentor review. NO-GO for manuscript promotion, submission, public figures/tables, new clinical signs, diagnostic '
        'claims or method-superiority claims until PI review of permissions, disclosure risk, clinical wording and contribution relative to FaceKit. '
        'Raw-component and real-image experiments require compatible inputs and explicit reprocessing authorization. '
        'Absence from a clinical summary does not establish novelty.\n')
    checklist(out)
    figure_plan(out)
    print('Private editorial memo, insertion drafts, claims ledger, figure plan and reporting checklist created')


def checklist(out):
    rows = [
        ('Design/question', 'Draft complete', 'Exploratory observational geometry/robustness study; no diagnosis or clinical-validation claim.'),
        ('Participants/flow/eligibility', 'Local verified; setting description PI review', '1845 recorded GMDB patients;872 in48 primary cohorts;3409 FF paths. Retain frozen flow/exclusions: current120,valid ages,frontal/QC,one photo/patient,syndrome n>=10. Do not infer unique FF people.'),
        ('Variables/measurement', 'Provenance complete; gaps explicit', 'Pinned FaceKit/model/S6; units/denominators; categorical age and recorded population labels; no absolute dimensions.'),
        ('Missingness', 'Documented; manuscript integration pending', 'Proxy inventory and full adjusted/age support tables identify missingness/exclusions. Unsupported estimates/display q missing. Crop/confidence/native resolution/exposure absent.'),
        ('Bias/study size', 'Sensitivity complete; identification unresolved', 'Source and case confounded; other conditions not normal controls; no retrospective power justification invented.'),
        ('Statistical methods', 'Draft complete', 'Separate BH families; signed effects; patient grouping/fixed comparators; seeds; n thresholds; post hoc choices; no biased concordance CI.'),
        ('Main/other results', 'Private evidence complete', 'Full grids and unfavorable results; photo opportunity counts; age/quality/reference/dependence/estimator checks; no independent replication interpretation.'),
        ('Discussion/generalizability', 'Draft complete', 'Different estimands, coverage, literature scope, missing normal controls; no exact clinical validation or novel-method claim.'),
        ('Code/environment', 'Local verification complete', 'Completion/provenance record commands/versions/hashes. No new remote CI or push.'),
        ('IRB/consent/data use', 'PI REVIEW REQUIRED', 'PI must supply protocol ID, approved scope and consent/waiver; original image reprocessing/publication rights unresolved. No approval inferred.'),
        ('FairFace/bibliography', 'Source review done; terms PI review', 'Margin0.25 provenance and attribution/use conditions retained; source dates/passages/access states documented; exact measure equivalence not established.'),
        ('Synthetic records', 'Separate; not rerun', 'Zero PDIDB observations in inference; lineage/independence/illustration terms need PI/provider review.'),
        ('Privacy/release/history', 'PRIVATE; PI REVIEW REQUIRED', 'No photos/individual rows/landmarks/credentials exported. Aggregate rare-cohort disclosure review remains. No upload, permission change, history rewrite or manifest promotion.'),
        ('Authorship/funding/conflicts', 'AUTHOR/PI INPUT REQUIRED', 'No fabricated contributions, affiliations, funding or conflict statements.'),
    ]
    frame = pd.DataFrame(rows, columns=['reporting_item', 'status', 'evidence_or_action'])
    write_aggregate(frame, out / 'P3_reporting_checklist.csv')
    text = ('# Reporting and permission checklist\n\nThis task-specific review is not certification of every STROBE item. '
            'The [official STROBE checklist collection](https://www.strobe-statement.org/checklists/) was accessed November 10, 2026; '
            'use the appropriate full design checklist for the final manuscript.\n\n')
    text += '\n'.join(f'- **{r.reporting_item} — {r.status}:** {r.evidence_or_action}' for r in frame.itertuples()) + '\n'
    (out / 'P3_reporting_checklist.md').write_text(text)


def figure_plan(out):
    details = {
        'P0_comparator_scatter': ('Main', 'x primary g;y other-GMDB g;panels case n10–19/20–39/40+;four q categories', '5760 cells,48syndromes;872cases,FF3409,GMDB1845'),
        'P0_comparator_categories': ('Supplement', 'x hypothesis count;y four categories', '5760 cells;separate5760 BH families'),
        'P0_discordant_cases': ('Main; negative findings', 'x g/original marginal intervals;y shortened syndrome,canonical ID,feature', '8 post hoc discordant/attenuating strong pairs;exact n in source'),
        'P0_negative_control_distribution': ('Main', 'x discoveries per5760 tests;y repeats;design panels', '100 repeats/design;48pseudo-case groups,872paths;disjoint reference half'),
        'P0_recurrent_feature_quality': ('Supplement', 'x pose proxy;y feature;color/text rho', 'GMDB1845,FF3409;source/group720-cell families'),
        'P0_acquisition_proxy_distributions': ('Supplement', 'x recorded signed pose;y recorded broad group;points median/bars IQR;source colors', 'GMDB1845,FF3409 overall;per-group n/missingness in source table'),
        'P1_source_correlation_heatmap': ('Supplement', 'axes S6 index1–120;color rho/source difference', '120features/source;GMDB1845,FF3409'),
        'P1_feature_dependency_network': ('Supplement', 'solid arrows algebraic;dashed shared normalization', 'original schematic;complete120-feature edges in CSV'),
        'P1_reference_cluster_signatures': ('Supplement', 'x g/original marginal intervals;y representatives', '4predeclared cohorts:CdLS64,Noonan34,WBS67,RTS24;FF3409;top5 |g| display selection'),
        'P1_predeclared_case_specifications': ('Supplement', 'x g/marginal intervals;y five frozen specifications', '6prespecified pairs;exact n by specification in evidence matrix'),
        'P1_normalization_sensitivity': ('Main; negative findings', 'x g/descriptive case-bootstrap range;y ratio;source colors', '2cases×3definitions×2comparators=12new BH cells;CdLS64,WBS67;FF3409 or remainingGMDB'),
        'P1_ratio_simulation': ('Supplement', 'x denominator%;y ratio%;lines numerator%', '93 numerical conditions;no patients/physical units'),
        'P1_conditional_bootstrap_coverage': ('Supplement', 'x coverage/MonteCarlo Wilson interval;y regime;nominal.95 line', '4regimes×300independent datasets;200draws;fixed reference500'),
        'P1_candidate_rule_sensitivity': ('Supplement', 'axes six rules;color/text shared top20', 'same patients;case n>=20;post hoc,no held-out performance'),
        'P1_geometry_landmark_schematic': ('Main or supplement', 'invented landmarks/denominators/equations', '2constructs;S6 formulas;no real coordinates'),
    }
    rows = []
    for path in sorted((out / 'figures').glob('*.json')):
        item = json.loads(path.read_text())
        require(path.stem in details, 'Figure lacks editorial/axis review')
        role, axes, denominator = details[path.stem]
        rows.append(dict(figure=path.stem, suggested_role=role, axes_and_legend=axes, denominator=denominator,
                         caption=item['caption_and_alt_text'], source_files='; '.join(item['source_sha256']),
                         metadata_json=str(path.relative_to(out)), status='PRIVATE; PI/disclosure review pending'))
    require(len(rows) == len(details), 'Expected figure set incomplete')
    frame = pd.DataFrame(rows)
    write_aggregate(frame, out / 'P3_figure_plan.csv')
    text = '# Private figure plan\n\nPNG drafts have editable SVGs and JSON source/code/output hashes. No individual vectors or real photographs. Main-figure choices need PI/disclosure review.\n\n'
    for r in frame.itertuples():
        text += (f'## {r.figure}\n\nRole: {r.suggested_role}. Denominator: {r.denominator}.\n\nAxes/legend: {r.axes_and_legend}.\n\n'
                 f'Caption/alt text: {r.caption}\n\nSource: {r.source_files}; provenance: {r.metadata_json}.\n\n')
    (out / 'P3_figure_plan.md').write_text(text)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    out = private_output(args.output_dir)
    protect_baseline(out)
    build(out)
    protect_baseline(out)


if __name__ == '__main__':
    main()

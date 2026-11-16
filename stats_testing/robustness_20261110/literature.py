import argparse
import json
from pathlib import Path

import pandas as pd

from .common import ROOT, DEFAULT_OUT, KEY, features, private_output, protect_baseline, require, write_aggregate, write_json

SOURCES = {
    'NBK1104': ('Cornelia de Lange consensus', 'Kline AD et al.', '2018-07-11',
                'https://pmc.ncbi.nlm.nih.gov/articles/PMC7136165/', 'Clinical diagnostic criteria, Box 1 and Figure 3',
                'verified_full_text_passage', 'Length/smoothness is a cardinal philtrum description; no normalized lip-peak width equivalence'),
    'NBK1124': ('Noonan Syndrome', 'Roberts AE', '2025-12-04', 'https://www.ncbi.nlm.nih.gov/sites/books/NBK1124/',
                'Suggestive Findings; Clinical Description / facial features', 'verified_full_text_passage',
                'Eye spacing is described across ages; appearance varies with age; pupil/bizygomatic ratio is not directly specified'),
    'NBK1249': ('Williams Syndrome', 'Morris CA', '2023-04-13', 'https://ncbi.nlm.nih.gov/sites/books/NBK1249/?report=classic',
                'Suggestive Findings / distinctive facies', 'verified_full_text_passage',
                'Periorbital and cheek appearance is described across ages; neither eyebrow area nor exact normalized outer-canthal span is established'),
    'NBK1526': ('Rubinstein-Taybi Syndrome', 'Stevens CA', '2023-11-09', 'https://www.ncbi.nlm.nih.gov/sites/books/NBK1526/',
                'Summary / craniofacial features', 'verified_full_text_passage',
                'Descriptive craniofacial findings are not a calibration of the projected cheek-to-oval index'),
    'NBK580718': ('Wiedemann-Steiner Syndrome', 'Sheppard SE; Quintero-Rivera F', '2022-05-26',
                  'https://www.ncbi.nlm.nih.gov/sites/books/NBK580718/', 'Summary; Table 2; characteristic dysmorphic features',
                  'verified_full_text_passage', 'Spacing-related appearance is documented; exact inner-canthal ratio remains a different construct'),
    'NBK1479': ('Sotos Syndrome', 'Ocansey S et al.', '2025-06-05', 'https://www.ncbi.nlm.nih.gov/sites/books/NBK1479/',
                'Suggestive Findings / characteristic facial appearance', 'verified_full_text_passage',
                'Forehead shape is described and most recognizable at ages 1–6; normalized forehead width is not directly validated'),
    'NBK62111': ('Kabuki Syndrome', 'Adam MP; Hannibal M', '2026-04-23', 'https://www.ncbi.nlm.nih.gov/sites/books/NBK62111/',
                 'Suggestive Findings; diagnostic facial criteria', 'verified_full_text_passage',
                 'Individual fissure length and eyelid eversion differ from total outer-canthal or inner-canthal span'),
    'NBK487886': ('ANKRD11-Related KBG Syndrome', 'Swols DM; Tekin M', '2026-06-11', 'https://www.ncbi.nlm.nih.gov/sites/books/NBK487886/',
                  'Summary; Suggestive Findings / characteristic facial appearance', 'verified_full_text_passage',
                  'Triangular outline does not directly measure either bigonial/bizygomatic or jaw/forehead ratio'),
    'NBK1396': ('Alpha-Mannosidosis', 'Ficicioglu C; Stepien KM', '2024-06-13', 'https://www.ncbi.nlm.nih.gov/sites/books/NBK1396/',
                'Clinical Description / facial findings', 'indirectly_reported_for_exact_candidate',
                'Metadata/full chapter retrieved; exact pupil-spacing passage not established in this review; no absence inference'),
    'NBK1274': ('Mucopolysaccharidosis Type II', 'Scarpa M; Lampe C', '2026-04-28', 'https://www.ncbi.nlm.nih.gov/sites/books/NBK1274/',
                'Suggestive Findings / symptomatic individual', 'verified_full_text_passage',
                'Type-II childhood descriptions cannot validate a pooled mucopolysaccharidoses cohort or a specific chin angle'),
    'NBK1530': ('Holoprosencephaly Overview', 'Tekendo-Ngongang C; Muenke M; Kruszka P', '2026-06-09',
                'https://www.ncbi.nlm.nih.gov/sites/books/NBK1530/', 'Clinical Characteristics / craniofacial spectrum',
                'verified_full_text_passage', 'Structural midline anomalies and clefts are not the photograph-dependent mouth-opening index'),
    'NBK131811': ('Coffin-Siris Syndrome', 'Schrier Vergano S et al.', '2026-11-05', 'https://www.ncbi.nlm.nih.gov/sites/books/NBK131811/',
                  'Clinical Description; separate 2012 clinical series below', 'verified_full_text_chapter_but_candidate_not_established',
                  'Current chapter retrieval alone does not settle the exact length interpretation; independent clinical reports have genotype/case-specific heterogeneity')}


def build(out):
    definitions = features()
    potential = []
    for row in definitions.itertuples():
        label = row.feature.replace('_', ' ')
        term, hpo, relation = label + ' geometry', '', 'anatomical construct only; clinical equivalence not established'
        if row.feature == 'philtrum_length':
            term, hpo = 'Long philtrum (potential related trait)', 'HP:0000343'
            relation = 'related absolute length trait; normalized midface denominator prevents direct equivalence'
        elif row.feature == 'philtrum_width':
            term = "Cupid’s-bow upper-lip-peak separation; not necessarily philtral-column width"
        elif row.feature.startswith('eb_thickness'):
            term = 'Eyebrow polygon area; not clinician-defined hair density/thickness'
        elif row.feature in ['inter_canthal_distance', 'inter_pupillary_distance', 'outer_canthal_distance']:
            term = 'Ocular landmark spacing; clinical spacing/anthropometry requires separate reference and calibration'
        potential.append(dict(feature=row.feature, potential_term=term, potential_hpo_id=hpo, relation=relation,
                              formula=row.formula, units=row.units, anatomical_region=row.anatomical_region,
                              mapping_source='current pinned S6; HPO philtrum ID verified in Kline et al. Box 1',
                              exact_clinical_equivalence=False))
    potential = pd.DataFrame(potential)
    write_aggregate(potential, out / 'aggregate_candidates/P2_feature_to_potential_clinical_terms.csv')
    prior = pd.read_csv(ROOT / 'publication/final_research_120/clinical_literature_audit_17.csv')
    current = pd.read_csv(out / 'evidence/P0_exhaustive_evidence_matrix.csv', dtype={'canonical_syndrome_id': 'string'}, low_memory=False)
    prior = prior.merge(current[['syndrome', 'canonical_syndrome_id', 'feature', 'ff_hedges_g', 'ff_mannwhitney_q_fdr_global']],
                        on=['syndrome', 'feature'], validate='one_to_one')
    require((prior.hedges_g - prior.ff_hedges_g).abs().lt(1e-10).all()
            and (prior.global_q - prior.ff_mannwhitney_q_fdr_global).abs().lt(1e-10).all(), 'Prior clinical claims disagree with active output')
    rows = []
    for row in prior.itertuples():
        book = row.clinical_source.rstrip('/').split('/')[-1]
        title, authors, date, url, section, access, scope = SOURCES[book]
        classification = 'related_but_not_directly_comparable'
        if row.clinical_support in ['unsupported_exact_trait', 'scope_and_construct_mismatch', 'not_verified']:
            classification = 'insufficient_evidence'
        extra = ''
        if row.feature == 'philtrum_length' and row.syndrome.startswith('Coffin-Siris'):
            extra = ('A 2012 15-case report mentions elongation in an individual; a 2018 DPF2 cohort instead reports short/broad appearance. '
                     'These scopes are heterogeneous and neither validates this pooled normalized measurement.')
            access = 'indirectly_reported_indexed_full_text_case_passages'
        rows.append(dict(canonical_syndrome_id=row.canonical_syndrome_id, syndrome=row.syndrome, feature=row.feature,
                         primary_n=row.n, primary_g=row.hedges_g, primary_global_q=row.global_q, formula=row.formula,
                         clinical_comparison_class=classification, related_trait_support=row.clinical_support,
                         exact_comparable_measure=False, source_title=title, source_authors=authors, source_date=date,
                         source_url=url, source_section=section, access_state=access, accessed='2026-11-10',
                         source_scope=scope, numerator_denominator_caveat=row.interpretation, additional_scope_note=extra,
                         absence_is_novelty_evidence=False))
    sources = pd.DataFrame(rows)
    write_aggregate(sources, out / 'P2_phenotype_source_matrix.csv')
    candidates = pd.read_csv(out / 'P1_ranked_candidate_evidence.csv', dtype={'canonical_syndrome_id': 'string'}, low_memory=False)
    candidates = candidates[KEY + ['syndrome']].drop_duplicates()
    seed_candidates = out / 'aggregate_candidates/P1_seed_candidate_lists.csv'
    if seed_candidates.exists():
        seeded = pd.read_csv(seed_candidates, dtype={'canonical_syndrome_id': 'string'})
        candidates = pd.concat([candidates, seeded[KEY + ['syndrome']]], ignore_index=True).drop_duplicates(KEY)
    candidates = candidates.merge(sources[KEY + ['clinical_comparison_class', 'access_state', 'source_url', 'source_scope']], on=KEY, how='left', validate='one_to_one')
    candidates['clinical_comparison_class'] = candidates.clinical_comparison_class.fillna('unmapped')
    candidates['access_state'] = candidates.access_state.fillna('not_verified')
    candidates['source_scope'] = candidates.source_scope.fillna('No candidate-specific clinical source assessment; absence does not establish novelty')
    write_aggregate(candidates, out / 'aggregate_candidates/P2_ranked_candidate_clinical_status.csv')
    notes = {book: dict(title=value[0], authors=value[1], date=value[2], url=value[3], section=value[4], access_state=value[5], scope=value[6])
             for book, value in SOURCES.items()}
    notes['additional_Coffin_Siris_2012'] = {'url': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC3402612/',
        'title': 'The Coffin-Siris syndrome: A Proposed Diagnostic Approach and Assessment of 15 Overlapping Cases',
        'access_state': 'indexed_full_text_case_passage; live page inaccessible', 'scope': 'individual case appearance in a heterogeneous clinical series; not exact normalized geometry'}
    notes['additional_DPF2_2018'] = {'url': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC5985265/',
        'title': 'Mutations in the BAF-Complex Subunit DPF2 Are Associated with Coffin-Siris Syndrome',
        'access_state': 'indexed_full_text_passage', 'scope': 'DPF2-related subset; short/broad philtrum appearance cannot be generalized to pooled CSS'}
    notes['Porras_2021'] = {'url': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC8172204/',
        'title': 'Objective differential diagnosis of Noonan and Williams–Beuren syndromes in diverse populations using quantitative facial phenotyping',
        'authors': 'Porras AR; Summar M; Linguraru MG', 'date': '2021-03-27',
        'access_state': 'abstract_only_this_review; live full text inaccessible', 'scope': 'Noonan/Williams facial metrics across recorded populations; exact normalization/comparability not verified here'}
    write_json(notes, out / 'PUBLIC_SOURCE_NOTES.json')
    text = ('# P2 phenotype correspondence and claim audit\n\n'
        'All 17 earlier numerical claims match the exact active S1 rows and S6 formulas. The reviewed denominator is 17 selected questions, '
        'not all 955 primary associations. All 120 features receive a potential anatomical-term mapping, but this is not 120 validated HPO equivalences. '
        'Only the long-philtrum HPO identifier was verified explicitly in the retrieved consensus passage; other IDs are left blank rather than invented.\n\n'
        'No reviewed pair has been validated on an exactly comparable normalized clinical measurement. This does not imply the associations are false or the traits absent. '
        'The categories distinguish related appearance from anatomical equivalence; no apparent contradiction is asserted when ages, scale, landmarks or denominators differ. '
        'Candidate-specific source assessment is limited to the 17 prior questions; other ranked pairs are explicitly unmapped/not verified. '
        'There is no literature-match percentage over all 955 discoveries.\n\n'
        'Full-text source sections were retrieved for the accessible GeneReviews chapters and Cornelia de Lange consensus; some live PMC pages failed browser checks. '
        'For those, indexed passages or abstract-only access are labeled separately. Source title/author/date metadata and scope are recorded in PUBLIC_SOURCE_NOTES.json. '
        'No absence in one summary is treated as absence of a trait or evidence of novelty. No new clinical sign is declared.\n\n'
        'An important refinement of the earlier review is Coffin–Siris: targeted retrieval found a long/smooth description in an individual from a 2012 clinical series, '
        'and a short/broad description in a DPF2 subset. This adds evidence of heterogeneous reported appearance, not confirmation or contradiction of the pooled FaceKit ratio. '
        'The new audit does not overwrite the earlier frozen appendix.\n\n'
        'FaceKit already supplies the geometry platform and studies phenotype correspondence, reference conditioning and synthetic consistency. '
        'The current work contributes explicit comparator/normalization sensitivity and transparent exploratory evidence tables, not an extractor or an established new phenotype. '
        'See the frozen publication/final_research_120/FACEKIT_OVERLAP_AUDIT.md for precise preprint sections/figures; no methodological novelty has been established.\n')
    (out / 'P2_literature_claim_audit.md').write_text(text)
    draft = ('# P2 clinical comparison manuscript draft\n\n'
        '## Cornelia de Lange\n\n'
        'The consensus description supports an elongated/smooth philtrum, while FaceKit philtrum_length measures subnasale-to-upper-lip distance relative to midface height. '
        'The positive primary and within-GMDB length contrasts are directionally related to that description, but the exact ratio is not clinically calibrated. '
        'philtrum_width instead spans the upper-lip peaks; its negative primary effect cannot be interpreted as clinically narrow philtral columns. '
        '[Kline et al., 2018, Box 1 and Fig. 3](https://pmc.ncbi.nlm.nih.gov/articles/PMC7136165/).\n\n'
        '## Noonan\n\n'
        'The clinical description includes increased eye spacing and changes in appearance with age. The observed interpupillary/bizygomatic association is related to that appearance, '
        'but relative pupil spacing does not establish absolute hypertelorism or a diagnostic threshold. The within-GMDB and adjusted results answer different comparator questions. '
        '[Roberts, GeneReviews, revised December 4, 2025, Suggestive Findings](https://www.ncbi.nlm.nih.gov/sites/books/NBK1124/).\n\n'
        '## Williams–Beuren\n\n'
        'Periorbital fullness is a documented appearance, but neither outer-canthal/bizygomatic distance nor eyebrow polygon area is interchangeable with that trait. '
        'The same data give a much smaller or reversed outer-span effect under height normalization, which changes the estimand rather than proving the original index wrong. '
        'Eyebrow area is geometry, not hair density. Absolute anthropometric findings cannot establish a contradiction without matched landmarks, scale, population and age scope. '
        '[Morris, GeneReviews, updated April 13, 2023, distinctive facies](https://ncbi.nlm.nih.gov/sites/books/NBK1249/?report=classic). '
        '[Porras et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8172204/) addresses a related quantitative question; exact measure comparability was not verified from full text in this review.\n\n'
        'Numerical inserts and each correction family are supplied separately in P3_manuscript_insertions.md and the aggregate evidence tables. '
        'These paragraphs are drafts for clinical/PI review, not diagnostic guidance.\n')
    (out / 'P2_clinical_comparison_manuscript_draft.md').write_text(draft)
    print(f'Clinical source audit complete: 120 potential constructs, 17 assessed prior questions, {len(candidates)} unique ranked pairs classified')


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

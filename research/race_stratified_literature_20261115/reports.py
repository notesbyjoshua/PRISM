import argparse
import json
import re
from pathlib import Path

import pandas as pd

from directions import BASELINE, OUTPUT, ROOT, common, markdown_table

NOTES = {
    'ARBOLEDA-THAM': 'Lower normalized brow arch opposes the arched-brow descriptor. Sparse hair does not resolve a geometric arch difference; age, landmark placement and brow-length denominator remain unverified.',
    'CLARK-BARAITSER': 'Pointed chin is annotated in only 2/8. The angle is at menton: a larger angle is blunter. Thin upper lip alone cannot predict total vermilion height, so the total-lip feature is not classified as a directional contradiction.',
    'FLOATING-HARBOR': 'A lower mouth-corner-drop metric does not establish clinically upturned corners. The signed landmark baseline differs from the clinical definition; frequency is unspecified.',
    'Hyperphosphatasia': 'Short philtrum is occasional in a broad genetically heterogeneous syndrome group. A higher normalized length is a relative proxy mismatch, not exclusion of Mabry syndrome.',
    'INTELLECTUAL DEVELOPMENTAL DISORDER WITH DYSMORPHIC FACIES': 'Downslant is annotated in 5/20 FBXO11 cases. The PRISM mean angle is positive; this is worth checking against photographs and age-matched controls, but no race-specific expectation was established.',
    'INTELLECTUAL DEVELOPMENTAL DISORDER WITH PERSISTENCE': 'The BCL11A downslant annotation is 4/8 and facial findings vary. The positive measured angle and positive relative shift merit image/landmark review; the metric is not clinically calibrated.',
    'INTELLECTUAL DEVELOPMENTAL DISORDER, AUTOSOMAL DOMINANT 21': 'Pointed chin is annotated in 1/46, which does not establish a typical phenotype. Larger menton angle means a blunter angle. The source paper describes 39 new individuals; the catalog denominator must not be mistaken for that new-case count.',
    'MOWAT-WILSON': 'Upward and downward slant are both reported: two 2023 cases explicitly have upward slant. Upper-lid/iris distance remains a gaze-dependent proxy for ptosis, with no quantified clinical threshold.',
    'NICOLAIDES-BARAITSER': 'Short fissures and narrow bridge are catalog findings; narrow bridge includes a 1/3 annotation. Fissure/bridge ratios use face width as denominator and do not measure absolute clinical length or breadth.',
    'NEURODEVELOPMENTAL DISORDER WITH DYSMORPHIC FACIES AND DISTAL LIMB ANOMALIES': 'A BPTF expansion-cohort case has synophrys. Higher PRISM medial eyebrow spacing is a possible proxy mismatch with that case, not a contradiction of a typical or universal phenotype. Hair continuity and landmark spacing are different constructs.',
    'SMITH-MAGENIS': 'GeneReviews describes close-spaced eyes, whereas the catalogs also list hypertelorism. Lower PRISM spacing is consistent with the former; this is source heterogeneity, not a confirmed contradiction.',
    'SOTOS': 'The narrow-fissure annotation is a single 1/1 report and is not established as typical. Tall chin is a recognized descriptor, but normalized chin height, age and facial width require review.',
    'WIEDEMANN-STEINER': 'Synophrys is supported by case series, but medial landmark spacing does not measure hair continuity. A higher distance is only a proxy mismatch.',
    'WILLIAMS-BEUREN': 'Occasional hypotelorism is not a universal expectation: a Chinese Han cohort reports hypertelorism in 24/29. Nose-length ratio still opposes the recognized short-nose descriptor; compare absolute geometry, age and face-width denominator before calling contradiction.',
    'WITTEVEEN-KOLK': 'Facial asymmetry is recorded without frequency. A lower 2D mirror-landmark discrepancy is not evidence that no patient has clinical asymmetry; acquisition/pose and deletion scope need review.',
    'ZTTK': 'Narrow forehead is annotated in 2/7 and as occasional. GeneReviews forehead prominence describes depth, not breadth; it cannot settle the transverse-width mismatch.',
}


def note(disease):
    return next((value for key, value in NOTES.items() if disease.startswith(key)), '')


def links(group, registry):
    found = []
    for row in group.itertuples():
        if row.source.startswith('Targeted_narrative_'):
            key = row.reference
            title = registry.loc[key, 'title'] if key in registry.index else key
            found.append(f'[{title}]({row.source_url})')
        else:
            for reference in str(row.reference).split(';'):
                for pmid in re.findall(r'PMID:(\d+)', reference):
                    found.append(f'[PMID {pmid}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/) (catalog citation; full text not necessarily retrieved)')
            if row.source.startswith('Orpha'):
                found.append('[Orphadata phenotype annotations](https://www.orphadata.com/data/xml/en_product4.xml)')
    return '; '.join(sorted(set(found))) or '[HPO phenotype annotations](https://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa)'


def build(output):
    associations = pd.read_csv(output / 'all_982_association_literature_audit.csv', float_precision='round_trip').fillna('')
    evidence = pd.read_csv(output / 'matched_annotation_evidence.csv').fillna('')
    registry = pd.read_csv(output / 'targeted_narrative_source_registry.csv').set_index('source_id').fillna('')
    diseases = pd.read_csv(output / 'disease_identity_and_source_coverage.csv').fillna('')
    summary = json.loads((output / 'annotation_screen_summary.json').read_text())
    associations['clinical_review_note'] = associations.disease.map(note)
    keys = ['disease', 'canonical_race', 'feature']
    associations['reviewed_source_links'] = [links(evidence.loc[(evidence[keys] == row[keys]).all(axis=1)], registry)
                                            if row.matched_evidence_sources else '' for _, row in associations.iterrows()]
    associations.to_csv(output / 'all_982_association_literature_audit.csv', index=False)
    missing = associations.loc[associations.literature_status.eq('not_identified_in_reviewed_annotation_catalogs')]
    mismatch = associations.loc[associations.literature_status.eq('potential_relative_direction_mismatch')]
    mixed = associations.loc[associations.literature_status.eq('heterogeneous_reported_trait_directions')]
    initial = associations.loc[associations.catalog_only_status.eq('potential_relative_direction_mismatch')]
    initial.to_csv(output / 'catalog_direction_flags_with_targeted_review.csv', index=False)
    candidates = associations.loc[associations.literature_status.isin(['not_identified_in_reviewed_annotation_catalogs', 'potential_relative_direction_mismatch'])]
    candidates.to_csv(output / 'not_identified_or_potentially_opposed_pairs.csv', index=False)
    for status, group in associations.groupby('literature_status'):
        group.to_csv(output / (status + '.csv'), index=False)
    lines = ['# Main-analysis comparison with clinical literature', '',
             f'Review date: 2026-11-15. This inventory covers all **982 significant syndrome/category/feature results** '
             f'(original global q<.05) across **62 syndromes, 111 features and 891 distinct syndrome–feature comparisons**. '
             'All 9,600 planned primary tests remain unchanged. Nonsignificant results are not treated as discoveries needing a novelty claim.', '',
             f'**{len(missing)} results have no matching trait in the reviewed source annotations; {len(mismatch)} have a '
             f'potential opposing relative direction; {len(mixed)} encounter reported traits in both directions.** '
             'No exact FaceKit-measurement contradiction was confirmed. This does not establish that all results agree with literature.', '',
             '**Scope of “not identified”:** an exhaustive row inventory and a structured screen of the versioned HPO/Orphadata '
             'corpus, supplemented by targeted GeneReviews/primary-paper checks. It is not an exhaustive search of every publication. '
             'An absent annotation does not mean a feature is absent clinically, unprecedented, or novel. '
             'The source registry identifies which syndromes received additional narrative checks; the others have catalog screening only. '
             'The complete CSV is the authoritative filterable inventory.', '',
             '## Sources and interpretation', '',
             '[HPO ontology v2026-09-01](https://github.com/obophenotype/human-phenotype-ontology/releases/tag/v2026-09-01), '
             '[HPO annotations dated 2026-09-02](https://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa), and '
             '[Orphadata phenotypes dated 2026-06-23](https://www.orphadata.com/data/xml/en_product4.xml) were downloaded, '
             'versioned and hashed. Orphadata nomenclature supplied exact names/aliases. HPO and Orphadata can reuse the same '
             'underlying cases and are not independent replications. Targeted narrative annotations preserve their own source URLs.', '',
             'Source access varied: many NCBI direct pages returned reCAPTCHA, so targeted verification often used indexed clinical '
             'passages rather than an entire retrieved chapter. A PubMed identifier retained from HPO is a catalog citation, '
             'not a claim that its full text was read. No patient images or individual data were submitted to web services.', '',
             'Feature formulas were reviewed individually. Clinical terms were joined to explicitly curated geometry families, '
             'including HPO descendant terms; composite metrics were not assigned a direction from only one component. '
             'None of these mappings has validated exact clinical equivalence. Positive g means a higher mean than the '
             'same recorded-category FairFace reference; it does not by itself imply clinical hypertelorism, upslant, thinness, '
             'or an absolute abnormality. Frequencies refer to source cohorts and are not predictions for every patient.', '',
             'The African/Asian/European labels are recorded dataset groups, not verified genetic ancestry. The Asian category '
             'includes harmonized subgroups. Most sources are not category-specific. Patient age, sex, genotype, source images, '
             'expression, gaze, landmarks and normalization denominators can affect apparent differences.', '',
             '## Classification counts', '', markdown_table(pd.DataFrame(summary['status_counts'].items(), columns=['Classification', 'Results'])), '',
             '“Measurement construct not directionally mapped” means the metric has no defensible one-trait clinical sign mapping. '
             '“Related trait not directionally comparable” means a related clinical annotation was found but cannot settle the metric. '
             '“Directionally related support” is proxy agreement, not validation. Unresolved broad diagnoses and missing parent-level '
             'coverage are listed separately.', '',
             '## Every potential relative-direction mismatch', '',
             'These are review candidates, not confirmed biological contradictions or proven errors. Each is significant in the original '
             'main analysis. The complete CSV includes raw case/reference means, sample counts, exact formulas and all original p/q values.', '']
    for row in mismatch.sort_values(keys).itertuples():
        lines += [f'### {row.disease.split("__")[0]} / {row.canonical_race} / `{row.feature}`', '',
                  f'g={row.hedges_g:+.4f}; q={row.q:.5g}; cases={row.case_n}; references={row.reference_n}. '
                  f'Case mean={row.case_mean:.7g}; reference mean={row.reference_mean:.7g}.', '',
                  f'Reported trait: **{row.reported_directional_traits}**. Source frequency: {row.reported_trait_frequencies or "unspecified"}. '
                  f'{row.comparability_note}', '', row.clinical_review_note or 'No additional syndrome-specific narrative resolution established.', '',
                  f'Evidence: {row.reviewed_source_links}', '']
    lines += ['## Flags changed by targeted source review', '',
              'Smith–Magenis eye spacing illustrates catalog versus GeneReviews disagreement. Mowat–Wilson slant has upward-slant case '
              'reports as well as downward-slant annotations. Williams–Beuren spacing has hypertelorism in a primary cohort despite '
              'occasional hypotelorism in HPO. These are retained as heterogeneous evidence. '
              'Total vermilion height was never assigned a sign from thin upper lip alone. See catalog_direction_flags_with_targeted_review.csv '
              'for every catalog direction flag and its final classification. The targeted BPTF case-table review also adds a '
              'synophrys-related proxy flag that was absent from the catalogs; it appears in the complete mismatch list above.', '',
              markdown_table(mixed[['disease', 'canonical_race', 'feature', 'hedges_g', 'q', 'reported_directional_traits']]), '',
              '## Exhaustive list without a matching source annotation', '',
              'Each line below is one syndrome/recorded-category/feature result. “Not identified” is bounded to the reviewed annotations. '
              'Use the feature mapping CSV for the exact clinical concepts searched and the source coverage table below to assess gaps.', '']
    columns = ['canonical_race', 'feature', 'hedges_g', 'q', 'case_n', 'feature_family']
    for disease, group in missing.sort_values(keys).groupby('disease', sort=True):
        lines += [f'### {disease.split("__")[0]}', '', markdown_table(group[columns]), '']
    lines += ['## Complete coverage for every significant syndrome', '',
              markdown_table(diseases[['display_name', 'disease_ids', 'annotation_rows', 'scope_unresolved', 'targeted_narrative_sources']]), '',
              'Distal arthrogryposis, cutis laxa, mucopolysaccharidoses and the ambiguous autism-susceptibility label were not '
              'substituted with one molecular subtype. Broad Joubert has no direct catalog annotations at the selected parent ID; '
              'a GeneReviews illustration supplies related context but does not resolve its two significant metrics.', '',
              '## Targeted narrative source registry', '']
    for source, row in registry.iterrows():
        lines += [f'- **{source}:** [{row.title}]({row.url}). Version/date: {row.source_version}. Access: '
                  f'{row.access_and_section}. {row.review_note}']
    lines += ['', '## Files and provenance', '',
              '- `all_982_association_literature_audit.csv`: every significant result, classification, notes and source links.',
              '- `not_identified_or_potentially_opposed_pairs.csv`: the requested candidate subset, with distinct uncertainty classes.',
              '- Per-classification CSVs retain measurement, scope and heterogeneous cases separately.',
              '- `matched_annotation_evidence.csv`: term-level matches, ontology mappings, negation, frequency, catalog references and URLs.',
              '- `feature_to_clinical_trait_mapping.csv`: all 120 formulas and their curated directional/related concepts.',
              '- `DIRECTION_REVERSALS_REPORT.md`: separate cross-category sign-change report.', '',
              f'Primary SHA-256: `{summary["primary_sha256"]}`. Results stay in ignored `local_review/` storage. '
              'No main-analysis statistics, demographic analyses, sensitivity analyses, patient data or manuscript claims were changed.']
    (output / 'LITERATURE_COMPARISON_REPORT.md').write_text('\n'.join(lines) + '\n')
    manifest = {'review_date': '2026-11-15', 'primary_sha256': common.sha(BASELINE / 'primary_all_tests.csv'),
                'source_downloads': json.loads((output / 'public_sources/download_manifest.json').read_text()),
                'code_sha256': {str(p.relative_to(ROOT)): common.sha(p) for p in Path(__file__).parent.glob('*.py')},
                'result_sha256': {p.name: common.sha(p) for p in output.iterdir() if p.suffix in {'.csv', '.md', '.json'} and p.name != 'review_manifest.json'},
                'clinical_absence_or_novelty_proven': False, 'statistics_modified': False}
    common.write_json(output / 'review_manifest.json', manifest)
    print(json.dumps({'missing_source_annotation': len(missing), 'potential_relative_mismatch': len(mismatch),
                      'heterogeneous_directions': len(mixed), 'candidate_subset': len(candidates)}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, default=OUTPUT)
    build(common.private_output(parser.parse_args().output_dir))

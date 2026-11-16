import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from curation import disease_mapping, rules
from directions import BASELINE, OUTPUT, common, load_primary, markdown_table
from sources import annotation_tables, nomenclature, ontology
from narrative import source_registry, supplementary_annotations


def feature_mapping(definitions, terms):
    rows = []
    for definition in definitions.itertuples():
        matches = [rule for rule in rules() if any(re.search(pattern, definition.feature) for pattern in rule['patterns'])]
        if len(matches) != 1:
            raise ValueError(f'Feature requires exactly one reviewed mapping: {definition.feature} ({len(matches)})')
        rule = matches[0]
        for key in list(rule['expected']) + list(rule['related']):
            if key not in terms or terms[key].get('obsolete'):
                raise ValueError('Invalid/obsolete curated HPO term: ' + key)
        rows.append({'feature': definition.feature, 'anatomical_region': definition.anatomical_region,
                     'formula': definition.formula, 'units': definition.units, 'family': rule['family'],
                     'relation': rule['relation'], 'clinical_comparability_note': rule['note'],
                     'expected_trait_directions': json.dumps(rule['expected'], sort_keys=True),
                     'expected_trait_labels': '|'.join(terms[key]['name'] for key in rule['expected']),
                     'related_nondirectional_terms': '|'.join(rule['related']),
                     'related_nondirectional_labels': '|'.join(terms[key]['name'] for key in rule['related']),
                     'validated_exact_clinical_equivalence': False})
    return pd.DataFrame(rows)


def association_review(primary, mappings, diseases, annotations, terms, ancestors):
    rule_map = mappings.set_index('feature').to_dict('index')
    disease_map = diseases.set_index('disease').to_dict('index')
    records_map = {key: annotations.loc[annotations.disease_id.isin(value['disease_ids'].split('|'))]
                   for key, value in disease_map.items()}
    evidence_rows, association_rows = [], []
    for row in primary.loc[primary.supported & primary.q.lt(.05)].itertuples():
        mapping, disease = rule_map[row.feature], disease_map[row.disease]
        expected = json.loads(mapping['expected_trait_directions'])
        related = set(filter(None, mapping['related_nondirectional_terms'].split('|')))
        records = records_map[row.disease]
        matched_expected, related_records, excluded = [], [], []
        for item in records.itertuples():
            term_ancestors = ancestors(item.term_id)
            match = [key for key in expected if key in term_ancestors]
            match_related = [key for key in related if key in term_ancestors]
            if not match and not match_related:
                continue
            if item.negated:
                excluded.append(item)
            elif match:
                matched_expected.extend((item, key, int(expected[key])) for key in match)
            else:
                related_records.append(item)
            for key in match + match_related:
                evidence_rows.append({'disease': row.disease, 'canonical_race': row.canonical_race, 'feature': row.feature,
                                      'disease_id': item.disease_id, 'annotation_disease_name': item.disease_name,
                                      'annotated_term_id': item.term_id,
                                      'annotated_term': terms.get(item.term_id, {}).get('name', item.term_id),
                                      'mapped_trait_id': key, 'mapped_trait': terms[key]['name'],
                                      'expected_relative_direction': expected.get(key, 0), 'negated': item.negated,
                                      'frequency': item.frequency, 'onset': item.onset, 'sex': item.sex,
                                      'frequency_label': terms.get(item.frequency, {}).get('name', item.frequency),
                                      'reference': item.reference, 'source': item.source,
                                      'source_url': getattr(item, 'source_url', '') or ('https://www.orphadata.com/data/xml/en_product4.xml' if item.source.startswith('Orpha')
                                      else 'https://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa'),
                                      'annotation_not_race_specific': True})
        directions = set(direction for _, _, direction in matched_expected)
        sign = int(np.sign(row.hedges_g))
        status = 'not_identified_in_reviewed_annotation_catalogs'
        if disease['scope_unresolved']:
            status = 'diagnosis_scope_unresolved'
        elif not len(records):
            status = 'no_exact_scope_annotation_coverage'
        elif len(directions) > 1:
            status = 'heterogeneous_reported_trait_directions'
        elif directions and sign in directions:
            status = 'directionally_related_support'
        elif directions:
            status = 'potential_relative_direction_mismatch'
        elif related_records:
            status = 'related_trait_not_directionally_comparable'
        elif excluded:
            status = 'only_negated_related_annotation'
        elif not expected:
            status = 'measurement_construct_not_directionally_mapped'
        absolute_slant = ''
        if mapping['family'] == 'palpebral_slant' and len(directions) == 1:
            expected_direction = next(iter(directions))
            absolute_slant = 'mean_angle_opposes_reported_slant' if np.sign(row.case_mean) == -expected_direction else 'mean_angle_not_opposite_or_zero'
        association_rows.append({**row._asdict(), 'literature_status': status,
                                 'recorded_group': row.canonical_race, 'feature_family': mapping['family'],
                                 'clinical_comparability': mapping['relation'],
                                 'comparability_note': mapping['clinical_comparability_note'],
                                 'mapped_disease_ids': disease['disease_ids'], 'diagnosis_scope_note': disease['scope_note'],
                                 'annotation_coverage_rows': len(records),
                                 'reported_directional_traits': '|'.join(sorted({terms[key]['name'] for _, key, _ in matched_expected})),
                                 'expected_relative_signs': '|'.join(map(str, sorted(directions))),
                                 'related_nondirectional_traits': '|'.join(sorted({terms.get(item.term_id, {}).get('name', item.term_id) for item in related_records})),
                                 'reported_trait_frequencies': '|'.join(sorted({terms.get(item.frequency, {}).get('name', str(item.frequency)) for item, _, _ in matched_expected if str(item.frequency)})),
                                 'matched_evidence_sources': '|'.join(sorted({item.source for item, _, _ in matched_expected} | {item.source for item in related_records})),
                                 'absolute_slant_check': absolute_slant,
                                 'confirmed_literature_contradiction': False,
                                 'exact_FaceKit_measurement_verified_in_literature': False,
                                 'novel_phenotype_claim': False,
                                 'source_scope': 'All annotations for exact curated disease IDs; not race-specific; broad Asian harmonization retained'})
    return pd.DataFrame(association_rows).drop(columns='Index'), pd.DataFrame(evidence_rows)


def build(output):
    primary, manifest = load_primary()
    sources = output / 'public_sources'
    for name, record in json.loads((sources / 'download_manifest.json').read_text()).items():
        if common.sha(sources / name) != record['sha256']:
            raise ValueError('Public source checksum differs: ' + name)
    terms, names, ancestors = ontology(sources / 'hp.obo')
    annotations = annotation_tables(sources)
    definitions = pd.read_csv(BASELINE / 'exact_120_feature_definitions.csv')
    mappings = feature_mapping(definitions, terms)
    significant = primary.loc[primary.supported & primary.q.lt(.05)]
    diseases = disease_mapping(significant.disease.unique(), annotations,
                               nomenclature(sources / 'orphanet_nomenclature.xml'))
    catalog, _ = association_review(primary, mappings, diseases, annotations, terms, ancestors)
    supplements = supplementary_annotations()
    if not supplements.term_id.isin(terms).all():
        raise ValueError('Invalid narrative HPO mapping')
    annotations = pd.concat([annotations, supplements], ignore_index=True).fillna('')
    registry = source_registry()
    diseases['targeted_narrative_sources'] = diseases.disease_ids.map(
        lambda ids: '|'.join(registry.loc[registry.disease_ids.map(lambda values: bool(set(values.split('|')) & set(ids.split('|')))), 'source_id']))
    mappings.to_csv(output / 'feature_to_clinical_trait_mapping.csv', index=False)
    diseases.to_csv(output / 'disease_identity_and_source_coverage.csv', index=False)
    associations, evidence = association_review(primary, mappings, diseases, annotations, terms, ancestors)
    associations['catalog_only_status'] = catalog.literature_status.to_numpy()
    associations = associations.merge(diseases[['disease', 'targeted_narrative_sources']], on='disease', validate='many_to_one')
    associations['targeted_narrative_review_performed'] = associations.targeted_narrative_sources.ne('')
    assert len(associations) == len(significant) == 982
    assert not associations[['disease', 'canonical_race', 'feature']].duplicated().any()
    associations.to_csv(output / 'all_982_association_literature_audit.csv', index=False)
    evidence.to_csv(output / 'matched_annotation_evidence.csv', index=False)
    registry.to_csv(output / 'targeted_narrative_source_registry.csv', index=False)
    supplements.to_csv(output / 'targeted_narrative_annotations.csv', index=False)
    for old in output.glob('*.csv'):
        if old.stem in set(catalog.literature_status) | set(associations.literature_status):
            old.unlink()
    statuses = associations.literature_status.value_counts().to_dict()
    for status in statuses:
        associations.loc[associations.literature_status.eq(status)].to_csv(output / f'{status}.csv', index=False)
    summary = {'primary_baseline': manifest['study_id'], 'primary_family': len(primary),
               'significant_associations_reviewed': len(associations), 'unique_syndrome_feature_comparisons': len(significant[['disease', 'feature']].drop_duplicates()),
               'syndromes': significant.disease.nunique(), 'distinct_features': significant.feature.nunique(),
               'status_counts': statuses, 'confirmed_exact_measurement_contradictions': 0,
               'scope': 'Exhaustive result-row inventory and structured HPO/Orphadata screen plus targeted narrative verification; not a systematic literature-wide absence search',
               'primary_sha256': common.sha(BASELINE / 'primary_all_tests.csv'),
               'HPO_annotation_version': '2026-09-02', 'HPO_ontology_version': '2026-09-01', 'Orphadata_version': '2026-06-23'}
    common.write_json(output / 'annotation_screen_summary.json', summary)
    print(json.dumps(summary, indent=2))
    return associations, evidence, mappings, diseases


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, default=OUTPUT)
    args = parser.parse_args()
    build(common.private_output(args.output_dir))

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

CODE = Path(__file__).resolve().parents[1] / 'research/race_stratified_literature_20261115'
sys.path.insert(0, str(CODE))
import audit
import curation
import directions
import narrative
import sources


def test_clinical_composites_do_not_inherit_single_component_direction():
    families = {rule['family']: rule for rule in curation.rules()}
    for name in ['total_vermilion', 'outer_canthal_span', 'philtrum_named_width', 'lower_lid_iris_distance', 'nasal_alar_height']:
        assert families[name]['expected'] == {}
    assert families['chin_angle']['expected'] == {'HP:0000307': -1}


def test_curated_disease_identities_reject_similar_wrong_diagnoses():
    alpha = curation.DIAGNOSIS_OVERRIDES['MANNOSIDOSIS, ALPHA B, LYSOSOMAL; MANSA'][0]
    assert 'OMIM:248500' in alpha and 'OMIM:248510' not in alpha
    assert 'OMIM:617991' in curation.DIAGNOSIS_OVERRIDES['CHUNG-JANSEN SYNDROME; CHUJANS'][0]
    assert 'ORPHA:528084' not in curation.DIAGNOSIS_OVERRIDES['FBXW7'][0]


def test_no_fuzzy_disease_join():
    annotations = pd.DataFrame([{'disease_id': 'OMIM:1', 'disease_name': 'Beta mannosidosis'}])
    mapped = curation.disease_mapping(['Alpha mannosidosis__x'], annotations, [])
    assert mapped.iloc[0].disease_ids == ''
    assert mapped.iloc[0].mapping_status == 'unresolved_identity'


def clinical_case(annotations, g=1, expected=None, scope=False):
    primary = pd.DataFrame([{'disease': 'test__x', 'canonical_race': 'Asian', 'feature': 'x',
                             'supported': True, 'q': .01, 'hedges_g': g, 'case_mean': 2.}])
    mappings = pd.DataFrame([{'feature': 'x', 'expected_trait_directions': json.dumps({'HP:1': 1} if expected is None else expected),
                              'related_nondirectional_terms': 'HP:3', 'family': 'test', 'relation': 'proxy',
                              'clinical_comparability_note': 'not exact'}])
    diseases = pd.DataFrame([{'disease': 'test__x', 'disease_ids': 'OMIM:1', 'scope_unresolved': scope, 'scope_note': ''}])
    terms = {f'HP:{i}': {'name': f'Trait {i}'} for i in [1, 2, 3]}
    records = pd.DataFrame([{'disease_id': 'OMIM:1', 'disease_name': 'test', 'term_id': term, 'negated': negated,
                             'frequency': '1/1', 'onset': '', 'sex': '', 'reference': 'PMID:1',
                             'source': 'HPO_annotations_test'} for term, negated in annotations],
                           columns=['disease_id', 'disease_name', 'term_id', 'negated', 'frequency', 'onset', 'sex', 'reference', 'source'])
    result, evidence = audit.association_review(primary, mappings, diseases, records, terms, lambda key: {key})
    return result.iloc[0], evidence


@pytest.mark.parametrize('g,status', [(1, 'directionally_related_support'), (-1, 'potential_relative_direction_mismatch')])
def test_matching_trait_and_relative_direction(g, status):
    result, evidence = clinical_case([('HP:1', False)], g)
    assert result.literature_status == status
    assert not result.confirmed_literature_contradiction
    assert len(evidence) == 1


def test_both_reported_directions_override_single_mismatch():
    result, _ = clinical_case([('HP:1', False), ('HP:2', False)], expected={'HP:1': 1, 'HP:2': -1})
    assert result.literature_status == 'heterogeneous_reported_trait_directions'


def test_negation_is_not_support_or_missing_evidence():
    result, evidence = clinical_case([('HP:1', True)])
    assert result.literature_status == 'only_negated_related_annotation'
    assert evidence.negated.all()


def test_unmapped_measurement_and_unresolved_diagnosis_are_separate():
    result, _ = clinical_case([('HP:2', False)], expected={})
    assert result.literature_status == 'measurement_construct_not_directionally_mapped'
    result, _ = clinical_case([('HP:1', False)], scope=True)
    assert result.literature_status == 'diagnosis_scope_unresolved'


def test_no_scope_coverage_does_not_become_a_missing_trait_claim():
    result, _ = clinical_case([])
    assert result.literature_status == 'no_exact_scope_annotation_coverage'


def test_ontology_ancestors_and_alt_ids(tmp_path):
    path = tmp_path / 'hp.obo'
    path.write_text('[Term]\nid: HP:1\nname: Parent\n\n[Term]\nid: HP:2\nname: Child\nalt_id: HP:9\nis_a: HP:1 ! Parent\n')
    _, _, ancestors = sources.ontology(path)
    assert ancestors('HP:9') == {'HP:1', 'HP:2'}
    assert 'HP:2' not in ancestors('HP:1')


def direction_frame():
    rows = []
    for race, g, q in [('African', 1., .01), ('Asian', -1., .01), ('European', -1e-14, .6)]:
        rows.append({'disease': 'same__id', 'feature': 'x', 'canonical_race': race, 'hedges_g': g,
                     'q': q, 'p': q / 10, 'supported': True, 'case_n': 5, 'reference_n': 50,
                     'case_mean': 1. + g, 'reference_mean': 1., 'anatomical_region': 'nose', 'formula': 'x', 'units': 'ratio'})
    return pd.DataFrame(rows)


def test_sign_reversal_tolerance_and_significance_tiers():
    result = directions.comparisons(direction_frame())
    assert len(result) == 3
    flips = result.loc[result.opposite_g_sign]
    assert len(flips) == 1 and flips.iloc[0].significance_tier == 'both_significant'
    assert flips.iloc[0].case_minus_reference_mean_1 == 1.
    assert flips.iloc[0].case_minus_reference_mean_2 == -1.


def test_one_significant_is_not_reported_as_both():
    frame = direction_frame()
    frame.loc[frame.canonical_race.eq('Asian'), 'q'] = .051
    flips = directions.comparisons(frame).query('opposite_g_sign')
    assert flips.iloc[0].significance_tier == 'one_significant'


def test_ineligible_and_different_diagnoses_never_join():
    frame = direction_frame()
    frame.loc[frame.canonical_race.eq('Asian'), 'supported'] = False
    frame.loc[frame.canonical_race.eq('European'), 'disease'] = 'different__id'
    assert directions.comparisons(frame).empty


def test_narrative_source_provenance_is_explicit():
    registry = narrative.source_registry()
    assert registry.source_id.is_unique and registry.url.str.startswith('https://').all()
    assert not registry.literature_absence_can_be_inferred.any()
    added = narrative.supplementary_annotations()
    assert not added.negated.any()
    assert added.source_url.str.startswith('https://').all()


def test_private_real_result_inventory_when_present():
    output = directions.OUTPUT
    if not (output / 'all_982_association_literature_audit.csv').exists():
        pytest.skip('Private real outputs are not distributed')
    main, _ = directions.load_primary()
    reviewed = pd.read_csv(output / 'all_982_association_literature_audit.csv', float_precision='round_trip')
    keys = ['disease', 'canonical_race', 'feature']
    expected = main.loc[main.supported & main.q.lt(.05)]
    assert len(reviewed) == len(expected) == 982 and not reviewed[keys].duplicated().any()
    merged = reviewed.merge(expected, on=keys, validate='one_to_one', suffixes=('_review', '_main'))
    for field in ['hedges_g', 'p', 'q', 'case_mean', 'reference_mean']:
        assert np.array_equal(merged[field + '_review'], merged[field + '_main'])
    assert not reviewed.confirmed_literature_contradiction.any()
    assert not reviewed.novel_phenotype_claim.any()

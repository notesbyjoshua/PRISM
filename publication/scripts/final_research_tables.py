import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'publication/supplementary'
OUT = ROOT / 'publication/final_research_120'
FEATURES = ['philtrum_length', 'inter_pupillary_distance', 'forehead_height']
LABELS = {'eb_lateral_thickness_mean': 'Bilateral lateral eyebrow thickness index',
          'forehead_height': 'Forehead height / bizygomatic width', 'chin_height': 'Chin height / bizygomatic width',
          'eb_lateral_thickness_l': 'Left lateral eyebrow thickness index',
          'eb_lateral_thickness_r': 'Right lateral eyebrow thickness index'}
APPENDIX = [
    ('WILLIAMS-BEUREN', 'eb_thickness_mean', 'NBK1249', 'Clinical Description / facial features', 'related_construct_only', 'Periorbital fullness does not establish eyebrow polygon area or hair density'),
    ('Cornelia de Lange', 'philtrum_length', 'NBK1104', 'Suggestive Findings / craniofacial appearance', 'direct_related_trait', 'Long/smooth philtrum supports length-related interpretation; midface denominator not clinically validated'),
    ('WILLIAMS-BEUREN', 'outer_canthal_distance', 'NBK1249', 'Clinical Description / facial features', 'related_construct_only', 'Periorbital fullness does not establish increased outer-canthal/bizygomatic ratio'),
    ('Cornelia de Lange', 'philtrum_width', 'NBK1104', 'Suggestive Findings / craniofacial appearance', 'unsupported_exact_trait', 'Long/smooth philtrum does not establish philtrum width'),
    ('Noonan syndrome__', 'inter_pupillary_distance', 'NBK1124', 'Suggestive Findings / characteristic facies', 'direct_related_trait', 'Widely spaced eyes supports spacing-related interpretation; normalized pupil ratio is a distinct measurement'),
    ('Mucopolysaccharidoses', 'chin_pointedness_angle', 'NBK1274', 'Clinical Description / facial features', 'scope_and_construct_mismatch', 'An MPS-II description cannot validate a pooled Mucopolysaccharidoses cohort or this jaw angle'),
    ('Noonan syndrome__', 'nose_bridge_width', 'NBK1124', 'Clinical Description / age-dependent facies', 'related_construct_only', 'Depressed bridge and wide nasal base are not lateral bridge width; age-dependent description'),
    ('KBG SYNDROME', 'jaw_bigonial_width', 'NBK487886', 'Summary / characteristic facial features', 'related_construct_only', 'Triangular face does not establish a particular bigonial/bizygomatic ratio'),
    ('WIEDEMANN-STEINER', 'inter_canthal_distance', 'NBK580718', 'Clinical Description / characteristic dysmorphic features', 'direct_related_trait', 'Widely spaced eyes supports spacing-related interpretation; inner-canthal normalization remains distinct'),
    ('Kabuki syndrome', 'outer_canthal_distance', 'NBK62111', 'Summary / clinical characteristics', 'related_construct_only', 'Long palpebral fissures do not directly establish total outer-canthal separation'),
    ('Holoprosencephaly', 'mouth_opening', 'NBK1530', 'Clinical Characteristics / craniofacial spectrum', 'unsupported_exact_trait', 'Structural midline anomalies do not validate expression-dependent mouth opening in a photograph'),
    ('Coffin-Siris syndrome', 'philtrum_length', 'NBK131811', 'Clinical Description / facial features', 'not_verified', 'Retrieved clinical source does not verify the exact long-philtrum claim; retain as unresolved clinical interpretation'),
    ('KBG SYNDROME', 'jaw_to_forehead_ratio', 'NBK487886', 'Summary / characteristic facial features', 'related_construct_only', 'Triangular face is related to outline, not proof of the exact jaw/forehead ratio'),
    ('SOTOS SYNDROME;', 'forehead_width_mid', 'NBK1479', 'Suggestive Findings / characteristic facial appearance', 'direct_related_trait', 'Broad prominent forehead supports a width-related interpretation; facial-width denominator remains distinct'),
    ('Kabuki syndrome', 'inter_canthal_distance', 'NBK62111', 'Summary / clinical characteristics', 'related_construct_only', 'Long palpebral fissures and eyelid eversion do not establish inner-canthal distance'),
    ('MANNOSIDOSIS, ALPHA', 'inter_pupillary_distance', 'NBK1396', 'Clinical Description / facial features', 'unsupported_exact_trait', 'Coarse facies and widely spaced teeth do not establish pupil spacing'),
    ('Rubinstein-Taybi', 'malar_bulge_mean', 'NBK1526', 'Clinical Description / craniofacial features', 'related_construct_only', 'Facial gestalt and grimacing smile do not clinically validate the projected cheek/oval index')]


def save(frame, name):
    frame.to_csv(OUT / name, index=False, float_format='%.17g', lineterminator='\n')


def load(name):
    return pd.read_csv(SOURCE / name)


def recurrent_and_regions(primary, definitions):
    rows = []
    for feature, group in primary.groupby('feature', sort=True):
        significant = group.loc[group.mannwhitney_q_fdr_global.lt(.05)]
        rows.append({'feature': feature, 'feature_name': LABELS.get(feature, feature),
                     'significant_syndromes': significant.disease.nunique(), 'tested_syndromes': group.disease.nunique(),
                     'positive': int(significant.hedges_g.gt(0).sum()), 'negative': int(significant.hedges_g.lt(0).sum()),
                     'significance': 'q < 0.05; global BH across all 5760 primary hypotheses'})
    recurrent = pd.DataFrame(rows).sort_values(['significant_syndromes', 'feature'], ascending=[False, True])
    recurrent = recurrent.merge(definitions[['feature', 'anatomical_region', 'formula', 'units']], on='feature', validate='one_to_one')
    merged = primary.merge(definitions[['feature', 'anatomical_region']], on='feature', validate='many_to_one')
    regions = merged.assign(significant=merged.mannwhitney_q_fdr_global.lt(.05)).groupby('anatomical_region', sort=True).agg(
        tested=('feature', 'size'), significant=('significant', 'sum'), features=('feature', 'nunique')).reset_index()
    regions['fraction_significant'] = regions.significant / regions.tested
    pd.testing.assert_frame_equal(regions, load('regional_summary_120.csv'), check_dtype=False, check_exact=False, rtol=1e-12)
    for region, significant, tested in [('chin / jaw', 87, 240), ('forehead', 60, 192), ('eyebrows', 197, 960), ('cheeks / malar region', 35, 384)]:
        row = regions.set_index('anatomical_region').loc[region]
        assert row.significant == significant and row.tested == tested
    return recurrent.head(5), regions


def demographic_examples():
    examples, group_details, source_files = [], [], []
    for layer in ['layer1_broad_ethnicity', 'layer2_subcategory']:
        a_name = f'{layer}_analysis_A_reference_kruskal_all_groups.csv'
        b_name = f'{layer}_analysis_B_matched_disease_vs_reference.csv'
        c_name = f'{layer}_analysis_C_all_interactions.csv'
        a, b, c = [load(name) for name in [a_name, b_name, c_name]]
        source_files += [a_name, b_name, c_name]
        for feature in FEATURES:
            row = a.set_index('feature').loc[feature]
            sizes, means, medians = [{pair.split(':')[0]: float(pair.split(':')[1]) for pair in row[k].split('|')}
                                     for k in ['group_sample_sizes', 'group_means', 'group_medians']]
            examples.append(dict(layer=layer, analysis='A', syndrome='reference only', feature=feature,
                                 group_1='all supported recorded groups', group_2='', family_hypotheses=len(a),
                                 q=row.fdr_q_value, omnibus_H=row.kruskal_H, uncertainty='Omnibus test; no effect CI in A output',
                                 source=a_name))
            for group in sorted(sizes):
                group_details.append(dict(layer=layer, feature=feature, group=group, n=sizes[group], mean=means[group], median=medians[group]))
        for analysis, frame, source in [('B', b, b_name), ('C', c, c_name)]:
            frame = frame.copy()
            if analysis == 'B':
                frame['minimum_disease_n'], frame['minimum_reference_n'] = frame.n_disease, frame.n_reference_matched
                sort = ['minimum_disease_n', 'minimum_reference_n', 'syndrome_name', 'group']
                cells = ['syndrome_name', 'group']
            else:
                frame['minimum_disease_n'] = frame[['n_disease_group_1', 'n_disease_group_2']].min(axis=1)
                frame['minimum_reference_n'] = frame[['n_reference_group_1', 'n_reference_group_2']].min(axis=1)
                sort = ['minimum_disease_n', 'minimum_reference_n', 'syndrome_name', 'group_1', 'group_2']
                cells = ['syndrome_name', 'group_1', 'group_2']
            supported = frame.loc[frame.minimum_disease_n.ge(10) & frame.minimum_reference_n.ge(100)]
            selected = supported.sort_values(sort, ascending=[False, False] + [True] * (len(sort) - 2))
            selected = selected.drop_duplicates(cells).drop_duplicates('syndrome_name').head(3)
            for feature, (_, cell) in zip(FEATURES, selected.iterrows()):
                subset = frame.loc[frame.feature.eq(feature)]
                for column in cells:
                    subset = subset.loc[subset[column].eq(cell[column])]
                assert len(subset) == 1
                row = subset.iloc[0]
                item = dict(layer=layer, analysis=analysis, syndrome=row.syndrome_name, feature=feature,
                            family_hypotheses=len(frame), q=row.global_fdr_q_value, source=source)
                if analysis == 'B':
                    item.update(group_1=row.group, group_2='', n_disease_1=row.n_disease, n_reference_1=row.n_reference_matched,
                                effect=row.hedges_g, ci_lower=row.hedges_g_ci_lower, ci_upper=row.hedges_g_ci_upper,
                                disease_mean_1=row.disease_mean, reference_mean_1=row.reference_mean,
                                mean_difference=row.mean_difference, effect_type='Hedges g; syndrome minus matched reference',
                                uncertainty='1000-draw marginal bootstrap CI; not multiplicity-adjusted')
                else:
                    item.update(group_1=row.group_1, group_2=row.group_2, n_disease_1=row.n_disease_group_1,
                                n_disease_2=row.n_disease_group_2, n_reference_1=row.n_reference_group_1,
                                n_reference_2=row.n_reference_group_2, effect=row.interaction_beta_group2_minus_group1,
                                ci_lower=row.interaction_ci_lower, ci_upper=row.interaction_ci_upper,
                                disease_mean_1=row.disease_mean_group_1, disease_mean_2=row.disease_mean_group_2,
                                reference_mean_1=row.reference_mean_group_1, reference_mean_2=row.reference_mean_group_2,
                                effect_type='Interaction beta: (disease-reference) group2 minus group1',
                                uncertainty='95% HC3 interval in feature units; sparse-cell exploratory model')
                examples.append(item)
    return pd.DataFrame(examples), pd.DataFrame(group_details), source_files


def adjusted_examples(primary):
    rows = []
    for row in load('within_gmdb_age_sex_adjusted.csv').itertuples(index=False):
        rows.append(dict(syndrome=row.syndrome, feature=row.feature, comparator='all other age-eligible GMDB patients',
                         adjustment='age and sex', beta=row.target_beta, ci_lower=row.ci_95_lower, ci_upper=row.ci_95_upper,
                         n_complete=row.n_complete, n_disease=row.n_target, n_comparator=row.n_other,
                         q=row.q_four_models, correction_family='four prespecified within-GMDB models', family_hypotheses=4,
                         source='within_gmdb_age_sex_adjusted.csv'))
    for row in load('three_predefined_reference_adjusted_models.csv').itertuples(index=False):
        rows.append(dict(syndrome=row.syndrome, feature=row.feature, comparator='frozen FairFace reference',
                         adjustment='recorded sex; broad group; yaw; pitch; roll', beta=row.adjusted_disease_coefficient,
                         ci_lower=row.ci_95_lower, ci_upper=row.ci_95_upper, n_complete=row.n_total, n_disease=row.n_disease,
                         n_comparator=row.n_reference, q=row.bh_q_across_three_predefined_models,
                         correction_family='three prespecified reference models', family_hypotheses=3,
                         source='three_predefined_reference_adjusted_models.csv'))
    full = load('full_adjusted_screen.csv')
    anchors = load('clinical_highlights_120.csv')
    for row in anchors.itertuples(index=False):
        syndrome = getattr(row, 'disease', None) or getattr(row, 'syndrome', None)
        match = full.loc[full.syndrome.str.split('__').str[0].eq(syndrome) & full.feature.eq(row.feature)]
        assert len(match) == 1
        item = match.iloc[0]
        rows.append(dict(syndrome=item.syndrome, feature=item.feature, comparator='all other age-eligible GMDB patients',
                         adjustment='age and sex', beta=item.target_beta, ci_lower=item.ci_95_lower, ci_upper=item.ci_95_upper,
                         n_complete=item.n_complete, n_disease=item.n_target, n_comparator=item.n_other, q=item.q_global_bh,
                         correction_family='full prespecified 49-target x 120-feature screen', family_hypotheses=5880,
                         source='full_adjusted_screen.csv'))
    return pd.DataFrame(rows)


def literature_audit(primary, definitions):
    rows = []
    for index, (label, feature, book, section, support, interpretation) in enumerate(APPENDIX, 1):
        match = primary.loc[primary.disease.str.contains(label, regex=False) & primary.feature.eq(feature)]
        assert len(match) == 1, (label, feature)
        row = match.iloc[0]
        rows.append(dict(appendix_entry=f'A{index}', syndrome=row.disease, feature=feature, n=row.n_disease,
                         hedges_g=row.hedges_g, ci_lower=row.hedges_g_ci_95_low, ci_upper=row.hedges_g_ci_95_high,
                         global_q=row.mannwhitney_q_fdr_global, significant=bool(row.mannwhitney_q_fdr_global < .05),
                         clinical_support=support, exact_ratio_validated=False,
                         retain_supported_clinical_example=bool(support == 'direct_related_trait' and row.mannwhitney_q_fdr_global < .05),
                         interpretation=interpretation, clinical_source=f'https://www.ncbi.nlm.nih.gov/books/{book}/',
                         clinical_source_section=section, accessed='2026-11-09'))
    return pd.DataFrame(rows).merge(definitions[['feature', 'formula', 'units']], on='feature', validate='many_to_one')


def build():
    OUT.mkdir(parents=True, exist_ok=True)
    primary = load('S1_primary_all_5760_hypotheses.csv')
    definitions = load('Supplementary_Table_S6_exact_feature_definitions.csv')
    assert len(primary) == 5760 and primary.disease.nunique() == 48 and primary.mannwhitney_q_fdr_global.lt(.05).sum() == 955
    top, regions = recurrent_and_regions(primary, definitions)
    demos, groups, files = demographic_examples()
    demos = demos.merge(definitions[['feature', 'formula', 'units', 'anatomical_region']], on='feature', validate='many_to_one')
    adjusted = adjusted_examples(primary).merge(definitions[['feature', 'formula', 'units']], on='feature', validate='many_to_one')
    literature = literature_audit(primary, definitions)
    within = load('within_gmdb_results.csv')
    a = set(map(tuple, primary.loc[primary.mannwhitney_q_fdr_global.lt(.05), ['disease', 'feature']].to_numpy()))
    b = set(map(tuple, within.loc[within.q_global_bh.lt(.05), ['syndrome', 'feature']].to_numpy()))
    assert [len(a), len(b), len(a & b), len(a - b), len(b - a)] == [955, 222, 173, 782, 49]
    contrasts = pd.DataFrame([dict(fairface=955, within_GMDB=222, shared=173, fairface_only=782, within_only=49,
                                   interpretation='Different estimands; not an empirical false-positive rate')])
    for frame, name in [(top, 'top_five_recurrent_features.csv'), (regions, 'regional_discovery_source.csv'),
                        (demos, 'demographic_examples.csv'), (groups, 'demographic_A_group_summaries.csv'),
                        (adjusted, 'adjusted_examples.csv'), (literature, 'clinical_literature_audit_17.csv'),
                        (contrasts, 'comparator_overlap.csv')]:
        save(frame, name)
    names = files + ['S1_primary_all_5760_hypotheses.csv', 'Supplementary_Table_S6_exact_feature_definitions.csv',
                     'regional_summary_120.csv', 'within_gmdb_results.csv', 'within_gmdb_age_sex_adjusted.csv',
                     'three_predefined_reference_adjusted_models.csv', 'full_adjusted_screen.csv', 'clinical_highlights_120.csv']
    manifest = {'source_sha256': {name: hashlib.sha256((SOURCE / name).read_bytes()).hexdigest() for name in sorted(names)},
                'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                'primary': {'features': 120, 'syndromes': 48, 'hypotheses': 5760, 'significant': 955},
                'selection_rules': 'SELECTION_RULES.md; support-first selection without significance filtering'}
    (OUT / 'TABLE_SOURCE_MANIFEST.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps({'status': 'PASS', 'demographic_examples': len(demos), 'adjusted_examples': len(adjusted), 'literature_entries': len(literature)}))


if __name__ == '__main__':
    build()

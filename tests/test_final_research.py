import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('targeted_tables', ROOT / 'publication/scripts/final_research_tables.py')
tables = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tables)


def test_final_research_package_has_verified_source_correspondence():
    from publication.scripts.verify_final_research import verify
    assert verify()['status'] == 'PASS'


def test_demographic_selection_keeps_nonsignificant_supported_examples():
    rows, groups, files = tables.demographic_examples()
    assert len(rows) == 18 and len(files) == 6
    assert rows.groupby(['layer', 'analysis']).size().eq(3).all()
    assert rows.loc[rows.analysis.ne('A'), 'n_disease_1'].ge(10).all()
    assert rows.loc[rows.analysis.eq('C'), 'n_disease_2'].ge(10).all()
    assert rows.q.gt(.05).any() and rows.q.lt(.05).any()
    a = rows.loc[rows.analysis.eq('A')]
    assert a.feature.tolist() == tables.FEATURES * 2
    assert groups.n.ge(100).all()


def test_regional_source_change_is_detected():
    primary = tables.load('S1_primary_all_5760_hypotheses.csv')
    definitions = tables.load('Supplementary_Table_S6_exact_feature_definitions.csv')
    changed = primary.copy()
    index = changed.index[changed.mannwhitney_q_fdr_global.lt(.05)][0]
    changed.loc[index, 'mannwhitney_q_fdr_global'] = .9
    with pytest.raises(AssertionError):
        tables.recurrent_and_regions(changed, definitions)


def test_clinical_trait_support_never_validates_exact_ratio():
    primary = tables.load('S1_primary_all_5760_hypotheses.csv')
    definitions = tables.load('Supplementary_Table_S6_exact_feature_definitions.csv')
    audit = tables.literature_audit(primary, definitions)
    assert len(audit) == 17 and audit.retain_supported_clinical_example.sum() == 4
    assert not audit.exact_ratio_validated.any()
    assert not audit.loc[audit.feature.eq('malar_bulge_mean'), 'retain_supported_clinical_example'].any()

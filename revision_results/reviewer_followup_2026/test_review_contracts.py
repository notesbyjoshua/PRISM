import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest

import followup as f


def test_uninitialized_facekit_is_not_parent_repository(tmp_path, monkeypatch):
    monkeypatch.setattr(f, 'ROOT', tmp_path)
    source = tmp_path/'publication/historical_facekit/extractor.py'
    source.parent.mkdir(parents=True)
    source.write_text('historical snapshot fixture')
    (tmp_path/'FaceKit').mkdir()
    with patch.object(f.subprocess, 'check_output', return_value='recorded-submodule-sha\n') as git:
        provenance = f.facekit_provenance()
    assert provenance['facekit_checkout'] is None
    assert provenance['facekit_checkout_status'] == 'uninitialized'
    assert provenance['facekit_recorded_gitlink'] == 'recorded-submodule-sha'
    assert len(git.call_args_list) == 1


def test_facekit_top_level_must_be_distinct(tmp_path, monkeypatch):
    monkeypatch.setattr(f, 'ROOT', tmp_path)
    (tmp_path/'FaceKit').mkdir()
    (tmp_path/'FaceKit/.git').touch()
    source = tmp_path/'publication/historical_facekit/extractor.py'
    source.parent.mkdir(parents=True)
    source.write_text('fixture')
    with patch.object(f.subprocess, 'check_output', side_effect=[str(tmp_path), 'recorded-sha']):
        provenance = f.facekit_provenance()
    assert provenance['facekit_checkout'] is None
    assert provenance['facekit_checkout_status'] == 'not a distinct FaceKit repository'


def test_actual_facekit_checkout_is_informational(tmp_path, monkeypatch):
    monkeypatch.setattr(f, 'ROOT', tmp_path)
    (tmp_path/'FaceKit').mkdir()
    (tmp_path/'FaceKit/.git').touch()
    source = tmp_path/'publication/historical_facekit/extractor.py'
    source.parent.mkdir(parents=True)
    source.write_text('fixture')
    with patch.object(f.subprocess, 'check_output', side_effect=[str(tmp_path/'FaceKit'), 'actual-sha', 'recorded-sha']):
        provenance = f.facekit_provenance()
    assert provenance['facekit_checkout'] == 'actual-sha'
    assert provenance['facekit_recorded_gitlink'] == 'recorded-sha'


def test_independent_mannwhitney_handles_ties_and_equal_constants():
    path = Path(__file__).resolve().parent/'independent_review.py'
    spec = importlib.util.spec_from_file_location('independent_reviewer', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    from scipy.stats import mannwhitneyu
    for x, y in [(np.ones(12), np.ones(20)), (np.array([1., 1., 2., 3., 3.]), np.array([1., 2., 2., 4., 4.]))]:
        u, p = module.ranked_test(x, y)
        expected = mannwhitneyu(x, y, method='asymptotic', alternative='two-sided', use_continuity=True)
        assert u == expected.statistic
        assert np.isclose(p, expected.pvalue)


def test_hpo_offline_snapshot_records_acquisition_not_false_download(tmp_path, monkeypatch):
    import verify_hpo as hpo
    ontology = tmp_path/'source.obo'
    ontology.write_text('format-version: 1.2\ndata-version: test-release\n[Term]\nid: HP:0000343\nname: Long philtrum\n')
    output = tmp_path/'output'
    output.mkdir()
    monkeypatch.setattr(hpo, 'OUT', output)
    monkeypatch.setattr('sys.argv', ['verify_hpo.py', '--ontology-file', str(ontology)])
    with patch.object(hpo.urllib.request, 'urlopen', side_effect=AssertionError('Offline replay must not access network')):
        hpo.main()
    record = json.loads((output/'hpo_source_provenance.json').read_text())
    assert record['acquisition'] == 'local_snapshot'
    assert record['retrieved_date'] is None
    assert record['verified_utc']
    assert (output/'local_only/hp.obo').read_bytes() == ontology.read_bytes()


def test_exact_thirty_percent_missingness_is_eligible():
    target = np.array([70, 69, 10])
    other = np.array([100, 100, 100])
    assert f.eligible_feature_counts(target, other, 100, 100).tolist() == [True, False, False]
    assert bool(f.eligible_feature_counts(np.array(10), np.array(10), 14, 14))
    assert not bool(f.eligible_feature_counts(np.array(10), np.array(10), 15, 15))

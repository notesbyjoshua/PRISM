import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_active_primary_matches_validated_manifest():
    config = json.loads((ROOT/'ACTIVE_ANALYSIS.json').read_text())
    assert config['status'] == 'VALIDATED_ACTIVE'
    assert config['feature_count'] == 120 and config['synthetic_rows_in_primary'] == 0
    results = pd.read_csv(ROOT/config['results_directory']/'stats_each_disease_vs_healthy.csv')
    assert len(results) == config['primary_hypotheses'] == 5760
    assert results.feature.nunique() == 120 and results.disease.nunique() == 48
    assert results.mannwhitney_q_fdr_global.lt(.05).sum() == config['primary_significant'] == 955
    assert results.n_healthy.eq(config['reference_images']).all()


def test_active_api_uses_current_measurements():
    pytest.importorskip('fastapi')
    pytest.importorskip('httpx')
    from fastapi.testclient import TestClient
    from backend.app.main import app
    config = json.loads((ROOT/'ACTIVE_ANALYSIS.json').read_text())
    if not (ROOT/config['disease_csv']).exists():
        pytest.skip('Authorized private current measurements are local-only')
    with TestClient(app) as client:
        summary = client.get('/api/summary').json()
        assert summary['feature_count'] == 120 and summary['disease_count'] == 48
        assert summary['data_source'] == 'analysis' and not summary['is_demo']
        results = pd.read_csv(ROOT/config['results_directory']/'stats_each_disease_vs_healthy.csv')
        first = results.iloc[0]
        response = client.get('/api/distributions', params={'disease': first.disease, 'feature': first.feature})
        assert response.status_code == 200
        values = response.json()
        assert len(values['healthy']) == config['reference_images'] == 3409
        assert len(values['disease_values']) == first.n_disease


def test_active_manuscript_and_figures_are_current():
    text = (ROOT/'manuscript/main.tex').read_text()
    assert '5,760' in text and '955' in text and '120' in text
    assert '6,125' not in text and '811' not in text
    figures = ROOT/'publication/figures'
    assert (figures/'fig2_clinical_indices.png').is_file()
    assert not (figures/'fig2_demographic_analyses.pdf').exists()
    assert not (figures/'fig3_selection_overlap.pdf').exists()

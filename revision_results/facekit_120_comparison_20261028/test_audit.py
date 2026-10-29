import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

from audit import HERE, ROOT, FROZEN, bh, g, independent_mw, sha

OUT = Path(os.environ.get('PRISM_FACEKIT_AUDIT', HERE/'local_only/run_20261028_final'))


@pytest.mark.parametrize('values', [[.001, .05, .2, .04], [.01, np.nan, .5, .01], [1., 1., 1.]])
def test_bh_matches_independent_library(values):
    values = np.asarray(values)
    valid = np.isfinite(values)
    result = bh(values)
    np.testing.assert_allclose(result[valid], multipletests(values[valid], method='fdr_bh')[1])
    assert np.isnan(result[~valid]).all()


@pytest.mark.parametrize('seed', [1, 42, 2026])
def test_rank_based_u_and_tie_corrected_p(seed):
    rng = np.random.default_rng(seed)
    x, y = rng.integers(0, 8, 20), rng.integers(0, 8, 31)
    np.testing.assert_allclose(independent_mw(x, y), mannwhitneyu(x, y, alternative='two-sided', method='asymptotic'))


def test_all_ties_and_effect_sign():
    assert independent_mw(np.ones(10), np.ones(12))[1] == 1
    x, y = np.arange(15.)+2, np.arange(15.)
    assert g(x, y) > 0 and g(x, y) == -g(y, x)


def test_crosswalk_covers_both_universes_and_preserves_merge():
    frame = pd.read_csv(OUT/'feature_crosswalk.csv').fillna('')
    assert len(frame) == 125 and frame.old_feature.nunique() == 125
    assert frame.loc[frame.new_feature.ne(''), 'new_feature'].nunique() == 120
    assert frame.loc[frame.old_feature.eq('hairline_height'), 'new_feature'].item() == 'forehead_height'
    assert not frame.mapping_status.isin(['unchanged_exact', 'mathematically_equivalent']).any()
    assert frame.manual_review_required.all()


def test_missing_120_values_are_not_zero_or_historical_copies():
    metrics = json.loads((OUT/'comparison_metrics.json').read_text())
    assert metrics['decision'] == 'INDETERMINATE'
    for key in ['direction_agreement_meaningful', 'effect_similarity_meaningful', 'significance_retention']:
        assert metrics[key]['status'] == 'not_run'
        assert metrics[key]['numerator'] is metrics[key]['denominator'] is metrics[key]['fraction'] is None
    highlights = pd.read_csv(OUT/'clinical_highlights_comparison.csv')
    assert len(highlights) == 5 and highlights.g_120.isna().all() and highlights.q_120_global.isna().all()


def test_baseline_and_patient_aligned_photo_sensitivity():
    result = json.loads((OUT/'baseline_audit.json').read_text())
    assert result['historical_hypotheses'] == result['independent_rank_checks'] == 6125
    assert result['historical_significant'] == 811
    assert result['alternative_selected_images_changed'] == 94
    assert max(result['max_absolute_error'].values()) < 1e-10


def test_frozen_manifest_integrity_and_no_modifications():
    manifest = json.loads((FROZEN/'publication_manifest.json').read_text())
    for name, digest in manifest['files_sha256'].items():
        assert sha(FROZEN/name) == digest
    run = json.loads((OUT/'reproducibility_manifest.json').read_text())
    for name, digest in run['source_sha256'].items():
        assert sha(ROOT/name) == digest


def test_all_generated_results_are_git_ignored():
    import subprocess
    paths = [str(p.relative_to(ROOT)) for p in OUT.rglob('*') if p.is_file()]
    result = subprocess.run(['git', 'check-ignore', '--stdin'], cwd=ROOT, input='\n'.join(paths), text=True, capture_output=True)
    assert result.returncode == 0 and set(result.stdout.splitlines()) == set(paths)

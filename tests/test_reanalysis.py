import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'stats_testing'))
from bootstrap_utils import bootstrap_g
from prepare_cohorts import deduplicate_reference, select_patients


class ReanalysisTests(unittest.TestCase):
    def test_reference_duplicates_and_conflicts(self):
        frame = pd.DataFrame({'file': ['train/1.jpg', 'train/1.jpg', 'train/2.jpg'], 'value': [1, 1, 2]})
        selected, manifest = deduplicate_reference(frame)
        self.assertEqual(len(selected), 2)
        self.assertEqual(manifest.selection_reason.eq('repeated_reference_path').sum(), 1)
        frame.loc[1, 'value'] = 3
        with self.assertRaises(ValueError):
            deduplicate_reference(frame)

    def test_age_filter_precedes_patient_selection(self):
        frame = pd.DataFrame({'patient_id': ['a', 'a', 'b'], 'image_id': ['1', '2', '3'],
                              'age_year': [20, 5, 0], 'age_month': [0, 0, 0], 'frontal_ok': [True]*3,
                              'pose_yaw': [0, 3, 0], 'pose_pitch': [0]*3, 'pose_roll': [0]*3,
                              'internal_syndrome_name': ['s']*3, 'image_dataset_ethnicity_category': ['European']*3,
                              'gender': ['male']*3})
        selected, manifest = select_patients(frame, '3-9')
        self.assertEqual(list(selected.image_id), ['2'])
        self.assertEqual(manifest.selection_reason.eq('outside_or_unresolved_age_3_9').sum(), 2)
        shuffled, _ = select_patients(frame.sample(frac=1, random_state=1), '3-9')
        self.assertEqual(list(shuffled.image_id), ['2'])

    def test_bootstrap_resamples_both_groups_and_matches_formula(self):
        x, y = np.array([1., 2., 4.]), np.array([0., 3., 5., 6.])
        rng = np.random.default_rng(4)
        xb = x[rng.integers(3, size=(20, 3))]
        yb = y[rng.integers(4, size=(20, 4))]
        pooled = (2 * xb.var(axis=1, ddof=1) + 3 * yb.var(axis=1, ddof=1)) / 5
        with np.errstate(divide='ignore', invalid='ignore'):
            expected = (1 - 3/19) * (xb.mean(axis=1) - yb.mean(axis=1)) / np.sqrt(pooled)
        actual = bootstrap_g(x, y, 20, np.random.default_rng(4))
        np.testing.assert_allclose(actual, expected[np.isfinite(expected)])
        self.assertEqual(len(bootstrap_g([1], y, 20, rng)), 0)


if __name__ == '__main__':
    unittest.main()

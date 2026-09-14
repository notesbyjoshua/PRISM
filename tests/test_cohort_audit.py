import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from stats_testing.cohort_audit import ROOT, age_diagnostics, clean_id, feature_counts, patient_summary, quality_status, read_mapping_functions


class CohortAuditTests(unittest.TestCase):
    def test_repeated_and_missing_ids_are_distinct(self):
        frame = pd.DataFrame({'patient_id': ['a', 'a', 'b', None, '', 'unknown']})
        result = patient_summary(frame)
        self.assertEqual(result['known_patients'], 2)
        self.assertEqual(result['missing_patient_id_rows'], 3)
        self.assertEqual(result['patients_with_repeats'], 1)
        self.assertEqual(result['excess_rows_among_known_patients'], 1)
        self.assertEqual(result['rows_from_repeated_patients'], 2)

    def test_absent_patient_column_does_not_invent_people(self):
        result = patient_summary(pd.DataFrame({'file': ['a', 'b']}))
        self.assertEqual(result['known_patients'], 0)
        self.assertEqual(result['missing_patient_id_rows'], 2)

    def test_identifiers_preserve_leading_zeros(self):
        ids = clean_id(pd.Series(['001', '1', ' 001 ', 'NA']))
        self.assertEqual(ids.nunique(), 2)
        self.assertTrue(pd.isna(ids.iloc[-1]))

    def test_malformed_quality_is_not_truthy(self):
        result = quality_status(pd.Series([True, False, 'True', ' false ', 'yes', 1, None, '']))
        self.assertEqual(list(result), ['pass', 'fail', 'pass', 'fail', 'malformed', 'malformed', 'missing', 'missing'])

    def test_zero_age_remains_unresolved(self):
        frame = pd.DataFrame({'age_year': [0, 0, 3, 9, 10, None, 4, -1],
                              'age_month': [0, 6, 0, 11, 0, 0, 12, 0]})
        result = age_diagnostics(frame)
        self.assertEqual(result.age_status.iloc[0], 'zero_zero_unresolved')
        self.assertTrue(pd.isna(result.conditional_age_years.iloc[0]))
        self.assertEqual(result.conditional_age_years.iloc[1], .5)
        self.assertEqual(list(result.conditional_age_band.iloc[2:5]), ['3_to_under_10', '3_to_under_10', '10_to_under_18'])
        self.assertEqual(result.age_status.iloc[6], 'invalid_components')
        self.assertEqual(result.age_status.iloc[7], 'invalid_components')

    def test_missing_and_infinite_features_are_not_valid(self):
        frame = pd.DataFrame({'feature': [1, '2', None, 'bad', np.inf]})
        result = feature_counts(frame, ['feature']).iloc[0]
        self.assertEqual(result.valid_n, 2)
        self.assertEqual(result.missing_n, 2)
        self.assertEqual(result.infinite_n, 1)
        self.assertAlmostEqual(result.invalid_fraction, .6)

    def test_empty_feature_group(self):
        result = feature_counts(pd.DataFrame({'feature': []}), ['feature']).iloc[0]
        self.assertEqual(result.valid_n, 0)
        self.assertTrue(pd.isna(result.invalid_fraction))

    def test_unknown_mapping_agrees_between_AB_and_C(self):
        ab = read_mapping_functions(ROOT / 'stats_testing/ancestry_analysis.py')
        c = read_mapping_functions(ROOT / 'stats_testing/ancestry_disease_interaction.py')
        self.assertTrue(pd.isna(ab['map_fairface_layer1']('unrecognized')))
        self.assertEqual(ab['map_fairface_layer2']('Middle Eastern'), 'Middle Eastern')
        self.assertTrue(pd.isna(ab['map_gmdb_layer2']('Unknown', 'East Asian')))
        self.assertTrue(pd.isna(c['map_gmdb_layer2']('Unknown', 'East Asian')))

    def test_loading_mappings_does_not_run_analysis(self):
        source = (ROOT / 'stats_testing/ancestry_analysis.py').read_text()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'mapping.py'
            path.write_text("raise RuntimeError('must not run')\n" + source)
            mappings = read_mapping_functions(path)
            self.assertEqual(mappings['map_gmdb_layer1']('European'), 'European')

    def test_conflicting_metadata_excludes_missing_ids(self):
        frame = pd.DataFrame({'patient_id': ['a', 'a', 'b', 'b', None, None],
                              'gender': ['male', 'female', 'male', None, 'male', 'female']})
        counts = frame.assign(patient_id=clean_id(frame.patient_id)).groupby('patient_id').gender.nunique()
        self.assertEqual(int(counts.gt(1).sum()), 1)


if __name__ == '__main__':
    unittest.main()

import importlib.util
import tempfile
import unittest
from pathlib import Path

import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/'publication/scripts/generate_feature_definitions.py'


def load_generator():
    spec=importlib.util.spec_from_file_location('current_definitions',SCRIPT)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class PublicationSourceTests(unittest.TestCase):
    def test_current_source_matches_the_published_schema(self):
        source=load_generator().current_source()
        self.assertIn('GeometricFeatureExtractor',source)
        schema=pd.read_csv(ROOT/'publication/feature_definitions/current_120_feature_formulas.csv')
        self.assertEqual(schema.feature.nunique(),120)
        self.assertNotIn('historical_mapping',schema)

    def test_corrupted_source_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'extractor.py'
            path.write_text('changed source')
            with self.assertRaisesRegex(ValueError,'differs from the pinned'):
                load_generator().current_source(source_path=path)

    def test_noncurrent_schema_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            schema=pd.read_csv(ROOT/'publication/feature_definitions/current_120_feature_formulas.csv')
            schema=pd.concat([schema,schema.iloc[:1]],ignore_index=True)
            path=Path(directory)/'schema.csv';schema.to_csv(path,index=False)
            with self.assertRaisesRegex(ValueError,'Exactly 120'):
                load_generator().current_source(schema_path=path)


if __name__=='__main__':
    unittest.main()

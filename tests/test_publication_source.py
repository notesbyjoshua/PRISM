import hashlib
import importlib.util
import os
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'publication/scripts/generate_feature_definitions.py'
COMMIT = 'e474628f3cb17a83edad59ac6e1d9e5438753032'


def load_generator():
    spec = importlib.util.spec_from_file_location('publication_definitions', SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class PublicationSourceTests(unittest.TestCase):
    def test_snapshot_works_without_application_submodule_or_network(self):
        with patch.dict(os.environ, {}, clear=True):
            module = load_generator()
            with patch.object(module.subprocess, 'check_output', side_effect=AssertionError('Git access not allowed')):
                source = module.historical_source()
            self.assertEqual(hashlib.sha256(source.encode()).hexdigest(), module.SOURCE_SHA256)
            self.assertEqual(module.COMMIT, COMMIT)

    def test_corrupted_snapshot_is_rejected(self):
        with patch.dict(os.environ, {}, clear=True), tempfile.TemporaryDirectory() as folder:
            module = load_generator()
            module.SNAPSHOT = Path(folder) / 'extractor.py'
            module.SNAPSHOT.write_text('changed historical source')
            with self.assertRaisesRegex(ValueError, 'differs from the pinned'):
                module.historical_source()

    def test_explicit_repository_uses_pinned_revision_and_checks_content(self):
        with patch.dict(os.environ, {'PUBLICATION_FACEKIT_REPOSITORY': '/explicit/source'}):
            module = load_generator()
            snapshot = module.SNAPSHOT.read_text()
            with patch.object(module.subprocess, 'check_output', return_value=snapshot) as git:
                self.assertEqual(module.historical_source(), snapshot)
                git.assert_called_once_with(['git', '-C', '/explicit/source', 'show', f'{COMMIT}:{module.SOURCE_PATH}'], text=True)
            with patch.object(module.subprocess, 'check_output', return_value='a newer extractor'):
                with self.assertRaises(ValueError):
                    module.historical_source()


if __name__ == '__main__':
    unittest.main()

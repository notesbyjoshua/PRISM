import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('extraction_validation',ROOT/'data/dataset_creation/validate_current_120.py')
validation=importlib.util.module_from_spec(spec)
spec.loader.exec_module(validation)


@pytest.mark.parametrize('writable',[False,True])
def test_completed_run_hashes_use_verified_read_only_execution_archive(tmp_path,monkeypatch,writable):
    archive=tmp_path/'local_review/archive';archive.mkdir(parents=True)
    source=archive/'source.csv';source.write_text('feature,value\nsynthetic,1\n')
    digest=hashlib.sha256(source.read_bytes()).hexdigest()
    source.chmod(0o600 if writable else 0o400)
    manifest={'files':[{'original_path':'data/source.csv','archived_path':'source.csv','sha256':digest}]}
    (archive/'ARCHIVE_MANIFEST.json').write_text(json.dumps(manifest))
    (tmp_path/'ACTIVE_ANALYSIS.json').write_text(json.dumps({'retired_input_manifest':'local_review/archive/ARCHIVE_MANIFEST.json'}))
    monkeypatch.setattr(validation,'ROOT',tmp_path)
    if writable:
        with pytest.raises(ValueError,match='checksum could not be verified'):
            validation.recorded_file(tmp_path/'data/source.csv',digest)
    else:
        assert validation.recorded_file(tmp_path/'data/source.csv',digest)==source

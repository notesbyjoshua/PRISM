import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('current_package', ROOT/'publication/scripts/verify_reproducible_outputs.py')
package = importlib.util.module_from_spec(spec)
spec.loader.exec_module(package)


@pytest.fixture
def aggregate(tmp_path):
    path = tmp_path/'stats_results/summary.csv'
    path.parent.mkdir()
    path.write_text('feature,mean\nx,1\n')
    return tmp_path, path, {'files': {'stats_results/summary.csv': package.sha(path)}}


def test_committed_aggregate_verification_needs_no_private_input(aggregate):
    root, path, manifest = aggregate
    assert package.verify_manifest(root, manifest) == 1
    assert package.check_public_csvs(root, package.artifact_paths(root)) == 1
    assert not (root/'local_review').exists()


@pytest.mark.parametrize('change', ['modified', 'missing', 'unlisted'])
def test_changed_aggregate_exports_fail_closed(aggregate, change):
    root, path, manifest = aggregate
    if change == 'modified':
        path.write_text('feature,mean\nx,2\n')
    elif change == 'missing':
        path.unlink()
    else:
        (path.parent/'extra.csv').write_text('feature,mean\ny,3\n')
    with pytest.raises(ValueError, match='checksum|inventory'):
        package.verify_manifest(root, manifest)


def test_person_measurement_columns_cannot_enter_aggregate_package(aggregate):
    root, path, manifest = aggregate
    path.write_text('patient_id,feature,value\nsynthetic,x,1\n')
    with pytest.raises(ValueError, match='Person-level'):
        package.check_public_csvs(root, package.artifact_paths(root))


def test_manifest_cannot_read_outside_aggregate_inventory(aggregate):
    root, path, manifest = aggregate
    manifest['files']['../restricted.csv'] = 'not a valid path'
    with pytest.raises(ValueError, match='inventory'):
        package.verify_manifest(root, manifest)

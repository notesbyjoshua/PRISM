import importlib.util
import io
import sqlite3
import zipfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('extract120', ROOT/'data/dataset_creation/extract_current_120.py')
extract = importlib.util.module_from_spec(spec)
spec.loader.exec_module(extract)


def png():
    buffer = io.BytesIO()
    Image.new('RGB', (8, 8), color=(100, 100, 100)).save(buffer, format='PNG')
    return buffer.getvalue()


def fake_pipeline(faces=1, pose=True, fail=False, wrong_schema=False):
    pipeline = extract.Pipeline.__new__(extract.Pipeline)
    pipeline.features = [f'f{i}' for i in range(120)]
    pipeline.mp = SimpleNamespace(Image=lambda **kwargs: kwargs['data'], ImageFormat=SimpleNamespace(SRGB='RGB'))
    point = SimpleNamespace(x=.5, y=.5)
    result = SimpleNamespace(face_landmarks=[[point]*478 for _ in range(faces)],
                             facial_transformation_matrixes=[np.eye(4) for _ in range(faces)])
    pipeline.screen = pipeline.detector = SimpleNamespace(detect=lambda image: result)
    def measure(*args, **kwargs):
        if fail:
            raise RuntimeError('synthetic extractor failure')
        values = {feature: float(index) for index, feature in enumerate(pipeline.features)}
        if wrong_schema:
            values['unexpected'] = values.pop('f0')
        values.update(pose_yaw=0., pose_pitch=0., pose_roll=0., frontal_ok=pose, derotated=True)
        return values
    pipeline.extractor = SimpleNamespace(extract_from_landmarks=measure)
    return pipeline


@pytest.mark.parametrize('faces,reason', [(0, 'no_face'), (2, 'multi_face_or_ambiguous')])
def test_face_failures_are_explicit(faces, reason):
    result = fake_pipeline(faces).process(png())
    assert result['reason'] == reason and result['measurements'] is None


def test_corrupt_image_and_extractor_failure():
    assert fake_pipeline().process(b'corrupt')['reason'] == 'decode_failed'
    result = fake_pipeline(fail=True).process(png())
    assert result['reason'] == 'extractor_error' and result['measurements'] is None


def test_exact_schema_and_pose_failed_measurement_retention():
    with pytest.raises(extract.ExtractionContractError, match='schema'):
        fake_pipeline(wrong_schema=True).process(png())
    result = fake_pipeline(pose=False).process(png())
    assert result['reason'] == 'pose_failed' and len(result['measurements']) == 120
    assert result['feature_finite'] and not result['quality_passed']
    accepted = fake_pipeline().process(png())
    assert accepted['quality_passed'] and accepted['reason'] == 'accepted'


def test_model_checksum_mismatch_is_fatal(tmp_path):
    model = tmp_path/'bad.task'
    model.write_bytes(b'wrong weights')
    with pytest.raises(extract.ExtractionContractError, match='checksum'):
        extract.fingerprint(model)


def labels_and_archive(tmp_path, duplicate=False, conflict=False, missing_column=False, corrupt_zip=False):
    train = pd.DataFrame({'file': ['train/1.jpg'], 'age': ['3-9'], 'race': ['White'], 'gender': ['Female']})
    if duplicate:
        train = pd.concat([train, train], ignore_index=True)
        if conflict:
            train.loc[1, 'race'] = 'Black'
    if missing_column:
        train = train.drop(columns='age')
    val = pd.DataFrame({'file': ['val/1.jpg'], 'age': ['3-9'], 'race': ['Black'], 'gender': ['Male']})
    train.to_csv(tmp_path/'train.csv', index=False)
    val.to_csv(tmp_path/'val.csv', index=False)
    archive_path = tmp_path/'margin025.zip'
    if corrupt_zip:
        archive_path.write_bytes(b'not a zip')
    else:
        with zipfile.ZipFile(archive_path, 'w') as archive:
            archive.writestr('train/1.jpg', png())
            archive.writestr('val/1.jpg', png())
    return archive_path, tmp_path/'train.csv', tmp_path/'val.csv'


def test_exact_duplicate_labels_and_distinct_split_basenames(tmp_path):
    files = labels_and_archive(tmp_path, duplicate=True)
    candidates, lookup, audit = extract.prepare_fairface(*files, tmp_path)
    assert len(candidates) == 2 and audit['exact_duplicate_label_rows'] == 1
    assert set(lookup) == {'train/1.jpg', 'val/1.jpg'}
    assert audit['cross_split_path_collisions'] == 0


@pytest.mark.parametrize('variant', ['conflict', 'missing_column', 'corrupt_zip'])
def test_invalid_labels_and_archives_stop(variant, tmp_path):
    options = {variant: True}
    if variant == 'conflict':
        options['duplicate'] = True
    files = labels_and_archive(tmp_path, **options)
    with pytest.raises((extract.ExtractionContractError, zipfile.BadZipFile)):
        extract.prepare_fairface(*files, tmp_path)


def test_missing_archive_image_is_counted_and_loader_fails(tmp_path):
    files = labels_and_archive(tmp_path)
    with zipfile.ZipFile(files[0], 'w') as archive:
        archive.writestr('train/1.jpg', png())
    candidates, lookup, audit = extract.prepare_fairface(*files, tmp_path)
    assert audit['candidate_paths_missing_in_archive'] == 1
    loaders = extract.source_records('fairface', candidates, files[0], lookup)
    next(loaders)[2]()
    with pytest.raises(KeyError):
        next(loaders)[2]()


def test_checkpoint_interrupt_resume_and_duplicate_guard(tmp_path):
    path = tmp_path/'checkpoint.sqlite3'
    checkpoint = extract.Checkpoint(path, {'model': 'a'}, False)
    for index in range(26):
        checkpoint.add(str(index), {'index': index})
    checkpoint.connection.close()
    resumed = extract.Checkpoint(path, {'model': 'a'}, True)
    assert len(resumed.keys()) == 25
    resumed.add('25', {'index': 25})
    with pytest.raises(sqlite3.IntegrityError):
        resumed.add('0', {})
    resumed.close()
    with pytest.raises(extract.ExtractionContractError, match='Resume rejected'):
        extract.Checkpoint(path, {'model': 'different'}, True)
    with pytest.raises(FileExistsError):
        extract.Checkpoint(path, {'model': 'a'}, False)


def test_sampling_is_stable_balanced_and_does_not_oversample():
    rows = []
    for race, sex, n in [('A', 'Female', 4), ('A', 'Male', 3), ('B', 'Female', 2), ('B', 'Male', 1), ('B', 'Unknown', 1)]:
        rows += [{'image_path_key': f'{race}/{sex}/{i}', 'race': race, 'gender': sex} for i in range(n)]
    frame = pd.DataFrame(rows)
    first, shortfalls = extract.sampling(frame, ['A', 'B'])
    second, _ = extract.sampling(frame.sample(frac=1, random_state=7), ['A', 'B'])
    for flag in ['sample_A', 'sample_B', 'sample_C']:
        assert set(first.loc[first[flag], 'image_path_key']) == set(second.loc[second[flag], 'image_path_key'])
    assert first.sample_A.sum() == first.sample_B.sum() == len(frame)
    assert first.loc[first.sample_C].groupby(['race', 'gender']).size().eq(1).all()
    assert not first.loc[first.gender.eq('Unknown'), 'sample_C'].any()
    assert shortfalls.shortfall.gt(0).all()


def test_missing_race_sex_cell_produces_explicit_zero_capacity():
    frame = pd.DataFrame([{'image_path_key': 'train/1.jpg', 'race': 'A', 'gender': 'Female'}])
    selected, capacity = extract.sampling(frame, ['A', 'B'])
    assert selected.sample_C.sum() == 0
    assert len(capacity.loc[capacity.design.eq('C')]) == 4
    assert capacity.loc[capacity.design.eq('C'), 'target'].eq(0).all()


def test_gmdb_candidates_need_only_verified_arrow_metadata(tmp_path,monkeypatch):
    import hashlib
    import json
    import pyarrow as pa
    backup=tmp_path/'backup';backup.mkdir()
    count=10980
    values={'image_id':[str(i) for i in range(count)],'patient_id':[str(i) for i in range(count)],
            'internal_syndrome_id':['synthetic']*count,'internal_syndrome_name':['Synthetic']*count,
            'age_year':[5]*count,'age_month':[0]*count,'gender':['Unknown']*count,
            'ethnicity_category':['Unknown']*count,'ethnicity_sub_category':['Unknown']*count}
    table=pa.Table.from_pydict(values);path=backup/'train.arrow'
    with pa.OSFile(str(path),'wb') as sink:
        with pa.ipc.new_stream(sink,table.schema) as writer:writer.write_table(table)
    (tmp_path/'backup_manifest.json').write_text(json.dumps({'files':[{'file':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}]}))
    def retired_csv(*args,**kwargs):
        raise AssertionError('Retired CSV read')
    monkeypatch.setattr(extract.pd,'read_csv',retired_csv)
    out=tmp_path/'output';out.mkdir()
    candidates=extract.gmdb_candidates(backup,out)
    assert len(candidates)==count and candidates.key.is_unique
    assert not any('historical' in column or '125' in column for column in candidates)

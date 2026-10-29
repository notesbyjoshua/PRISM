import argparse
import hashlib
import io
import json
import os
import platform
import socket
import sqlite3
import subprocess
import sys
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath

os.environ['HF_HUB_OFFLINE'] = '1'
os.environ['HF_DATASETS_OFFLINE'] = '1'
os.environ['HF_HUB_DISABLE_TELEMETRY'] = '1'
os.umask(0o077)
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'FaceKit/src'))
sys.path.insert(0, str(ROOT/'stats_testing'))

import numpy as np
import pandas as pd
import pyarrow as pa
from PIL import Image, UnidentifiedImageError

COMMIT = '5764bdd4af9360b55a046a16579fb5a890a8182f'
MODEL_SHA = '64184e229b263107bc2b804c6625db1341ff2bb731874b0bcc2fe6544e0bc9ff'
SEED = 20261028
PRIVATE = ROOT/'local_review/feature_extraction_120'


class ExtractionContractError(ValueError):
    pass


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(8*1024*1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path, value):
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    os.replace(temporary, path)


def atomic_csv(path, frame):
    temporary = path.with_suffix(path.suffix+'.tmp')
    frame.to_csv(temporary, index=False)
    os.replace(temporary, path)


def now():
    return datetime.now(timezone.utc).isoformat()


def plain(value):
    if isinstance(value, (list, dict, str, bool)) or value is None:
        return value
    if pd.isna(value):
        return None
    return value.item() if isinstance(value, np.generic) else value


def image_key(value):
    path = PurePosixPath(str(value).replace('\\', '/'))
    if len(path.parts) < 2 or path.parts[-2] not in {'train', 'val'} or path.name in {'', '.', '..'}:
        raise ExtractionContractError('Invalid FairFace split/image path')
    return '/'.join(path.parts[-2:])


def fingerprint(model):
    from facekit.core.geometric.extractor import GeometricFeatureExtractor
    from facekit.api import get_registered
    commit = subprocess.check_output(['git', '-C', str(ROOT/'FaceKit'), 'rev-parse', 'HEAD'], text=True).strip()
    if commit != COMMIT or sha(model) != MODEL_SHA:
        raise ExtractionContractError('FaceKit commit or landmark-model checksum differs from the pinned pipeline')
    if subprocess.check_output(['git', '-C', str(ROOT/'FaceKit'), 'status', '--porcelain'], text=True).strip():
        raise ExtractionContractError('FaceKit working tree must be clean')
    features = GeometricFeatureExtractor._init_feature_columns()
    if len(features) != 120 or len(set(features)) != 120 or get_registered():
        raise ExtractionContractError('Expected exactly 120 base morphology features and no plugins')
    import mediapipe as mp
    import PIL
    import cv2
    signature = {'facekit_commit': commit, 'model_sha256': sha(model), 'model_asset': 'MediaPipe Face Landmarker float16 version 1',
        'extractor_sha256': sha(ROOT/'FaceKit/src/facekit/core/geometric/extractor.py'),
        'canonical_depth_sha256': sha(ROOT/'FaceKit/src/facekit/core/geometric/canonical_depth.py'),
        'script_sha256': sha(Path(__file__)), 'features': features, 'seed': SEED,
        'thresholds_degrees': {'yaw': 15, 'pitch': 15, 'roll': 10},
        'versions': {'python': platform.python_version(), 'numpy': np.__version__, 'pandas': pd.__version__,
                     'pyarrow': pa.__version__, 'mediapipe': mp.__version__, 'Pillow': PIL.__version__, 'opencv': cv2.__version__},
        'preprocessing': 'Pillow full decode to RGB; stored orientation preserved; no EXIF transpose, crop, resize, enhancement or pixel export',
        'face_selection': 'same model with num_faces=2 screens ambiguity; separate num_faces=1 detector supplies FaceKit-default float pixel landmarks and matrix',
        'geometry': '478 x/y float pixel landmarks; SVD pure rotation plus fixed canonical depth; derotated must be true',
        'unmeasured_quality': ['occlusion', 'expression', 'clinical validity', 'identity uniqueness', 'undetected additional faces'],
        'paper_execution_pin': 'unverified; no claim that the paper used these exact execution assets'}
    return signature


class Pipeline:
    def __init__(self, model, features):
        import mediapipe as mp
        from mediapipe.tasks import python
        from mediapipe.tasks.python import vision
        from facekit.core.geometric.extractor import GeometricFeatureExtractor
        self.mp, self.features = mp, features
        self.extractor = GeometricFeatureExtractor(model_path=model, yaw_thresh=15, pitch_thresh=15, roll_thresh=10)
        options = {'base_options': python.BaseOptions(model_asset_path=str(model)),
                   'running_mode': vision.RunningMode.IMAGE, 'output_facial_transformation_matrixes': True,
                   'output_face_blendshapes': False}
        self.screen = vision.FaceLandmarker.create_from_options(vision.FaceLandmarkerOptions(**options, num_faces=2))
        self.detector = vision.FaceLandmarker.create_from_options(vision.FaceLandmarkerOptions(**options, num_faces=1))

    def close(self):
        self.screen.close()
        self.detector.close()

    def process(self, data):
        result = {'reason': '', 'measurements': None, 'quality_passed': False, 'feature_finite': False}
        try:
            with Image.open(io.BytesIO(data)) as image:
                result.update(width=image.width, height=image.height, exif_orientation=int(image.getexif().get(274, 1)))
                rgb = np.array(image.convert('RGB'))
        except (UnidentifiedImageError, OSError, ValueError, Image.DecompressionBombError):
            result['reason'] = 'decode_failed'
            return result
        image = self.mp.Image(image_format=self.mp.ImageFormat.SRGB, data=np.ascontiguousarray(rgb))
        try:
            screened = self.screen.detect(image)
            count = len(screened.face_landmarks)
            result['screen_detected_faces'] = count
            if count != 1:
                result['reason'] = 'no_face' if count == 0 else 'multi_face_or_ambiguous'
                return result
            detected = self.detector.detect(image)
            if len(detected.face_landmarks) != 1 or len(detected.facial_transformation_matrixes) != 1:
                result['reason'] = 'multi_face_or_ambiguous'
                return result
            lm = np.array([[point.x*rgb.shape[1], point.y*rgb.shape[0]] for point in detected.face_landmarks[0]])
            matrix = np.asarray(detected.facial_transformation_matrixes[0])
            if lm.shape != (478, 2) or matrix.shape != (4, 4) or not np.isfinite(lm).all() or not np.isfinite(matrix).all():
                result['reason'] = 'invalid_landmarks_or_pose'
                return result
            values = self.extractor.extract_from_landmarks(lm, matrix, frontal_check=False)
            actual = [name for name in values if name not in {'pose_yaw', 'pose_pitch', 'pose_roll', 'frontal_ok', 'derotated'}]
            if actual != self.features:
                raise ExtractionContractError('Actual morphology feature schema differs from the validated 120-feature order')
            measurements = [float(values[name]) for name in self.features]
            result.update({key: plain(values[key]) for key in ['pose_yaw', 'pose_pitch', 'pose_roll', 'frontal_ok', 'derotated']})
            result['landmark_count'] = 478
            if not np.isfinite(measurements).all():
                result['reason'] = 'nonfinite_features'
                return result
            result['feature_finite'] = True
            result['measurements'] = measurements
            if not values['derotated']:
                result['reason'] = 'pose_correction_unavailable'
            elif not values['frontal_ok']:
                result['reason'] = 'pose_failed'
            else:
                result.update(reason='accepted', quality_passed=True)
            return result
        except ExtractionContractError:
            raise
        except (RuntimeError, ValueError, np.linalg.LinAlgError, ZeroDivisionError, IndexError, TypeError) as error:
            result.update(reason='extractor_error', error_class=type(error).__name__)
            return result


class Checkpoint:
    def __init__(self, path, signature, resume):
        exists = path.exists()
        if exists and not resume:
            raise FileExistsError('Checkpoint already exists; pass --resume for an identical run')
        self.connection = sqlite3.connect(path)
        self.connection.execute('PRAGMA journal_mode=WAL')
        self.connection.execute('PRAGMA synchronous=FULL')
        self.connection.execute('CREATE TABLE IF NOT EXISTS settings (signature TEXT NOT NULL)')
        self.connection.execute('CREATE TABLE IF NOT EXISTS records (image_key TEXT PRIMARY KEY, record TEXT NOT NULL)')
        encoded = json.dumps(signature, sort_keys=True, allow_nan=False)
        previous = self.connection.execute('SELECT signature FROM settings').fetchone()
        if previous and previous[0] != encoded:
            self.connection.close()
            raise ExtractionContractError('Resume rejected: model, code, environment, source or schema changed')
        if not previous:
            self.connection.execute('INSERT INTO settings VALUES (?)', (encoded,))
        self.connection.commit()
        self.pending = 0

    def keys(self):
        return {row[0] for row in self.connection.execute('SELECT image_key FROM records')}

    def add(self, key, record):
        self.connection.execute('INSERT INTO records VALUES (?, ?)', (key, json.dumps(record, allow_nan=False)))
        self.pending += 1
        if self.pending >= 25:
            self.commit()

    def commit(self):
        self.connection.commit()
        self.pending = 0

    def records(self):
        return [json.loads(row[0]) for row in self.connection.execute('SELECT record FROM records ORDER BY image_key')]

    def close(self):
        self.commit()
        self.connection.execute('PRAGMA wal_checkpoint(TRUNCATE)')
        self.connection.close()


def prepare_fairface(archive_path, train_path, val_path, out):
    if train_path.resolve() == val_path.resolve():
        raise ExtractionContractError('Train and validation labels must be separate files')
    parts = []
    required = ['file', 'age', 'race', 'gender']
    for split, path in [('train', train_path), ('val', val_path)]:
        frame = pd.read_csv(path)
        if not set(required).issubset(frame):
            raise ExtractionContractError('Label file is missing required columns')
        if frame[['file', 'age', 'race']].isna().any().any():
            raise ExtractionContractError('Label file contains missing required path, age or race values')
        for field in required:
            frame[field] = frame[field].astype('string').str.strip()
        frame['image_path_key'] = frame.file.map(image_key)
        if not frame.image_path_key.str.startswith(split+'/').all():
            raise ExtractionContractError('Train/validation path collision or split-label inconsistency')
        frame['split'] = split
        parts.append(frame)
    full = pd.concat(parts, ignore_index=True)
    conflicts = full.groupby('image_path_key')[['race', 'gender', 'age', 'split']].nunique(dropna=False).gt(1).any(axis=1)
    if conflicts.any():
        raise ExtractionContractError('Conflicting duplicate FairFace labels')
    duplicate_count = int(full.image_path_key.duplicated().sum())
    unique = full.drop_duplicates('image_path_key').copy()
    candidates = unique.loc[unique.age.eq('3-9'), ['image_path_key', 'split', 'race', 'gender', 'age']].copy()
    candidates['gender'] = candidates.gender.fillna('Unknown')
    candidates['dataset'] = 'FairFace'
    candidates['image_variant'] = 'margin025'
    candidates['key'] = candidates.image_path_key
    if candidates.empty:
        raise ExtractionContractError('No age-3-9 FairFace candidates')
    with zipfile.ZipFile(archive_path) as archive:
        bad_member = archive.testzip()
        if bad_member:
            raise ExtractionContractError('ZIP CRC integrity check failed')
        lookup = {}
        for member in archive.infolist():
            if member.is_dir():
                continue
            try:
                key = image_key(member.filename)
            except ExtractionContractError:
                continue
            if key in lookup:
                raise ExtractionContractError('Duplicate image key in FairFace archive')
            lookup[key] = member.filename
    audit = {'all_label_rows': len(full), 'unique_label_paths': len(unique), 'exact_duplicate_label_rows': duplicate_count,
             'conflicting_label_paths': 0, 'cross_split_path_collisions': 0, 'age_3_9_candidates': len(candidates),
             'race_categories': sorted(unique.race.dropna().unique()), 'archive_indexed_image_paths': len(lookup),
             'all_label_paths_missing_in_archive': int((~unique.image_path_key.isin(lookup)).sum()),
             'candidate_paths_missing_in_archive': int((~candidates.key.isin(lookup)).sum()), 'zip_crc_check': 'PASS',
             'identity_uniqueness': 'Unique image paths do not establish unique people', 'image_variant': 'margin025',
             'license': 'CC BY 4.0, as stated by https://github.com/joojs/fairface; no redistribution performed'}
    atomic_json(out/'input_audit.json', audit)
    counts = full.groupby(['split', 'age', 'race', 'gender'], dropna=False).size().rename('label_rows').reset_index()
    atomic_csv(out/'label_counts_before_filtering.csv', counts)
    atomic_csv(out/'candidate_manifest.csv', candidates)
    return candidates, lookup, audit


def gmdb_candidates(backup, out):
    source = backup.parent
    manifest = json.loads((source/'backup_manifest.json').read_text())
    for entry in manifest['files']:
        if sha(backup/entry['file']) != entry['sha256']:
            raise ExtractionContractError('Private GMDB backup checksum mismatch')
    columns = ['image_id', 'patient_id', 'internal_syndrome_id', 'internal_syndrome_name', 'age_year', 'age_month',
               'gender', 'ethnicity_category', 'ethnicity_sub_category']
    frames = []
    for path in sorted(backup.glob('*.arrow')):
        with pa.memory_map(str(path), 'r') as stream:
            for batch in pa.ipc.open_stream(stream):
                frame = batch.select(columns).to_pandas()
                frame['shard'] = path.name
                frames.append(frame)
    candidates = pd.concat(frames, ignore_index=True)
    candidates['source_image_id'] = candidates.image_id.astype(str)
    candidates['key'] = 'gmdb:'+candidates.source_image_id
    candidates['dataset'] = 'GMDB_restored'
    candidates['cache_split'] = 'train'
    if candidates.key.duplicated().any() or len(candidates) != 10980:
        raise ExtractionContractError('Unexpected GMDB cache key universe')
    atomic_csv(out/'candidate_manifest.csv', candidates)
    atomic_json(out/'input_audit.json', {'candidates': len(candidates), 'private_backup_verified': True,
        'selection': 'All cache records; no historical measurement inputs or comparison membership flags',
        'permission': 'User explicitly requested GMDB extraction in this turn; private local processing only'} )
    return candidates


def source_records(dataset, candidates, source, lookup=None):
    metadata = {row['key']: {name: plain(value) for name, value in row.items()} for row in candidates.to_dict('records')}
    if dataset == 'fairface':
        with zipfile.ZipFile(source) as archive:
            for key in sorted(metadata):
                yield key, metadata[key], (lambda k=key: archive.read(lookup[k]))
    else:
        for path in sorted(source.glob('*.arrow')):
            with pa.memory_map(str(path), 'r') as stream:
                for batch in pa.ipc.open_stream(stream):
                    ids = batch.column(batch.schema.get_field_index('image_id'))
                    images = batch.column(batch.schema.get_field_index('image_restored'))
                    for index in range(batch.num_rows):
                        key = 'gmdb:'+str(ids[index].as_py())
                        yield key, metadata[key], (lambda cell=images[index]: cell.as_py()['bytes'])


def sampling(accepted, races, seed=SEED):
    accepted = accepted.sort_values('image_path_key').copy()
    accepted['sample_A'] = False
    accepted['sample_B'] = True
    accepted['sample_C'] = False
    rng = np.random.default_rng(seed)
    shortfalls = []
    for race in sorted(races):
        indexes = accepted.index[accepted.race.eq(race)].to_numpy()
        chosen = rng.permutation(indexes)[:500]
        accepted.loc[chosen, 'sample_A'] = True
        shortfalls.append({'design': 'A', 'race': race, 'gender': 'all', 'capacity': len(indexes),
                           'target': 500, 'selected': len(chosen), 'shortfall': max(0, 500-len(chosen))})
    cells = [(race, sex) for race in sorted(races) for sex in ['Female', 'Male']]
    capacities = [int((accepted.race.eq(race) & accepted.gender.eq(sex)).sum()) for race, sex in cells]
    target = min([250]+capacities)
    for (race, sex), capacity in zip(cells, capacities):
        indexes = accepted.index[accepted.race.eq(race) & accepted.gender.eq(sex)].to_numpy()
        chosen = rng.permutation(indexes)[:target]
        accepted.loc[chosen, 'sample_C'] = True
        shortfalls.append({'design': 'C', 'race': race, 'gender': sex, 'capacity': capacity,
                           'target': target, 'selected': len(chosen), 'shortfall': max(0, 250-len(chosen))})
    return accepted, pd.DataFrame(shortfalls)


def export(dataset, out, checkpoint, candidates, signature, completed, input_audit):
    records = checkpoint.records()
    features = signature['features']
    rows = []
    for item in records:
        row = dict(item['metadata'])
        qc = item['qc']
        row.update({key: value for key, value in qc.items() if key != 'measurements'})
        row.update(dict(zip(features, qc['measurements'] or [None]*120)))
        row['image_sha256'] = item['image_sha256']
        row['facekit_commit'] = COMMIT
        row['model_sha256'] = MODEL_SHA
        rows.append(row)
    frame = pd.DataFrame(rows)
    if frame.empty:
        return
    assert frame.key.is_unique and set(frame.key).issubset(set(candidates.key))
    assert len(features) == 120
    has_values = frame[features].notna().all(axis=1)
    assert np.isfinite(frame.loc[has_values, features].to_numpy(float)).all()
    assert frame.loc[frame.quality_passed, features].notna().all().all()
    status = 'FULL' if completed and len(frame) == len(candidates) else 'PARTIAL'
    accepted = frame.loc[frame.quality_passed].copy()
    if dataset == 'fairface':
        accepted, shortfalls = sampling(accepted, input_audit['race_categories'])
        atomic_csv(out/'sampling_capacities_and_shortfalls.csv', shortfalls)
        for flag, name in [('sample_A', 'balanced'), ('sample_B', 'all_eligible'), ('sample_C', 'race_sex')]:
            columns = ['image_path_key', 'split', 'race', 'gender', 'age', 'image_sha256']
            atomic_csv(out/f'fairface_120_{name}_manifest.csv', accepted.loc[accepted[flag], columns])
        frame = frame.merge(accepted[['key', 'sample_A', 'sample_B', 'sample_C']], how='left', on='key', validate='one_to_one')
        for flag in ['sample_A', 'sample_B', 'sample_C']:
            frame[flag] = frame[flag].fillna(False).astype(bool)
        grouping = ['race', 'gender', 'split']
        counts = frame.groupby(grouping, dropna=False).agg(attempted=('key', 'size'), quality_passed=('quality_passed', 'sum'),
            selected_A=('sample_A', 'sum'), selected_B=('sample_B', 'sum'), selected_C=('sample_C', 'sum')).reset_index()
        capacity = candidates.groupby(grouping, dropna=False).size().rename('candidates').reset_index()
        counts = capacity.merge(counts, how='left', on=grouping).fillna(0)
        failure = frame.loc[~frame.quality_passed].groupby(grouping+['reason'], dropna=False).size().rename('failed_n').reset_index()
        atomic_csv(out/'cohort_counts.csv', counts)
        atomic_csv(out/'quality_failures_by_group.csv', failure)
        path = out/'fairface_120_all_eligible_measurements.parquet'
        temporary = path.with_suffix('.parquet.tmp')
        accepted.to_parquet(temporary, index=False)
        os.replace(temporary, path)
        atomic_csv(out/'fairface_120_extraction_status.csv', frame.drop(columns=features))
        sizes = {flag: int(accepted[flag].sum()) for flag in ['sample_A', 'sample_B', 'sample_C']}
    else:
        for subset, name in [(frame, 'all_records'), (accepted, 'quality_passed')]:
            path = out/f'gmdb_120_{name}_measurements.parquet'
            temporary = path.with_suffix('.parquet.tmp')
            subset.to_parquet(temporary, index=False)
            os.replace(temporary, path)
        atomic_csv(out/'gmdb_120_extraction_status.csv', frame.drop(columns=features))
        counts = frame.groupby(['reason'], dropna=False).size().rename('n').reset_index()
        atomic_csv(out/'cohort_counts.csv', counts)
        sizes = {'all_records_with_120_finite': int(has_values.sum())}
    failure_counts = {str(key): int(value) for key, value in frame.reason.value_counts().items()}
    summary = {'status': status, 'candidates': len(candidates), 'attempted': len(frame), 'pending': len(candidates)-len(frame),
               'quality_passed': len(accepted), 'reason_counts': failure_counts, 'samples': sizes,
               'finite_morphology_columns': 120, 'completed_at_utc': now(), 'input_audit': input_audit,
               'analysis': 'NOT RUN; extraction and QC bookkeeping only'}
    atomic_json(out/'extraction_summary.json', summary)
    text = [f'# {dataset.upper()} 120-feature extraction QA', '', f'**Execution status: {status}.** No PRISM statistical analyses were run.', '',
        f'Candidates: {len(candidates):,}; attempted: {len(frame):,}; quality-passed: {len(accepted):,}; pending: {len(candidates)-len(frame):,}.', '',
        '## Quality accounting', '', '```json', json.dumps(failure_counts, indent=2), '```', '',
        '## Cohort / retention counts', '', '```json', json.dumps(sizes, indent=2), '```', '',
        'Accepted rows have exactly 120 named, finite morphology measurements. Image keys, source labels and quality/pose metadata are separate fields. Duplicate image keys are prohibited. All accepted rows use full matrix-based canonicalization; unknown/roll-only corrections are rejected.', '',
        'Images were decoded to RGB in memory without crop, resize, EXIF rotation or restoration. A two-face detector screens ambiguity, and an independent one-face detector supplies the standard FaceKit geometry input. This cannot rule out undetected faces, occlusion or expression artifacts. No such unmeasured quality outcomes were invented.', '',
        'The 15/15/10-degree gate uses current FaceKit frontal_ok semantics. For GMDB, finite measurements from pose-failed records remain in all_records with explicit QC flags; quality_passed is a separate subset. Patient selection belongs to the separate current statistical preparation. No patient images, landmarks, embeddings or binary pixels are exported.', '',
        'FairFace A uses up to 500 per source race category. B includes all quality-passed images. C uses the same count in every race x Female/Male cell: min(250, smallest observed cell capacity). Unknown sex is preserved in A/B but excluded from C. Shortfalls from the nominal 250-per-cell goal are explicit. A/B/C overlap and are not independent cohorts. Labels remain source labels, not verified ancestry, exact age or clinical health.', '',
        'Restart requires --resume with identical inputs, code, model, versions and feature schema. SQLite primary keys prevent duplicate rows; transactions checkpoint every 25 records. Atomic exports may be rebuilt within that same run. A partial scan must never be used as a final balanced reference.', '']
    (out/'EXTRACTION_QA.md').write_text('\n'.join(text))
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('dataset', choices=['gmdb', 'fairface'])
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--archive', type=Path)
    parser.add_argument('--train-labels', type=Path)
    parser.add_argument('--val-labels', type=Path)
    parser.add_argument('--gmdb-backup', type=Path)
    parser.add_argument('--model', type=Path, default=Path.home()/'.cache/facekit/face_landmarker.task')
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--pilot-only', action='store_true')
    args = parser.parse_args()
    out = args.output.resolve()
    if PRIVATE not in out.parents:
        raise ExtractionContractError('Extraction outputs must stay in ignored local_review/feature_extraction_120')
    if out.exists() and not args.resume:
        raise FileExistsError('Choose a new output directory or use --resume')
    out.mkdir(parents=True, mode=0o700, exist_ok=args.resume)
    if subprocess.run(['git', 'check-ignore', '-q', str(out/'provenance.json')], cwd=ROOT).returncode:
        raise ExtractionContractError('Output destination is not Git-ignored')
    signature = fingerprint(args.model)
    features = signature['features']
    if args.dataset == 'fairface':
        if not all([args.archive, args.train_labels, args.val_labels]):
            raise ExtractionContractError('FairFace archive and both official labels are required')
        source_paths = [args.archive, args.train_labels, args.val_labels]
        candidates, lookup, input_audit = prepare_fairface(args.archive, args.train_labels, args.val_labels, out)
        source = args.archive
        pilot_groups = ['race', 'gender', 'split']
    else:
        if not args.gmdb_backup:
            raise ExtractionContractError('Verified private GMDB backup is required')
        source_paths = sorted(args.gmdb_backup.glob('*.arrow'))
        source_paths += [args.gmdb_backup/'dataset_info.json']
        candidates = gmdb_candidates(args.gmdb_backup, out)
        input_audit = json.loads((out/'input_audit.json').read_text())
        lookup, source = None, args.gmdb_backup
        pilot_groups = ['ethnicity_category', 'gender', 'cache_split']
    signature['input_sha256'] = {str(path.resolve()): sha(path) for path in source_paths}
    checkpoint = Checkpoint(out/'measurements_checkpoint.sqlite3', signature, args.resume)
    done = checkpoint.keys()
    if not done.issubset(set(candidates.key)):
        raise ExtractionContractError('Checkpoint contains keys outside the frozen candidate universe')
    provenance = dict(signature, status='RUNNING', started_or_resumed_at_utc=now(), dataset=args.dataset,
        model_path=str(args.model), input_paths=[str(path.resolve()) for path in source_paths], output_directory=str(out),
        privacy='Private local processing authorized by user; no upload, statistical analysis or image export',
        sources={'FairFace': 'https://github.com/joojs/fairface', 'variant': 'margin025',
                 'GMDB': 'GestaltGANwglabshare/gmdb_enhanced_v1.1.0_new; verified restored-image Arrow backup'})
    atomic_json(out/'provenance.json', provenance)
    schema_rows = []
    definitions = pd.read_csv(ROOT/'publication/feature_definitions/current_120_feature_formulas.csv')
    if definitions.feature.tolist() != features or not definitions.source_commit.eq(COMMIT).all() or not definitions.extractor_sha256.eq(signature['extractor_sha256']).all():
        raise ExtractionContractError('Current feature definitions differ from the pinned extraction schema')
    schema_rows = definitions.to_dict('records')
    atomic_csv(out/'feature_schema.csv', pd.DataFrame(schema_rows))
    pipeline = Pipeline(args.model, features)
    pilot = candidates.groupby(pilot_groups, dropna=False, group_keys=False).head(1)
    pilot_keys = set(pilot.key)
    pilot_checks = []
    started, processed = time.monotonic(), len(done)
    try:
        for phase in ['pilot', 'full']:
            if phase == 'full' and args.pilot_only:
                break
            for key, metadata, loader in source_records(args.dataset, candidates, source, lookup):
                if key in done or (phase == 'pilot' and key not in pilot_keys):
                    continue
                try:
                    data = loader()
                    if not data:
                        raise ValueError('No encoded image data')
                    qc = pipeline.process(data)
                    if phase == 'pilot':
                        repeat = pipeline.process(data)
                        equal = qc['reason'] == repeat['reason']
                        if qc['measurements'] is not None:
                            equal = equal and np.allclose(qc['measurements'], repeat['measurements'], rtol=1e-10, atol=1e-10)
                        pilot_checks.append({'key': key, 'deterministic': bool(equal), 'reason': qc['reason'],
                                             'max_absolute_feature_change': float(np.max(np.abs(np.asarray(qc['measurements'])-repeat['measurements']))) if qc['measurements'] is not None else None})
                        if not equal:
                            raise ExtractionContractError('Repeated image extraction exceeded deterministic tolerance')
                    digest = hashlib.sha256(data).hexdigest()
                except (KeyError, zipfile.BadZipFile, OSError) as error:
                    qc = {'reason': 'missing_in_archive' if isinstance(error, KeyError) else 'decode_failed',
                          'quality_passed': False, 'feature_finite': False, 'measurements': None, 'error_class': type(error).__name__}
                    digest = None
                checkpoint.add(key, {'metadata': metadata, 'qc': qc, 'image_sha256': digest})
                done.add(key)
                processed += 1
                if processed % 100 == 0:
                    checkpoint.commit()
                    atomic_json(out/'progress.json', {'completed': len(done), 'candidates': len(candidates),
                        'elapsed_this_session_seconds': round(time.monotonic()-started, 1), 'phase': phase, 'updated_at_utc': now()})
                    print(f'{args.dataset}: {len(done)}/{len(candidates)} records checkpointed', flush=True)
            checkpoint.commit()
            if phase == 'pilot':
                if pilot_checks:
                    atomic_csv(out/'private_pilot_determinism.csv', pd.DataFrame(pilot_checks))
                    if not any(item['reason'] == 'accepted' for item in pilot_checks):
                        raise ExtractionContractError('Pilot yielded no fully accepted images')
                elif not (out/'private_pilot_determinism.csv').exists():
                    raise ExtractionContractError('No pilot determinism evidence for resumed run')
                print(f'{args.dataset}: pilot passed; cached {len(pilot_keys)} demographic/selection cells', flush=True)
        summary = export(args.dataset, out, checkpoint, candidates, signature, not args.pilot_only, input_audit)
        for path, expected in signature['input_sha256'].items():
            if sha(Path(path)) != expected:
                raise ExtractionContractError('An input changed during extraction')
        if sha(args.model) != MODEL_SHA or sha(Path(__file__)) != signature['script_sha256']:
            raise ExtractionContractError('Model or extraction script changed during execution')
        provenance.update(status=summary['status'], completed_at_utc=now(), counts=summary,
                          source_hashes_unchanged=True, model_hash_unchanged=True)
        atomic_json(out/'provenance.json', provenance)
        print(json.dumps({key: summary[key] for key in ['status', 'candidates', 'attempted', 'quality_passed', 'reason_counts', 'samples']}, indent=2), flush=True)
    except BaseException as error:
        checkpoint.commit()
        provenance.update(status='INTERRUPTED' if isinstance(error, KeyboardInterrupt) else 'FAILED',
                          error_class=type(error).__name__, stopped_at_utc=now(), checkpoint_rows=len(checkpoint.keys()))
        atomic_json(out/'provenance.json', provenance)
        raise
    finally:
        pipeline.close()
        checkpoint.close()
    command = f'/Library/Frameworks/Python.framework/Versions/3.11/bin/python3 data/dataset_creation/extract_current_120.py {args.dataset} --output {out} --resume'
    if args.dataset == 'fairface':
        command += f' --archive {args.archive} --train-labels {args.train_labels} --val-labels {args.val_labels}'
    else:
        command += f' --gmdb-backup {args.gmdb_backup}'
    (out/'RUN_INSTRUCTIONS.md').write_text('# Offline extraction-only reproduction\n\nUse the same pinned files, model and environment. For a new run, choose a new ignored output directory and omit --resume. For restart:\n\n```sh\n'+command+'\n```\n\nThe extraction process makes no network requests and never downloads model weights. Inputs must exist. No credentials, original photos, thumbnails, landmarks or embeddings are exported.\n')


if __name__ == '__main__':
    def deny_network(event, args):
        if event in {'socket.connect', 'socket.getaddrinfo', 'http.client.connect', 'urllib.Request'}:
            raise RuntimeError('Extraction is offline; network requests are disabled')
    sys.addaudithook(deny_network)
    main()

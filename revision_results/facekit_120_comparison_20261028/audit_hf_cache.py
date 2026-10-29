import argparse
import ast
import hashlib
import io
import json
import os
import platform
import re
import stat
import subprocess
import sys
from collections import Counter
from pathlib import Path

for key in ['HF_DATASETS_OFFLINE', 'HF_HUB_OFFLINE', 'HF_HUB_DISABLE_TELEMETRY', 'DO_NOT_TRACK']:
    os.environ[key] = '1'
os.umask(0o077)


def offline(event, args):
    if event in {'socket.connect', 'socket.getaddrinfo', 'http.client.connect', 'urllib.Request'}:
        raise RuntimeError('Network access is disabled for this audit')


sys.addaudithook(offline)

import numpy as np
import pandas as pd
import pyarrow as pa
from PIL import Image, UnidentifiedImageError, __version__ as pillow_version

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'stats_testing'))
from prepare_cohorts import select_patients


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(8*1024*1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def normalized(value, field=''):
    if isinstance(value, (list, tuple, np.ndarray)):
        return json.dumps(sorted(str(item).strip() for item in value), ensure_ascii=False)
    if value is None or pd.isna(value):
        return None
    if field == 'hpo_terms' and isinstance(value, str):
        try:
            parsed = ast.literal_eval(value)
        except (ValueError, SyntaxError):
            parsed = None
        if isinstance(parsed, (list, tuple)):
            return normalized(parsed)
    if isinstance(value, (bool, np.bool_)):
        return str(bool(value)).lower()
    if isinstance(value, (int, float, np.number)):
        return format(value, '.15g')
    value = str(value).strip()
    if value.lower() in {'', 'nan', 'none', 'null', '<na>'}:
        return None
    if field in {'image_id', 'patient_id', 'internal_syndrome_id', 'person_type', 'syndrome_id', 'OMIM', 'age_year', 'age_month', 'distinctiveness'}:
        if re.fullmatch(r'-?\d+(\.\d+)?', value):
            return format(float(value), '.15g')
    return value


def compare_metadata(cache, historical, source, fields, prefix=''):
    rows, mismatches = [], []
    left = cache.set_index('_image_key')
    right = historical.set_index('_image_key')
    common = left.index.intersection(right.index, sort=False)
    for field in fields:
        oldfield = prefix+field
        if oldfield not in right:
            continue
        a = left.loc[common, field].map(lambda value: normalized(value, field))
        b = right.loc[common, oldfield].map(lambda value: normalized(value, field))
        equal = a.eq(b) | (a.isna() & b.isna())
        rows.append({'historical_source': source, 'cache_field': field, 'historical_field': oldfield,
                     'matched_image_rows': len(common), 'equal_normalized_n': int(equal.sum()),
                     'different_n': int((~equal).sum()), 'cache_missing_n': int(a.isna().sum()),
                     'historical_missing_n': int(b.isna().sum())})
        for key in common[~equal.to_numpy()]:
            mismatches.append({'historical_source': source, 'image_id': key, 'field': field,
                               'historical_field': oldfield, 'kind': 'normalized_metadata_difference'})
    return rows, mismatches


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=HERE/'local_only/hf_cache_audit_20261028')
    args = parser.parse_args()
    out = args.output.resolve()
    if HERE/'local_only' not in out.parents:
        raise ValueError('Private audit destination must stay in local_only')
    backup_manifest = json.loads((out/'backup_manifest.json').read_text())
    backup = out/'private_backup'
    results = out/'audit_results'
    results.mkdir(mode=0o700, exist_ok=False)
    source = Path(backup_manifest['source_cache'])
    before = {item['file']: sha(source/item['file']) for item in backup_manifest['files']}
    for item in backup_manifest['files']:
        assert before[item['file']] == item['sha256'] == sha(backup/item['file'])
    historical_paths = [ROOT/'data'/name for name in ['master_dataset_disease.csv', 'metadata_with_images.csv', 'phenotypes_all_disease.csv']]
    historical_before = {str(path.relative_to(ROOT)): sha(path) for path in historical_paths}
    info = json.loads((backup/'dataset_info.json').read_text())
    shards = sorted(backup.glob('*train-*.arrow'))
    assert len(shards) == len(info['splits']['train']['shard_lengths']) == 2
    frames, image_rows, schema_rows, summaries = [], [], [], []
    schema = None
    for shard_index, path in enumerate(shards):
        count, batches, byte_count, null_structs, embedded, path_values = 0, 0, 0, 0, 0, 0
        formats, modes, dimensions, image_errors = Counter(), Counter(), Counter(), Counter()
        with pa.memory_map(str(path), 'r') as handle:
            reader = pa.ipc.open_stream(handle)
            if schema is None:
                schema = reader.schema
                (results/'arrow_schema.txt').write_text(str(schema)+'\n')
                for field in schema:
                    schema_rows.append({'field': field.name, 'arrow_type': str(field.type), 'nullable': field.nullable,
                                        'hf_feature': json.dumps(info['features'].get(field.name, {}))})
            assert schema.equals(reader.schema, check_metadata=True), 'Shard schemas differ'
            for batch in reader:
                batch.validate(full=True)
                batches += 1
                metadata = batch.select([name for name in batch.schema.names if name != 'image_restored']).to_pandas()
                frames.append(metadata)
                column = batch.column(batch.schema.get_field_index('image_restored'))
                for row_index, item in enumerate(column.to_pylist()):
                    image_key = normalized(metadata.iloc[row_index]['image_id'], 'image_id')
                    data = item.get('bytes') if item else None
                    path_value = item.get('path') if item else None
                    null_structs += item is None
                    path_values += path_value is not None
                    embedded += bool(data)
                    byte_count += len(data) if data else 0
                    valid, reason, width, height, image_format, mode = False, '', None, None, None, None
                    if data:
                        try:
                            with Image.open(io.BytesIO(data)) as image:
                                image_format, mode = image.format, image.mode
                                width, height = image.size
                                image.verify()
                            with Image.open(io.BytesIO(data)) as image:
                                image.load()
                            valid = True
                            formats[image_format] += 1
                            modes[mode] += 1
                            dimensions[f'{width}x{height}'] += 1
                        except (UnidentifiedImageError, OSError, ValueError, Image.DecompressionBombError) as error:
                            reason = type(error).__name__
                            image_errors[reason] += 1
                    else:
                        reason = 'missing_embedded_bytes'
                        image_errors[reason] += 1
                    image_rows.append({'image_id': image_key, 'shard': path.name, 'shard_row': count+row_index,
                                       'embedded_bytes_n': len(data) if data else 0, 'path_value_present': path_value is not None,
                                       'byte_sha256': hashlib.sha256(data).hexdigest() if data else '',
                                       'decoded_ok': valid, 'width': width, 'height': height,
                                       'image_format': image_format, 'mode': mode, 'error_class': reason})
                count += batch.num_rows
        expected = info['splits']['train']['shard_lengths'][shard_index]
        assert count == expected, 'Shard count does not match manifest'
        summaries.append({'shard': path.name, 'rows': count, 'record_batches': batches,
                          'all_record_batches_validated': True, 'expected_rows': expected, 'null_image_structs': null_structs,
                          'embedded_images': embedded, 'non_null_paths': path_values, 'embedded_image_bytes': byte_count,
                          'decoded_ok': count-sum(image_errors.values()), 'image_error_counts': dict(image_errors),
                          'formats': dict(formats), 'modes': dict(modes), 'dimensions': dict(dimensions)})
        print(f'Shard {shard_index+1}: {count} records checked; {count-sum(image_errors.values())} images decode.', flush=True)
    cached = pd.concat(frames, ignore_index=True)
    cached['_image_key'] = cached.image_id.map(lambda value: normalized(value, 'image_id'))
    assert cached._image_key.notna().all() and cached._image_key.is_unique, 'Cache image keys missing or duplicated'
    images = pd.DataFrame(image_rows)
    images.to_csv(results/'private_image_integrity_manifest.csv', index=False)
    pd.DataFrame(schema_rows).to_csv(results/'schema.csv', index=False)
    metadata_rows, differences, coverage = [], [], []
    historical_frames = {}
    fields = [field for field in schema.names if field not in {'image_restored', 'image_id'}]
    for path in historical_paths:
        historical = pd.read_csv(path, dtype={'image_id': 'string', 'patient_id': 'string'}, low_memory=False)
        historical['_image_key'] = historical.image_id.map(lambda value: normalized(value, 'image_id'))
        assert historical._image_key.notna().all() and historical._image_key.is_unique
        historical_frames[path.name] = historical
        cache_keys, old_keys = set(cached._image_key), set(historical._image_key)
        common = cache_keys & old_keys
        coverage.append({'historical_source': path.name, 'historical_rows': len(historical), 'cache_rows': len(cached),
                         'overlapping_image_ids': len(common), 'historical_ids_missing_from_cache': len(old_keys-cache_keys),
                         'cache_ids_absent_from_historical': len(cache_keys-old_keys),
                         'same_order_if_same_set': cached._image_key.tolist() == historical._image_key.tolist(),
                         'overlap_with_decodable_image_n': int(images.loc[images.image_id.isin(common), 'decoded_ok'].sum())})
        direct, discordant = compare_metadata(cached, historical, path.name, fields)
        metadata_rows.extend(direct)
        differences.extend(discordant)
        prefixed, discordant = compare_metadata(cached, historical, path.name, fields, 'image_dataset_')
        metadata_rows.extend(prefixed)
        differences.extend(discordant)
        for key in sorted(cache_keys-old_keys):
            differences.append({'historical_source': path.name, 'image_id': key, 'field': '', 'historical_field': '', 'kind': 'cache_only_image_id'})
        for key in sorted(old_keys-cache_keys):
            differences.append({'historical_source': path.name, 'image_id': key, 'field': '', 'historical_field': '', 'kind': 'historical_only_image_id'})
    pd.DataFrame(coverage).to_csv(results/'historical_record_coverage.csv', index=False)
    pd.DataFrame(metadata_rows).to_csv(results/'metadata_comparison_counts.csv', index=False)
    pd.DataFrame(differences, columns=['historical_source', 'image_id', 'field', 'historical_field', 'kind']).to_csv(results/'private_discrepancy_keys.csv', index=False)
    disease = historical_frames['master_dataset_disease.csv']
    selected, _ = select_patients(disease.drop(columns='_image_key'), '3-9', 'pose')
    selected_keys = selected.image_id.map(lambda value: normalized(value, 'image_id'))
    joined = images.set_index('image_id').reindex(selected_keys)
    selected_summary = {'historical_selected_rows': len(selected), 'matched_cache_images': int(joined.decoded_ok.notna().sum()),
                        'decodable_selected_images': int(joined.decoded_ok.fillna(False).sum()),
                        'unmatched_selected_images': int(joined.decoded_ok.isna().sum())}
    patients = cached.patient_id.map(lambda value: normalized(value, 'patient_id'))
    duplicate_bytes = images.loc[images.byte_sha256.ne(''), 'byte_sha256'].value_counts()
    summary = {'status': 'OFFLINE_INPUT_AUDIT_COMPLETE', 'date': '2026-10-28', 'dataset': 'GestaltGANwglabshare/gmdb_enhanced_v1.1.0_new',
               'shards': summaries, 'total_rows': len(cached), 'expected_rows': info['splits']['train']['num_examples'],
               'schema_field_count': len(schema), 'matching_shard_schemas': True,
               'cached_unique_image_ids': cached._image_key.nunique(), 'cached_unique_known_patient_ids': patients.nunique(),
               'cached_missing_patient_id_rows': int(patients.isna().sum()),
               'cached_unique_syndrome_ids': int(cached.internal_syndrome_id.nunique()),
               'all_images_embedded': bool(images.embedded_bytes_n.gt(0).all()),
               'all_images_decode': bool(images.decoded_ok.all()), 'decoded_images': int(images.decoded_ok.sum()),
               'byte_duplicate_groups': int(duplicate_bytes.gt(1).sum()),
               'duplicate_byte_excess_rows': int((duplicate_bytes[duplicate_bytes.gt(1)]-1).sum()),
               'historical_coverage': coverage, 'selected_primary_coverage': selected_summary,
               'metadata_fields_with_differences': sum(row['different_n']>0 for row in metadata_rows),
               'metadata_field_comparisons': len(metadata_rows),
               'fairface_status': 'Not present in these GMDB shards; no FairFace images recovered by this audit',
               'historical_image_byte_identity': 'Unverified: historical source image hashes unavailable; ID/metadata coverage is not proof of identical image bytes or preprocessing',
               'pilot_status': 'NOT RUN; current request says no new analysis and separate image-reprocessing approval remains unconfirmed',
               'analysis_status': 'No FaceKit extraction, disease tests, comparisons of effects or new statistical analyses run'}
    assert len(cached) == info['splits']['train']['num_examples']
    dump(results/'audit_summary.json', summary)
    validation = []
    for item in backup_manifest['files']:
        source_path, backup_path = source/item['file'], backup/item['file']
        same = sha(source_path) == sha(backup_path) == before[item['file']] == item['sha256']
        assert same and source_path.stat().st_mtime_ns == item['source_mtime_ns']
        assert stat.S_IMODE(backup_path.stat().st_mode) == 0o600
        validation.append({'check': 'source_and_backup_integrity', 'file': item['file'], 'passed': same})
    for name, digest in historical_before.items():
        assert sha(ROOT/name) == digest
        validation.append({'check': 'historical_csv_unchanged', 'file': name, 'passed': True})
    assert stat.S_IMODE(out.stat().st_mode) == stat.S_IMODE(backup.stat().st_mode) == 0o700
    ignored = subprocess.run(['git', 'check-ignore', '--stdin'], cwd=ROOT,
        input='\n'.join(str(path.relative_to(ROOT)) for path in out.rglob('*') if path.is_file()),
        text=True, capture_output=True, check=True)
    expected_files = {str(path.relative_to(ROOT)) for path in out.rglob('*') if path.is_file()}
    assert set(ignored.stdout.splitlines()) == expected_files
    validation += [{'check': 'private_backup_permissions', 'file': 'private_backup', 'passed': True},
                   {'check': 'all_outputs_git_ignored', 'file': 'audit and private backup', 'passed': True},
                   {'check': 'shard_schema_and_expected_rows', 'file': 'both shards', 'passed': True}]
    pd.DataFrame(validation).to_csv(results/'validation_checks.csv', index=False)
    dump(results/'reproducibility_manifest.json', {'python': platform.python_version(), 'pyarrow': pa.__version__,
        'pandas': pd.__version__, 'numpy': np.__version__, 'Pillow': pillow_version, 'network': 'blocked by Python audit hook; HF offline flags enabled',
        'cache_before_and_after_sha256': before, 'historical_csv_before_and_after_sha256': historical_before,
        'audit_script_sha256': sha(Path(__file__)), 'user_helper_sha256': sha(ROOT/'data/dataset_creation/reloadHFCache.py'),
        'facekit_commit': subprocess.check_output(['git', '-C', str(ROOT/'FaceKit'), 'rev-parse', 'HEAD'], text=True).strip(),
        'image_handling': 'Embedded bytes verified and decoded in memory for integrity only; no images exported, transformed, displayed or extracted into FaceKit'} )
    lines = ['# Offline GMDB Arrow cache audit', '', '**Status: inputs recovered and inspected; FaceKit pilot NOT RUN.**', '',
        'The author requested no new analysis. This audit only validates storage, image decoding, record coverage and metadata. The original cache and historical PRISM files were preserved. A private independent backup was verified before decoding.', '',
        '## Shards and embedded images', '', '| Shard | Records | Embedded images | Decoded images |', '| --- | ---: | ---: | ---: |']
    lines += [f'| {row["shard"]} | {row["rows"]} | {row["embedded_images"]} | {row["decoded_ok"]} |' for row in summaries]
    lines += ['', f'Combined records: {len(cached):,}; expected: {info["splits"]["train"]["num_examples"]:,}. Both Arrow streams were fully traversed and every record batch validated. Shard schemas match, including Hugging Face metadata.', '',
        f'The schema has {len(schema)} fields. `image_restored` is `struct<bytes: binary, path: string>`, and its Hugging Face feature type is Image. All embedded: {summary["all_images_embedded"]}; all decode: {summary["all_images_decode"]}. This verifies encoded-image integrity, not clinical image quality or original-photo provenance.', '',
        'Per-shard formats, dimensions, modes, missingness and embedded byte totals are in audit_summary.json; the full field/type list is in schema.csv and arrow_schema.txt.', '',
        '## Historical record coverage', '', '| Historical file | Historical rows | Matched image IDs | Missing from cache | Cache-only IDs |', '| --- | ---: | ---: | ---: | ---: |']
    lines += [f'| {row["historical_source"]} | {row["historical_rows"]} | {row["overlapping_image_ids"]} | {row["historical_ids_missing_from_cache"]} | {row["cache_ids_absent_from_historical"]} |' for row in coverage]
    lines += ['', f'Historical selected childhood images: {selected_summary["historical_selected_rows"]:,}; cache matches: {selected_summary["matched_cache_images"]:,}; decodable matches: {selected_summary["decodable_selected_images"]:,}. The frozen selection function was reused solely to locate the existing selected image IDs; selection was not changed.', '',
        f'Cache distinct known patient IDs: {summary["cached_unique_known_patient_ids"]:,}; missing patient IDs: {summary["cached_missing_patient_id_rows"]}. Distinct syndrome IDs: {summary["cached_unique_syndrome_ids"]}. Exact-byte duplicate groups: {summary["byte_duplicate_groups"]}; excess duplicate-byte records: {summary["duplicate_byte_excess_rows"]}. No records were dropped.', '',
        '## Metadata checks', '', f'{summary["metadata_field_comparisons"]} field/source comparisons were made on matching image IDs; {summary["metadata_fields_with_differences"]} have at least one normalized difference. See metadata_comparison_counts.csv for numerators and denominators. Both direct and image_dataset_-prefixed historical fields are checked where available.', '',
        'Normalization treats equivalent numeric representations as equal, trims surrounding whitespace, standardizes null encodings, and compares HPO lists independent of order. It does not infer or relabel ancestry, diagnoses or patient identities. No metadata values or patient names are included in this aggregate report. Mismatch keys and image hashes remain in private files.', '',
        '## Limits and next step', '', 'These shards recover restored GMDB images, not verified original photographs. Matching image IDs and metadata do not establish that their encoded pixels exactly match the historical extraction inputs; historical input hashes and crop/restoration provenance are still needed. FairFace photographs are not contained in this GMDB cache.', '',
        'The earlier feasibility report is superseded specifically on GMDB image availability. Full 120-versus-125 analysis remains INDETERMINATE. The requested pilot is held: no FaceKit features or new statistical results were generated, and separate image-reprocessing approval remains unconfirmed. Confirm approval and pilot scope before a later extraction-only run.', '',
        '## Privacy and reproducibility', '', 'All work was offline, with network connections blocked in the audit process. No patient images were displayed or uploaded. Source and backup SHA-256 hashes and original modification times were checked before and after. Backup directory permissions are 0700 and backup files 0600; all generated files are Git-ignored. The user helper was reviewed and preserved without edits. Direct Arrow reading avoids remote dataset loaders and preserves the cache.', '',
        'Run with the existing system Python (pyarrow/datasets are already installed there):', '', '```sh',
        'HF_DATASETS_OFFLINE=1 HF_HUB_OFFLINE=1 HF_HUB_DISABLE_TELEMETRY=1 \\',
        '  /Library/Frameworks/Python.framework/Versions/3.11/bin/python3 \\',
        '  revision_results/facekit_120_comparison_20261028/audit_hf_cache.py \\',
        '  --output revision_results/facekit_120_comparison_20261028/local_only/hf_cache_audit_20261028',
        '```', '', 'The results directory must not already exist; keep the verified backup and source cache unchanged. The script validates an existing backup manifest and refuses to overwrite its audit results.', '']
    (results/'README.md').write_text('\n'.join(lines))
    print(json.dumps({'rows': len(cached), 'all_images_embedded': summary['all_images_embedded'],
        'all_images_decode': summary['all_images_decode'], 'historical_coverage': coverage,
        'selected_primary_coverage': selected_summary, 'metadata_fields_with_differences': summary['metadata_fields_with_differences'],
        'pilot_status': summary['pilot_status']}, indent=2), flush=True)


if __name__ == '__main__':
    main()

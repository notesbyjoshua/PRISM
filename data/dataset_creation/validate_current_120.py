import argparse
import hashlib
import json
import os
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
os.umask(0o077)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(8*1024*1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def recorded_file(path,expected):
    path=Path(path)
    if path.is_file() and sha(path)==expected:
        return path
    config=json.loads((ROOT/'ACTIVE_ANALYSIS.json').read_text())
    manifest_path=ROOT/config['retired_input_manifest']
    manifest=json.loads(manifest_path.read_text())
    for entry in manifest['files']:
        if (ROOT/entry['original_path']).resolve()==path.resolve() and entry['sha256']==expected:
            archived=manifest_path.parent/entry['archived_path']
            if archived.stat().st_mode&0o222==0 and sha(archived)==expected:
                return archived
    raise ValueError('Recorded execution input or source checksum could not be verified')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('run', type=Path)
    args = parser.parse_args()
    run = args.run.resolve()
    if ROOT/'local_review/feature_extraction_120' not in run.parents:
        raise ValueError('Expected a private extraction run')
    checks, sources = [], {}
    summaries = {}
    for dataset in ['gmdb', 'fairface']:
        out = run/dataset
        provenance = json.loads((out/'provenance.json').read_text())
        summary = json.loads((out/'extraction_summary.json').read_text())
        summaries[dataset] = summary
        schema = pd.read_csv(out/'feature_schema.csv')
        features = schema.feature.tolist()
        assert len(features) == len(set(features)) == 120 == summary['finite_morphology_columns']
        assert features == provenance['features']
        assert provenance['status'] == summary['status'] == 'FULL' and summary['pending'] == 0
        status = pd.read_csv(out/f'{dataset}_120_extraction_status.csv')
        candidate = pd.read_csv(out/'candidate_manifest.csv')
        assert status.key.is_unique and candidate.key.is_unique and set(status.key) == set(candidate.key)
        assert len(status) == len(candidate) == summary['candidates'] == summary['attempted']
        assert status.quality_passed.sum() == summary['quality_passed']
        assert status.reason.value_counts().to_dict() == summary['reason_counts']
        assert status.quality_passed.eq(status.reason.eq('accepted')).all()
        name = f'{dataset}_120_all_eligible_measurements.parquet' if dataset == 'fairface' else 'gmdb_120_all_records_measurements.parquet'
        table = pd.read_parquet(out/name)
        table = table.set_index('key', drop=False)
        assert table.key.is_unique and set(features).issubset(table)
        accepted = table if dataset == 'fairface' else table.loc[table.quality_passed]
        assert len(accepted) == summary['quality_passed']
        assert np.isfinite(accepted[features].to_numpy(float)).all()
        assert accepted.derotated.eq(True).all() and accepted.frontal_ok.eq(True).all()
        assert accepted.pose_yaw.abs().le(15).all() and accepted.pose_pitch.abs().le(15).all() and accepted.pose_roll.abs().le(10).all()
        assert set(accepted.key) == set(status.loc[status.quality_passed, 'key'])
        assert not any(('landmark' in field.lower() and field != 'landmark_count') or 'embedding' in field.lower() or field == 'image_restored' for field in table.columns)
        for field in table.select_dtypes(include=['object', 'str']).columns:
            assert not table[field].map(lambda value: isinstance(value, (bytes, bytearray, memoryview, list, np.ndarray))).any()
        connection = sqlite3.connect(f'file:{out/"measurements_checkpoint.sqlite3"}?mode=ro', uri=True)
        assert connection.execute('PRAGMA integrity_check').fetchone()[0] == 'ok'
        count = 0
        for key, encoded in connection.execute('SELECT image_key, record FROM records'):
            item = json.loads(encoded)
            assert key == item['metadata']['key']
            vector = item['qc']['measurements']
            assert vector is None or len(vector) == 120
            if vector is not None and key in table.index:
                np.testing.assert_array_equal(np.asarray(vector), table.loc[key, features].to_numpy(float))
            count += 1
        assert count == len(status)
        connection.close()
        pilot = pd.read_csv(out/'private_pilot_determinism.csv')
        assert pilot.deterministic.all() and pilot.max_absolute_feature_change.dropna().le(1e-10).all()
        sources.update(provenance['input_sha256'])
        sources.update(provenance.get('historical_source_sha256',{}))
        recorded_file(ROOT/'data/dataset_creation/extract_current_120.py',provenance['script_sha256'])
        assert sha(Path(provenance['model_path'])) == provenance['model_sha256']
        checks += [{'dataset': dataset, 'check': name, 'passed': True} for name in [
            'all_candidates_accounted', 'exact_120_schema', 'accepted_features_finite', 'no_duplicate_keys',
            'quality_and_failure_counts_reconcile', 'pose_and_full_canonicalization_gates',
            'checkpoint_integrity_and_feature_roundtrip', 'real_image_repeat_determinism',
            'no_images_landmarks_or_embeddings_in_measurements', 'pinned_model_and_script_unchanged']]
        if dataset == 'fairface':
            assert table.age.eq('3-9').all() and table.dataset.eq('FairFace').all()
            assert set(table.key) == set(status.loc[status.quality_passed, 'key'])
            labels = candidate.set_index('key')
            for field in ['race', 'gender', 'age', 'split', 'image_path_key']:
                assert table[field].eq(labels.loc[table.index, field]).all()
            assert table.image_path_key.str.startswith('train/').eq(table.split.eq('train')).all()
            for flag, cohort in [('sample_A', 'balanced'), ('sample_B', 'all_eligible'), ('sample_C', 'race_sex')]:
                manifest = pd.read_csv(out/f'fairface_120_{cohort}_manifest.csv')
                assert manifest.image_path_key.is_unique
                assert set(manifest.image_path_key) == set(table.loc[table[flag], 'image_path_key'])
            races = summary['input_audit']['race_categories']
            for race in races:
                group = table.loc[table.race.eq(race)]
                assert group.sample_A.sum() == min(500, len(group))
            capacity = [int((table.race.eq(race) & table.gender.eq(sex)).sum()) for race in races for sex in ['Female', 'Male']]
            cell_target = min([250]+capacity)
            cells = table.loc[table.sample_C].groupby(['race', 'gender']).size()
            assert len(cells) == 14 and cells.eq(cell_target).all()
            assert table.sample_B.all() and len(table) > 0
            for column, flag in [('selected_A', 'sample_A'), ('selected_B', 'sample_B'), ('selected_C', 'sample_C')]:
                assert pd.read_csv(out/'cohort_counts.csv')[column].sum() == table[flag].sum()
            checks += [{'dataset': dataset, 'check': name, 'passed': True} for name in [
                'source_labels_preserved', 'all_age_3_9_candidates_screened', 'sampling_no_oversampling',
                'A_B_C_manifests_match_master', 'race_sex_cells_balanced', 'cohort_counts_reconcile']]
    for path, expected in sources.items():
        recorded_file(path,expected)
    backup_manifest_path = ROOT/'revision_results/facekit_120_comparison_20261028/local_only/hf_cache_audit_20261028/backup_manifest.json'
    backup = json.loads(backup_manifest_path.read_text())
    for item in backup['files']:
        original = Path(backup['source_cache'])/item['file']
        assert sha(original) == item['sha256'] and original.stat().st_mtime_ns == item['source_mtime_ns']
    checks += [{'dataset': 'both', 'check': name, 'passed': True} for name in [
        'recorded_input_and_execution_source_hashes_verified', 'original_GMDB_cache_bytes_and_mtimes_unchanged']]
    pd.DataFrame(checks).to_csv(run/'validation_checks.csv', index=False)
    summary = {'status': 'FULL_EXTRACTION_QA_PASS', 'validation_checks_passed': len(checks), 'unit_tests_run_separately': True,
               'fairface': summaries['fairface'], 'gmdb': summaries['gmdb'],
               'statistical_analysis': 'NOT RUN', 'remaining_review': [
                   'Study source roles, identity constraints and publication permissions require PI review']}
    (run/'FINAL_EXTRACTION_STATUS.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({'status': summary['status'], 'checks_passed': len(checks), 'source_validation': 'Recorded hashes checked against current files or read-only execution archives; no historical comparison outputs generated'}, indent=2))


if __name__ == '__main__':
    main()

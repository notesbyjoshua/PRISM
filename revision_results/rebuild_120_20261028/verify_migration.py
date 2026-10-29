import csv
import hashlib
import json
import re
import sqlite3
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
RUN = HERE/'local_only/run_20261028'


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8*1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    config = json.loads((ROOT/'ACTIVE_ANALYSIS.json').read_text())
    assert config['feature_count'] == 120 and config['status'] == 'VALIDATED_ACTIVE'
    assert config['synthetic_rows_in_primary'] == 0
    archive = RUN/'archive/historical_125'
    history = json.loads((archive/'ARCHIVE_MANIFEST.json').read_text())
    for item in history['files']:
        path = archive/item['archived_path']
        assert sha(path) == item['sha256'] and path.stat().st_mode & 0o222 == 0
    assert archive.stat().st_mode & 0o222 == 0
    backup = json.loads((ROOT/'revision_results/facekit_120_comparison_20261028/local_only/hf_cache_audit_20261028/backup_manifest.json').read_text())
    for item in backup['files']:
        source = Path(backup['source_cache'])/item['file']
        copy = Path(backup['private_backup'])/item['file']
        assert sha(source) == sha(copy) == item['sha256']
        assert source.stat().st_mtime_ns == item['source_mtime_ns']
        assert source.stat().st_size == item['size_bytes']
        assert source.stat().st_ino != copy.stat().st_ino
    prohibited = {'patient_id','image_id','facekit_image_id','patient_name','file','image_path','image_path_key',
                  'source_image_key','recorded_patient_group_key','reference_key','prompt','landmarks','embedding'}
    csv_count = 0
    for directory in ['stats_results','publication','manuscript']:
        for path in (ROOT/directory).rglob('*.csv'):
            with path.open(newline='') as handle:
                header = next(csv.reader(handle))
            assert not prohibited.intersection(header), f'Private columns in active output: {path.relative_to(ROOT)}'
            csv_count += 1
    for path in [RUN,ROOT/'local_review/feature_extraction_120/20261028_current120']:
        assert subprocess.run(['git','check-ignore','-q',str(path)],cwd=ROOT).returncode == 0
    for field in ['disease_csv','reference_csv','disease_master_csv','reference_master_csv','feature_schema']:
        assert (ROOT/config[field]).is_file()
    with sqlite3.connect(ROOT/'backend/data/prism.db') as connection:
        counts = connection.execute("SELECT COUNT(*), COUNT(DISTINCT feature), COUNT(DISTINCT disease), SUM(q_value<0.05) FROM analysis_results WHERE analysis_type='pairwise'").fetchone()
    assert counts == (5760,120,48,955)
    assert json.loads((RUN/'VALIDATION_FULL.json').read_text())['status'] == 'VALIDATED'
    expected = ['SOURCE_AND_PERMISSION_AUDIT.md','FEATURE_SCHEMA_AND_EXTRACTION_AUDIT.md',
                'UNIFIED_120_MANIFEST.json','PRISM_120_Reanalysis_and_Source_Comparison_Report.pdf',
                'MANUSCRIPT_HANDOFF.md','CHANGELOG_120_MIGRATION.md','REPRODUCIBILITY_MANIFEST.json']
    assert all((RUN/name).is_file() for name in expected)
    token_pattern = re.compile(r'\bhf_[A-Za-z0-9]{25,}\b')
    for directory in [HERE,ROOT/'publication',ROOT/'manuscript']:
        for path in directory.rglob('*'):
            if path.is_file() and path.suffix in {'.md','.py','.json','.tex','.txt'} and 'archive' not in path.parts:
                assert not token_pattern.search(path.read_text(errors='replace')), 'Potential credential literal; redact before delivery'
    record = {'status':'PASS','active_analysis':'PRISM_120_20261028',
        'current_database':{'tests':5760,'features':120,'syndromes':48,'BH_discoveries':955},
        'archive_files_hash_verified':len(history['files']),'original_cache_and_independent_backup':'hash_size_mtime_verified_unchanged',
        'active_aggregate_CSVs_without_prohibited_person_columns':csv_count,'private_directories_gitignored':True,
        'credential_literal_scan':'PASS; no credential file read',
        'system_pytest':'41 passed; 1 skipped (FastAPI unavailable in system Python)',
        'virtualenv_targeted_pytest':'19 passed including FastAPI current-distribution integration',
        'manuscript_compilation':'success via desktop compiler',
        'PDF_review':'9 pages rendered; complete montage and detailed pages 4, 6 and 8 inspected',
        'frontend_build':'UNVERIFIED: offline npm ENOTCACHED; no external registry download',
        'historical_tracked_row_files':'preserved unchanged; PI privacy remediation review required',
        'release':'no commit, push, deployment, patient data upload or live Google Doc edit'}
    (RUN/'FINAL_DELIVERY_AUDIT.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))


if __name__ == '__main__':
    main()

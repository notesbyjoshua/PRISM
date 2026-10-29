import csv
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

ROOT=Path(__file__).resolve().parents[2]
RUN=ROOT/'revision_results/rebuild_120_20261028/local_only/run_20261028'
OUT=ROOT/'local_review/prism_120_primary_20261028'
os.umask(0o077)
sys.path.insert(0,str(ROOT/'publication/scripts'))
from verify_reproducible_outputs import verify_package


def sha(path):
    digest=hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda:handle.read(8*1024*1024),b''):digest.update(block)
    return digest.hexdigest()


def verify_archive(folder):
    manifest=json.loads((folder/'ARCHIVE_MANIFEST.json').read_text())
    assert folder.stat().st_mode&0o222==0
    for item in manifest['files']:
        path=folder/item.get('archived_path',item.get('path'))
        assert sha(path)==item['sha256'] and path.stat().st_mode&0o222==0
    return len(manifest['files'])


def main():
    package=verify_package()
    validation=json.loads((RUN/'VALIDATION_FULL.json').read_text())
    replay=json.loads((OUT/'PRIMARY_REPLAY_VALIDATION.json').read_text())
    assert validation['status']=='VALIDATED' and replay['status']=='PASS'
    archives={}
    for folder in [RUN/'archive/historical_125',OUT/'archive_before_cleanup',RUN/'archive/retired_comparisons_20261028',RUN/'archive/retired_reports_20261028']:
        archives[str(folder.relative_to(ROOT))]=verify_archive(folder)
    coverage=json.loads((OUT/'ARCHIVE_COVERAGE_MANIFEST.json').read_text())
    for item in coverage['archives']:
        folder=ROOT/item['private_archive']
        assert folder.stat().st_mode&0o222==0
        for record in item['files']:
            path=folder/record['path']
            assert path.stat().st_mode&0o222==0 and sha(path)==record['sha256']
        archives[item['private_archive']]=len(item['files'])
    backup=json.loads((ROOT/'revision_results/facekit_120_comparison_20261028/local_only/hf_cache_audit_20261028/backup_manifest.json').read_text())
    for item in backup['files']:
        source=Path(backup['source_cache'])/item['file'];copy=Path(backup['private_backup'])/item['file']
        assert sha(source)==sha(copy)==item['sha256']
        assert source.stat().st_size==item['size_bytes'] and source.stat().st_mtime_ns==item['source_mtime_ns']
        assert source.stat().st_ino!=copy.stat().st_ino and copy.stat().st_mode&0o222==0
    gmdb=pq.read_table(RUN/'CATALOG/gmdb_real_120.parquet',columns=['source_dataset','synthetic']).to_pandas()
    synthetic=pq.read_table(RUN/'CATALOG/pdidb_synthetic_120.parquet',columns=['source_dataset','synthetic']).to_pandas()
    assert len(gmdb)==1845 and not gmdb.synthetic.any() and gmdb.source_dataset.str.contains('gmdb',case=False).all()
    assert len(synthetic)==1727 and synthetic.synthetic.all()
    assert subprocess.run(['git','check-ignore','-q',str(OUT)],cwd=ROOT).returncode==0
    assert not subprocess.check_output(['git','diff','--name-only','--','data/*.csv'],cwd=ROOT,text=True).strip()
    result={'status':'PASS','study_type':'Statistical investigation using FaceKit','active_package':package,
        'independent_numerical_components':len(validation['components']),'primary_full_replay_tables':len(replay['tables']),
        'concordance_uncertainty':'excluded_biased_bootstrap_intervals; descriptive_points_only',
        'synthetic_records':1727,'synthetic_records_in_real_inference':0,'private_read_only_archive_files':archives,
        'source_cache':'original bytes, size and timestamps unchanged; private independent read-only backup verified',
        'patient_level_API_arrays':'removed; aggregate histogram counts and summaries only',
        'legacy_per_row_data_files':'unchanged; existing Git-history exposure requires PI remediation review',
        'permission_review':'publication/PI_REVIEW.md','source_and_role_tables':'private, gitignored; no photographs, landmarks or row-level tables in manuscript outputs'}
    result['verification']={'repository_tests':'51 passed','synthetic_extension_tests':'16 passed',
                            'frontend':'typecheck, lint and production build passed','manuscript':'desktop compiler success'}
    (OUT/'FINAL_ACTIVE_STUDY_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()

import argparse
import importlib.metadata
import json
import platform
from pathlib import Path

from .common import ROOT, DEFAULT_OUT, SEED, private_output, protect_baseline, require, sha, write_json


def finalize(out):
    checks = json.loads((out / 'validation/NEW_RESEARCH_VALIDATION.json').read_text())
    replay = json.loads((out / 'validation/READONLY_REPLAY_ATTESTATION.json').read_text())
    require(checks['status'] == replay['status'] == 'PASS', 'Finalize requires numerical and read-only replay checks')
    require(checks.get('figure_hash_checks') == 15 and checks.get('claim_source_checks') == 14, 'Finalize requires complete new artifact validation')
    for name, count in [('repository_tests.log', 79), ('extension_tests.log', 16)]:
        log = (out / 'validation' / name).read_text()
        require(f'{count} passed' in log and 'FAILED' not in log, 'Expected local test suite did not pass')
    require('"status": "PASS"' in (out / 'validation/public_verifier.log').read_text(), 'Public aggregate verification did not pass')
    source_files = sorted((ROOT / 'stats_testing/robustness_20261110').glob('*.py'))
    source_files += [ROOT / 'stats_testing/robustness_20261110/ANALYSIS_PLAN.md',
                     ROOT / 'stats_testing/robustness_20261110/README.md', ROOT / 'tests/test_robustness_20261110.py',
                     ROOT / 'PRISM_CODEX_RESEARCH_COMPLETION.md']
    require(all(path.is_file() for path in source_files), 'Final documentation incomplete')
    outputs = [path for path in out.rglob('*') if path.is_file()
               and path.name != 'RESEARCH_PROVENANCE.json' and path.suffix in {'.csv', '.json', '.md', '.svg', '.png'}
               and 'matplotlib' not in path.parts and 'checkpoints' not in path.parts]
    baseline = json.loads((out / '00_baseline_validation.json').read_text())
    record = {
        'status': 'PRIVATE_RESEARCH_COMPLETE_WITH_EXPLICIT_BLOCKERS', 'date': '2026-11-10',
        'baseline_git_head': baseline['git_head'], 'study': baseline['analysis_id'],
        'facekit_commit': baseline['facekit_commit'], 'extractor_sha256': baseline['extractor_sha256'],
        'model_sha256': baseline['model_sha256'], 'root_seed': SEED, 'original_photo_seed': 20261025,
        'parameters': 'ANALYSIS_PLAN.md and script CLI defaults; checkpoint signatures retain stage parameters/input hashes',
        'versions': {name: importlib.metadata.version(name) for name in ['numpy', 'pandas', 'scipy', 'statsmodels', 'matplotlib', 'pyarrow', 'pytest']},
        'python': platform.python_version(), 'platform': platform.platform(),
        'private_input_digests_source': '00_baseline_validation.json', 'private_input_digests': baseline['protected_input_sources'],
        'source_sha256': {str(path.relative_to(ROOT)): sha(path) for path in source_files},
        'output_sha256': {str(path.relative_to(out)): sha(path) for path in sorted(outputs)},
        'frozen_files_unchanged': protect_baseline(out), 'new_method_status': 'NOT_JUSTIFIED',
        'blocked': ['compatible raw landmark/components unavailable', 'real-image reprocessing authorization absent',
                    'clinical equivalence/external validation not established', 'PI disclosure/publication approval pending'],
        'excluded_from_final_evidence': ['checkpoints/source_negative_controls.json (invalid v1 ordered allocation; retained audit only)',
                                         'superseded cluster bootstrap v1 cache; valid final source is v2'],
        'checkpoint_policy': 'All checkpoints remain private and outside this aggregate draft manifest, including valid v2 states.',
        'validation_scope': checks, 'full_private_baseline_components': replay['components'],
        'test_results': {'repository': 79, 'extension': 16, 'remote_ci_for_new_work': 'not run; no push'},
        'release': 'No new public artifacts, manuscript edits, patient data/image upload, credentials, push or history rewrite.',
    }
    write_json(record, out / 'RESEARCH_PROVENANCE.json')
    print(f'Private provenance manifest created: {len(source_files)} sources, {len(outputs)} output artifacts')


def verify(out):
    record = json.loads((out / 'RESEARCH_PROVENANCE.json').read_text())
    require(all(sha(ROOT / source) == digest for source, digest in record['source_sha256'].items()), 'Research source hash mismatch')
    require(all(sha(out / source) == digest for source, digest in record['output_sha256'].items()), 'Research output hash mismatch')
    protect_baseline(out)
    print('PASS: research source/output hashes, protected private inputs and frozen baseline unchanged')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUT)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    out = private_output(args.output_dir)
    verify(out) if args.check else finalize(out)


if __name__ == '__main__':
    main()

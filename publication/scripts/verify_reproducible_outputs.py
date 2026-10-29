import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DIRECTORIES = ['publication/figures', 'publication/supplementary', 'publication/feature_definitions',
               'publication/feature_correlation', 'stats_results', 'manuscript/figures', 'manuscript/supplementary']
PERSON_COLUMNS = {'patient_id', 'image_id', 'facekit_image_id', 'patient_name', 'file', 'image_path',
                  'image_path_key', 'source_image_key', 'recorded_patient_group_key', 'reference_key',
                  'prompt', 'landmarks', 'embedding'}


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8*1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def artifact_paths(root):
    paths = [root/name for name in ['ACTIVE_ANALYSIS.json', 'manuscript/main.tex'] if (root/name).is_file()]
    for directory in DIRECTORIES:
        paths.extend(path for path in (root/directory).rglob('*') if path.is_file() and path.suffix in {'.csv', '.png', '.svg'})
    return sorted(path.relative_to(root).as_posix() for path in paths)


def check_public_csvs(root, paths):
    count = 0
    for name in paths:
        if name.endswith('.csv'):
            columns = pd.read_csv(root/name, nrows=0).columns
            if PERSON_COLUMNS.intersection(columns):
                raise ValueError(f'Person-level fields are prohibited in aggregate exports: {name}')
            count += 1
    return count


def verify_manifest(root, manifest):
    expected = manifest['files']
    if set(artifact_paths(root)) != set(expected):
        raise ValueError('Aggregate artifact inventory differs from the reviewed manifest')
    for name, digest in expected.items():
        path = root/name
        if path.is_symlink() or not path.resolve().is_relative_to(root.resolve()):
            raise ValueError('Aggregate manifest paths must stay inside the checkout')
        if sha(path) != digest:
            raise ValueError(f'Aggregate artifact checksum differs: {name}')
    return len(expected)


def bh(values):
    values = np.asarray(values, float)
    order = np.argsort(values, kind='stable')
    adjusted = np.empty(len(values))
    adjusted[order] = np.minimum(1, np.minimum.accumulate((values[order]*len(values)/np.arange(1, len(values)+1))[::-1])[::-1])
    return adjusted


def verify_package(root=ROOT, manifest_required=True):
    config = json.loads((root/'ACTIVE_ANALYSIS.json').read_text())
    if config['status'] != 'VALIDATED_ACTIVE' or config['feature_count'] != 120 or config['synthetic_rows_in_primary'] != 0:
        raise ValueError('A validated real-source current 120-feature analysis is required')
    expected_tests = config['eligible_syndromes']*config['feature_count']
    if expected_tests != config['primary_hypotheses']:
        raise ValueError('Primary family denominator disagrees with the active configuration')
    schema_path = root/'publication/supplementary/Supplementary_Table_S6_exact_feature_definitions.csv'
    schema = pd.read_csv(schema_path)
    features = schema.feature.tolist()
    if len(features) != 120 or len(set(features)) != 120 or schema.order.tolist() != list(range(120)):
        raise ValueError('The aggregate package must define exactly 120 ordered distinct features')
    commit = subprocess.check_output(['git', '-C', str(root/'FaceKit'), 'rev-parse', 'HEAD'], text=True).strip()
    if commit != config['facekit_commit'] or not schema.source_commit.eq(commit).all():
        raise ValueError('FaceKit checkout and current feature definitions have different source pins')
    extractor = root/'FaceKit/src/facekit/core/geometric/extractor.py'
    if not schema.extractor_sha256.eq(sha(extractor)).all():
        raise ValueError('The current extractor differs from the measured source')
    sys.path.insert(0, str(root/'FaceKit/src'))
    from facekit.core.geometric.extractor import GeometricFeatureExtractor
    if features != GeometricFeatureExtractor._init_feature_columns():
        raise ValueError('Published feature order differs from the current pinned extractor')
    primary_path = root/config['results_directory']/'stats_each_disease_vs_healthy.csv'
    primary = pd.read_csv(primary_path)
    if len(primary) != expected_tests or primary.duplicated(['disease', 'feature']).any():
        raise ValueError('Primary hypothesis count or uniqueness failed')
    if set(primary.feature) != set(features) or primary.disease.nunique() != config['eligible_syndromes']:
        raise ValueError('Primary feature or syndrome coverage failed')
    if not primary.groupby('disease').feature.nunique().eq(120).all():
        raise ValueError('Each eligible syndrome must retain the complete 120-feature family')
    if not primary.n_healthy.eq(config['reference_images']).all() or not primary.n_disease.ge(10).all():
        raise ValueError('Primary sample counts differ from the selected design')
    if not np.isfinite(primary[['mannwhitney_p', 'mannwhitney_q_fdr_global', 'hedges_g']]).all().all():
        raise ValueError('Main primary tests must have finite effects and probabilities')
    if not primary.mannwhitney_p.between(0, 1).all():
        raise ValueError('Primary p-values must be probabilities')
    np.testing.assert_allclose(primary.mannwhitney_q_fdr_global, bh(primary.mannwhitney_p), rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(primary.hedges_g, primary.cohens_d*(1-3/(4*(primary.n_disease+primary.n_healthy)-9)), rtol=1e-10, atol=1e-12)
    significant = int(primary.mannwhitney_q_fdr_global.lt(.05).sum())
    if significant != config['primary_significant']:
        raise ValueError('Primary significance count disagrees with the active configuration')
    for path in [root/'publication/supplementary/S1_primary_all_5760_hypotheses.csv',
                 root/'manuscript/supplementary/S1_primary_all_5760_hypotheses.csv']:
        expected=primary.rename(columns=lambda c:c.replace('healthy','reference')).copy()
        expected['comparison']=expected.comparison.str.replace('healthy','FairFace reference',regex=False)
        pd.testing.assert_frame_equal(expected, pd.read_csv(path), check_exact=False, rtol=1e-12, atol=1e-14)
    for path in [root/'publication/feature_definitions/current_120_feature_formulas.csv',
                 root/'manuscript/supplementary/Supplementary_Table_S6_exact_feature_definitions.csv']:
        if sha(path) != sha(schema_path):
            raise ValueError('Manuscript/publication feature definitions disagree')
    for figure in ['fig1_cross_syndrome_heatmap', 'fig2_clinical_indices', 'fig3_reference_demographics', 'fig4_photo_selection']:
        for extension in ['png', 'svg']:
            first, second = [root/f'{folder}/figures/{figure}.{extension}' for folder in ['publication', 'manuscript']]
            if sha(first) != sha(second):
                raise ValueError('Manuscript and publication figures disagree')
    for filename,planned in [('within_gmdb_results.csv',5760),('within_gmdb_current_eligible_sensitivity.csv',5760),
                             ('full_adjusted_screen.csv',5880),('targeted_pairwise.csv',22320),('age_stratified.csv',11760)]:
        frame=pd.read_csv(root/'stats_results/extensions_120'/filename)
        if len(frame)!=planned or not set(frame.feature).issubset(features):
            raise ValueError('An active extension has the wrong feature schema or planned family size')
        p=frame.correction_p if 'correction_p' in frame else frame.p_value
        expected=bh(p)
        supported=frame.status.eq('ok') if 'status' in frame else frame.estimable
        np.testing.assert_allclose(frame.loc[supported,'q_global_bh'],expected[supported],rtol=1e-10,atol=1e-12)
    for folder in ['stats_results/extensions_120','publication/supplementary','manuscript/supplementary']:
        concordance=pd.read_csv(root/folder/'concordance_summary.csv')
        if any('bootstrap' in column or '_ci_' in column for column in concordance):
            raise ValueError('Concordance bootstrap uncertainty must not appear in active exports')
        if not concordance.uncertainty_status.str.contains('descriptive_point_estimates_only').all():
            raise ValueError('Concordance must be explicitly descriptive')
    text=(root/'manuscript/main.tex').read_text().lower()
    if any(term in text for term in ['125-feature','historical','compatibility','rebuilt','[0.731','[0.725']):
        raise ValueError('Retired analysis content has entered the active manuscript')
    if 'statistical investigation' not in text or 'facekit' not in text:
        raise ValueError('PRISM must be framed as a statistical investigation built on FaceKit')
    for path in [root/'publication/historical_facekit',root/'stats_results/sensitivities/reference_fixed_historical']:
        if path.exists():
            raise ValueError('A retired analysis directory remains in the active package')
    paths = artifact_paths(root)
    count = check_public_csvs(root, paths)
    if manifest_required:
        manifest = json.loads((root/'publication/ARTIFACT_MANIFEST.json').read_text())
        if manifest['analysis_id'] != config['analysis_id'] or manifest['feature_count'] != 120:
            raise ValueError('Artifact manifest identifies a different analysis')
        verify_manifest(root, manifest)
    return {'status': 'PASS', 'analysis_id': config['analysis_id'], 'features': 120,
            'syndromes': config['eligible_syndromes'], 'primary_hypotheses': len(primary),
            'global_BH_discoveries': significant, 'aggregate_CSVs': count, 'aggregate_artifacts': len(paths),
            'scope': 'Committed aggregate integrity, current schema/source pins, primary BH and effect arithmetic; no private inputs or image reprocessing'}


def main():
    parser = argparse.ArgumentParser(description='Verify the current aggregate package without restricted local inputs')
    parser.add_argument('--write-manifest', action='store_true', help='Freeze reviewed aggregate artifact hashes after intentional regeneration')
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    result = verify_package(manifest_required=not args.write_manifest)
    if args.write_manifest:
        manifest = {'analysis_id': result['analysis_id'], 'feature_count': 120,
                    'scope': 'Reviewed code-derived aggregate working exports; excludes private data and reports',
                    'files': {name: sha(ROOT/name) for name in artifact_paths(ROOT)}}
        (ROOT/'publication/ARTIFACT_MANIFEST.json').write_text(json.dumps(manifest, indent=2)+'\n')
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        text = '# Current 120-feature aggregate verification\n\n'+''.join(f'- {key}: {value}\n' for key, value in result.items())
        text += '\nRestricted extraction and patient-level numerical replication run locally. This check does not establish image redistribution rights or publication approval.\n'
        args.report.write_text(text)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()

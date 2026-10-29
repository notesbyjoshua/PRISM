import argparse
import ast
import hashlib
import importlib.metadata
import importlib.util
import json
import platform
import subprocess
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, norm, rankdata

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CURRENT = '5764bdd4af9360b55a046a16579fb5a890a8182f'
HISTORICAL = 'e474628f3cb17a83edad59ac6e1d9e5438753032'
FROZEN = ROOT / 'revision_results/final_20260914'
SOURCE = HERE / 'local_only/sources' / f'FaceKit-{CURRENT}'
sys.path.insert(0, str(ROOT / 'stats_testing'))
from prepare_cohorts import deduplicate_reference, select_patients


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def bh(values):
    values = np.asarray(values, float)
    result = np.full(len(values), np.nan)
    valid = np.flatnonzero(np.isfinite(values))
    order = valid[np.argsort(values[valid], kind='stable')]
    result[order] = np.minimum(1, np.minimum.accumulate((values[order] * len(order) / np.arange(1, len(order)+1))[::-1])[::-1])
    return result


def independent_mw(x, y):
    nx, ny = len(x), len(y)
    ranks = rankdata(np.r_[x, y])
    u = ranks[:nx].sum() - nx*(nx+1)/2
    _, ties = np.unique(np.r_[x, y], return_counts=True)
    n = nx + ny
    variance = nx*ny/12 * (n+1 - ((ties**3-ties).sum()/(n*(n-1))))
    if variance == 0:
        return float(u), 1.
    z = (max(u, nx*ny-u) - nx*ny/2 - .5) / np.sqrt(variance)
    return float(u), min(1., float(2*norm.sf(z)))


def g(x, y):
    nx, ny = len(x), len(y)
    pooled = np.sqrt(((nx-1)*np.var(x, ddof=1)+(ny-1)*np.var(y, ddof=1))/(nx+ny-2))
    return (np.mean(x)-np.mean(y))/pooled * (1-3/(4*(nx+ny)-9)) if pooled else np.nan


def load_formula_probe():
    spec = importlib.util.spec_from_file_location('formula_tools', ROOT / 'publication/scripts/generate_feature_definitions.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def source_probe(path, current):
    module = load_formula_probe()
    source, records = path.read_text(), {}
    def record(key, value, function, line, end_line, expression, values):
        local = module.safe_locals(values)
        records[key] = {'function': function, 'line': line, 'expression': module.instantiate(expression, local)}
        return value
    tree = module.Recorder(source).visit(ast.parse(source))
    ast.fix_missing_locations(tree)
    namespace = {'_record': record, '__name__': 'audited_extractor'}
    exec(compile(tree, str(path), 'exec'), namespace)
    cls = namespace['GeometricFeatureExtractor']
    instance = cls.__new__(cls)
    rng = np.random.default_rng(20261028)
    landmarks = rng.uniform(0, 1000, (478, 2))
    canonical = instance._canonicalize(landmarks)
    features = {}
    for method in module.METHODS:
        features.update(getattr(instance, method)(canonical[0], canonical[1]))
    definitions = module.assignment_definitions(ast.parse(source))
    for feature, item in records.items():
        support = module.supporting_definitions(item['function'], item['line'], item['expression'], definitions)
        item['support'] = '; '.join(f'L{line} {name} = {value}' for name, line, value in support)
        item['formula'] = item['expression'] + ('; supporting definitions: '+item['support'] if support else '')
    return namespace, instance, records, features


def compatibility(out):
    sys.path.insert(0, str(SOURCE / 'src'))
    oldpath = ROOT / 'publication/historical_facekit/extractor.py'
    newpath = SOURCE / 'src/facekit/core/geometric/extractor.py'
    oldns, old, oldrecords, oldfeatures = source_probe(oldpath, False)
    newns, new, newrecords, newfeatures = source_probe(newpath, True)
    schema = pd.read_csv(FROZEN / 'exact_feature_names.csv').feature.tolist()
    assert list(oldfeatures) == schema and len(newfeatures) == 120
    assert np.isfinite(list(oldfeatures.values())).all() and np.isfinite(list(newfeatures.values())).all()
    definitions = pd.read_csv(ROOT / 'publication/feature_definitions/historical_feature_formulas.csv').set_index('feature')
    rows, formulas = [], []
    for feature in schema:
        match = feature if feature in newfeatures else ('forehead_height' if feature == 'hairline_height' else '')
        olddef = definitions.loc[feature]
        reason = 'Full pose correction with fixed canonical depth replaces roll-only alignment; full measurement equivalence is unverified.'
        if match != feature and match:
            reason += ' Historical hairline_height is anatomically merged into forehead_height; this is a many-to-one relation.'
        if not match:
            reason = 'Historical feature is absent from the current base schema; no direct comparable replacement verified.'
        rows.append({'old_feature': feature, 'new_feature': match,
                     'mapping_status': 'changed_definition' if match else 'removed',
                     'anatomical_region': olddef.anatomical_region, 'old_formula': olddef.full_formula_description,
                     'new_formula': newrecords[match]['formula'] if match else '',
                     'landmark_dependency_change': 'See exact indexed formulas and complete source excerpts; full-pose path also depends on fixed depth and transformation matrix.' if match else 'not applicable',
                     'normalization_change': 'See old/new formulas; distances and areas may inherit changed canonical coordinates.' if match else 'not applicable',
                     'alignment_change': 'roll-only -> SVD pure-rotation frontalization using fixed canonical depth; fallback roll-only',
                     'units_or_sign_change': 'angles in degrees; other outputs unitless; no sign harmonization applied',
                     'reason': reason, 'manual_review_required': True})
        formulas.append({'old_feature': feature, 'new_feature': match, 'old_source': f'{HISTORICAL}:extractor.py:L{oldrecords[feature]["line"]}',
                         'new_source': f'{CURRENT}:extractor.py:L{newrecords[match]["line"]}' if match else '',
                         'old_formula': oldrecords[feature]['formula'], 'new_formula': newrecords[match]['formula'] if match else ''})
    crosswalk = pd.DataFrame(rows)
    crosswalk.to_csv(out/'feature_crosswalk.csv', index=False)
    pd.DataFrame(formulas).to_csv(out/'formula_details.csv', index=False)
    pd.DataFrame({'order': np.arange(120), 'feature': list(newfeatures)}).to_csv(out/'current_feature_order.csv', index=False)
    for name, path in [('historical', oldpath), ('current', newpath)]:
        (out/f'{name}_extractor_source.txt').write_text(path.read_text())
    changes = []
    for name in load_formula_probe().METHODS + ['_canonicalize', '_compute']:
        a = next(n for n in ast.walk(ast.parse(oldpath.read_text())) if isinstance(n, ast.FunctionDef) and n.name == name)
        b = next(n for n in ast.walk(ast.parse(newpath.read_text())) if isinstance(n, ast.FunctionDef) and n.name == name)
        changes.append({'function': name, 'ast_changed': ast.dump(a, include_attributes=False) != ast.dump(b, include_attributes=False),
                        'historical_start': a.lineno, 'current_start': b.lineno})
    pd.DataFrame(changes).to_csv(out/'source_method_changes.csv', index=False)
    fallback = []
    for seed in range(10):
        lm = np.random.default_rng(seed).uniform(0, 1000, (478, 2))
        a, sa = old._canonicalize(lm)
        b, sb, flag = new._canonicalize(lm)
        assert np.allclose(a, b) and not flag
        af, bf = {}, {}
        for method in load_formula_probe().METHODS:
            af.update(getattr(old, method)(a, sa))
            bf.update(getattr(new, method)(b, sb))
        fallback.append({f: bool(np.isclose(af[f], bf[f], rtol=1e-10, atol=1e-10)) for f in bf})
    stable = [f for f in newfeatures if all(row[f] for row in fallback)]
    pd.DataFrame({'feature': list(newfeatures), 'equal_on_10_synthetic_roll_only_probes': [f in stable for f in newfeatures],
                  'interpretation': 'diagnostic probe only; not mathematical equivalence or image validation'}).to_csv(out/'synthetic_fallback_diagnostics.csv', index=False)
    return {'old_features': 125, 'new_features': 120, 'mapping_status_counts': dict(Counter(crosswalk.mapping_status)),
            'strict_old_n': 0, 'strict_old_denominator': 125, 'strict_new_n': 0, 'strict_new_denominator': 120,
            'strict_new_fraction': 0, 'strict_old_fraction': 0, 'anatomically_related_old_n': int(crosswalk.new_feature.ne('').sum()),
            'anatomically_related_new_n': crosswalk.loc[crosswalk.new_feature.ne(''), 'new_feature'].nunique(),
            'absent_old_names': sorted(set(oldfeatures)-set(newfeatures)), 'added_names': sorted(set(newfeatures)-set(oldfeatures)),
            'synthetic_fallback_equal_count': len(stable), 'finite_synthetic_new_count': 120,
            'interpretation': 'Strict coverage concerns the complete default measurement, including changed alignment. No source-verified full-pipeline equivalences were certified. Zero verified does not mean zero anatomical overlap.'}


def baseline(out):
    disease = pd.read_csv(ROOT/'data/master_dataset_disease.csv', dtype={'patient_id': 'string', 'image_id': 'string'}, low_memory=False)
    reference = pd.read_csv(ROOT/'data/master_dataset_healthy.csv', low_memory=False)
    primary, manifest = select_patients(disease, '3-9', 'pose')
    alternate, _ = select_patients(disease, '3-9', 'image_id')
    ref, _ = deduplicate_reference(reference)
    frozen = pd.read_csv(FROZEN/'age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz')
    rows = []
    for disease_name, group in frozen.groupby('disease', sort=True):
        subset = primary.loc[primary.disease.eq(disease_name)]
        for item in group.itertuples():
            x = pd.to_numeric(subset[item.feature], errors='coerce').to_numpy(float)
            y = pd.to_numeric(ref[item.feature], errors='coerce').to_numpy(float)
            x, y = x[np.isfinite(x)], y[np.isfinite(y)]
            u, p = mannwhitneyu(x, y, alternative='two-sided')
            ui, pi = independent_mw(x, y)
            effect = g(x, y)
            rows.append({'disease': disease_name, 'feature': item.feature, 'n_disease': len(x), 'n_healthy': len(y),
                         'u_recomputed': float(u), 'u_frozen': item.mannwhitney_u,
                         'g_recomputed': effect, 'p_recomputed': p, 'g_frozen': item.hedges_g,
                         'p_frozen': item.mannwhitney_p, 'q_frozen': item.mannwhitney_q_fdr_global,
                         'n_match': len(x) == item.n_disease and len(y) == item.n_healthy,
                         'independent_u_error': abs(ui-u), 'independent_p_error': abs(pi-p)})
    results = pd.DataFrame(rows)
    results['q_recomputed'] = bh(results.p_recomputed)
    results.to_csv(out/'historical_125_point_estimate_checks.csv.gz', index=False)
    errors = {name: float((results[f'{name}_recomputed']-results[f'{name}_frozen']).abs().max()) for name in ['u', 'g', 'p', 'q']}
    assert all(v < 1e-10 for v in errors.values()) and results.n_match.all()
    assert len(primary) == 1850 and len(ref) == 942 and ref.age.eq('3-9').all()
    assert len(results) == 6125 and results.disease.nunique() == 49 and (results.q_recomputed < .05).sum() == 811
    assert results.independent_u_error.max() == 0 and results.independent_p_error.max() < 1e-10
    selection_matches = []
    for name, frame in [('primary', primary), ('reference', ref)]:
        selected = frame.patient_id.astype(str).tolist() if name == 'primary' else frame.file.astype(str).tolist()
        selection_matches.append({'dataset': name, 'count': len(selected),
                                  'ordered_selection_sha256': hashlib.sha256(json.dumps(selected).encode()).hexdigest()})
    dump(out/'private_selection_hashes.json', selection_matches)
    checks = {'primary_patients': len(primary), 'reference_unique_images': len(ref), 'retained_syndrome_labels': primary.disease.nunique(),
              'eligible_syndromes': results.disease.nunique(), 'historical_hypotheses': len(results),
              'historical_significant': int((results.q_recomputed < .05).sum()), 'max_absolute_error': errors,
              'n_counts_match_all': bool(results.n_match.all()), 'independent_rank_checks': len(results),
              'max_independent_p_error': float(results.independent_p_error.max()), 'bootstrap_intervals': 'not rerun',
              'image_extraction': 'not rerun', 'baseline_status': 'PASS_POINT_ESTIMATES_ONLY',
              'disease_master_rows': len(disease), 'reference_master_rows': len(reference),
              'alternative_selected_images_changed': int(primary.set_index('patient_id').image_id.ne(alternate.set_index('patient_id').image_id).sum())}
    flow = pd.DataFrame([{'stage': 'GMDB input rows', 'n': len(disease)}, {'stage': 'Selected age 3-9 patients', 'n': len(primary)},
                         {'stage': 'Selected syndrome labels', 'n': primary.disease.nunique()}, {'stage': 'Eligible syndrome labels', 'n': 49},
                         {'stage': 'FairFace input rows', 'n': len(reference)}, {'stage': 'Distinct reference images', 'n': len(ref)}])
    flow.to_csv(out/'cohort_flow.csv', index=False)
    manifest.selection_reason.value_counts().rename_axis('reason').reset_index(name='n').to_csv(out/'cohort_exclusions.csv', index=False)
    missing = []
    for label, frame in [('GMDB_selected', primary), ('FairFace_selected', ref)]:
        features = pd.read_csv(FROZEN/'exact_feature_names.csv').feature.tolist()
        missing.append({'source': label, 'n_rows': len(frame), 'missing_feature_cells': int(frame[features].isna().sum().sum()),
                        'total_feature_cells': len(frame)*125, 'missing_gender_rows': int(frame.gender.isna().sum())})
    pd.DataFrame(missing).to_csv(out/'historical_missingness.csv', index=False)
    available = []
    for label, frame, column in [('GMDB_selected', primary, 'facekit_relative_path'), ('FairFace_selected', ref, 'file')]:
        found = 0
        for value in frame[column].dropna().astype(str):
            path = Path(value).expanduser()
            candidates = [path] if path.is_absolute() else [ROOT/path, ROOT/'data'/path]
            found += any(p.is_file() for p in candidates)
        available.append({'source': label, 'selected_n': len(frame), 'referenced_files_resolvable_n': int(found),
                          'scope': 'Recorded absolute locations or repository/data-relative locations; no external-folder search'})
    pd.DataFrame(available).to_csv(out/'input_access_summary.csv', index=False)
    highlights = [('Cornelia de Lange', 'philtrum_length'), ('Noonan', 'inter_pupillary_distance'),
                  ('WILLIAMS-BEUREN', 'eb_thickness_mean'), ('WILLIAMS-BEUREN', 'outer_canthal_distance'),
                  ('Rubinstein-Taybi', 'malar_bulge_mean')]
    clinical = []
    for syndrome, feature in highlights:
        item = frozen.loc[frozen.disease.str.contains(syndrome, case=False, regex=False) & frozen.feature.eq(feature)].iloc[0]
        clinical.append({'syndrome': syndrome, 'old_feature': feature, 'new_feature': feature,
                         'g_125': item.hedges_g, 'q_125_global': item.mannwhitney_q_fdr_global,
                         'g_120': None, 'q_120_global': None, 'status': 'manual_review_changed_definition; 120 not run',
                         'reason': 'Local highlighted formulas retain their definitions on common canonical coordinates, but global alignment changes. No comparable image-level rerun.'})
    pd.DataFrame(clinical).to_csv(out/'clinical_highlights_comparison.csv', index=False)
    dump(out/'baseline_audit.json', checks)
    return checks


def privacy_and_inputs(out):
    tracked = subprocess.check_output(['git', 'ls-files', 'data'], cwd=ROOT, text=True).splitlines()
    audit = []
    for file in tracked:
        path = ROOT/file
        if path.suffix != '.csv':
            continue
        columns = list(pd.read_csv(path, nrows=0))
        categories = []
        for category, terms in {'identifiers': ['patient_id', 'image_id', 'patient_name'], 'locators': ['file', 'path'],
                                'demographics': ['age', 'gender', 'ethnicity', 'race'], 'genetic_metadata': ['hgvs', 'gene', 'hpo']}.items():
            if any(any(term in column.lower() for term in terms) for column in columns):
                categories.append(category)
        audit.append({'repository_file': file, 'sha256': sha(path), 'column_count': len(columns),
                      'risk_categories': '; '.join(categories), 'human_review_required': True,
                      'action': 'Review redistribution permissions; no file deletion or Git history change performed'})
    pd.DataFrame(audit).to_csv(out/'tracked_data_privacy_audit.csv', index=False)
    originals, landmarks, assets = [], [], []
    for folder in [ROOT/'data', ROOT/'local_review', ROOT/'FaceKit', SOURCE]:
        for path in folder.rglob('*'):
            if not path.is_file() or any(part in {'node_modules', '.venv', '.git', '__pycache__', 'test_deps', 'vendor'} for part in path.parts):
                continue
            if path.suffix.lower() in {'.jpg', '.jpeg', '.webp'}:
                originals.append(path)
            if path.suffix.lower() in {'.jsonl', '.npy', '.npz'}:
                landmarks.append(path)
            if path.suffix.lower() in {'.task', '.onnx', '.pth'}:
                assets.append(path)
    dump(out/'private_candidate_inventory.json', {'image_candidates': [str(p) for p in originals],
                                                  'landmark_candidates': [str(p) for p in landmarks],
                                                  'model_assets': [{'path': str(p), 'sha256': sha(p)} for p in assets]})
    summary = {'tracked_data_csvs_flagged': len(audit), 'original_image_candidate_files': len(originals),
               'landmark_array_or_jsonl_candidates': len(landmarks), 'repository_model_asset_files': len(assets),
               'authorization': 'Image reprocessing permission unconfirmed by author on 2026-10-28; extraction blocked.',
               'scope': 'data/, local_review/, FaceKit/ and pinned source. PNGs are reporting/test assets; no claim about other local folders.'}
    dump(out/'input_inventory_summary.json', summary)
    model = Path.home()/'.cache/facekit/face_landmarker.task'
    model_info = {'asset': 'face_landmarker.task', 'cached_asset_available': model.exists(),
                  'sha256': sha(model) if model.exists() else None,
                  'use_status': 'not executed; historical use of this asset is unverified',
                  'declared_url': 'https://storage.googleapis.com/mediapipe-models/face_landmarker/face_landmarker/float16/1/face_landmarker.task'}
    dump(out/'model_asset_audit.json', model_info)
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=HERE/'local_only/run_20261028')
    args = parser.parse_args()
    out = args.output.resolve()
    if HERE/'local_only' not in out.parents:
        raise ValueError('Outputs must remain inside this workflow local_only directory')
    out.mkdir(parents=True, exist_ok=False)
    source_paths = [ROOT/'data/master_dataset_disease.csv', ROOT/'data/master_dataset_healthy.csv',
                    ROOT/'publication/historical_facekit/extractor.py', ROOT/'FaceKit/src/facekit/core/geometric/extractor.py',
                    SOURCE/'src/facekit/core/geometric/extractor.py', SOURCE/'src/facekit/core/geometric/canonical_depth.py',
                    SOURCE/'pyproject.toml', FROZEN/'publication_manifest.json',
                    FROZEN/'age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz']
    source_paths += list((ROOT/'stats_testing').glob('*.py'))
    source_paths += [p for p in FROZEN.rglob('*') if p.is_file()]
    source_paths += list(HERE.glob('*.py'))
    before = {str(p.relative_to(ROOT)): sha(p) for p in source_paths}
    mapping = compatibility(out)
    control = baseline(out)
    inventory = privacy_and_inputs(out)
    metrics = {'decision': 'INDETERMINATE', 'date': '2026-10-28', 'mapping': mapping,
               'direction_agreement_meaningful': {'status': 'not_run', 'numerator': None, 'denominator': None, 'fraction': None, 'wilson_95': None},
               'effect_similarity_meaningful': {'status': 'not_run', 'numerator': None, 'denominator': None, 'fraction': None},
               'significance_retention': {'status': 'not_run', 'numerator': None, 'denominator': None, 'fraction': None},
               'threshold': .90, 'quality_gate': 'BLOCKED: authorized original image/landmark inputs unavailable',
               'clinical_gate': 'MANUAL_REVIEW: all five highlighted measurements have changed full-pipeline definitions',
               'cross_analysis': 'not_run', 'baseline': control, 'inventory': inventory,
               'missing_metric_policy': 'Undefined fractions are null, never 0/0 or 0% agreement.'}
    dump(out/'comparison_metrics.json', metrics)
    pd.DataFrame([{'scope': 'all verified comparable image-level hypotheses', 'status': 'not_run',
                   'n_discordant': None, 'denominator': None, 'reason': 'No authorized 120-feature extraction; strict full-pipeline mapping not certified'}]).to_csv(out/'discordance_summary.csv', index=False)
    (out/'analysis_120').mkdir()
    (out/'analysis_120/README.md').write_text('# 120-feature analysis: NOT RUN\n\nAuthorized image/landmark inputs and a paper-linked code/model pin are required. No statistical values were inferred from historical columns.\n')
    for path in source_paths:
        assert before[str(path.relative_to(ROOT))] == sha(path), 'Source or frozen input changed during audit'
    versions = {}
    for name in ['numpy', 'pandas', 'scipy', 'statsmodels', 'mediapipe', 'opencv-contrib-python', 'pytest']:
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = 'not installed'
    manifest = {'prism_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                'historical_facekit': HISTORICAL, 'current_facekit_audited': CURRENT,
                'application_facekit': subprocess.check_output(['git', '-C', str(ROOT/'FaceKit'), 'rev-parse', 'HEAD'], text=True).strip(),
                'paper': 'https://arxiv.org/html/2610.09288v1', 'paper_commit_link': 'unverified; candidate source only',
                'python': platform.python_version(), 'platform': platform.platform(), 'versions': versions,
                'seed_synthetic_schema_probe': 20261028, 'source_sha256': before,
                'frozen_files_unchanged': True, 'source_acquisition': json.loads((HERE/'local_only/source_acquisition.json').read_text()),
                'author_input_permission': 'Unconfirmed; extraction not authorized', 'execution_scope': 'source-only feature probes and historical point-estimate reconstruction'}
    dump(out/'reproducibility_manifest.json', manifest)
    (out/'environment_freeze.txt').write_text(subprocess.check_output([sys.executable, '-m', 'pip', 'freeze', '--all'], text=True))
    print(json.dumps({'decision': metrics['decision'], 'baseline': control, 'mapping': mapping, 'inventory': inventory}, indent=2))


if __name__ == '__main__':
    main()

import argparse
import hashlib
import json
import os
import re
import shutil
import stat
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

os.umask(0o077)
ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
EXTRACTION = ROOT/'local_review/feature_extraction_120/20261028_current120'
sys.path.insert(0, str(ROOT/'stats_testing'))
from prepare_cohorts import select_patients, deduplicate_reference
from cohort_audit import read_mapping_functions


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8*1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, default=str)+'\n')


def archive(run):
    out = run/'archive/historical_125'
    if out.exists():
        manifest = json.loads((out/'ARCHIVE_MANIFEST.json').read_text())
        assert all(sha(out/item['archived_path']) == item['sha256'] for item in manifest['files'])
        return
    out.mkdir(parents=True, mode=0o700)
    sources = [ROOT/path for path in ['stats_results', 'pds_results_v2', 'publication', 'manuscript',
        'revision_results/final_20260914', 'local_review/final_reanalysis_20260914', 'stats_testing']]
    files = []
    for folder in sources:
        files += [path for path in folder.rglob('*') if path.is_file() and '__pycache__' not in path.parts and path.name != '.DS_Store']
    files += list((ROOT/'data').glob('*.csv'))
    files += [ROOT/'data/dataset_creation/fairface_facekit_pipeline.py', ROOT/'data/dataset_creation/format_fairface.py']
    followup = ROOT/'revision_results/reviewer_followup_2026'
    files += [path for path in followup.iterdir() if path.is_file() and path.suffix in {'.py', '.md', '.csv', '.json', '.txt', '.pdf'}]
    if (followup/'figures').exists():
        files += [p for p in (followup/'figures').rglob('*') if p.is_file()]
    rows = []
    for path in sorted(set(files)):
        relative = path.relative_to(ROOT)
        target = out/relative
        target.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
        digest, before = sha(path), path.stat()
        shutil.copy2(path, target)
        assert sha(target) == digest == sha(path)
        rows.append({'source_path': str(relative), 'archived_path': str(relative), 'size_bytes': before.st_size,
                     'mtime_ns': before.st_mtime_ns, 'sha256': digest})
    manifest = {'status': 'verified_read_only_copy', 'git_head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'historical_extractor_commit': 'e474628f3cb17a83edad59ac6e1d9e5438753032',
        'scope': 'Historical master CSVs, original extraction/cohort/statistical code, frozen final package, full local final reanalysis, current historical statistics/PDS/publication/manuscript outputs, and reviewer follow-up aggregate outputs/source. Redundant clone/test-cache copies remain untouched in their existing locations.',
        'files': rows}
    write_json(out/'ARCHIVE_MANIFEST.json', manifest)
    for path in out.rglob('*'):
        path.chmod(0o500 if path.is_dir() else 0o400)
    out.chmod(0o500)
    print(f'Read-only historical archive verified: {len(rows)} files', flush=True)


def harmonize(frame, reference=False):
    maps = read_mapping_functions(ROOT/'stats_testing/ancestry_analysis.py')
    frame = frame.copy()
    frame['sex_label'] = frame.gender.astype('string').str.strip().str.lower().map({'male': 'Male', 'm': 'Male', 'female': 'Female', 'f': 'Female'})
    if reference:
        frame['broad_group'] = frame.race.map(maps['map_fairface_layer1'])
        frame['detailed_group'] = frame.race.map(maps['map_fairface_layer2'])
    else:
        frame['broad_group'] = frame.image_dataset_ethnicity_category.map(maps['map_gmdb_layer1'])
        frame['detailed_group'] = [maps['map_gmdb_layer2'](a, b) for a,b in zip(frame.image_dataset_ethnicity_category, frame.image_dataset_ethnicity_sub_category)]
    return frame


def save_cohort(out, name, frame, features, reference=False):
    frame.to_parquet(out/f'{name}.parquet', index=False)
    columns = ['image_id', 'frontal_ok', 'pose_yaw', 'pose_pitch', 'pose_roll']+features
    phenotype = frame[columns].copy()
    phenotype.insert(0, 'disease', 'healthy' if reference else frame.disease)
    phenotype.to_csv(out/f'{name}_phenotypes.csv', index=False)
    metadata = [field for field in ['image_id', 'patient_id', 'file', 'disease', 'internal_syndrome_name',
        'image_dataset_ethnicity_category', 'image_dataset_ethnicity_sub_category', 'gender', 'race', 'age',
        'age_year', 'age_month', 'conditional_age_years', 'broad_group', 'detailed_group', 'sex_label',
        'frontal_ok', 'pose_yaw', 'pose_pitch', 'pose_roll'] if field in frame]
    master = frame[metadata].copy()
    for column in features:
        master[column] = frame[column].to_numpy()
    master.to_csv(out/f'{name}_master.csv', index=False)


def prepare(run):
    out = run/'COHORTS'
    out.mkdir(exist_ok=True)
    features = pd.read_csv(EXTRACTION/'gmdb/feature_schema.csv').feature.tolist()
    assert len(features) == 120 == len(set(features))
    for source in ['gmdb', 'fairface']:
        evidence = json.loads((EXTRACTION/source/'provenance.json').read_text())
        assert evidence['status'] == 'FULL' and evidence['features'] == features
        assert evidence['facekit_commit'] == '5764bdd4af9360b55a046a16579fb5a890a8182f'
    new = pd.read_parquet(EXTRACTION/'gmdb/gmdb_120_all_records_measurements.parquet')
    metadata = pd.read_csv(ROOT/'data/metadata_with_images.csv', dtype={'image_id':'string','patient_id':'string'}, low_memory=False)
    old = pd.read_csv(ROOT/'data/master_dataset_disease.csv', dtype={'image_id':'string','patient_id':'string'}, low_memory=False)
    labels = pd.read_csv(ROOT/'data/syndrome_folder_map.csv', dtype={'internal_syndrome_id':'string'})
    assert labels.internal_syndrome_id.is_unique
    mapping = labels.set_index('internal_syndrome_id').facekit_cohort.to_dict()
    raw = metadata[['image_id', 'patient_id', 'internal_syndrome_id', 'internal_syndrome_name', 'age_year', 'age_month',
        'gender', 'image_dataset_ethnicity_category', 'image_dataset_ethnicity_sub_category']].copy()
    raw['disease'] = raw.internal_syndrome_id.astype('string').map(mapping)
    assert raw.disease.notna().all()
    new['image_id'] = new.image_id.astype('string')
    raw = raw.merge(new[['image_id','quality_passed','reason','frontal_ok','pose_yaw','pose_pitch','pose_roll']+features],
        on='image_id', validate='one_to_one', how='left')
    assert len(raw) == 10980 and raw.reason.notna().all()
    raw['frontal_ok'] = raw.quality_passed
    primary, manifest = select_patients(raw, '3-9', 'pose')
    alternate, alt_manifest = select_patients(raw, '3-9', 'image_id')
    primary, alternate = harmonize(primary), harmonize(alternate)
    raw = harmonize(raw)
    assert primary.patient_id.is_unique and set(primary.patient_id) == set(alternate.patient_id)
    old_primary, _ = select_patients(old, '3-9', 'pose')
    old_alternate, _ = select_patients(old, '3-9', 'image_id')
    fair = pd.read_parquet(EXTRACTION/'fairface/fairface_120_all_eligible_measurements.parquet')
    fair['file'] = fair.image_path_key
    fair['image_id'] = fair.image_path_key
    fair = harmonize(fair, True)
    references = {'balanced': fair.loc[fair.sample_A].copy(), 'all_eligible': fair.copy(), 'race_sex': fair.loc[fair.sample_C].copy()}
    historical_new = pd.read_parquet(EXTRACTION/'fairface/fairface_120_historical942_measured_with_qc.parquet')
    historical_new['file'] = historical_new.image_path_key
    historical_new['image_id'] = historical_new.image_path_key
    references['fixed_historical'] = harmonize(historical_new.loc[historical_new.quality_passed].copy(), True)
    historic_ref, _ = deduplicate_reference(pd.read_csv(ROOT/'data/master_dataset_healthy.csv'))
    old_features = pd.read_csv(ROOT/'revision_results/final_20260914/exact_feature_names.csv').feature.tolist()
    old_primary_ids = set(old_primary.image_id)
    paired_new = harmonize(raw.loc[raw.image_id.isin(old_primary_ids) & raw.quality_passed].copy())
    paired_old = old_primary.loc[old_primary.image_id.isin(paired_new.image_id)].copy()
    old_ref_key = historic_ref.file.astype('string').str.replace('\\','/',regex=False).str.split('/').str[-2:].str.join('/')
    paired_old_ref = historic_ref.loc[old_ref_key.isin(references['fixed_historical'].image_path_key)].copy()
    paired_old_ref['image_id'] = old_ref_key.loc[paired_old_ref.index]
    assert set(paired_old.image_id) == set(paired_new.image_id) and len(paired_old_ref) == len(references['fixed_historical'])
    for name, frame in [('primary', primary), ('alternate', alternate), ('fixed_common_new', paired_new)]:
        save_cohort(out, name, frame, features)
    save_cohort(out, 'fixed_common_old', paired_old, old_features)
    save_cohort(out, 'fixed_common_reference_old', paired_old_ref, old_features, True)
    for name, ref in references.items():
        assert ref.file.is_unique and ref.age.eq('3-9').all()
        save_cohort(out, 'reference_'+name, ref, features, True)
    raw.to_parquet(out/'all_gmdb_candidates.parquet', index=False)
    manifest.to_csv(out/'private_primary_selection_manifest.csv', index=False)
    alt_manifest.to_csv(out/'private_alternate_selection_manifest.csv', index=False)
    targets = sorted(primary.disease.value_counts().loc[lambda counts: counts >= 10].index)
    original_targets = sorted(pd.read_csv(ROOT/'revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz', usecols=['disease']).disease.unique())
    counts = primary.disease.value_counts().rename_axis('disease').reset_index(name='n_patients')
    counts['eligible'] = counts.n_patients.ge(10)
    counts.to_csv(out/'syndrome_eligibility.csv', index=False)
    selected = old_primary[['patient_id','image_id']].merge(primary[['patient_id','image_id']], on='patient_id', how='outer', suffixes=('_old','_new'), indicator=True)
    selected['photo_changed'] = selected.image_id_old.ne(selected.image_id_new) & selected['_merge'].eq('both')
    selected.to_csv(out/'private_old_new_photo_selection.csv', index=False)
    plan = {'status':'cohort_rules_frozen_before_inference','date':'2026-10-28','primary_patient_n':len(primary),
        'alternate_patient_n':len(alternate),'eligible_syndromes':len(targets),'targets':targets,
        'historical_targets':original_targets,'newly_eligible':sorted(set(targets)-set(original_targets)),
        'no_longer_eligible':sorted(set(original_targets)-set(targets)), 'primary_planned_hypotheses':len(targets)*120,
        'reference_n':{name:len(frame) for name,frame in references.items()},'fixed_common_gmdb_n':len(paired_new),
        'fixed_common_reference_n':len(paired_old_ref),'primary_photo_changes_in_shared_patients':int(selected.photo_changed.sum()),
        'lost_historical_patients':int(selected['_merge'].eq('left_only').sum()),
        'new_current_patients':int(selected['_merge'].eq('right_only').sum()),
        'alternative_photo_changes':int(primary.set_index('patient_id').image_id.ne(alternate.set_index('patient_id').image_id).sum()),
        'primary_rule':'Original age/ID/label/group maps; all cached real GMDB candidates; current QC; one eligible image per patient using current normalized pose score then lexical image ID',
        'fixed_photo_benchmark':'Intersection of historical chosen photos and current QC; recompute BOTH old 125 and new 120 on exactly those images, with the same 941 current-QC historical references',
        'missingness_policy':'Main tests n>=10, <=30% missing in either group; nonestimable cells retained separately; tested-family BH matches historical finite-p policy',
        'bootstrap_settings':{'main':500,'rank':200,'main_seed':2026,'demographic_B':1000,'demographic_seed':42,
            'within_gmdb':1000,'photo_replicates':100,'followup_seed':20261024,'extension':1000,'concordance':1000,'extension_seed':20261026},
        'sources_excluded_from_primary':['PDIDB'], 'reference_designs_overlap':True,
        'authorizations':'Author-provided mentor/IRB scope attestation; no independent review of underlying IRB paperwork or publication permission',
        'analysis_execution_scope':'No synthetic pooling; no public patient-level exports; no live Google Doc edits'}
    write_json(run/'ANALYSIS_PLAN.json',plan)
    pd.DataFrame({'feature':features}).to_csv(out/'exact_feature_names.csv',index=False)
    dictionary = pd.read_csv(run/'archive/historical_125/publication/supplementary/Supplementary_Table_S6_exact_feature_definitions.csv')
    dictionary = dictionary.set_index('feature').loc[features].reset_index()
    dictionary['definition_status']='Current formula pointers in extracted feature_schema.csv; historical anatomical labels retained for grouping only'
    dictionary.to_csv(out/'feature_anatomy_grouping.csv',index=False)
    print(json.dumps({k:plan[k] for k in ['primary_patient_n','eligible_syndromes','primary_planned_hypotheses','reference_n','primary_photo_changes_in_shared_patients','lost_historical_patients','new_current_patients']},indent=2),flush=True)


def inventory(run):
    roots = ['stats_testing','publication','manuscript','backend','frontend','revision_results/final_20260914','revision_results/reviewer_followup_2026']
    pattern = re.compile(r'\b125\b|6,125|\b6125\b|\b811\b|\b379\b|72/125|89/125')
    rows=[]
    for name in roots:
        for path in (ROOT/name).rglob('*'):
            if not path.is_file() or path.suffix not in {'.py','.md','.tex','.ts','.tsx','.json'} or any(p in {'node_modules','.next','local_only','__pycache__'} for p in path.parts):
                continue
            for line,text in enumerate(path.read_text(errors='replace').splitlines(),1):
                if pattern.search(text):
                    relative=str(path.relative_to(ROOT))
                    historical = relative.startswith(('revision_results/final_20260914','publication/historical_facekit','publication/feature_definitions'))
                    rows.append({'path':relative,'line':line,'excerpt':text[:300],
                                 'role':'historical_preserve' if historical else 'review_for_active_migration',
                                 'migration_status':'not_changed_before_validation'})
    pd.DataFrame(rows).to_csv(run/'MIGRATION_INVENTORY.csv',index=False)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--run',type=Path,default=HERE/'local_only/run_20261028')
    parser.add_argument('--dry-run',action='store_true')
    args=parser.parse_args();run=args.run.resolve()
    if HERE/'local_only' not in run.parents:
        raise ValueError('Private run required')
    if args.dry_run:
        print('Dry run: inspect validated 120 inputs, create read-only historical archive, freeze real-patient cohorts, inventory migration paths; no inference or active migration')
        return
    run.mkdir(parents=True,exist_ok=True,mode=0o700)
    archive(run)
    if not (run/'ANALYSIS_PLAN.json').exists():
        prepare(run)
    inventory(run)


if __name__=='__main__':
    main()

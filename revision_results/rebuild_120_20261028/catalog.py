import hashlib
import json
import os
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
RUN=HERE/'local_only/run_20261028'
EXTRACTION=ROOT/'local_review/feature_extraction_120/20261028_current120'
os.umask(0o077)


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(8*1024*1024),b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    out=RUN/'CATALOG';out.mkdir(exist_ok=True,mode=0o700)
    features=pd.read_csv(RUN/'COHORTS/exact_feature_names.csv').feature.tolist()
    source_provenance=json.loads((EXTRACTION/'gmdb/provenance.json').read_text())
    gmdb=pd.read_parquet(EXTRACTION/'gmdb/gmdb_120_all_records_measurements.parquet')
    fair_rows=[]
    connection=sqlite3.connect(f'file:{EXTRACTION/"fairface/measurements_checkpoint.sqlite3"}?mode=ro',uri=True)
    for (encoded,) in connection.execute('SELECT record FROM records'):
        item=json.loads(encoded);row=dict(item['metadata']);qc=item['qc']
        row.update({k:v for k,v in qc.items() if k!='measurements'})
        row.update(dict(zip(features,qc['measurements'] or [None]*120)))
        row['image_sha256']=item['image_sha256'];fair_rows.append(row)
    connection.close()
    fair=pd.DataFrame(fair_rows)
    pdidb=pd.read_parquet(RUN/'PDIDB/pdidb_synthetic_120_all_records.parquet')
    pdidb_revision=json.loads((RUN/'PDIDB/provenance.json').read_text())['dataset_revision']
    fair_revision=json.loads((EXTRACTION/'inputs/fairface/zip_download_provenance.json').read_text())['sha256']
    frames=[]
    for name,table,revision in [('GMDB',gmdb,'cache_fingerprint:a94fb856b30eaeeea670b9f06db903afc44e29db'),
        ('FairFace',fair,'official_margin025_zip_sha256:'+fair_revision),('PDIDB',pdidb,pdidb_revision)]:
        base=pd.DataFrame({'source_dataset':[name]*len(table),'synthetic':[name=='PDIDB']*len(table),
            'dataset_revision':[revision]*len(table),'extractor_commit':[source_provenance['facekit_commit']]*len(table),
            'model_hash':[source_provenance['model_sha256']]*len(table),
            'feature_schema_version':['facekit120_5764bdd4']*len(table),
            'image_quality_status':table.reason.to_numpy(),'quality_passed':table.quality_passed.to_numpy(),
            'source_image_key':table.key.to_numpy(),'image_sha256':table.image_sha256.to_numpy()})
        base['recorded_patient_group_key']=('GMDB:'+table.patient_id.astype('string')).to_numpy() if name=='GMDB' else None
        base['syndrome_label']=table.historical_disease_label.to_numpy() if name=='GMDB' else table.disease.to_numpy() if name=='PDIDB' else None
        base['age_category']=table.age.to_numpy() if name=='FairFace' else None
        base['recorded_age_year']=table.age_year.to_numpy() if name=='GMDB' else None
        base['recorded_age_month']=table.age_month.to_numpy() if name=='GMDB' else None
        base['recorded_gender_label']=table.gender.astype('string').to_numpy() if name!='PDIDB' else None
        base['fairface_source_race_label']=table.race.to_numpy() if name=='FairFace' else None
        base['gmdb_ethnicity_category']=table.ethnicity_category.to_numpy() if name=='GMDB' else None
        base['source_independence_status']='recorded_patient_ids_only' if name=='GMDB' else 'distinct_image_paths_not_verified_people' if name=='FairFace' else 'unknown_generation_lineage'
        frames.append(pd.concat([base.reset_index(drop=True),table[features].reset_index(drop=True)],axis=1))
    unified=pd.concat(frames,ignore_index=True)
    assert len(unified)==10980+11764+1727 and unified[['source_dataset','source_image_key']].duplicated().sum()==0
    assert unified.loc[unified.source_dataset.eq('PDIDB'),'synthetic'].all()
    unified.to_parquet(out/'unified_120_private_catalog.parquet',index=False)
    roles={}
    primary=pd.read_parquet(RUN/'COHORTS/primary.parquet')
    primary['source_dataset']='GMDB';primary['synthetic']=False
    roles['gmdb_real_120']=primary
    for name,cohort in [('fairface_balanced_120','balanced'),('fairface_all_eligible_120','all_eligible'),('fairface_race_sex_120','race_sex')]:
        frame=pd.read_parquet(RUN/f'COHORTS/reference_{cohort}.parquet')
        frame['source_dataset']='FairFace';frame['synthetic']=False
        roles[name]=frame
    roles['pdidb_synthetic_120']=pdidb
    entries=[]
    for name,frame in roles.items():
        path=out/f'{name}.parquet';frame.to_parquet(path,index=False)
        synthetic=name=='pdidb_synthetic_120'
        entries.append({'table':name,'path':str(path),'records':len(frame),
            'known_recorded_subjects':int(frame.patient_id.nunique()) if name=='gmdb_real_120' else None,
            'synthetic':synthetic,'primary_role':'disease_input' if name=='gmdb_real_120' else 'excluded_synthetic' if synthetic else 'reference_alternative',
            'sha256':sha(path),'missing_feature_cells':int(frame[features].isna().sum().sum())})
    assert not any(frame.synthetic.any() for name,frame in roles.items() if name!='pdidb_synthetic_120')
    raw=pd.read_parquet(RUN/'COHORTS/all_gmdb_candidates.parquet')
    raw.to_csv(RUN/'COHORTS/all_gmdb_candidates_master.csv',index=False)
    manifest={'feature_count':120,'features':features,'catalog_records':len(unified),'catalog_sha256':sha(out/'unified_120_private_catalog.parquet'),
        'tables':entries,'source_roles_enforced':True,'synthetic_rows_in_primary':0,
        'catalog_is_not_an_analytic_pool':True,'source_patient_independence_PDIDB':'unresolved',
        'FairFace_subject_count':'unknown; record counts are distinct image paths',
        'publication_permissions':'GMDB images prohibited; PDIDB illustration redistribution unverified; no images exported',
        'credentials':'No authentication material stored in catalog or manifest'}
    (RUN/'UNIFIED_120_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    pd.DataFrame([{k:v for k,v in row.items() if k!='path'} for row in entries]).to_csv(RUN/'AGGREGATE_SOURCE_MANIFEST.csv',index=False)
    print(json.dumps({'catalog_records':len(unified),'role_counts':{name:len(frame) for name,frame in roles.items()},'synthetic_rows_in_primary':0},indent=2))


if __name__=='__main__':
    main()

import argparse
import hashlib
import io
import json
import os
import shutil
import sqlite3
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from PIL import Image, UnidentifiedImageError

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'data/dataset_creation'))
import extract_current_120 as e


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--run',type=Path,default=HERE/'local_only/run_20261028')
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args();run=args.run.resolve()
    if HERE/'local_only' not in run.parents:
        raise ValueError('Private run required')
    audit=json.loads((HERE/'local_only/pdidb_download_manifest.json').read_text())
    if audit['status']!='downloaded_with_existing_authenticated_session':
        raise ValueError('Provider access is not verified')
    out=run/'PDIDB';out.mkdir(parents=True,exist_ok=True,mode=0o700)
    raw=Path(audit['local_cache_path'])
    backup=out/'private_source_backup.parquet'
    if not backup.exists():
        shutil.copy2(raw,backup)
        assert e.sha(backup)==audit['sha256']==e.sha(raw)
        backup.chmod(0o400)
    assert e.sha(backup)==audit['sha256']
    model=Path.home()/'.cache/facekit/face_landmarker.task'
    signature=e.fingerprint(model)
    signature.update(source_dataset='PDIDB',dataset_revision=audit['revision'],input_sha256=audit['sha256'],
        synthetic=True,synthetic_provenance='Author supplied synthetic designation; v1/v2 model tags; generation method and real-patient lineage not documented',
        source_patient_independence='unknown',publication_permission='unverified; no images exported or used in report')
    checkpoint=e.Checkpoint(out/'measurements_checkpoint.sqlite3',signature,args.resume)
    done=checkpoint.keys();features=signature['features']
    metadata=pq.read_table(backup,columns=['image_id','disease','disease_abbr','model','verification','age_range','gender','race']).to_pandas()
    assert len(metadata)==1727 and metadata.image_id.is_unique
    pilot=set(metadata.groupby(['model','verification','disease'],dropna=False).head(1).image_id)
    pipeline=e.Pipeline(model,features)
    diagnostics=[]
    try:
        for phase in ['pilot','full']:
            cursor=0
            for batch in pq.ParquetFile(backup).iter_batches(batch_size=100):
                batch.validate(full=True)
                for record in batch.to_pylist():
                    identifier=record['image_id'];key='PDIDB:'+identifier
                    cursor+=1
                    if key in done or (phase=='pilot' and identifier not in pilot):
                        continue
                    data=record['image'].get('bytes') if record['image'] else None
                    qc=pipeline.process(data or b'')
                    if phase=='pilot':
                        repeated=pipeline.process(data or b'')
                        equal=qc['reason']==repeated['reason']
                        delta=None
                        if qc['measurements'] is not None:
                            delta=float(np.max(np.abs(np.asarray(qc['measurements'])-repeated['measurements'])))
                            equal=equal and np.allclose(qc['measurements'],repeated['measurements'],rtol=1e-10,atol=1e-10)
                        if not equal:
                            raise ValueError('PDIDB deterministic replay failed')
                        diagnostics.append({'image_key':key,'deterministic':bool(equal),'max_absolute_feature_change':delta,'reason':qc['reason']})
                    fields={name:e.plain(record[name]) for name in ['image_id','disease','disease_abbr','model','verification','age_range','gender','race']}
                    fields.update(key=key,source_dataset='PDIDB',synthetic=True,dataset_revision=audit['revision'],
                        recorded_patient_group_key=None,source_linkage_status='unresolved; no source-patient field supplied',
                        analytic_role='synthetic_separate_only; illustration_rights_unverified')
                    digest=hashlib.sha256(data).hexdigest() if data else None
                    checkpoint.add(key,{'metadata':fields,'qc':qc,'image_sha256':digest})
                    done.add(key)
                    if len(done)%100==0:
                        print(f'PDIDB: {len(done)}/1727 checkpointed',flush=True)
            checkpoint.commit()
            if phase=='pilot' and diagnostics:
                pd.DataFrame(diagnostics).to_csv(out/'private_pilot_determinism.csv',index=False)
                assert any(item['reason']=='accepted' for item in diagnostics)
        rows=[]
        for item in checkpoint.records():
            row=dict(item['metadata']);qc=item['qc']
            row.update({k:v for k,v in qc.items() if k!='measurements'})
            row.update(dict(zip(features,qc['measurements'] or [None]*120)))
            row.update(image_sha256=item['image_sha256'],extractor_commit=e.COMMIT,model_hash=e.MODEL_SHA)
            rows.append(row)
        table=pd.DataFrame(rows)
        assert len(table)==1727 and table.key.is_unique
        valid=table.loc[table.quality_passed]
        assert np.isfinite(valid[features].to_numpy(float)).all()
        table.to_parquet(out/'pdidb_synthetic_120_all_records.parquet',index=False)
        valid.to_parquet(out/'pdidb_synthetic_120_quality_passed.parquet',index=False)
        table.drop(columns=features).to_csv(out/'private_extraction_status.csv',index=False)
        images=table.image_sha256.dropna().value_counts()
        real=pd.read_parquet(ROOT/'local_review/feature_extraction_120/20261028_current120/gmdb/gmdb_120_all_records_measurements.parquet',columns=['image_sha256'])
        overlap=int(table.image_sha256.isin(set(real.image_sha256.dropna())).sum())
        old=pd.read_csv(ROOT/'data/syndrome_folder_map.csv')
        def normalized(value):
            return ''.join(c.lower() for c in str(value) if c.isalnum())
        lookup={}
        for row in old.itertuples():
            for name in [row.syndrome_name,row.internal_syndrome_name]:
                lookup.setdefault(normalized(name),set()).add(row.facekit_cohort)
        cross=[]
        for name,count in table.disease.value_counts().items():
            matches=lookup.get(normalized(name),set())
            cross.append({'pdidb_disease':name,'records':int(count),'historical_label':next(iter(matches)) if len(matches)==1 else None,
                'mapping_status':'normalized_exact_label_match' if len(matches)==1 else 'unresolved',
                'independence_status':'unresolved; label match is not source-patient linkage'})
        pd.DataFrame(cross).to_csv(out/'taxonomy_crosswalk.csv',index=False)
        summary={'status':'FULL_EXTRACTION','revision':audit['revision'],'records':len(table),'quality_passed':len(valid),
            'reason_counts':{str(k):int(v) for k,v in table.reason.value_counts().items()},
            'same_pipeline_as_real_sources':True,'feature_count':120,'duplicate_byte_groups':int(images.gt(1).sum()),
            'exact_encoded_byte_overlap_with_GMDB_rows':overlap,
            'overlap_interpretation':'Exact-byte check only; no finding establishes privacy or source independence',
            'source_patient_ids_available':False,'independent_subject_count':None,
            'demographic_missingness':{key:int(table[key].isna().sum()) for key in ['age_range','gender','race']},
            'model_tags':table.model.value_counts().to_dict(),'verification_tags':table.verification.value_counts().to_dict(),
            'generation_provenance':'Author designates this dataset synthetic; schema supplies v1/v2 tags but no generator identity or real-patient lineage',
            'publication_permission':'NOT ESTABLISHED; gated access is not redistribution permission',
            'primary_inferential_rows':0,'optional_source_agreement_analysis':'NOT RUN: no verified source linkage or independence plan'}
        e.atomic_json(out/'extraction_summary.json',summary)
        e.atomic_json(out/'provenance.json',dict(signature,status='FULL',source_backup_sha256=e.sha(backup),output_rows=len(table),credentials='Existing SDK-authenticated session used only for upstream download; no credentials copied into workflow files'))
        assert e.sha(raw)==audit['sha256']==e.sha(backup)
        print(json.dumps(summary,indent=2),flush=True)
    finally:
        pipeline.close();checkpoint.close()


if __name__=='__main__':
    def offline(event,args):
        if event in {'socket.connect','socket.getaddrinfo','http.client.connect','urllib.Request'}:
            raise RuntimeError('Extraction is offline')
    sys.addaudithook(offline)
    main()

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
os.umask(0o077)
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from scipy.cluster.hierarchy import linkage,fcluster
from scipy.spatial.distance import squareform

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
RUN=HERE/'local_only/run_20261028'
sys.path.insert(0,str(ROOT/'revision_results/reviewer_followup_2026'))
import followup as f
import extended_analysis as ext


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def contrast(frame,targets,features,comparator,iterations):
    rows=[]
    for index,syndrome in enumerate(targets):
        target,other=f.comparator_frames(frame,syndrome,targets,comparator)
        x,y=target[features].to_numpy(float),other[features].to_numpy(float)
        nx,ny=np.isfinite(x).sum(0),np.isfinite(y).sum(0)
        g=f.effect(x,y);valid=f.eligible_feature_counts(nx,ny,len(x),len(y))&np.isfinite(g)
        test=mannwhitneyu(x,y,axis=0,alternative='two-sided',method='asymptotic',use_continuity=True,nan_policy='omit')
        draws=f.bootstrap_effect(x,y,iterations,[f.SEED,index])
        low,high=np.nanquantile(draws,[.025,.975],axis=0)
        for j,feature in enumerate(features):
            rows.append({'syndrome':syndrome,'feature':feature,'comparator':comparator,
                'target_total':len(x),'comparator_total':len(y),'n_target':nx[j],'n_comparator':ny[j],
                'target_missing_fraction':1-nx[j]/len(x),'comparator_missing_fraction':1-ny[j]/len(y),
                'estimable':bool(valid[j]),'hedges_g':g[j] if valid[j] else np.nan,
                'ci_95_lower':low[j] if valid[j] else np.nan,'ci_95_upper':high[j] if valid[j] else np.nan,
                'bootstrap_iterations':iterations,'bootstrap_finite':int(np.isfinite(draws[:,j]).sum()),
                'mannwhitney_u':test.statistic[j],'p_value':test.pvalue[j] if valid[j] else 1.})
        if (index+1)%10==0 or index==len(targets)-1:
            print(f'Within GMDB {comparator}: {index+1}/{len(targets)}',flush=True)
    result=pd.DataFrame(rows);result['q_global_bh']=f.bh(result.p_value)
    result['significant']=result.estimable&result.q_global_bh.lt(.05)
    assert len(result)==len(targets)*len(features)
    return result


def random_photos(raw,reference,primary,targets,features,replicates):
    age=f.age_diagnostics(raw).conditional_age_years
    eligible=raw.loc[raw.frontal_ok.eq(True)&age.ge(3)&age.lt(10)].copy()
    eligible['patient_id']=f.clean_id(eligible.patient_id)
    assert eligible.patient_id.notna().all()
    groups=list(eligible.groupby('patient_id',sort=True).indices.values())
    assert len(groups)==len(primary)
    values=eligible[features].to_numpy(float);reference=reference[features].to_numpy(float)
    labels=eligible.disease.to_numpy();rng=np.random.default_rng(f.SEED+1)
    effects,qvalues,replicate_rows=[],[],[]
    for rep in range(replicates):
        chosen=np.array([indices[rng.integers(len(indices))] for indices in groups])
        samples,diseases=values[chosen],labels[chosen]
        gs,ps=[],[]
        for syndrome in targets:
            x=samples[diseases==syndrome]
            assert f.eligible_feature_counts(np.isfinite(x).sum(0),np.isfinite(reference).sum(0),len(x),len(reference)).all()
            gs.extend(f.effect(x,reference))
            ps.extend(mannwhitneyu(x,reference,axis=0,method='asymptotic',alternative='two-sided',nan_policy='omit').pvalue)
        qs=f.bh(ps);effects.append(gs);qvalues.append(qs)
        replicate_rows.append({'replicate':rep+1,'hypotheses':len(qs),'significant':int((qs<.05).sum())})
        if (rep+1)%20==0:
            print(f'Random photographs: {rep+1}/{replicates}',flush=True)
    gs,qs=np.array(effects),np.array(qvalues);rows=[]
    for index,(syndrome,feature) in enumerate((s,v) for s in targets for v in features):
        frequency=(qs[:,index]<.05).mean()
        rows.append({'syndrome':syndrome,'feature':feature,'replicates':replicates,
            'median_g':np.median(gs[:,index]),'g_q025':np.quantile(gs[:,index],.025),'g_q975':np.quantile(gs[:,index],.975),
            'min_g':gs[:,index].min(),'max_g':gs[:,index].max(),'positive_g_frequency':(gs[:,index]>0).mean(),
            'significant_frequency':frequency,'frequency_mc_standard_error':np.sqrt(frequency*(1-frequency)/replicates)})
    return pd.DataFrame(rows),pd.DataFrame(replicate_rows)


def reference_redundancy(reference,features,mapping,entries,adjusted,out):
    corr=reference[features].corr(method='spearman')
    rho=corr.to_numpy();assert np.isfinite(rho).all()
    distance=np.clip(1-np.abs(rho),0,1);np.fill_diagonal(distance,0)
    clusters=fcluster(linkage(squareform(distance,checks=True),method='average'),.3,criterion='distance')
    groups=pd.DataFrame({'feature':features,'cluster_id':clusters}).merge(mapping[['feature','anatomical_region','feature_family']],on='feature',validate='one_to_one')
    groups['cluster_size']=groups.groupby('cluster_id').feature.transform('size')
    medians={}
    for cluster,group in groups.groupby('cluster_id'):
        index=[features.index(feature) for feature in group.feature]
        upper=np.abs(rho[np.ix_(index,index)])[np.triu_indices(len(index),1)]
        medians[cluster]=float(np.median(upper)) if len(upper) else np.nan
    groups['median_absolute_within_cluster_rho']=groups.cluster_id.map(medians)
    corr.rename_axis('feature').to_csv(out/'feature_correlation_matrix.csv')
    groups.to_csv(out/'feature_clusters.csv',index=False)
    upper=np.abs(rho[np.triu_indices(len(features),1)])
    pd.DataFrame([{'absolute_rho_threshold':value,'pairs':int((upper>=value).sum()),'total_pairs':len(upper),
        'source':'new_120_balanced_FairFace_reference'} for value in [.7,.8,.9,.95]]).to_csv(out/'redundancy_thresholds.csv',index=False)
    joined=entries.loc[entries.comparator.eq('all')].merge(groups[['feature','cluster_id']],on='feature',validate='many_to_one')
    joined=joined.merge(adjusted[['syndrome','feature','status','significant']].rename(columns={'status':'adjusted_status','significant':'adjusted_significant'}),on=['syndrome','feature'],validate='one_to_one')
    rows=[]
    for cluster,group in joined.groupby('cluster_id'):
        subset=groups.loc[groups.cluster_id.eq(cluster)]
        rows.append({'cluster_id':cluster,'features':len(subset),'hypotheses':len(group),
            'fairface_discoveries':int(group.fairface_significant.sum()),'within_discoveries':int(group.within_significant.sum()),
            'adjusted_supported':int(group.adjusted_status.eq('ok').sum()),'adjusted_discoveries':int(group.adjusted_significant.sum()),
            'anatomical_regions':';'.join(sorted(subset.anatomical_region.unique())),
            'median_absolute_attenuation':group.absolute_attenuation.median()})
    return pd.DataFrame(rows)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--phase',choices=['within','models','concordance','photos','all'],default='all')
    parser.add_argument('--dry-run',action='store_true');args=parser.parse_args()
    plan=json.loads((RUN/'ANALYSIS_PLAN.json').read_text());targets=plan['targets'];fixed=plan['historical_targets']
    out=RUN/'RESULTS_120/extensions';out.mkdir(parents=True,exist_ok=True,mode=0o700)
    cohorts=RUN/'COHORTS'
    features=pd.read_csv(cohorts/'exact_feature_names.csv').feature.tolist()
    primary=pd.read_csv(cohorts/'primary_master.csv',dtype={'patient_id':'string','image_id':'string'})
    reference=pd.read_csv(cohorts/'reference_balanced_master.csv')
    raw=pd.read_csv(cohorts/'all_gmdb_candidates_master.csv',dtype={'patient_id':'string','image_id':'string'})
    f.OUT=out
    if args.dry_run:
        print(json.dumps({'within_targets':len(targets),'historical_extension_targets':len(fixed),
            'adjusted_planned':len(fixed)*120,'pairwise_planned':len(ext.planned_pairs(fixed))*120,
            'age_strata_planned':len(fixed)*120*2,'bootstrap':1000,'concordance_bootstrap':1000,'random_photos':100},indent=2));return
    signature={'primary':sha(cohorts/'primary_master.csv'),'reference':sha(cohorts/'reference_balanced_master.csv'),
               'script':sha(Path(__file__)),'followup_source':sha(Path(f.__file__)),'extension_source':sha(Path(ext.__file__))}
    def stage(name,fn):
        state=out/(name+'_STATUS.json')
        if state.exists():
            old=json.loads(state.read_text())
            assert old['signature']==signature,'Changed inputs or code; choose a new run'
            if old['status']=='COMPLETE':
                assert all(sha(out/file)==digest for file,digest in old['outputs'].items())
                print(name+': completed checkpoint verified',flush=True);return
        before=set(out.glob('*.csv'))
        state.write_text(json.dumps({'status':'RUNNING','signature':signature},indent=2))
        try:
            fn()
            files=set(out.glob('*.csv'))-before
            state.write_text(json.dumps({'status':'COMPLETE','signature':signature,'outputs':{p.name:sha(p) for p in files}},indent=2)+'\n')
        except BaseException as error:
            state.write_text(json.dumps({'status':'FAILED','signature':signature,'error_class':type(error).__name__},indent=2)+'\n');raise
    def within():
        contrast(primary,targets,features,'all',1000).to_csv(out/'within_gmdb_results.csv',index=False)
        contrast(primary,targets,features,'eligible49',1000).to_csv(out/'within_gmdb_current_eligible_sensitivity.csv',index=False)
        f.adjusted(primary,features)
    def models():
        ext.full_adjusted(primary,fixed,features).to_csv(out/'full_adjusted_screen.csv',index=False)
        pairs=ext.pairwise(primary,fixed,features,1000)
        pairs['fdr_family']='historical_fixed_anchor_pairs_120_features';pairs.to_csv(out/'targeted_pairwise.csv',index=False)
        age=ext.age_stratified(primary,fixed,features,1000);age.to_csv(out/'age_stratified.csv',index=False)
        ext.age_agreement(age).to_csv(out/'age_stratum_agreement.csv',index=False)
        import importlib.util
        spec=importlib.util.spec_from_file_location('candidate_covariates',RUN/'EXECUTED_CODE/reference_covariate_models.py')
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        result=module.adjusted_models(primary,reference)
        result['estimand']='Adjusted syndrome-minus-reference mean difference in current normalized feature units'
        result.to_csv(out/'three_predefined_reference_adjusted_models.csv',index=False)
    def concordance():
        result=pd.read_csv(RUN/'RESULTS_120/primary_balanced/stats_each_disease_vs_healthy.csv')
        mapping=pd.read_csv(cohorts/'feature_anatomy_grouping.csv')
        existing={'all':pd.read_csv(out/'within_gmdb_results.csv'),
            'eligible49':pd.read_csv(out/'within_gmdb_current_eligible_sensitivity.csv')}
        summary,entries=ext.concordance(primary,reference,result,targets,features,mapping,existing,1000)
        summary.to_csv(out/'concordance_summary.csv',index=False);entries.to_csv(out/'concordance_entries.csv',index=False)
        adjusted=pd.read_csv(out/'full_adjusted_screen.csv')
        reference_redundancy(reference,features,mapping,entries,adjusted,out).to_csv(out/'redundancy_cluster_discoveries.csv',index=False)
    def photos():
        result,replicates=random_photos(raw,reference,primary,targets,features,100)
        result.to_csv(out/'random_photo_sensitivity.csv',index=False);replicates.to_csv(out/'random_photo_replicates.csv',index=False)
    for name,fn in [('within',within),('models',models),('concordance',concordance),('photos',photos)]:
        if args.phase in [name,'all']:
            stage(name,fn)
    (out/'FAMILY_DEFINITIONS.json').write_text(json.dumps({'primary_current_targets':len(targets),
        'within_each_family':len(targets)*120,'within_restricted_comparator_count':len(targets),
        'legacy_internal_label_eligible49_means_current_target_set':True,
        'adjusted_fixed_historical_targets':len(fixed),'adjusted_planned':len(fixed)*120,
        'targeted_pairs':len(ext.planned_pairs(fixed)),'pairwise_planned':len(ext.planned_pairs(fixed))*120,
        'age_strata_planned':len(fixed)*120*2,'unsupported_extension_cells':'correction p=1, estimate/q missing; retained in planned denominator',
        'four_highlighted_within_GMDB_models':'separate original 4-model family',
        'three_highlighted_reference_models':'separate original 3-model family'},indent=2)+'\n')


if __name__=='__main__':
    main()

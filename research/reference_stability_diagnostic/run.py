import argparse
import importlib.metadata
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(HERE))
from research.race_stratified_sensitivity import core
from engine import Bank,digest,ladder,stability

OLD=ROOT/'local_review/race_stratified_sensitivity'
PUBLIC=ROOT/'research/race_stratified_sensitivity_results_20261115/summary.json'
KEYS=core.KEYS


def read_original(bank):
    artifacts=json.loads((OLD/'artifact_sha256.json').read_text())
    paths=[OLD/'artifact_sha256.json',PUBLIC]
    for name in ['analysis_plan.json','reference_resampling_aggregate.csv','reference_resampling_repetitions.csv',
                 'checkpoints/reference-resampling/state.json']:
        path=OLD/name
        if core.baseline.sha(path)!=artifacts[name]:
            raise ValueError('BLOCKED: original sensitivity artifact changed: '+name)
        paths.append(path)
    oldplan=json.loads((OLD/'analysis_plan.json').read_text())
    if oldplan['features']!=bank.features or oldplan['planned_hypotheses']!=bank.total:
        raise ValueError('BLOCKED: original sensitivity schema differs')
    for name,h in oldplan['code_sha256'].items():
        if core.baseline.sha(ROOT/name)!=h:
            raise ValueError('BLOCKED: original execution code differs')
    state=json.loads((OLD/'checkpoints/reference-resampling/state.json').read_text())
    gs,qs=[],[]
    stop=0
    for chunk in state['chunks']:
        path=OLD/'checkpoints/reference-resampling'/chunk['file']
        if chunk['start']!=stop or core.baseline.sha(path)!=chunk['sha256']:
            raise ValueError('BLOCKED: original checkpoint corrupt/noncontiguous')
        if core.baseline.sha(path)!=artifacts[str(path.relative_to(OLD))]:
            raise ValueError('BLOCKED: checkpoint differs from final artifact map')
        with np.load(path) as saved:
            gs.append(saved['g'].copy())
            qs.append(saved['q'].copy())
        paths.append(path)
        stop=chunk['stop']
    if stop!=state['completed'] or stop!=1000:
        raise ValueError('BLOCKED: all original 1,000 draws required')
    g,q=np.vstack(gs),np.vstack(qs)
    if g.shape!=(1000,9600) or q.shape!=g.shape:
        raise ValueError('BLOCKED: original draw dimensions wrong')
    public=json.loads(PUBLIC.read_text())
    agg=pd.read_csv(OLD/'reference_resampling_aggregate.csv',float_precision='round_trip')
    if len(agg)!=9600 or agg[KEYS].duplicated().any():
        raise ValueError('BLOCKED: sensitivity keys wrong')
    index=pd.MultiIndex.from_frame(bank.primary[KEYS])
    agg=agg.set_index(KEYS).loc[index].reset_index()
    freq,sign,stable=stability(g,q,bank.primary.hedges_g.to_numpy())
    np.testing.assert_allclose(freq,agg.significant_frequency,atol=1e-14)
    np.testing.assert_allclose(sign,agg.sign_agreement_frequency,atol=1e-14)
    np.testing.assert_array_equal(stable,agg.descriptive_stable)
    np.testing.assert_allclose(np.quantile(g,[.025,.5,.975],axis=0),agg[['g_resampling_q025','median_g','g_resampling_q975']].to_numpy().T,rtol=1e-12)
    expected={'African':(135,76),'Asian':(252,193),'European':(595,580)}
    for race,(discoveries,stables) in expected.items():
        mask=bank.primary.canonical_race.eq(race).to_numpy()&bank.primary.q.lt(.05).to_numpy()
        if mask.sum()!=discoveries or stable[mask].sum()!=stables:
            raise ValueError('BLOCKED: primary/stability totals differ')
        row=next(r for r in public['reference_stability_by_race'] if r['canonical_race']==race)
        if row['original_discoveries']!=discoveries or row['descriptive_stable']!=stables:
            raise ValueError('BLOCKED: public/private stability mismatch')
    return g,q,oldplan,agg,paths


def selections(bank,sizes,seed,stream,repetition,replacement=False):
    return {race:np.random.default_rng(np.random.SeedSequence([seed,stream,repetition,i])).choice(len(data['y']),sizes[race],replace=replacement)
            for i,(race,data) in enumerate(bank.races.items())}


def draw_metadata(bank,selections):
    hashes=[]
    counts=[]
    for race,data in bank.races.items():
        w=np.bincount(selections[race],minlength=len(data['y']))
        hashes.append(digest(w))
        counts.append(np.count_nonzero(w))
    return np.array(hashes,dtype='S64'),np.asarray(counts)


def verify_direct(bank,selections,p,g):
    checked=0
    for race,data in bank.races.items():
        y=data['y'][selections[race]]
        for cell in sorted(set([0,len(data['cases'])-1])):
            x=data['cases'][cell]
            u0,pp,gg=core.tests_direct(x,y)
            where=data['where'].reshape(len(data['cases']),-1)[cell]
            np.testing.assert_allclose(p[where],pp,rtol=1e-9,atol=1e-14)
            np.testing.assert_allclose(g[where],gg,rtol=1e-10,atol=1e-12)
            for feature in [0,25,97,119]:
                direct=mannwhitneyu(x[:,feature],y[:,feature],alternative='two-sided',method='auto')
                np.testing.assert_allclose(p[where[feature]],direct.pvalue,rtol=1e-9,atol=1e-14)
                checked+=1
    np.testing.assert_allclose(core.bh(p),multipletests(p,method='fdr_bh')[1],rtol=1e-12,atol=1e-15)
    return checked


def summarize(bank,g,q,label,scheme,sizes,ref_means=None,ref_sds=None,digests=None,unique_counts=None):
    primary=bank.primary
    freq,sign,stable=stability(g,q,primary.hedges_g.to_numpy())
    frame=primary.copy()
    frame['design']=label
    frame['sampling_scheme']=scheme
    frame['repetitions']=len(g)
    frame['evaluated_reference_n']=frame.canonical_race.map(sizes)
    frame['significant_frequency']=freq
    frame['sign_agreement_frequency']=sign
    frame['descriptive_stable']=stable
    frame['g_q025'],frame['median_g'],frame['g_q975']=np.nanquantile(g,[.025,.5,.975],axis=0)
    frame['g_interval_width']=frame.g_q975-frame.g_q025
    frame['g_empirical_sd']=np.nanstd(g,axis=0,ddof=1) if len(g)>1 else np.nan
    frame['median_q']=np.median(q,axis=0)
    frame['minimum_q']=q.min(0)
    frame['maximum_q']=q.max(0)
    frame['pass_q_recurrence']=freq>=.8
    frame['pass_sign_agreement']=sign>=.95
    frame['ever_opposite_sign']=(g*primary.hedges_g.to_numpy()[None,:]<0).any(0)
    frame['ever_nonmatching_sign']=(np.sign(g)!=np.sign(primary.hedges_g.to_numpy())[None,:]).any(0)
    frame['g_interval_covers_zero']=frame.g_q025.le(0)&frame.g_q975.ge(0)
    frame['undefined_g_draws']=(~np.isfinite(g)).sum(0)
    frame['original_q_near_005']=frame.q.ge(.025)&frame.q.lt(.05)
    summaries=[]
    for i,(race,data) in enumerate(bank.races.items()):
        mask=primary.canonical_race.eq(race).to_numpy()
        original=mask&primary.q.lt(.05).to_numpy()
        subset_varies=scheme!='deterministic_full_pool' and not (scheme in ['without_replacement_subset','original_without_replacement_subset'] and sizes[race]==len(data['y']))
        frame.loc[mask,'reference_selection_varies']=subset_varies
        if ref_means is not None:
            for field,values in [('reference_mean_empirical_sd',np.std(ref_means[:,i],axis=0,ddof=1) if len(g)>1 else np.full(120,np.nan)),
                                 ('reference_SD_empirical_sd',np.std(ref_sds[:,i],axis=0,ddof=1) if len(g)>1 else np.full(120,np.nan))]:
                frame.loc[mask,field]=frame.loc[mask,'feature'].map(dict(zip(bank.features,values)))
        summaries.append({'design':label,'canonical_race':race,'sampling_scheme':scheme,'reference_n':sizes[race],
            'pool_n':len(data['y']),'repetitions':len(g),'reference_selection_varies':subset_varies,
            'original_discoveries':int(original.sum()),'stable_original':int(stable[original].sum()) if subset_varies else np.nan,
            'stable_original_percent':100*stable[original].mean() if subset_varies else np.nan,
            'pass_q_original':int((freq[original]>=.8).sum()) if subset_varies else np.nan,
            'pass_sign_original':int((sign[original]>=.95).sum()) if subset_varies else np.nan,
            'median_sign_agreement_original':float(np.median(sign[original])) if subset_varies else np.nan,
            'median_g_interval_width_original':float(frame.loc[original,'g_interval_width'].median()) if subset_varies else np.nan,
            'mean_g_interval_width_original':float(frame.loc[original,'g_interval_width'].mean()) if subset_varies else np.nan,
            'all_hypotheses':int(mask.sum()),'any_sign_reversal_original':int(frame.loc[original,'ever_opposite_sign'].sum()),
            'any_sign_reversal_all':int(frame.loc[mask,'ever_opposite_sign'].sum()),
            'mean_discoveries':float((q[:,mask]<.05).sum(1).mean()),'mean_q_recurrence_original':float(freq[original].mean()),
            'median_reference_mean_sd':float(frame.loc[mask,'reference_mean_empirical_sd'].median()) if ref_means is not None and len(g)>1 else np.nan,
            'truly_distinct_draws':int(len(np.unique(digests[:,i]))) if digests is not None else 1,
            'mean_distinct_reference_photographs':float(unique_counts[:,i].mean()) if unique_counts is not None else sizes[race],
            'replacement_duplicate_draws':int((unique_counts[:,i]<sizes[race]).sum()) if scheme=='with_replacement_reference_bootstrap' else 0})
        if not subset_varies:
            frame['descriptive_stable']=frame['descriptive_stable'].astype('Float64')
            frame.loc[mask,'descriptive_stable']=pd.NA
        if scheme=='deterministic_full_pool':
            for column in ['significant_frequency','sign_agreement_frequency','pass_q_recurrence','pass_sign_agreement']:
                frame[column]=frame[column].astype('Float64')
                frame.loc[mask,column]=pd.NA
    return frame,pd.DataFrame(summaries)


def main(args):
    output=args.output_dir.resolve()
    if not output.name.startswith('reference_stability_diagnostic') or output.is_relative_to(OLD) or output.is_relative_to(core.baseline.PRIVATE):
        raise ValueError('New isolated reference_stability_diagnostic output directory required')
    loaded=core.load_baseline()
    bank=Bank(loaded,core)
    old_g,old_q,oldplan,oldagg,oldpaths=read_original(bank)
    pools={r:len(d['y']) for r,d in bank.races.items()}
    configs=ladder(bank.sampling['quotas'],pools)
    paths=bank.sources+oldpaths+[ROOT/'ACTIVE_ANALYSIS.json',ROOT/'SENSITIVITY_ANALYSIS.json']
    hashes={str(p):core.baseline.sha(p) for p in set(paths)}
    codes={str(p.relative_to(ROOT)):core.baseline.sha(p) for p in [HERE/'engine.py',HERE/'run.py']}
    plan={'study_id':'PRISM_REFERENCE_STABILITY_DIAGNOSTIC_20261116','baseline_id':'PRISM_RACE_STRATIFIED_20261115',
          'requested_repetitions':args.repetitions,'execution_status':'FINAL' if args.repetitions>=1000 else 'PILOT',
          'seed':args.seed,'pools':pools,'original_quotas':bank.sampling['quotas'],'ladder':configs,
          'paired_ladder_rule':'double African/Asian quotas together, cap at pools; European fixed; stage 0 reuses original draws',
          'bootstrap':'full observed reference pools, with replacement, shared per category; cases fixed; separate conditional model-based scheme',
          'stability_rule':'q<.05 in >=80% and same finite nonzero original-g sign in >=95% of draws',
          'planned_family':9600,'q_near_threshold_rule':'original q in [0.025,0.05)',
          'case_bands':['5-9','10-19','20+'],'g_bands':['<0.8','0.8-1.2','>=1.2'],
          'q_bands':['<0.001','0.001-0.01','0.01-0.025','0.025-0.05'],
          'seeding':'SeedSequence([seed, stream, repetition, category index]); B bootstrap stream=1, C stream=100+stage',
          'source_hashes':hashes,'execution_code_hashes':codes,'fixed_cases':True,'feature_order':bank.features,
          'invalid_g':'NaN remains undefined and cannot count as stable; g=0 cannot count as a directional success',
          'full_subset_rule':'one distinct full subset; no subset-stability score; European g fixed in original/C',
          'interpretation':'Selection distributions, not clinical confidence intervals, independent replication, or resolution of source confounding'}
    signature=__import__('hashlib').sha256(json.dumps(plan,sort_keys=True).encode()).hexdigest()
    if args.dry_run:
        print(json.dumps({'status':'FEASIBLE','plan':plan,'original_arrays_verified':True},indent=2),flush=True)
        return
    output=core.baseline.private_output(output)
    planpath=output/'analysis_plan.json'
    if planpath.exists() and json.loads(planpath.read_text())!=plan:
        raise ValueError('Input/code/plan changed; use a fresh output directory')
    core.write_json(planpath,plan)
    protected=core.baseline.protected_hashes()
    start=time.monotonic()
    checks=[]
    original={r:np.flatnonzero(d['ref'].file.isin(core.baseline.sample_reference(bank.reference,bank.sampling['quotas'],bank.sampling['seed']).file).to_numpy())
              for r,d in bank.races.items()}
    p,g,_,_=bank.evaluate(original)
    np.testing.assert_allclose(p,bank.primary.p,rtol=1e-9,atol=1e-14)
    np.testing.assert_allclose(g,bank.primary.hedges_g,rtol=1e-10,atol=1e-12)
    np.testing.assert_allclose(core.bh(p),bank.primary.q,rtol=1e-9,atol=1e-14)
    checks.append(verify_direct(bank,original,p,g))
    summaries=[]
    tables=[]
    if args.stage in ['audit','A','all'] or not (output/'original_reference_per_feature.csv').exists():
        means=[];sds=[];digests=[];unique=[]
        for rep in range(1000):
            draw={r:core.rng_for(oldplan['seed'],3,rep,i).choice(len(d['y']),bank.sampling['quotas'][r],replace=False)
                  for i,(r,d) in enumerate(bank.races.items())}
            mm=[];ss=[]
            for r,d in bank.races.items():
                yy=d['y'][draw[r]]
                mm.append(yy.mean(0));ss.append(yy.std(0,ddof=1))
            if rep==0:
                pp,gg,_,_=bank.evaluate(draw)
                np.testing.assert_allclose(core.bh(pp),old_q[0],rtol=1e-9,atol=1e-13)
                np.testing.assert_allclose(gg,old_g[0],rtol=1e-10,atol=1e-12)
                checks.append(verify_direct(bank,draw,pp,gg))
            dd,nn=draw_metadata(bank,draw)
            means.append(mm);sds.append(ss);digests.append(dd);unique.append(nn)
        frame,summary=summarize(bank,old_g,old_q,'A_original','original_without_replacement_subset',bank.sampling['quotas'],
                                np.asarray(means),np.asarray(sds),np.asarray(digests),np.asarray(unique))
        frame.to_csv(output/'original_reference_per_feature.csv',index=False)
        summary.to_csv(output/'original_reference_summary.csv',index=False)
        failures=frame.loc[frame.q.lt(.05)&frame.canonical_race.isin(['African','Asian'])].copy()
        failures['failure_reason']=np.select([~failures.pass_q_recurrence.astype(bool)&failures.pass_sign_agreement.astype(bool),
                                              failures.pass_q_recurrence.astype(bool)&~failures.pass_sign_agreement.astype(bool),
                                              ~failures.pass_q_recurrence.astype(bool)&~failures.pass_sign_agreement.astype(bool)],
                                             ['significance_only','sign_only','both'],default='stable')
        failures.to_csv(output/'original_387_African_Asian_discovery_diagnosis.csv',index=False)
        unstable=failures.loc[failures.failure_reason.ne('stable')]
        if len(unstable)!=118:
            raise ValueError('BLOCKED: expected exactly 118 African/Asian failures')
        unstable.to_csv(output/'original_118_failure_reasons.csv',index=False)
        mechanism=failures.groupby(['canonical_race','failure_reason']).size().unstack(fill_value=0).reindex(
            index=['African','Asian'],columns=['significance_only','sign_only','both','stable'],fill_value=0)
        mechanism.to_csv(output/'failure_mechanism_counts.csv')
        print('Original 118 failures verified:',mechanism.to_dict(),flush=True)
    if args.stage in ['B','C','all']:
        full={r:np.arange(len(d['y'])) for r,d in bank.races.items()}
        fullp,fullg,fullmeans,fullsds=bank.evaluate(full)
        checks.append(verify_direct(bank,full,fullp,fullg))
        fullq=core.bh(fullp)
        comparison=bank.primary.copy()
        comparison['B_reference_n']=comparison.canonical_race.map(pools)
        comparison['B_p']=fullp;comparison['B_q_global_9600']=fullq;comparison['B_g']=fullg
        for i,(r,d) in enumerate(bank.races.items()):
            mask=comparison.canonical_race.eq(r)
            comparison.loc[mask,'B_reference_mean']=comparison.loc[mask,'feature'].map(dict(zip(bank.features,fullmeans[i])))
        comparison['B_case_minus_reference_mean']=comparison.case_mean-comparison.B_reference_mean
        comparison['g_change_B_minus_A']=comparison.B_g-comparison.hedges_g
        comparison['abs_g_change']=comparison.g_change_B_minus_A.abs()
        comparison['sign_reversal_A_B']=comparison.B_g*comparison.hedges_g<0
        comparison['original_significant']=comparison.q.lt(.05)
        comparison['B_significant']=comparison.B_q_global_9600.lt(.05)
        comparison.to_csv(output/'A_vs_B_primary_effects.csv',index=False)
        deterministic,deterministic_summary=summarize(bank,fullg[None,:],fullq[None,:],'C_full_cap','deterministic_full_pool',pools)
        deterministic.to_csv(output/'full_reference_per_feature.csv',index=False)
        deterministic_summary.to_csv(output/'full_reference_summary.csv',index=False)
        if args.stage in ['B','all']:
            stages=[('B_reference_bootstrap',pools,True,1)]
        else:
            stages=[]
        if args.stage in ['C','all']:
            stages.extend((f'C_stage_{i}',size,False,100+i) for i,size in enumerate(configs[1:-1],start=1))
        for label,sizes,replacement,stream in stages:
            stage_signature=__import__('hashlib').sha256((signature+label).encode()).hexdigest()
            cp=core.Checkpoints(output,label,stage_signature,args.resume)
            for first in range(cp.state['completed'],args.repetitions,25):
                gg=[];qq=[];mm=[];ss=[];dd=[];nn=[];discoveries=[]
                for rep in range(first,min(first+25,args.repetitions)):
                    draw=selections(bank,sizes,args.seed,stream,rep,replacement)
                    pp,g,mean,sd=bank.evaluate(draw)
                    q=core.bh(pp)
                    if rep==0:
                        checks.append(verify_direct(bank,draw,pp,g))
                    gs,counts_unique=draw_metadata(bank,draw)
                    gg.append(g);qq.append(q);mm.append(mean);ss.append(sd);dd.append(gs);nn.append(counts_unique)
                    discoveries.append(int((q<.05).sum()))
                cp.append(first+len(gg),g=np.asarray(gg),q=np.asarray(qq),means=np.asarray(mm),sds=np.asarray(ss),
                          digests=np.asarray(dd),unique=np.asarray(nn),discoveries=np.asarray(discoveries))
                print(f'{label}: {cp.state["completed"]}/{args.repetitions}; elapsed {time.monotonic()-start:.1f}s',flush=True)
            arrays={name:np.vstack(cp.arrays(name)) for name in ['g','q','means','sds','digests','unique']}
            scheme='with_replacement_reference_bootstrap' if replacement else 'without_replacement_subset'
            frame,summary=summarize(bank,arrays['g'],arrays['q'],label,scheme,sizes,arrays['means'],arrays['sds'],arrays['digests'],arrays['unique'])
            if replacement:
                freq,sg,stable=stability(arrays['g'],arrays['q'],fullg)
                frame['sign_agreement_with_B']=sg
                frame['stable_with_B_sign']=stable
            frame.to_csv(output/(label+'_per_feature.csv'),index=False)
            summary.to_csv(output/(label+'_summary.csv'),index=False)
            del arrays
    all_summary=[];all_tables=[]
    for path in [output/'original_reference_summary.csv']+sorted(output.glob('C_stage_*_summary.csv'))+[output/'full_reference_summary.csv',output/'B_reference_bootstrap_summary.csv']:
        if path.exists():
            all_summary.append(pd.read_csv(path))
    for path in [output/'original_reference_per_feature.csv']+sorted(output.glob('C_stage_*_per_feature.csv'))+[output/'full_reference_per_feature.csv',output/'B_reference_bootstrap_per_feature.csv']:
        if path.exists():
            all_tables.append(pd.read_csv(path))
    if all_summary:
        pd.concat(all_summary,ignore_index=True).to_csv(output/'reference_size_ladder_summary.csv',index=False)
        pd.concat(all_tables,ignore_index=True).to_csv(output/'reference_size_ladder_per_feature.csv',index=False)
    if protected!=core.baseline.protected_hashes():
        raise ValueError('Protected artifact changed')
    for name,h in hashes.items():
        if core.baseline.sha(Path(name))!=h:
            raise ValueError('Original input/checkpoint changed')
    verification={'status':'PASS','primary_discoveries_verified':982,'old_stable_totals_verified':[76,193,580],
                  'full_family_per_draw':9600,'direct_real_SciPy_checks_this_invocation':sum(checks),
                  'independent_statsmodels_BH_verified':True,'original_checkpoints_reused':True,'protected_artifacts_unchanged':True}
    core.write_json(output/'verification_report.json',verification)
    manifest={'plan_signature':signature,'plan':plan,'command':vars(args)|{'output_dir':str(output)},'elapsed_seconds':time.monotonic()-start,
              'git_revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
              'software':{x:importlib.metadata.version(x) for x in ['numpy','pandas','scipy','statsmodels','matplotlib']},
              'python':sys.version,'protected_hashes':protected,'source_audit':bank.audit,'original_arrays_available':True,
              'recorded_crosswalk':oldplan['crosswalk'],'all_raw_rows_remain_in_memory':True}
    core.write_json(output/'source_and_execution_manifest.json',manifest)
    from reporting import report
    report(output)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--stage',choices=['audit','A','B','C','all'],default='all')
    parser.add_argument('--repetitions',type=int,default=1000)
    parser.add_argument('--seed',type=int,default=20261116)
    parser.add_argument('--output-dir',type=Path,default=ROOT/'local_review/reference_stability_diagnostic')
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--dry-run',action='store_true')
    args=parser.parse_args()
    if args.repetitions<2:
        parser.error('At least two repetitions required; <1000 is explicitly PILOT')
    main(args)

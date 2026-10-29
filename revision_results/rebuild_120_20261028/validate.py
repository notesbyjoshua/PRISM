import argparse
import ast
import hashlib
import json
import os
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata,norm,mannwhitneyu,kruskal

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
RUN=HERE/'local_only/run_20261028'
os.umask(0o077)


def bh(values):
    values=np.asarray(values,float);q=np.full(len(values),np.nan)
    selected=np.flatnonzero(np.isfinite(values));order=selected[np.argsort(values[selected],kind='stable')]
    q[order]=np.minimum(1,np.minimum.accumulate((values[order]*len(order)/np.arange(1,len(order)+1))[::-1])[::-1])
    return q


def independent(x,y):
    x,y=np.asarray(x,float),np.asarray(y,float)
    n,m=len(x),len(y);all_values=np.r_[x,y]
    u=rankdata(all_values)[:n].sum()-n*(n+1)/2
    _,ties=np.unique(all_values,return_counts=True)
    variance=n*m/12*((n+m+1)-np.sum(ties**3-ties)/((n+m)*(n+m-1)))
    p=min(1.,float(2*norm.sf((max(u,n*m-u)-n*m/2-.5)/np.sqrt(variance)))) if variance else 1.
    pooled=((n-1)*x.var(ddof=1)+(m-1)*y.var(ddof=1))/(n+m-2)
    g=(1-3/(4*(n+m)-9))*(x.mean()-y.mean())/np.sqrt(pooled) if pooled>0 else np.nan
    return u,p,g


def close(a,b,name,rtol=1e-8,atol=1e-10):
    np.testing.assert_allclose(a,b,rtol=rtol,atol=atol,equal_nan=True,err_msg=name)


def main_tests(results,disease,reference):
    maximum={'u':0.,'p':0.,'g':0.}
    for label,frame in results.groupby('disease',sort=False):
        target=disease.loc[disease.disease.eq(label)]
        for row in frame.itertuples():
            x=target[row.feature].dropna().to_numpy(float);y=reference[row.feature].dropna().to_numpy(float)
            u,p,g=independent(x,y)
            assert len(x)==row.n_disease and len(y)==row.n_healthy
            for name,value,expected in [('u',u,row.mannwhitney_u),('p',p,row.mannwhitney_p),('g',g,row.hedges_g)]:
                close(value,expected,'main '+name)
                maximum[name]=max(maximum[name],abs(value-expected))
    close(bh(results.mannwhitney_p),results.mannwhitney_q_fdr_global,'main global BH')
    for _,frame in results.groupby('disease'):
        close(bh(frame.mannwhitney_p),frame.mannwhitney_q_fdr_within_disease,'main within-syndrome BH')
    return maximum


def matrix_hc3(y,design):
    inverse=np.linalg.inv(design.T@design)
    coefficient=inverse@design.T@y
    residual=y-design@coefficient
    leverage=np.einsum('ij,jk,ik->i',design,inverse,design)
    meat=(design.T*(residual/(1-leverage))**2)@design
    covariance=inverse@meat@inverse
    se=np.sqrt(np.diag(covariance))
    return coefficient,se


def validate_demographics(disease,reference,features,rows):
    root=RUN/'RESULTS_120/demographics'
    for layer,column in [('layer1_broad_ethnicity','broad_group'),('layer2_subcategory','detailed_group')]:
        out=root/layer
        h={group:frame for group,frame in reference.dropna(subset=[column]).groupby(column)}
        d={(name,group):frame for (name,group),frame in disease.dropna(subset=[column]).groupby(['internal_syndrome_name',column])}
        a=pd.read_csv(out/'analysis_A_healthy_kruskal_all_groups.csv')
        for item in a.itertuples():
            groups=[frame[item.feature].to_numpy(float) for frame in h.values() if len(frame)>=20]
            test=kruskal(*groups);close(test.statistic,item.kruskal_H,'A H');close(test.pvalue,item.p_value,'A p')
        close(bh(a.p_value),a.fdr_q_value,'A KW BH')
        pair=pd.read_csv(out/'analysis_A_healthy_pairwise_groups.csv')
        for (g1,g2),frame in pair.groupby(['group_1','group_2']):
            for item in frame.itertuples():
                u,p,g=independent(h[g1][item.feature].to_numpy(),h[g2][item.feature].to_numpy())
                close([u,p,g],[item.mann_whitney_u,item.p_value,item.hedges_g],'A pair')
            close(bh(frame.p_value),frame.within_pair_fdr_q_value,'A per-pair BH')
        b=pd.read_csv(out/'analysis_B_matched_disease_vs_healthy.csv')
        for (syndrome,group),frame in b.groupby(['syndrome_name','group']):
            target,control=d[syndrome,group],h[group]
            selected=frame.feature.tolist()
            x,y=target[selected].to_numpy(float),control[selected].to_numpy(float)
            test=mannwhitneyu(x,y,axis=0,alternative='two-sided')
            close(test.statistic,frame.mann_whitney_u,'B U');close(test.pvalue,frame.p_value,'B p')
            n,m=len(x),len(y);pooled=((n-1)*x.var(0,ddof=1)+(m-1)*y.var(0,ddof=1))/(n+m-2)
            effect=(1-3/(4*(n+m)-9))*(x.mean(0)-y.mean(0))/np.sqrt(pooled)
            close(effect,frame.hedges_g,'B g');close(bh(frame.p_value),frame.within_group_fdr_q_value,'B within family')
        close(bh(b.p_value),b.global_fdr_q_value,'B global BH')
        c=pd.read_csv(out/'analysis_C_all_interactions.csv')
        for (syndrome,g1,g2),frame in c.groupby(['syndrome_name','group_1','group_2']):
            columns=frame.feature.tolist()
            cells=[d[syndrome,g1][columns].to_numpy(float),d[syndrome,g2][columns].to_numpy(float),
                   h[g1][columns].to_numpy(float),h[g2][columns].to_numpy(float)]
            coefficient=cells[1].mean(0)-cells[3].mean(0)-cells[0].mean(0)+cells[2].mean(0)
            se=np.sqrt(sum(cell.var(0,ddof=1)/(len(cell)-1) for cell in cells))
            p=2*norm.sf(np.abs(coefficient/se))
            close(coefficient,frame.interaction_beta_group2_minus_group1,'C independent four-cell beta',atol=1e-8)
            close(se,frame.interaction_standard_error,'C closed-form HC3',atol=1e-8)
            close(p,frame.interaction_p_value,'C independent p',atol=1e-8)
            close(coefficient-norm.ppf(.975)*se,frame.interaction_ci_lower,'C interval lower',atol=1e-8)
            close(coefficient+norm.ppf(.975)*se,frame.interaction_ci_upper,'C interval upper',atol=1e-8)
        close(bh(c.interaction_p_value),c.global_fdr_q_value,'C global BH')
        for _,frame in c.groupby('syndrome_name'):
            close(bh(frame.interaction_p_value),frame.within_syndrome_fdr_q_value,'C within syndrome BH')
        rows.append({'component':layer,'status':'PASS','A_KW_tests':len(a),'A_pair_tests':len(pair),'B_tests':len(b),'C_tests':len(c),
                     'independent_C_method':'saturated four-cell means and sum(sample variance/(n-1)) for HC3'})


def selected_cis(disease,reference,features,results):
    source=RUN/'archive/historical_125/stats_testing/statistical_analysis.py'
    tree=ast.parse(source.read_text());names={'cohens_d','hedges_g','cliffs_delta_fast','bootstrap_effect_ci'}
    nodes=[node for node in tree.body if isinstance(node,ast.FunctionDef) and node.name in names]
    scope={'np':np,'pd':pd};exec(compile(ast.Module(body=nodes,type_ignores=[]),'<frozen-bootstrap>','exec'),scope)
    rng=np.random.default_rng(2026)
    for _ in features:
        for _ in range(500):
            rng.integers(0,len(disease),len(disease));rng.integers(0,len(reference),len(reference))
    for _ in range(200):
        rng.integers(0,len(disease),len(disease));rng.integers(0,len(reference),len(reference))
    first=sorted(results.disease.unique())[0]
    target=disease.loc[disease.disease.eq(first)]
    for feature in features[:3]:
        actual=scope['bootstrap_effect_ci'](target[feature].to_numpy(),reference[feature].to_numpy(),500,rng)
        row=results.loc[results.disease.eq(first)&results.feature.eq(feature)].iloc[0]
        expected=[row.hedges_g_ci_95_low,row.hedges_g_ci_95_high,row.cliffs_delta_ci_95_low,row.cliffs_delta_ci_95_high,row.effect_direction_bootstrap_stability]
        close(actual,expected,'actual-data CIs replayed with frozen literal code',atol=1e-10)
    return {'selected_actual_CI_replays':3,'seed':2026,'method':'Exact historical draw stream, including all combined-feature and rank-bootstrap draws before first syndrome'}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--scope',choices=['primary','full'],default='full');args=parser.parse_args()
    cohorts=RUN/'COHORTS';results=RUN/'RESULTS_120';rows=[]
    plan=json.loads((RUN/'ANALYSIS_PLAN.json').read_text())
    disease=pd.read_csv(cohorts/'primary_master.csv',dtype={'patient_id':'string','image_id':'string'})
    reference=pd.read_csv(cohorts/'reference_balanced_master.csv')
    features=pd.read_csv(cohorts/'exact_feature_names.csv').feature.tolist()
    primary=pd.read_csv(results/'primary_balanced/stats_each_disease_vs_healthy.csv')
    assert len(primary)==len(plan['targets'])*120==5760 and primary.feature.nunique()==120 and primary.disease.nunique()==48
    assert disease.patient_id.is_unique and reference.image_id.is_unique
    manifest=json.loads((RUN/'UNIFIED_120_MANIFEST.json').read_text())
    assert manifest['synthetic_rows_in_primary']==0
    error=main_tests(primary,disease,reference)
    rows.append({'component':'primary_main','status':'PASS','tests':len(primary),'max_errors':error})
    validate_demographics(disease,reference,features,rows)
    ci=selected_cis(disease,reference,features,primary)
    rows.append({'component':'literal_actual_bootstrap_CI_replay','status':'PASS',**ci})
    if args.scope=='full':
        jobs=[('reference_all_eligible','primary','reference_all_eligible'),('reference_race_sex','primary','reference_race_sex'),
              ('alternative_photo','alternate','reference_balanced')]
        for name,d,r in jobs:
            table=pd.read_csv(results/name/'stats_each_disease_vs_healthy.csv')
            x=pd.read_csv(cohorts/f'{d}_phenotypes.csv');y=pd.read_csv(cohorts/f'{r}_phenotypes.csv')
            errors=main_tests(table,x,y)
            rows.append({'component':name,'status':'PASS','tests':len(table),'max_errors':errors})
        extension=results/'extensions'
        for filename,comparator in [('within_gmdb_results.csv','all'),('within_gmdb_current_eligible_sensitivity.csv','restricted')]:
            table=pd.read_csv(extension/filename)
            for syndrome,frame in table.groupby('syndrome'):
                x=disease.loc[disease.disease.eq(syndrome)]
                y=disease.loc[~disease.disease.eq(syndrome)]
                if comparator=='restricted':y=y.loc[y.disease.isin(plan['targets'])]
                for item in frame.itertuples():
                    u,p,g=independent(x[item.feature].to_numpy(),y[item.feature].to_numpy())
                    close([u,p,g],[item.mannwhitney_u,item.p_value,item.hedges_g],filename)
            close(bh(table.p_value),table.q_global_bh,'within BH')
            rows.append({'component':filename,'status':'PASS','tests':len(table)})
        for filename,denominator in [('full_adjusted_screen.csv',5880),('targeted_pairwise.csv',22320),('age_stratified.csv',11760)]:
            table=pd.read_csv(extension/filename);assert len(table)==denominator
            valid=table.status.eq('ok')&table.p_value.notna()
            expected=bh(np.where(valid,table.p_value,1.));expected[~valid]=np.nan
            close(expected,table.q_global_bh,filename+' planned-family BH')
            if filename=='full_adjusted_screen.csv':
                for syndrome,frame in table.loc[valid].groupby('syndrome'):
                    item=frame.sort_values('feature').iloc[0]
                    complete=disease.loc[disease.sex_label.isin(['Male','Female'])&disease.conditional_age_years.notna()&disease[item.feature].notna()]
                    design=np.column_stack([np.ones(len(complete)),complete.disease.eq(syndrome),complete.conditional_age_years-6,complete.sex_label.eq('Female')])
                    coefficient,se=matrix_hc3(complete[item.feature].to_numpy(),design)
                    close([coefficient[1],se[1]],[item.target_beta,item.hc3_se],'independent adjusted HC3 matrix',atol=1e-8)
            elif filename=='targeted_pairwise.csv':
                for _,frame in table.loc[valid].groupby('pair_index'):
                    item=frame.sort_values('feature').iloc[0]
                    a=disease.loc[disease.disease.eq(item.syndrome_a),item.feature].to_numpy()
                    b=disease.loc[disease.disease.eq(item.syndrome_b),item.feature].to_numpy()
                    u,p,g=independent(a,b);close([u,p,g],[item.mannwhitney_u,item.p_value,item.hedges_g],'pairwise independent check')
            rows.append({'component':filename,'status':'PASS','planned_tests':denominator,'supported':int(valid.sum())})
        highlighted=pd.read_csv(extension/'within_gmdb_age_sex_adjusted.csv')
        for item in highlighted.itertuples():
            complete=disease.loc[disease.sex_label.isin(['Male','Female'])&disease.conditional_age_years.notna()&disease[item.feature].notna()]
            design=np.column_stack([np.ones(len(complete)),complete.disease.eq(item.syndrome),complete.conditional_age_years-6,complete.sex_label.eq('Female')])
            coefficient,se=matrix_hc3(complete[item.feature].to_numpy(),design)
            close([coefficient[1],se[1],2*norm.sf(abs(coefficient[1]/se[1]))],[item.target_beta,item.hc3_se,item.p_value],'four highlighted models')
        close(bh(highlighted.p_value),highlighted.q_four_models,'four-model BH')
        covariates=pd.read_csv(extension/'three_predefined_reference_adjusted_models.csv')
        for item in covariates.itertuples():
            case=disease.loc[disease.disease.str.replace(r'__[0-9a-f]{8}$','',regex=True).eq(item.syndrome)].copy()
            control=reference.copy();case['target']=1.;control['target']=0.
            frame=pd.concat([case,control],ignore_index=True)
            frame=frame.dropna(subset=[item.feature,'sex_label','broad_group','pose_yaw','pose_pitch','pose_roll'])
            frame=frame.loc[frame.sex_label.isin(['Female','Male'])&frame.broad_group.isin(['European','African','Asian','Others'])]
            groups=sorted(frame.broad_group.unique())
            design=np.column_stack([np.ones(len(frame)),frame.target,frame.sex_label.eq('Female'),
                *[frame.broad_group.eq(group) for group in groups[1:]],frame.pose_yaw,frame.pose_pitch,frame.pose_roll]).astype(float)
            coefficient,se=matrix_hc3(frame[item.feature].to_numpy(),design)
            close([coefficient[1],se[1],2*norm.sf(abs(coefficient[1]/se[1]))],
                [item.adjusted_disease_coefficient,item.hc3_standard_error,item.raw_p],'three reference covariate models',atol=1e-8)
        close(bh(covariates.raw_p),covariates.bh_q_across_three_predefined_models,'three-model BH')
        rows.append({'component':'all_seven_predefined_adjusted_models','status':'PASS','models':7})
        within=pd.read_csv(extension/'within_gmdb_results.csv')
        for index,syndrome in enumerate(plan['targets']):
            if syndrome not in set(highlighted.syndrome):continue
            feature=highlighted.loc[highlighted.syndrome.eq(syndrome),'feature'].iloc[0]
            a=disease.loc[disease.disease.eq(syndrome),feature].to_numpy(float)
            b=disease.loc[~disease.disease.eq(syndrome),feature].to_numpy(float)
            rng=np.random.default_rng([20261024,index])
            aa=a[rng.integers(len(a),size=(1000,len(a)))];bb=b[rng.integers(len(b),size=(1000,len(b)))]
            pooled=((len(a)-1)*aa.var(1,ddof=1)+(len(b)-1)*bb.var(1,ddof=1))/(len(a)+len(b)-2)
            draws=(1-3/(4*(len(a)+len(b))-9))*(aa.mean(1)-bb.mean(1))/np.sqrt(pooled)
            row=within.loc[within.syndrome.eq(syndrome)&within.feature.eq(feature)].iloc[0]
            close(np.quantile(draws,[.025,.975]),[row.ci_95_lower,row.ci_95_upper],'literal within-GMDB CI replay',atol=1e-8)
        rows.append({'component':'four_actual_within_GMDB_CI_replays','status':'PASS','bootstrap':1000})
        concord=pd.read_csv(extension/'concordance_summary.csv')
        assert not any('bootstrap' in c or '_ci_' in c for c in concord.columns)
        assert concord.uncertainty_status.str.contains('descriptive_point_estimates_only').all()
        entries=pd.read_csv(extension/'concordance_entries.csv')
        for mode,group in entries.groupby('comparator'):
            record=concord.loc[concord.comparator.eq(mode)&concord.level.eq('overall')].iloc[0]
            close((np.sign(group.fairface_g)==np.sign(group.within_g)).mean(),record.direction_agreement,'independent concordance direction')
            close(np.corrcoef(group.fairface_g,group.within_g)[0,1],record.pearson_r,'independent concordance Pearson')
            close(np.corrcoef(rankdata(group.fairface_g),rankdata(group.within_g))[0,1],record.spearman_rho,'independent concordance Spearman')
            close(np.median(group.fairface_g.abs()-group.within_g.abs()),record.median_absolute_attenuation,'independent attenuation')
        rows.append({'component':'descriptive_concordance_and_sampling','status':'PASS','uncertainty':'excluded_biased_bootstrap_intervals'})
        random=pd.read_csv(extension/'random_photo_replicates.csv')
        assert len(random)==100 and random.hypotheses.eq(5760).all()
    result={'status':'VALIDATED' if args.scope=='full' else 'PRIMARY_VALIDATED','components':rows,
        'primary_tests':len(primary),'primary_global_BH_significant':int(primary.mannwhitney_q_fdr_global.lt(.05).sum()),
        'primary_patients':len(disease),'reference_images':len(reference),'feature_count':120,'synthetic_primary_rows':0}
    (RUN/f'VALIDATION_{args.scope.upper()}.json').write_text(json.dumps(result,indent=2)+'\n')
    pd.DataFrame([{k:v for k,v in row.items() if not isinstance(v,dict)} for row in rows]).to_csv(RUN/f'validation_{args.scope}_summary.csv',index=False)
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    main()

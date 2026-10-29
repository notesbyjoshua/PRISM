import json
import os
import shutil
import textwrap
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
RUN=HERE/'local_only/run_20261028'
os.umask(0o077)
sys.path.insert(0,str(ROOT))
from publication.scripts.current_120 import publication_table,copy_tables


def clean(label):
    label=label.rsplit('__',1)[0]
    if len(label)>45 and ';' in label:
        label=label.rsplit(';',1)[-1].strip()
    return textwrap.shorten(label,width=48,placeholder='...')


def tex(value):
    return str(value).replace('\\',r'\textbackslash{}').replace('_',r'\_').replace('&',r'\&').replace('%',r'\%').replace('#',r'\#')


def main():
    validation=json.loads((RUN/'VALIDATION_FULL.json').read_text())
    assert validation['status']=='VALIDATED'
    plan=json.loads((RUN/'ANALYSIS_PLAN.json').read_text())
    results=RUN/'RESULTS_120';primary=pd.read_csv(results/'primary_balanced/stats_each_disease_vs_healthy.csv')
    mapping=pd.read_csv(RUN/'COHORTS/feature_anatomy_grouping.csv')[['feature','anatomical_region','feature_family']]
    publication=RUN/'ACTIVE_PUBLICATION';publication.mkdir(exist_ok=True)
    figures=publication/'figures';figures.mkdir(exist_ok=True)
    tables=publication/'supplementary';tables.mkdir(exist_ok=True)
    for path in tables.glob('*.csv'):path.unlink()
    report=RUN/'ACTIVE_STUDY_REPORT';report.mkdir(exist_ok=True)
    report_figures=report/'figures';report_figures.mkdir(exist_ok=True)
    manuscript=RUN/'MANUSCRIPT_120';manuscript.mkdir(exist_ok=True)
    significant=primary.mannwhitney_q_fdr_global.lt(.05)
    numbers={'primary_patients':1845,'reference_images':3409,'eligible_syndromes':48,'morphological_features':120,
        'primary_hypotheses':len(primary),'primary_significant':int(significant.sum()),
        'primary_significant_percent':float(100*significant.mean()),
        'large_significant':int((significant&primary.hedges_g.abs().ge(.8)).sum()),
        'positive_significant':int((significant&primary.hedges_g.gt(0)).sum()),
        'negative_significant':int((significant&primary.hedges_g.lt(0)).sum())}
    patient=pd.read_csv(RUN/'COHORTS/primary_master.csv',usecols=['disease'])
    numbers['retained_syndrome_labels']=patient.disease.nunique()
    numbers['patients_in_eligible_syndromes']=int(patient.disease.isin(plan['targets']).sum())
    clinical=[]
    highlights=[('Cornelia de Lange','philtrum_length','Philtrum / midface height'),('Noonan','inter_pupillary_distance','Interpupillary / bizygomatic width'),
        ('WILLIAMS-BEUREN','eb_thickness_mean','Eyebrow geometric area index'),('WILLIAMS-BEUREN','outer_canthal_distance','Outer-canthal / bizygomatic width'),
        ('Rubinstein-Taybi','malar_bulge_mean','Malar geometric index')]
    for syndrome,feature,label in highlights:
        row=primary.loc[primary.disease.str.contains(syndrome,case=False,regex=False)&primary.feature.eq(feature)].iloc[0]
        clinical.append({'syndrome':clean(row.disease),'feature':feature,'display_label':label,'n_disease':int(row.n_disease),
            'n_reference':int(row.n_healthy),'hedges_g':row.hedges_g,'ci_low':row.hedges_g_ci_95_low,'ci_high':row.hedges_g_ci_95_high,
            'q_global':row.mannwhitney_q_fdr_global})
    clinical=pd.DataFrame(clinical);clinical.to_csv(tables/'clinical_highlights_120.csv',index=False)
    primary.to_csv(tables/'S1_primary_all_5760_hypotheses.csv',index=False)
    schema=pd.read_csv(ROOT/'local_review/feature_extraction_120/20261028_current120/fairface/feature_schema.csv').merge(mapping,on='feature',validate='one_to_one')
    schema.to_csv(tables/'Supplementary_Table_S6_exact_feature_definitions.csv',index=False)
    region=primary.merge(mapping,on='feature',validate='many_to_one').assign(significant=significant.to_numpy()).groupby('anatomical_region').agg(
        tested=('feature','size'),significant=('significant','sum'),features=('feature','nunique')).reset_index()
    region['fraction_significant']=region.significant/region.tested;region.to_csv(tables/'regional_summary_120.csv',index=False)
    demographic=[]
    for layer in ['layer1_broad_ethnicity','layer2_subcategory']:
        folder=results/'demographics'/layer
        for analysis,file,p,q in [('A_omnibus','analysis_A_healthy_kruskal_all_groups.csv','p_value','fdr_q_value'),
            ('B_matched','analysis_B_matched_disease_vs_healthy.csv','p_value','global_fdr_q_value'),
            ('C_interaction','analysis_C_all_interactions.csv','interaction_p_value','global_fdr_q_value')]:
            frame=pd.read_csv(folder/file)
            demographic.append({'layer':layer,'analysis':analysis,'hypotheses':len(frame),'significant':int(frame[q].lt(.05).sum()),
                'family':'separate global BH within this layer and analysis'})
            shutil.copyfile(folder/file,tables/f'{layer}_{file}')
    demographic=pd.DataFrame(demographic);demographic.to_csv(tables/'demographic_family_summary.csv',index=False)
    sensitivity=[]
    for name in ['primary_balanced','reference_all_eligible','reference_race_sex','alternative_photo']:
        frame=pd.read_csv(results/name/'stats_each_disease_vs_healthy.csv')
        sensitivity.append({'specification':name,'features':frame.feature.nunique(),'syndromes':frame.disease.nunique(),'hypotheses':len(frame),
            'significant':int(frame.mannwhitney_q_fdr_global.lt(.05).sum()),'reference_n':int(frame.n_healthy.iloc[0])})
    sensitivity=pd.DataFrame(sensitivity);sensitivity.to_csv(report/'reference_and_photo_sensitivity.csv',index=False)
    alternate=pd.read_csv(results/'alternative_photo/stats_each_disease_vs_healthy.csv')
    transition=primary.merge(alternate,on=['disease','feature'],suffixes=('_primary','_alternative'),validate='one_to_one')
    a,b=transition.mannwhitney_q_fdr_global_primary.lt(.05),transition.mannwhitney_q_fdr_global_alternative.lt(.05)
    numbers['photo_significant_both']=int((a&b).sum());numbers['photo_primary_only']=int((a&~b).sum())
    numbers['photo_alternative_only']=int((~a&b).sum());numbers['photo_alternative_significant']=int(b.sum())
    numbers['alternative_photo_changes']=plan['alternative_photo_changes']
    pdidb=json.loads((RUN/'PDIDB/extraction_summary.json').read_text())
    numbers.update(pdidb_records=pdidb['records'],pdidb_quality_passed=pdidb['quality_passed'],demographic_families=demographic.to_dict('records'))
    extensions=[]
    for filename in ['within_gmdb_results.csv','within_gmdb_current_eligible_sensitivity.csv','full_adjusted_screen.csv','targeted_pairwise.csv','age_stratified.csv']:
        frame=pd.read_csv(results/'extensions'/filename)
        q='q_global_bh';supported=frame.estimable if 'estimable' in frame else frame.status.eq('ok')
        extensions.append({'analysis':filename,'planned':len(frame),'supported':int(supported.sum()),'significant':int((supported&frame[q].lt(.05)).sum())})
        shutil.copyfile(results/'extensions'/filename,tables/filename)
    numbers['extension_families']=extensions
    for filename in ['within_gmdb_age_sex_adjusted.csv','three_predefined_reference_adjusted_models.csv','concordance_summary.csv','concordance_entries.csv','random_photo_sensitivity.csv','random_photo_replicates.csv','age_stratum_agreement.csv','redundancy_cluster_discoveries.csv']:
        publication_table(pd.read_csv(results/'extensions'/filename),concordance=filename=='concordance_summary.csv').to_csv(tables/filename,index=False)
    for path in list(tables.glob('*.csv')):
        frame=publication_table(pd.read_csv(path),concordance=path.name=='concordance_summary.csv')
        name=path.name.replace('healthy','reference')
        if name!=path.name:path.unlink()
        frame.to_csv(tables/name,index=False)
    (report/'REPORT_NUMBERS.json').write_text(json.dumps(numbers,indent=2)+'\n')
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'font.family':'DejaVu Sans'})
    def save(fig,name,active=True):
        destination=figures if active else report_figures
        fig.savefig(destination/(name+'.png'),dpi=300,bbox_inches='tight')
        fig.savefig(destination/(name+'.svg'),bbox_inches='tight')
        plt.close(fig)
    counts=primary.loc[significant].groupby('feature').size().sort_values(ascending=False,kind='stable')
    selected=counts.head(12).index.tolist()
    matrix=primary.pivot(index='disease',columns='feature',values='hedges_g').reindex(columns=selected)
    fig,ax=plt.subplots(figsize=(8.1,11.5),constrained_layout=True)
    image=ax.imshow(matrix.to_numpy(),cmap='RdBu_r',vmin=-2,vmax=2,aspect='auto')
    ax.set_yticks(range(len(matrix)),[clean(value) for value in matrix.index],fontsize=7.5)
    ax.set_xticks(range(len(selected)),[value.replace('_',' ') for value in selected],rotation=65,ha='right',fontsize=7.5)
    ax.set_title('Current 120-feature analysis: 48 eligible syndromes\n12 features with most global-BH discoveries; full table supplied',loc='left',fontsize=10)
    fig.colorbar(image,ax=ax,shrink=.5,label="Syndrome-minus-reference Hedges' g")
    save(fig,'fig1_cross_syndrome_heatmap')
    matrix.to_csv(tables/'figure1_heatmap_source.csv')
    fig,ax=plt.subplots(figsize=(7.3,3.8),constrained_layout=True)
    y=np.arange(len(clinical));values=clinical.hedges_g.to_numpy()
    ax.errorbar(values,y,xerr=[values-clinical.ci_low,clinical.ci_high-values],fmt='o',color='#176b80',capsize=4)
    labels=[clean(row.syndrome)+'\n'+row.display_label for row in clinical.itertuples()]
    ax.set_yticks(y,labels,fontsize=8);ax.invert_yaxis();ax.axvline(0,color='#777777',lw=.8)
    ax.set_xlabel("Hedges' g with marginal 95% bootstrap intervals (500 draws)")
    ax.set_title('Five predefined examples; reference n = 3,409',loc='left')
    save(fig,'fig2_clinical_indices')
    fig,ax=plt.subplots(figsize=(6.6,3.5),constrained_layout=True)
    kw=demographic.loc[demographic.analysis.eq('A_omnibus')]
    ax.bar(['Broad groups','Detailed groups'],kw.significant,color=['#176b80','#669f85'])
    for i,value in enumerate(kw.significant):ax.text(i,value+1,str(value),ha='center')
    ax.set_ylim(0,120);ax.set_ylabel('BH-significant features / 120')
    ax.set_title('Reference-only omnibus tests; separate 120-test families',loc='left')
    save(fig,'fig3_reference_demographics')
    fig,ax=plt.subplots(figsize=(5.5,3.6),constrained_layout=True)
    values=np.array([[int((~a&~b).sum()),int((~a&b).sum())],[int((a&~b).sum()),int((a&b).sum())]])
    ax.imshow(values,cmap='Blues')
    for i in range(2):
        for j in range(2):ax.text(j,i,str(values[i,j]),ha='center',va='center',color='white' if values[i,j]>values.max()/2 else 'black',fontsize=13)
    ax.set_xticks([0,1],['Alternative: no discovery','Alternative: discovery']);ax.set_yticks([0,1],['Primary: no discovery','Primary: discovery'])
    ax.set_title('Photograph selection: 5,760 common hypotheses',loc='left')
    save(fig,'fig4_photo_selection')
    fig,ax=plt.subplots(figsize=(6.5,3.6),constrained_layout=True)
    ss=sensitivity.loc[sensitivity.specification.isin(['primary_balanced','reference_all_eligible','reference_race_sex'])]
    labels=['Balanced\nn=3,409','All eligible\nn=5,864','Recorded category x sex\nn=2,114']
    ax.bar(labels,ss.significant,color='#176b80');ax.set_ylim(0,1100);ax.set_ylabel('Global-BH discoveries / 5,760')
    for i,value in enumerate(ss.significant):ax.text(i,value+12,str(value),ha='center')
    ax.set_title('Overlapping reference specifications; same current cases',loc='left')
    save(fig,'reference_sensitivity',False)
    sign_rows='\n'.join(f'{tex(row.syndrome)} & {tex(row.display_label)} & {row.hedges_g:.3f} & [{row.ci_low:.3f}, {row.ci_high:.3f}] & {row.q_global:.2g} \\\\' for row in clinical.itertuples())
    dem_rows='\n'.join(f'{tex(row.layer.replace("layer1_broad_ethnicity","Broad").replace("layer2_subcategory","Detailed"))} & {tex(row.analysis)} & {row.hypotheses:,} & {row.significant:,} \\\\' for row in demographic.itertuples())
    source=(HERE/'manuscript_120.tex').read_text()
    source=source.replace('CLINICAL_ROWS',sign_rows).replace('DEMOGRAPHIC_ROWS',dem_rows)
    (manuscript/'main.tex').write_text(source)
    shutil.copytree(figures,manuscript/'figures',dirs_exist_ok=True)
    (publication/'README.md').write_text('# Active 120-feature publication materials\n\nGenerated from the validated 48-syndrome, 5,760-test current analysis. Primary reference: 3,409 FairFace image paths. No synthetic observations in primary inference. Feature definitions and figures are current; older source material is archived separately. No patient photographs or per-patient tables are included.\n')
    (report/'FIGURE_SOURCES.json').write_text(json.dumps({'active_figures':{'fig1_cross_syndrome_heatmap':'primary_balanced/stats_each_disease_vs_healthy.csv; top 12 by discovery frequency',
        'fig2_clinical_indices':'supplementary/clinical_highlights_120.csv','fig3_reference_demographics':'demographic_family_summary.csv',
        'fig4_photo_selection':'primary vs alternative_photo on all 5,760 matched hypotheses'},
        'reference_sensitivity':'three overlapping FaceKit reference specifications'},indent=2)+'\n')
    print(json.dumps({'status':'private_current_materials_created','primary':numbers['primary_significant'],'hypotheses':len(primary),'clinical_examples':5},indent=2))


if __name__=='__main__':
    main()

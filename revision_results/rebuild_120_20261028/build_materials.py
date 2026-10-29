import json
import os
import shutil
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr,spearmanr

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
RUN=HERE/'local_only/run_20261028'
os.umask(0o077)


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
    report=RUN/'REPORT';report.mkdir(exist_ok=True)
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
    for name in ['primary_balanced','reference_all_eligible','reference_race_sex','reference_fixed_historical','alternative_photo','fixed_common_new','fixed_common_old']:
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
    old=pd.read_csv(results/'fixed_common_old/stats_each_disease_vs_healthy.csv')
    new=pd.read_csv(results/'fixed_common_new/stats_each_disease_vs_healthy.csv')
    paired=old.merge(new,on=['disease','feature'],suffixes=('_old','_new'),validate='one_to_one')
    meaningful=paired.hedges_g_old.abs().ge(.2)
    comparable=paired.loc[meaningful]
    same=np.sign(comparable.hedges_g_old)==np.sign(comparable.hedges_g_new)
    difference=(comparable.hedges_g_new-comparable.hedges_g_old).abs()
    tolerance=np.maximum(.2,.25*comparable.hedges_g_old.abs())
    old_sig=paired.mannwhitney_q_fdr_global_old.lt(.05);new_sig=paired.mannwhitney_q_fdr_global_new.lt(.05)
    comparison={'comparison_type':'matched anatomical names only; complete measurement definitions changed; not certified interchangeable',
        'paired_hypotheses':len(paired),'meaningful_old_effects':int(meaningful.sum()),'direction_agreement_n':int(same.sum()),
        'direction_agreement_fraction':float(same.mean()),'effect_tolerance_n':int(difference.le(tolerance).sum()),
        'effect_tolerance_fraction':float(difference.le(tolerance).mean()),'old_matched_significant_n':int(old_sig.sum()),
        'retained_significant_n':int((old_sig&new_sig).sum()),'retained_fraction':float((old_sig&new_sig).sum()/old_sig.sum()),
        'pearson':float(pearsonr(paired.hedges_g_old,paired.hedges_g_new).statistic),
        'spearman':float(spearmanr(paired.hedges_g_old,paired.hedges_g_new).statistic),
        'median_absolute_change':float((paired.hedges_g_new-paired.hedges_g_old).abs().median()),
        'p90_absolute_change':float((paired.hedges_g_new-paired.hedges_g_old).abs().quantile(.9)),
        'strict_unconditional_full_pipeline_equivalences_certified':0,
        'promotion_rule':'User authorized a validated new 120-feature primary analysis; earlier 90% interchangeability rule superseded',
        'fixed_photo_gmdb_n':plan['fixed_common_gmdb_n'],'fixed_photo_reference_n':plan['fixed_common_reference_n']}
    (report/'historical_named_feature_comparison.json').write_text(json.dumps(comparison,indent=2)+'\n')
    paired.to_csv(report/'paired_old_new_common_record_results.csv.gz',index=False)
    pdidb=json.loads((RUN/'PDIDB/extraction_summary.json').read_text())
    numbers.update(pdidb_records=pdidb['records'],pdidb_quality_passed=pdidb['quality_passed'],demographic_families=demographic.to_dict('records'))
    extensions=[]
    for filename in ['within_gmdb_results.csv','within_gmdb_current_eligible_sensitivity.csv','full_adjusted_screen.csv','targeted_pairwise.csv','age_stratified.csv']:
        frame=pd.read_csv(results/'extensions'/filename)
        q='q_global_bh';supported=frame.estimable if 'estimable' in frame else frame.status.eq('ok')
        extensions.append({'analysis':filename,'planned':len(frame),'supported':int(supported.sum()),'significant':int((supported&frame[q].lt(.05)).sum())})
        shutil.copyfile(results/'extensions'/filename,tables/filename)
    numbers['extension_families']=extensions
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
    fig,ax=plt.subplots(figsize=(6.5,4.2),constrained_layout=True)
    ax.scatter(paired.hedges_g_old,paired.hedges_g_new,s=5,alpha=.25,color='#176b80',rasterized=True)
    limits=[min(paired.hedges_g_old.min(),paired.hedges_g_new.min()),max(paired.hedges_g_old.max(),paired.hedges_g_new.max())]
    ax.plot(limits,limits,color='#888888',lw=1);ax.set_xlabel('Archived measurement g on matched records');ax.set_ylabel('Current measurement g on matched records')
    ax.set_title('5,760 matched-name tests; definitions may differ\n1,844 GMDB records and 941 reference paths',loc='left')
    save(fig,'historical_current_named_effects',False)
    fig,ax=plt.subplots(figsize=(6.5,3.6),constrained_layout=True)
    ss=sensitivity.loc[sensitivity.specification.isin(['primary_balanced','reference_all_eligible','reference_race_sex','reference_fixed_historical'])]
    labels=['Balanced\nn=3,409','All eligible\nn=5,864','Race x sex\nn=2,114','Historical paths\nn=941']
    ax.bar(labels,ss.significant,color='#176b80');ax.set_ylim(0,1100);ax.set_ylabel('Global-BH discoveries / 5,760')
    for i,value in enumerate(ss.significant):ax.text(i,value+12,str(value),ha='center')
    ax.set_title('Overlapping reference specifications; same current cases',loc='left')
    save(fig,'reference_sensitivity',False)
    sign_rows='\n'.join(f'{tex(row.syndrome)} & {tex(row.feature.replace("_"," "))} & {row.hedges_g:.3f} & [{row.ci_low:.3f}, {row.ci_high:.3f}] & {row.q_global:.2g} \\\\' for row in clinical.itertuples())
    dem_rows='\n'.join(f'{tex(row.layer.replace("layer1_broad_ethnicity","Broad").replace("layer2_subcategory","Detailed"))} & {tex(row.analysis)} & {row.hypotheses:,} & {row.significant:,} \\\\' for row in demographic.itertuples())
    source=r'''\documentclass[11pt]{article}
\usepackage[margin=1in]{geometry}
\usepackage{graphicx,booktabs,longtable}
\usepackage[hidelinks]{hyperref}
\title{Quantitative Facial Phenotyping Across Rare Syndromes Using a Rebuilt 120-Feature PRISM Workflow}
\author{Blinded for review}
\date{}
\begin{document}
\maketitle
\begin{abstract}
We rebuilt PRISM using 120 current FaceKit measurements from authorized restored GMDB images and an expanded FairFace reference. The age-eligible disease pool contained 1,845 unique recorded patients; 48 syndrome groups met the unchanged minimum size for testing against 3,409 reference images. Among 5,760 two-sided syndrome--feature tests, 955 met the global Benjamini--Hochberg threshold, including 482 with $|g|\geq0.8$. Reference, photograph, demographic and within-GMDB analyses were recomputed separately. Synthetic PDIDB records were excluded from all real-patient inference. These are exploratory geometric associations, not diagnostic-accuracy or clinical-validation results.
\end{abstract}
\section{Introduction}
Facial phenotype measurements can provide reproducible descriptions of variation across rare-disease image cohorts, while remaining sensitive to landmark detection, normalization and acquisition conditions. FaceKit derives interpretable geometric measurements from facial landmarks \cite{facekit2026}. PRISM examines cohort-level associations using GMDB \cite{gestaltmatcher} and FairFace \cite{fairface}. The present workflow measures current features from source images and recomputes every reported statistical family; it does not translate earlier feature tables or p-values.
\section{Materials and Methods}
\subsection{Sources, authorization and analytical roles}
Restricted restored GMDB images were processed locally under the author-provided mentor/IRB research authorization. Real patient images are not included in publication figures or released with this work. FairFace supplied reference images with recorded category labels rather than clinically screened healthy subjects. The official padding-0.25 archive and both label files were verified, including ZIP CRC integrity and complete path correspondence. All 11,764 unique age-category 3--9 candidates were screened before sampling. A balanced-by-source-category reference included up to 500 accepted images per category without replacement: 3,409 in total, including 409 Middle Eastern images and 500 in each other category. Alternatives used all 5,864 accepted paths or 151 in each of 14 category-by-recorded-sex cells, totaling 2,114. These overlapping specifications are not independent confirmations.
\subsection{Current feature extraction and cohort selection}
One pinned FaceKit revision, \texttt{5764bdd4af9360b55a046a16579fb5a890a8182f}, and one checksummed MediaPipe Face Landmarker float16 version-1 asset were used for all sources. Images were decoded to RGB without additional crop, resize, restoration or EXIF rotation. A two-face pass screened detected ambiguity; the standard one-face detector supplied 478 floating-point pixel coordinates and a facial transformation matrix. Current canonicalization uses an SVD-derived pure rotation and fixed canonical depth. Accepted measurements required full pose correction, exactly 120 finite outputs, and absolute yaw/pitch/roll no greater than $15^\circ/15^\circ/10^\circ$. Angular outputs are degrees; other outputs are unitless indices with feature-specific denominators. The exact paper-execution commit is not independently established. Supplementary Table S6 supplies current formula pointers and hashes.
GMDB ages use recorded year plus month/12, requiring nonnegative integer components, month 0--11 and $3\leq\mathrm{age}<10$; missing, invalid and zero-zero records are excluded. One eligible image per recorded patient was selected by minimum $\max(|yaw|/15,|pitch|/15,|roll|/10)$, then lexical image ID. Alternative selection uses the first lexical eligible image. The current pool has 1,845 patients. Six historical patients are no longer eligible and one newly measurable patient enters; one shared patient has a changed primary photograph. FBXW7 has nine retained patients and is below the unchanged primary minimum of ten. FairFace age remains the supplied category, not an inferred exact age. Distinct reference paths do not establish distinct people.
\subsection{Statistical specification}
The primary family comprises 48 eligible syndromes times 120 measurements: 5,760 two-sided Mann--Whitney tests. Effects are syndrome-minus-reference Hedges' $g$, with the historical small-sample correction. Feature eligibility requires at least ten finite observations per side and no more than 30\% missingness. One global BH correction spans the primary tested family. Marginal percentile intervals use 500 independent group-row bootstrap resamples, seed 2026; ranking stability uses 200 draws. A rank-count implementation accelerates Cliff's delta bootstraps while preserving the original interleaved random draws exactly, verified against frozen literal code.
Analyses A, B and C retain the frozen broad/detailed label maps and separate correction scopes. A uses reference-only Kruskal--Wallis tests with a 120-feature family per layer, and pairwise correction within each reference-group pair. B compares syndrome and reference observations within recorded groups, with 1,000 bootstrap draws and global BH per layer. C uses HC3 disease-by-group interaction coefficients with global BH per layer. The original computational support minimum is two disease and twenty reference observations in required cells; sparse estimates remain exploratory, not causal demographic effects.
Within-GMDB analyses compare each current eligible syndrome with all other age-eligible GMDB patients and with other current eligible syndromes, using separate 5,760-test families. The approved adjusted extension retains 49 historical target slots (5,880 planned models), 186 fixed anchor pairs (22,320 hypotheses) and two age strata (11,760 hypotheses). Unsupported cells retain correction p=1 and missing estimates/q-values. Four highlighted within-GMDB and three reference-covariate models keep separate four- and three-model families. Joint concordance bootstraps use 1,000 shared patient-row/reference-path resamples. Correlation clustering is descriptive and does not reduce BH denominators.
\section{Results}
\subsection{Primary associations}
The 48 eligible syndrome groups produced 5,760 primary hypotheses; 955 were global-BH significant (16.6\%), including 482 with $|g|\geq0.8$. These counts describe correlated measurements, not independent anatomical discoveries. The full table supplies negative and nonsignificant findings as well as supported associations. Figure~\ref{fig:clinical} presents the five predefined examples in their current measurement units.
\begin{table}[htbp]\centering\small
\begin{tabular}{p{2.8cm}p{3.8cm}rp{2.7cm}r}\toprule Syndrome & Measurement & $g$ & 95\% CI & Global $q$\\\midrule
CLINICAL_ROWS
\bottomrule\end{tabular}
\caption{Current predefined examples. The intervals are marginal bootstrap intervals; significance uses the complete 5,760-test primary family. Indices are normalized geometry rather than direct clinical diagnoses.}
\end{table}
\begin{figure}[htbp]\centering\includegraphics[width=\linewidth]{figures/fig2_clinical_indices.png}\caption{Current syndrome-minus-reference effects and marginal 95\% bootstrap intervals. No patient photographs are shown.}\label{fig:clinical}\end{figure}
\subsection{Demographic and sensitivity analyses}
\begin{table}[htbp]\centering
\begin{tabular}{llrr}\toprule Layer & Analysis & Hypotheses & BH discoveries\\\midrule
DEMOGRAPHIC_ROWS
\bottomrule\end{tabular}
\caption{Separate demographic correction families. Counts are not additive and overlapping observations do not provide independent replication.}\end{table}
Reference alternatives yield 971 discoveries with all eligible reference paths, 945 with category-by-sex balance, and 919 with the 941 currently passing historical reference paths. These specifications reuse cases and overlapping controls. The alternate eligible-photograph selection yields 996 discoveries; common and discordant hypotheses are supplied separately. The full within-GMDB, adjusted, pairwise and age-stratum tables retain declared planned denominators and unsupported cells.
\subsection{Synthetic source catalog}
PDIDB revision \texttt{d0907a56525dc3b25712d360f9ce3031150eb535} contains 1,727 author-designated synthetic records; 1,563 pass current QC. Its schema provides v1/v2 model tags and verification annotations, but no source-patient identifiers or substantive generation/license description. Demographic fields are missing for all records. PDIDB contributes zero observations to primary or demographic real-patient inference. Illustration redistribution and source independence remain unresolved; no PDIDB images are reproduced here.
\section{Discussion}
The rebuilt analysis supports quantitative cohort associations under a fixed, auditable current feature definition. The expanded reference and image-quality rules change the study population, so the new headline total should not be interpreted as a simple software-version effect. Matched-record-key and reference-size benchmarks are supplied as methodological sensitivity analyses. They do not establish interchangeability of changed feature formulas, exact historical pixel identity, or clinical validity. Within-GMDB comparisons address source heterogeneity but still share acquisition and annotation limitations.
\section{Limitations}
Disease status remains confounded with GMDB-versus-FairFace source. Landmark errors, restoration, image selection, residual pose, expression and occlusion can affect measured geometry. Ratios can change through their denominators; the indices are not calibrated three-dimensional anthropometry. FairFace labels are annotations rather than verified ancestry, exact ages or clinical health assessments. Recorded GMDB IDs support the implemented one-patient rule but do not establish undocumented biological independence; reference identities are unavailable. HC3 inference and marginal bootstrap intervals may be unreliable in sparse subgroups. The frozen concordance bootstrap has substantial bias: observed direction agreement is 0.713, while its percentile interval is [0.731, 0.775], excluding the observed point; Pearson intervals show the same problem. These intervals are retained with bias flags and should not be interpreted as reliable uncertainty confirmation. Correlated features and overlapping comparators do not establish realized FDR. Synthetic origin, generation lineage and publication rights for PDIDB require further provider/PI review.
\section{Conclusion}
PRISM now provides a validated 120-feature computational analysis with explicit source roles, quality accounting and versioned outputs. The findings are exploratory geometric associations; clinical interpretation and diagnostic utility require separate validation.
\section{Data, Code and Ethics}
Real GMDB images, patient identifiers and row-level measurements remain in restricted local storage. Aggregate artifacts require author privacy review before release. Historical analyses remain in a checksummed read-only archive; the active results and figure generators use current 120-feature outputs. The live Google Doc was not edited. [AUTHOR REVIEW REQUIRED: verify IRB identifier, exact permitted study scope, funding, contributions and final data/code release statement before submission.]
\begin{thebibliography}{9}
\bibitem{facekit2026} H. Chen et al., ``FaceKit: a Toolkit for Interpretable Facial Phenotyping, Synthetic Image Generation and Privacy Analysis in Rare Diseases,'' arXiv:2610.09288v1, 2026.
\bibitem{gestaltmatcher} T.-C. Hsieh et al., ``GestaltMatcher facilitates rare disease matching using facial phenotype descriptors,'' Nature Genetics, vol. 54, pp. 349--357, 2022.
\bibitem{fairface} K. Karkkainen and J. Joo, ``FairFace: Face Attribute Dataset for Balanced Race, Gender, and Age for Bias Measurement and Mitigation,'' WACV, pp. 1548--1558, 2021.
\end{thebibliography}
\end{document}
'''
    source=source.replace('CLINICAL_ROWS',sign_rows).replace('DEMOGRAPHIC_ROWS',dem_rows)
    (manuscript/'main.tex').write_text(source)
    shutil.copytree(figures,manuscript/'figures',dirs_exist_ok=True)
    (publication/'README.md').write_text('# Active 120-feature publication materials\n\nGenerated from the validated 48-syndrome, 5,760-test current analysis. Primary reference: 3,409 FairFace image paths. No synthetic observations in primary inference. Feature definitions and figures are current; older source material is archived separately. No patient photographs or per-patient tables are included.\n')
    (report/'FIGURE_SOURCES.json').write_text(json.dumps({'active_figures':{'fig1_cross_syndrome_heatmap':'primary_balanced/stats_each_disease_vs_healthy.csv; top 12 by discovery frequency',
        'fig2_clinical_indices':'supplementary/clinical_highlights_120.csv','fig3_reference_demographics':'demographic_family_summary.csv',
        'fig4_photo_selection':'primary vs alternative_photo on all 5,760 matched hypotheses'},
        'historical_named_effects':'fixed_common_old vs fixed_common_new; names not certified equivalent',
        'reference_sensitivity':'four overlapping new-extractor reference specifications'},indent=2)+'\n')
    print(json.dumps({'status':'private_current_materials_created','primary':numbers['primary_significant'],'hypotheses':len(primary),'clinical_examples':5},indent=2))


if __name__=='__main__':
    main()

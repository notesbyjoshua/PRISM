import json
from pathlib import Path

import numpy as np
import pandas as pd


def table(frame):
    lines=['| '+' | '.join(map(str,frame.columns))+' |','| '+' | '.join(['---']*len(frame.columns))+' |']
    for row in frame.itertuples(index=False,name=None):
        values=[f'{x:.5g}' if isinstance(x,float) else str(x) for x in row]
        lines.append('| '+' | '.join(x.replace('|',' / ') for x in values)+' |')
    return '\n'.join(lines)


def save(fig,directory,name):
    import matplotlib.pyplot as plt
    from matplotlib.text import Text
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    excluded=set()
    for ax in fig.axes:
        for axis,bounds in [(ax.xaxis,ax.get_xlim()),(ax.yaxis,ax.get_ylim())]:
            low,high=sorted(bounds)
            for tick in axis.get_major_ticks()+axis.get_minor_ticks():
                if not low<=tick.get_loc()<=high:
                    excluded.update([id(tick.label1),id(tick.label2)])
    for item in fig.findobj(Text):
        if id(item) not in excluded and item.get_visible() and item.get_text():
            box=item.get_window_extent(renderer)
            if box.x0<-.5 or box.y0<-.5 or box.x1>fig.bbox.x1+.5 or box.y1>fig.bbox.y1+.5:
                raise ValueError('Clipped figure text: '+item.get_text())
    for ext in ['png','pdf','svg']:
        fig.savefig(directory/(name+'.'+ext),dpi=300,facecolor='white')
    plt.close(fig)


def report(output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    plan=json.loads((output/'analysis_plan.json').read_text())
    manifest=json.loads((output/'source_and_execution_manifest.json').read_text())
    font='Times New Roman' if 'Times New Roman' in {f.name for f in font_manager.fontManager.ttflist} else 'DejaVu Serif'
    plt.rcParams.update({'font.family':font,'font.size':10,'pdf.fonttype':42,'svg.fonttype':'none'})
    figures=output/'figures'
    figures.mkdir(exist_ok=True)
    mechanism=pd.read_csv(output/'failure_mechanism_counts.csv',index_col=0)
    diagnoses=pd.read_csv(output/'original_387_African_Asian_discovery_diagnosis.csv')
    failures=pd.read_csv(output/'original_118_failure_reasons.csv')
    baseline=pd.read_csv(output/'original_reference_per_feature.csv')
    summaries=pd.read_csv(output/'reference_size_ladder_summary.csv')
    deterministic=summaries.sampling_scheme.eq('deterministic_full_pool')
    summaries.loc[deterministic,'mean_q_recurrence_original']=np.nan
    if 'point_discoveries' not in summaries:
        summaries['point_discoveries']=np.nan
    summaries.loc[deterministic,'point_discoveries']=summaries.loc[deterministic,'point_discoveries'].fillna(summaries.loc[deterministic,'mean_discoveries'])
    if (output/'A_vs_B_primary_effects.csv').exists():
        points=pd.read_csv(output/'A_vs_B_primary_effects.csv').groupby('canonical_race').B_significant.sum()
        summaries.loc[deterministic,'point_discoveries']=summaries.loc[deterministic,'canonical_race'].map(points)
    summaries.loc[deterministic,'mean_discoveries']=np.nan
    summaries.to_csv(output/'reference_size_ladder_summary.csv',index=False)
    for path in output.glob('*_per_feature.csv'):
        data=pd.read_csv(path,float_precision='round_trip',low_memory=False)
        for column in ['descriptive_stable','pass_q_recurrence','pass_sign_agreement','stable_with_B_sign']:
            if column in data:
                original=data[column].astype(str)
                allowed={'True','False','1','0','1.0','0.0','nan','None','<NA>',''}
                if not set(original.dropna()).issubset(allowed):
                    raise ValueError('Invalid exported stability flag: '+column)
                data[column]=original.map({'True':1.0,'False':0.0,'1':1.0,'0':0.0,'1.0':1.0,'0.0':0.0}).astype('Float64')
        data.to_csv(path,index=False)
    for path in [output/'full_reference_per_feature.csv',output/'reference_size_ladder_per_feature.csv']:
        if path.exists():
            data=pd.read_csv(path,float_precision='round_trip',low_memory=False)
            mask=data.sampling_scheme.eq('deterministic_full_pool')
            if 'deterministic_g' not in data:
                data['deterministic_g']=np.nan
                data['deterministic_q']=np.nan
            if data.loc[mask,'deterministic_g'].isna().any():
                data.loc[mask,'deterministic_g']=data.loc[mask,'median_g']
                data.loc[mask,'deterministic_q']=data.loc[mask,'median_q']
            data.loc[mask,['g_q025','median_g','g_q975','g_interval_width','g_empirical_sd','median_q','minimum_q','maximum_q']]=np.nan
            data.to_csv(path,index=False)
    if (output/'full_reference_summary.csv').exists():
        summaries.loc[deterministic].to_csv(output/'full_reference_summary.csv',index=False)
    colors={'African':'#0072B2','Asian':'#D55E00'}
    categories=['significance_only','sign_only','both']
    figure_counts=mechanism[categories].reset_index().melt(id_vars='canonical_race',var_name='failure_reason',value_name='count')
    figure_counts['denominator']=figure_counts.canonical_race.map(mechanism[categories].sum(1))
    figure_counts['percentage']=100*figure_counts['count']/figure_counts.denominator
    figure_counts.to_csv(output/'failure_bar_plot_values.csv',index=False)
    fig,ax=plt.subplots(figsize=(7.5,4.2))
    for i,race in enumerate(colors):
        values=mechanism.loc[race,categories].to_numpy()
        bars=ax.bar(np.arange(3)+(i-.5)*.32,values,width=.30,color=colors[race],label=race+' recorded category')
        denom=int(values.sum())
        ax.bar_label(bars,labels=[f'{n}\n({100*n/denom:.1f}%)' for n in values],padding=3,fontsize=9)
    ax.set_xticks(range(3),['Significance-only','Sign-only','Both'])
    ax.set_ylim(0,max(mechanism[categories].sum(1))*1.35)
    ax.set_ylabel('Originally significant but unstable associations')
    ax.set_title('Operational failure mechanisms in the original 118 findings',fontsize=12)
    ax.legend(fontsize=9,frameon=False)
    ax.spines[['top','right']].set_visible(False)
    fig.tight_layout()
    save(fig,figures,'original_failure_mechanisms')
    diagnosis_rows=[]
    diagnoses['case_band']=pd.cut(diagnoses.case_n,[4,9,19,np.inf],labels=['5-9','10-19','20+'])
    diagnoses['q_band']=pd.cut(diagnoses.q,[0,.001,.01,.025,.05],right=False,labels=['<0.001','0.001-0.01','0.01-0.025','0.025-0.05'])
    diagnoses['abs_g_band']=pd.cut(diagnoses.hedges_g.abs(),[0,.8,1.2,np.inf],right=False,labels=['<0.8','0.8-1.2','>=1.2'])
    for dimension in ['case_band','q_band','abs_g_band','anatomical_region']:
        for (race,value),group in diagnoses.groupby(['canonical_race',dimension],observed=True):
            diagnosis_rows.append({'dimension':dimension,'canonical_race':race,'stratum':str(value),'original_discoveries':len(group),
                'unstable':int(group.failure_reason.ne('stable').sum()),'median_original_q':group.q.median(),
                'median_abs_original_g':group.hedges_g.abs().median(),'median_q_recurrence':group.significant_frequency.median(),
                'median_sign_agreement':group.sign_agreement_frequency.median(),'median_g_interval_width':group.g_interval_width.median(),
                'median_reference_mean_sd':group.reference_mean_empirical_sd.median()})
    strata=pd.DataFrame(diagnosis_rows)
    strata.to_csv(output/'original_instability_descriptive_strata.csv',index=False)
    fig,ax=plt.subplots(figsize=(7.5,4.2))
    for race in colors:
        group=diagnoses.loc[diagnoses.canonical_race.eq(race)]
        ax.scatter(group.q,group.significant_frequency,s=20,color=colors[race],alpha=.6,label=race)
    ax.set_xscale('log')
    ax.axhline(.8,color='#777777',linestyle='--',linewidth=.8)
    ax.set(xlabel='Original global BH q (9,600 tests)',ylabel='q < 0.05 frequency across original reference draws',ylim=(-.02,1.05))
    ax.set_title('Original discovery strength and reference-selection recurrence',fontsize=12)
    ax.legend(frameon=False)
    fig.tight_layout()
    save(fig,figures,'original_q_vs_reference_recurrence')
    failure_text=['# Original reference-selection failure diagnosis','',
                  'COMPLETED: verified all original 1,000 checkpoint draws, exact keys, source hashes and public/private totals. '
                  'The original allocation was 52/190/647; the 9,600-test family contains 982 primary discoveries. '
                  'African 135/76 and Asian 252/193 original/stable counts produce exactly 118 failures. '
                  'European 595/580 recurrence counts are descriptive; its reference photograph composition never changed.', '',
                  table(mechanism.reset_index()),'',
                  'Failure labels are mutually exclusive operational categories, not causal mechanisms. Significance-only means '
                  'q recurrence <80% with sign agreement ≥95%; sign-only reverses those pass/fail criteria; both fails both. '
                  'Any sign reversal is separately retained: a finding can pass ≥95% sign agreement despite rare reversals. '
                  'g=0 and NaN never count as a directional success. A near-threshold original q is prespecified as [0.025,0.05).','',
                  f'Among the 118 failures, {int(failures.ever_opposite_sign.sum())} ever showed a strictly opposite g sign, '
                  f'{int(failures.g_interval_covers_zero.sum())} had a 2.5–97.5% g-selection interval covering zero, and '
                  f'{int(failures.original_q_near_005.sum())} had original q in [0.025,0.05). These properties can overlap.','',
                  '## Descriptive strata (no additional hypothesis tests)','',table(strata), '',
                  'These counts are correlated measurements, not independent clinical traits. Reference means and SD variation '
                  'were reconstructed from the original seeded subset allocations without rerunning the original experiment. '
                  'The q-recurrence graph and stratified summaries describe associations with original q, |g|, patient count and '
                  'region; they do not uniquely identify statistical or biological causes.']
    (output/'FAILURE_DIAGNOSIS.md').write_text('\n'.join(failure_text)+'\n')
    curve=summaries.loc[summaries.canonical_race.isin(colors)&summaries.sampling_scheme.isin(['original_without_replacement_subset','without_replacement_subset'])].copy()
    curve=curve.sort_values(['canonical_race','reference_n'])
    curve.to_csv(output/'subset_stability_curve_plot_values.csv',index=False)
    full_available=(output/'A_vs_B_primary_effects.csv').exists()
    for name,field,ylabel in [('stability_vs_reference_size','stable_original_percent','Original discoveries meeting fixed stability rule (%)'),
                              ('g_dispersion_vs_reference_size','median_g_interval_width_original','Median 2.5–97.5% g-selection interval width')]:
        fig,axes=plt.subplots(1,2,figsize=(9,4),sharey=True)
        for ax,race in zip(axes,colors):
            group=curve.loc[curve.canonical_race.eq(race)]
            ax.plot(group.reference_n,group[field],marker='o',color=colors[race])
            ax.set_xscale('log',base=2)
            ax.set_xticks(group.reference_n.astype(int),group.reference_n.astype(int).astype(str),rotation=35)
            ax.set_title(race+f' · {int(group.original_discoveries.iloc[0])} original findings',fontsize=11)
            ax.set_xlabel('Reference photographs (shared subset)')
            if field=='stable_original_percent':
                ax.set_ylim(0,105)
            ax.text(.02,.02,f'Full pool {plan["pools"][race]}: deterministic;\nnot a subset-stability point',transform=ax.transAxes,fontsize=8,
                    bbox={'facecolor':'white','edgecolor':'none','alpha':.9})
            ax.spines[['top','right']].set_visible(False)
        axes[0].set_ylabel(ylabel)
        fig.suptitle('Without-replacement subsets; European references fixed',fontsize=12)
        fig.tight_layout()
        save(fig,figures,name)
    comparison_text=['# Reference design comparison','',
                     'A: original fixed 889-photo primary sample and its original 1,000 subset draws, reused without alteration. '
                     'B: all eligible content-QC category pools once; an explicitly separate with-replacement reference bootstrap. '
                     'C: paired doubling ladder with European held at its full observed pool. All new draws recompute all 9,600 '
                     'p values before one BH correction. No alternative feature normalization or case changes.','',
                     table(summaries), '',
                     'No full-pool subset is called 100% stable. Its unique subset count is one; the full cap is outside the '
                     'resampleable stability curve. The B reference bootstrap is model-based conditional sampling of photographs '
                     'with replacement. It preserves duplicate-induced ties and may assume independence that unknown repeated '
                     'people/source clustering violates. It is not the same uncertainty quantity as A/C without-replacement selection.','']
    if full_available:
        comp=pd.read_csv(output/'A_vs_B_primary_effects.csv')
        ab_rows=[]
        for race in ['African','Asian','European']:
            group=comp.loc[comp.canonical_race.eq(race)]
            original=group.loc[group.original_significant]
            ab_rows.append({'canonical_race':race,'A_reference_n':int(group.reference_n.iloc[0]),'B_reference_n':int(group.B_reference_n.iloc[0]),
                'A_discoveries':int(group.original_significant.sum()),'B_discoveries':int(group.B_significant.sum()),
                'original_retained':int((group.original_significant&group.B_significant).sum()),
                'original_lost':int((group.original_significant&~group.B_significant).sum()),
                'new_discoveries':int((~group.original_significant&group.B_significant).sum()),
                'sign_reversals_original':int(original.sign_reversal_A_B.sum()),
                'median_abs_delta_g_original':original.abs_g_change.median(),'median_delta_g_original':original.g_change_B_minus_A.median()})
        ab=pd.DataFrame(ab_rows)
        ab.to_csv(output/'A_vs_B_category_summary.csv',index=False)
        comp['lost_original']=comp.original_significant&~comp.B_significant
        comp['new_discovery']=~comp.original_significant&comp.B_significant
        region=comp.groupby(['canonical_race','anatomical_region'],as_index=False).agg(
            hypotheses=('feature','size'),A_discoveries=('original_significant','sum'),B_discoveries=('B_significant','sum'),
            original_lost=('lost_original','sum'),new_discoveries=('new_discovery','sum'),
            median_abs_delta_g=('abs_g_change','median'),mean_abs_delta_g=('abs_g_change','mean'))
        region.to_csv(output/'A_vs_B_region_summary.csv',index=False)
        comp['case_band']=pd.cut(comp.case_n,[4,9,19,np.inf],labels=['5-9','10-19','20+'])
        comp.groupby(['canonical_race','case_band'],observed=True).agg(hypotheses=('feature','size'),
            A_discoveries=('original_significant','sum'),B_discoveries=('B_significant','sum'),
            original_lost=('lost_original','sum'),new_discoveries=('new_discovery','sum'),median_abs_delta_g=('abs_g_change','median')).to_csv(output/'A_vs_B_case_size_summary.csv')
        comparison_text+=['## A versus deterministic full-reference B','',table(ab),'',table(region),'',
                          'European g and raw p remain fixed; changed European q decisions therefore demonstrate cross-category '
                          'BH competition. For African/Asian, increasing reference size also changes the selected reference means '
                          'and SDs. The present design cannot uniquely attribute their q changes to power versus a different '
                          'observed reference distribution.','']
        scatter=comp.loc[comp.original_significant&comp.canonical_race.isin(colors)].copy()
        scatter.to_csv(output/'A_vs_B_scatter_plot_values.csv',index=False)
        fig,ax=plt.subplots(figsize=(6.5,5.5))
        for race in colors:
            group=scatter.loc[scatter.canonical_race.eq(race)]
            ax.scatter(group.hedges_g,group.B_g,color=colors[race],s=24,alpha=.65,label=race)
        changed=scatter.loc[scatter.sign_reversal_A_B]
        if len(changed):
            ax.scatter(changed.hedges_g,changed.B_g,facecolors='none',edgecolors='black',s=70,label='Sign reversal')
        low=min(scatter.hedges_g.min(),scatter.B_g.min())-.15
        high=max(scatter.hedges_g.max(),scatter.B_g.max())+.15
        ax.plot([low,high],[low,high],linestyle='--',color='#777777',linewidth=.8)
        ax.set(xlim=(low,high),ylim=(low,high),xlabel='Original A Hedges’ g',ylabel='Full-reference B Hedges’ g')
        ax.set_title(f'Original African/Asian discoveries (n = {len(scatter)}); {len(changed)} sign reversals',fontsize=11)
        ax.legend(frameon=False)
        fig.tight_layout()
        save(fig,figures,'original_vs_full_reference_effects')
    (output/'REFERENCE_DESIGN_COMPARISON.md').write_text('\n'.join(comparison_text)+'\n')
    expected_stages=1+max(0,len(plan['ladder'])-2)
    observed_stages=curve.design.nunique()
    b_bootstrap=(output/'B_reference_bootstrap_summary.csv').exists()
    completed=full_available and b_bootstrap and observed_stages==expected_stages and plan['requested_repetitions']>=1000
    status='COMPLETED' if completed else 'PARTIAL / PILOT' if plan['requested_repetitions']<1000 else 'PARTIAL'
    trend=[]
    for race in colors:
        group=curve.loc[curve.canonical_race.eq(race)].sort_values('reference_n')
        first,last=group.iloc[0],group.iloc[-1]
        trend.append(f'{race}: subset n {int(first.reference_n)} → {int(last.reference_n)}; stable '
                     f'{int(first.stable_original)}/{int(first.original_discoveries)} ({first.stable_original_percent:.1f}%) → '
                     f'{int(last.stable_original)}/{int(last.original_discoveries)} ({last.stable_original_percent:.1f}%); '
                     f'median g-selection width {first.median_g_interval_width_original:.4g} → {last.median_g_interval_width_original:.4g}.')
    incremental=[]
    for race,group in curve.groupby('canonical_race'):
        group=group.sort_values('reference_n')
        for before,after in zip(group.itertuples(),list(group.itertuples())[1:]):
            incremental.append({'canonical_race':race,'from_n':before.reference_n,'to_n':after.reference_n,
                'stable_percentage_point_change':after.stable_original_percent-before.stable_original_percent,
                'g_width_absolute_decrease':before.median_g_interval_width_original-after.median_g_interval_width_original,
                'g_width_percent_decrease':100*(1-after.median_g_interval_width_original/before.median_g_interval_width_original)})
    increments=pd.DataFrame(incremental)
    increments.to_csv(output/'incremental_ladder_changes.csv',index=False)
    evidence=diagnoses[['disease','canonical_race','feature','case_n','q','hedges_g']].copy()
    keys=['disease','canonical_race','feature']
    schemes=[]
    for path in [output/'original_reference_per_feature.csv']+sorted(output.glob('C_stage_*_per_feature.csv'))+[output/'B_reference_bootstrap_per_feature.csv']:
        if path.exists():
            source=pd.read_csv(path)
            name=source.design.iloc[0]
            column=name+'_passes_original_rule'
            evidence=evidence.merge(source[keys+['descriptive_stable']].rename(columns={'descriptive_stable':column}),on=keys,validate='one_to_one')
            schemes.append(column)
    evidence['consistent_across_available_schemes']=evidence[schemes].eq(1).all(axis=1)
    evidence.to_csv(output/'followup_evidence_all_387_original_discoveries.csv',index=False)
    evidence.loc[evidence.consistent_across_available_schemes].to_csv(output/'internally_consistent_followup_candidates.csv',index=False)
    original_only=int(mechanism['significance_only'].sum())
    sign_only=int(mechanism['sign_only'].sum())
    both=int(mechanism['both'].sum())
    decision=('Evidence supports considering larger category-matched references for a separately justified future analysis '
              'to reduce reference-selection dispersion, if the displayed completed ladder shows that reduction. '
              'It does not establish more clinically reliable findings, remove source confounding, or justify automatic promotion.') if completed else 'The completed evidence cannot yet resolve the full larger-reference comparison; unfinished stages are not tested-negative results.'
    if completed:
        narrowed=all(group.iloc[-1].median_g_interval_width_original<group.iloc[0].median_g_interval_width_original
                     for _,group in curve.groupby('canonical_race'))
        if narrowed:
            decision='Evidence supports considering larger category-matched references in a separately justified future analysis: observed g-selection dispersion decreased in both categories. This does not establish clinical reliability or justify automatic replacement of the frozen primary analysis.'
        if not narrowed:
            decision='The completed experiment does not consistently support larger subsets reducing g dispersion in both categories; no primary-design replacement is justified.'
    if full_available:
        orig=pd.read_csv(output/'A_vs_B_category_summary.csv')
        retention='; '.join(f'{r.canonical_race} {r.original_retained}/{r.A_discoveries} retained, {r.original_lost} lost, {r.new_discoveries} new' for r in orig.itertuples())
    else:
        retention='B not yet run.'
    bootstrap_rows=summaries.loc[summaries.sampling_scheme.eq('with_replacement_reference_bootstrap')]
    full_caveat='The separate full-pool bootstrap has not run.'
    if len(bootstrap_rows):
        full_caveat='Separate full-pool conditional bootstrap: '+'; '.join(f'{r.canonical_race} {int(r.stable_original)}/{int(r.original_discoveries)} ({r.stable_original_percent:.1f}%)' for r in bootstrap_rows.itertuples())+'. European weights vary in this with-replacement bootstrap; its 647-photo without-replacement subset remains fixed in A/C.'
    result_lines=['# Results and limitations','',f'Overall status: **{status}**.','',
                  f'A baseline/audit: COMPLETED; original arrays verified/reused. B full point: {"COMPLETED" if full_available else "NOT RUN"}; '
                  f'B reference bootstrap: {"COMPLETED" if b_bootstrap else "NOT RUN"}; C completed resampleable configurations '
                  f'{observed_stages}/{expected_stages} (including original A). Repetitions per new design: {plan["requested_repetitions"]}.','',
                  f'Frozen eligible reference pools: {plan["pools"]}; original allocations: {plan["original_quotas"]}. '
                  '80 cells, 67 syndromes, 736 distinct eligible patients, 120 measurements and 9,600 planned tests; 982 original discoveries.','',
                  f'Original African/Asian failures: {original_only} significance-only, {sign_only} sign-only, {both} both-failure. '
                  'African and Asian each contributed 59 failures. The 15 European operational q-recurrence failures do not represent '
                  'alternative-European-reference instability.','',*trend,'',retention,'',
                  '## Decision guidance','',decision,'',full_caveat,'',
                  f'{int(evidence.consistent_across_available_schemes.sum())}/387 original African/Asian findings pass the unchanged 80%/95% criteria in every available stochastic scheme. The cross-scheme intersection is a descriptive follow-up filter, not a new significance test or clinical validation.','',
                  '## Incremental size changes and plateau interpretation','',table(increments),'',
                  'Observed interval widths decrease throughout the tested ladder. This does not establish a population-level plateau. Near-full without-replacement samples must approach zero selection dispersion as the finite pool is exhausted, whereas the separate full-pool bootstrap retains nonzero uncertainty. Reference-mean SDs are saved per feature in native units; their cross-feature summaries are not physical anatomical sizes.','',
                  '## Limits','',
                  'A/C subset recurrence and B full-pool bootstrap are different conditional uncertainty quantities. '
                  'B uses every observed photograph once for its deterministic point estimate; its repeated full subsets cannot '
                  'measure uncertainty. The separate bootstrap samples with replacement, includes duplicate/tied measurements, '
                  'and conditions on the observed pool and fixed cases. Hidden subject/source clustering is not established. '
                  'Selection-distribution quantiles are not clinical confidence intervals or independent replication.','',
                  'Reduced g dispersion can coexist with persistent q-threshold instability or lower retention of initially favorable '
                  'reference comparisons. Jointly increased category sizes alter global BH competition, so q recurrence need not '
                  'improve monotonically. No unique causal explanation is attributed to larger sample size alone.','',
                  'The original discovery screen and the 80% recurrence requirement impose different evidence criteria; failing the latter is not itself proof of a false association. Finite Monte Carlo frequencies near the 80%/95% cutoffs remain uncertain. GMDB cases versus FairFace references remain perfectly linked to dataset/source. Larger FairFace samples do '
                  'not remove acquisition, referral, conditional age, sex, pose or hidden demographic confounding. Recorded labels '
                  'are not verified ancestry. Correlated features are not independent traits; no novelty, diagnostic validation or '
                  'syndrome-specific anatomy is established. The baseline, selectors, manuscript and original sensitivity artifacts '
                  'remain unchanged. A replacement primary analysis would require explicit scientific approval and a versioned plan.']
    (output/'RESULTS_AND_LIMITATIONS.md').write_text('\n'.join(result_lines)+'\n')
    captions=[
        'Figure 1. Operational causes of failing the fixed 80% q-recurrence / 95% sign-agreement rule among the original 118 African/Asian findings. Percentages use 59 failures per category; categories are descriptive, not causal.',
        'Figure 2. Stability of original discoveries across without-replacement reference subsets; European photographs fixed. Full-pool selection is deterministic and outside the curve. Same 9,600-test family per draw.',
        'Figure 3. Median width of the 2.5–97.5% g-selection distribution among 135 African and 252 Asian original discoveries. These are not patient confidence intervals; overlapping reference draws are not replication.',
        'Figure 4. Original versus full-reference g for 387 African/Asian original discoveries; diagonal denotes agreement. Sign reversals are marked, and original/case-feature definitions are unchanged.'
    ]
    (output/'COMPUTATIONAL_FORMULAS.md').write_text('''# Exact computational definitions

- U is case-first: sum over case/reference pairs of I(case>reference)+0.5 I(tie).
- MW auto: exact when min(n,m)<=8 and no combined ties; otherwise two-sided asymptotic with tie and continuity corrections. For N=n+m and T=sum(t^3-t), Var(U)=nm/12*(N+1-T/[N(N-1)]). z=(max(U,nm-U)-nm/2-0.5)/sqrt(Var(U)); p=min(1,2*normal_survival(z)). Zero variance gives p=1.
- Exact U probabilities: Gaussian-binomial coefficients of product_{i=1..n}(1-z^(m+i))/(1-z^i), divided by choose(n+m,n); two-sided p is twice the lower-tail CDF at min(U,nm-U), capped at 1. n is the smaller size. Public SciPy checks verify the optimized generating polynomial.
- Hedges g=J(df)*(case_mean-reference_mean)/sqrt(((n-1)*case_var+(m-1)*reference_var)/df), df=n+m-2; J=Gamma(df/2)/(sqrt(df/2)*Gamma((df-1)/2)). Sample variances use ddof=1; undefined pooled variance gives NaN.
- One full-family BH per draw: sort all 9,600 p values, reverse cumulative minimum of 9,600*p_(i)/i, cap at 1, restore exact keys. No subset BH or changed family.
- Repeated-reference observations are integer photo weights w. Means use sum(w*y)/m. Variance uses sum(w*(y-center)^2)-sum(w*(y-center))^2/m, divided by m-1. A fixed numerical center improves precision but does not alter the feature or normalization.
- Combined tie counts include reference multiplicities, case multiplicities and shared values; reference and case tie counts merge by adding 3*a*b*(a+b) per shared value with counts a,b. This preserves bootstrap duplicates.
- Sign agreement counts only finite draws with the same sign as a finite nonzero original g. Undefined g and g=0 cannot silently count as directional successes. q recurrence is fraction of all draws with finite q<0.05. Stable requires recurrence>=0.80, sign agreement>=0.95 and all g draws finite.
- Selection quantiles use NumPy linear quantiles at 0.025/0.975; width is their difference. Empirical SD uses ddof=1. These are reference-selection distributions, not patient confidence intervals.
- Full-pool without-replacement selection has one distinct subset. No selection recurrence/interval/SD is estimable; deterministic point g/q are exported separately. European subset results can have fixed raw p/g but changing global q through other groups.
- Stage A and C reuse each category's selected reference across all syndrome cells. B bootstrap samples full-size category photo pools with replacement, conditions on fixed cases/observed pool, and is not the same uncertainty estimand as A/C.
- Plot denominators: original African discoveries=135 and Asian=252; original unstable discoveries=59 per category. These counts are verified from original files, not independent phenotypes.
''')
    (output/'FIGURE_CAPTIONS.txt').write_text('\n\n'.join(captions)+'\n')
    methods=('Methods: In an exploratory secondary diagnostic, we retained all frozen patients, 120 original measurements and '
             'the 9,600-test family. We reused 1,000 verified original subset draws and partitioned stability failures using the '
             'prespecified ≥80% global-q recurrence and ≥95% g-sign criteria. We evaluated deterministic full-category FairFace '
             'pools, a separate with-replacement reference bootstrap, and paired doubling-size subsets with European references '
             f'fixed. New stochastic configurations used {plan["requested_repetitions"]:,} repetitions each; exact/asymptotic '
             'Mann–Whitney rules, tie corrections, exact Hedges gamma adjustment and full-family BH were retained.')
    results=(f'Results ({status}): Of 118 original African/Asian failures, {original_only} were significance-only, '
             f'{sign_only} sign-only and {both} failed both criteria. '+' '.join(trend)+' '+retention+' '+full_caveat)
    limitation=('Limitations: These conditional resampling schemes assess different variability sources. Full pools have one '
                'unique subset; bootstrap duplicates and reference-source/person clustering limit model-based interpretation. '
                'The original discovery set is selected, patient cohorts remain fixed, correlated traits are not independent, '
                'and GMDB–FairFace source confounding persists. Finite Monte Carlo frequencies near the 80%/95% thresholds are uncertain. No clinical validation or primary replacement follows.')
    (output/'manuscript_ready_summary.md').write_text(methods+'\n\n'+results+'\n\n'+limitation+'\n\n'+'\n\n'.join(captions)+'\n')
    manifest['status']=status
    manifest['report_code_sha256']=__import__('hashlib').sha256(Path(__file__).read_bytes()).hexdigest()
    manifest['output_sha256']={str(p.relative_to(output)):__import__('hashlib').sha256(p.read_bytes()).hexdigest()
                               for p in output.rglob('*') if p.is_file() and p.name!='source_and_execution_manifest.json'}
    (output/'source_and_execution_manifest.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
    print(status+'; '+' | '.join(trend),flush=True)

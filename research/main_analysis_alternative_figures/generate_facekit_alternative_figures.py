import argparse
import importlib.metadata
import json
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.special import gammaln

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from revision_results.race_stratified_20261115 import common
from research.main_analysis_panels.generate_facekit_panels import select, readable, COLORS, RACES, CRITERION

FEATURES = ['eye_fissure_height_mean', 'nose_tip_width', 'philtrum_length', 'chin_width']
FEATURE_LABELS = ['Eye fissure height', 'Nasal tip width', 'Philtrum length', 'Chin width']
REGION_MAP = {'forehead':'Forehead', 'eyebrows':'Eyebrows', 'eyes / orbital dimensions':'Eyes / orbital',
              'ocular position / iris / gaze':'Eyes / orbital', 'nose':'Nose', 'philtrum / mouth / lips':'Mouth / lips',
              'chin / jaw':'Chin / jaw', 'cheeks / malar region':'Face shape / cheeks',
              'global face shape / proportions':'Face shape / cheeks'}
REGIONS = list(dict.fromkeys(REGION_MAP.values()))
SYNDROME_COLORS = ['#332288','#88CCEE','#117733','#CC6677','#AA4499']
FOREST_CAPTION = (
    'Selected syndrome–feature comparisons for five syndromes chosen by eligible patient count among syndromes with '
    'at least two eligible recorded categories. The four features were fixed in advance. Points denote original Hedges’ g '
    'relative to same-category FairFace; intervals are conditional case-only 95% bootstrap percentile intervals '
    '(2,000 replications by default), holding the exact original references fixed. Filled markers denote original '
    'global BH q < 0.05 across 9,600 tests; hollow markers are nonsignificant. Interval exclusion of zero and original '
    'BH significance are different procedures and may disagree. Intervals omit reference sampling uncertainty, discovery '
    'selection, multiplicity, acquisition confounding and potential dependencies. They do not establish independent '
    'clinical validation. Labels are recorded demographic categories, not genetically verified ancestry. Missing '
    'categories are ineligible, not zero. Measurements are normalized 2D geometry, not absolute clinical dimensions.'
)
REGION_CAPTION = (
    'Percentage of significant original primary comparisons among supported FaceKit tests in seven fixed anatomical '
    'groups for five syndromes selected by eligible patient counts and their eligible recorded demographic categories. '
    'Significance is the original single global 9,600-test BH q < 0.05; there is no region-wise correction. Labels show '
    'the number of prespecified features in each group; the exported table records supported denominators for every bar. '
    'Correlated measurements are not independent phenotypic discoveries. Missing bars indicate ineligible comparisons; '
    'zero-height bars denote eligible groups with no discoveries. Colors identify syndromes; category labels are recorded '
    'annotations, not genetically verified ancestry. Sample size and GMDB/FairFace acquisition differences limit '
    'comparisons; discovery percentages do not measure diagnostic usefulness or establish causal category differences.'
)


def effect(x, y):
    x, y = np.asarray(x), np.asarray(y)
    df = len(x)+len(y)-2
    pooled = ((len(x)-1)*x.var(axis=0,ddof=1)+(len(y)-1)*y.var(axis=0,ddof=1))/df
    correction = np.exp(gammaln(df/2)-.5*np.log(df/2)-gammaln((df-1)/2))
    return correction*(x.mean(axis=0)-y.mean(axis=0))/np.sqrt(pooled)


def bootstrap(x, y, repetitions, seed):
    if repetitions<2 or min(len(x),len(y))<2 or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Finite real measurements and at least two observations/resamples required')
    rng = np.random.default_rng(seed)
    draws = x[rng.integers(0,len(x),size=(repetitions,len(x)))]
    df = len(x)+len(y)-2
    correction = np.exp(gammaln(df/2)-.5*np.log(df/2)-gammaln((df-1)/2))
    pooled = ((len(x)-1)*draws.var(axis=1,ddof=1)+(len(y)-1)*y.var(axis=0,ddof=1))/df
    if (pooled<=0).any():
        raise ValueError('Undefined bootstrap pooled variance')
    gs = correction*(draws.mean(axis=1)-y.mean(axis=0))/np.sqrt(pooled)
    if not np.isfinite(gs).all():
        raise ValueError('Undefined bootstrap effect; do not invent an interval')
    return np.quantile(gs,[.025,.975],axis=0)


def load():
    directory = common.PRIVATE
    paths = [directory/'statistics_manifest.json',directory/'primary_all_tests.csv',
             ROOT/'publication/feature_definitions/current_120_feature_formulas.csv']
    if any(not p.is_file() for p in paths):
        raise FileNotFoundError('BLOCKED: locked primary inputs unavailable')
    manifest = json.loads(paths[0].read_text())
    if manifest['study_id']!='PRISM_RACE_STRATIFIED_20261115' or common.sha(paths[1])!=manifest['output_sha256'][paths[1].name]:
        raise ValueError('BLOCKED: primary study or SHA256 mismatch')
    frame = pd.read_csv(paths[1],float_precision='round_trip')
    keys = ['disease','canonical_race','feature']
    definitions = pd.read_csv(paths[2])
    if len(frame)!=9600 or frame[keys].duplicated().any() or len(definitions)!=120:
        raise ValueError('BLOCKED: complete 9,600-test primary table required')
    if not frame.supported.isin([True,False]).all() or not frame.canonical_race.isin(RACES).all():
        raise ValueError('BLOCKED: invalid support/category flags')
    if frame.q.isna().any() or not frame.q.between(0,1).all():
        raise ValueError('BLOCKED: invalid original q values')
    schema = definitions.set_index('feature')
    if set(frame.feature)!=set(schema.index) or not set(FEATURES).issubset(schema.index):
        raise ValueError('BLOCKED: 120-feature schema mismatch')
    for col in ['formula','anatomical_region']:
        if not frame[col].eq(frame.feature.map(schema[col])).all():
            raise ValueError('BLOCKED: primary feature metadata mismatch')
    if set(frame.anatomical_region)!=set(REGION_MAP) or len(set(REGION_MAP.values()))!=7:
        raise ValueError('BLOCKED: nine source regions must partition into seven fixed display groups')
    frame['display_region'] = frame.anatomical_region.map(REGION_MAP)
    frame['original_significant'] = frame.supported & frame.q.lt(.05)
    cells = frame.groupby(['disease','canonical_race'],as_index=False).agg(
        case_n=('case_n','first'), n_variants=('case_n','nunique'), reference_n=('reference_n','first'),
        reference_variants=('reference_n','nunique'), tests=('feature','size'), distinct=('feature','nunique'))
    if (len(cells)!=80 or cells.disease.nunique()!=67 or not cells.tests.eq(120).all()
            or not cells.distinct.eq(120).all() or not cells.n_variants.eq(1).all()
            or not cells.reference_variants.eq(1).all() or cells.case_n.sum()!=736
            or not cells.case_n.ge(5).all() or frame.original_significant.sum()!=982):
        raise ValueError('BLOCKED: frozen cohort/family dimensions do not match')
    if not cells.reference_n.eq(cells.canonical_race.map(manifest['completed']['quotas'])).all():
        raise ValueError('BLOCKED: exact reference allocation differs from manifest')
    selected = select(cells)
    denominators = schema.reset_index().assign(display_region=lambda x:x.anatomical_region.map(REGION_MAP)).groupby('display_region').size().to_dict()
    return frame,cells,selected,denominators,manifest,paths


def save(fig,output,stem):
    fig.canvas.draw()
    from matplotlib.text import Text
    renderer = fig.canvas.get_renderer()
    for text in fig.findobj(Text):
        if text.get_visible() and text.get_text():
            box = text.get_window_extent(renderer)
            if box.x0<-.5 or box.y0<-.5 or box.x1>fig.bbox.x1+.5 or box.y1>fig.bbox.y1+.5:
                raise ValueError('Figure text clipped: '+text.get_text())
    for extension in ['png','pdf','svg']:
        fig.savefig(output/(stem+'.'+extension),dpi=300,facecolor='white')
    plt.close(fig)


def forest(primary,cells,selected,manifest,output,args):
    disease,all_reference,features,definitions,audit,sources = common.load_sources()
    for path in sources:
        name = str(path.relative_to(ROOT))
        if common.sha(path)!=manifest['source_sha256'][name]:
            raise ValueError('BLOCKED: baseline source provenance mismatch: '+name)
    selection_path = common.PRIVATE/'PRIVATE_reference_selection.csv'
    if common.sha(selection_path)!=manifest['output_sha256'][selection_path.name]:
        raise ValueError('BLOCKED: reference selection SHA256 mismatch')
    sampling = {'quotas':manifest['completed']['quotas'], 'seed':manifest['parameters']['seed']}
    if not cells.reference_n.eq(cells.canonical_race.map(sampling['quotas'])).all():
        raise ValueError('BLOCKED: reference counts differ from manifest allocation')
    reference = common.sample_reference(all_reference,sampling['quotas'],sampling['seed'])
    stored = pd.read_csv(selection_path,dtype={'image_id':'string'})
    for col in ['image_id','file','image_sha256','canonical_race']:
        if reference[col].astype(str).tolist()!=stored[col].astype(str).tolist():
            raise ValueError('BLOCKED: original FairFace sample cannot be reproduced')
    rows = []
    for i,entry in enumerate(selected.itertuples()):
        for j,race in enumerate(RACES):
            cell = cells.loc[cells.disease.eq(entry.disease)&cells.canonical_race.eq(race)]
            if not len(cell):
                continue
            case = disease.loc[disease.disease.eq(entry.disease)&disease.canonical_race.eq(race)].sort_values('patient_id')
            ref = reference.loc[reference.canonical_race.eq(race)]
            if not case.patient_id.is_unique or len(case)!=cell.case_n.iloc[0] or len(ref)!=cell.reference_n.iloc[0]:
                raise ValueError('BLOCKED: raw cell sizes/unique patient records do not match')
            x,y = case[FEATURES].to_numpy(float),ref[FEATURES].to_numpy(float)
            estimates = effect(x,y)
            low,high = bootstrap(x,y,args.bootstrap_reps,np.random.SeedSequence([args.seed,i,j]))
            for k,feature in enumerate(FEATURES):
                frozen = primary.loc[primary.disease.eq(entry.disease)&primary.canonical_race.eq(race)&primary.feature.eq(feature)].iloc[0]
                np.testing.assert_allclose([estimates[k],x[:,k].mean(),y[:,k].mean()],
                                           [frozen.hedges_g,frozen.case_mean,frozen.reference_mean],rtol=1e-11,atol=1e-13)
                rows.append({'disease':entry.disease,'canonical_race':race,'feature':feature,'case_n':len(x),'reference_n':len(y),
                             'g':frozen.hedges_g,'case_mean':frozen.case_mean,'reference_mean':frozen.reference_mean,
                             'original_global_q':frozen.q,'original_significant':bool(frozen.original_significant),
                             'bootstrap_lower':low[k],'bootstrap_upper':high[k],'bootstrap_reps':args.bootstrap_reps,
                             'interval_excludes_zero':bool(low[k]>0 or high[k]<0),'interval_method':'conditional case-only percentile; fixed original reference'})
    result = pd.DataFrame(rows)
    if len(result)!=4*cells.disease.isin(selected.disease).sum():
        raise ValueError('BLOCKED: missing forest comparisons')
    result.to_csv(output/'forest_plot_data.csv',index=False)
    fig,axes = plt.subplots(3,2,figsize=(10.6,9.4))
    fig.subplots_adjust(left=.14,right=.97,bottom=.07,top=.91,wspace=.66,hspace=.74)
    bounds = [min(-2.5,result.bootstrap_lower.min()-.15),max(2.5,result.bootstrap_upper.max()+.15)]
    for i,entry in enumerate(selected.itertuples()):
        ax = axes.flat[i]
        group = result.loc[result.disease.eq(entry.disease)]
        for j,race in enumerate(RACES):
            for k,feature in enumerate(FEATURES):
                part = group.loc[group.canonical_race.eq(race)&group.feature.eq(feature)]
                if not len(part):
                    continue
                row = part.iloc[0]
                yy = k+(j-1)*.20
                ax.hlines(yy,row.bootstrap_lower,row.bootstrap_upper,color=COLORS[race],linewidth=1.1)
                ax.plot(row.g,yy,marker='o',markersize=4.5,markeredgewidth=1,
                        markerfacecolor=COLORS[race] if row.original_significant else 'white',markeredgecolor=COLORS[race])
        counts = group.groupby('canonical_race').case_n.first()
        header = readable(entry.disease)+f'\nEligible total n = {entry.total_eligible_patients}\n'
        header += ' | '.join(f'{r} n={int(counts[r])}' for r in RACES if r in counts)
        ax.set_title(header,fontsize=9,pad=10,linespacing=1.2)
        ax.set_yticks(range(4),FEATURE_LABELS,fontsize=8)
        ax.set_ylim(3.5,-.5)
        ax.set_xlim(bounds)
        ax.axvline(0,linestyle='--',color='#777777',linewidth=.7)
        ax.tick_params(axis='x',labelsize=8)
        ax.spines[['top','right']].set_visible(False)
        ax.set_xlabel('Hedges’ g',fontsize=9)
    key = axes.flat[5]
    key.axis('off')
    handles = [Line2D([],[],color=COLORS[r],marker='o',label=r) for r in RACES]
    handles += [Line2D([],[],color='#333333',marker='o',linestyle='',label='Filled: original global q < 0.05'),
                Line2D([],[],color='#333333',marker='o',markerfacecolor='white',linestyle='',label='Hollow: original global q ≥ 0.05')]
    key.legend(handles=handles,loc='upper left',bbox_to_anchor=(0,1.04),frameon=False,fontsize=8,
               title='Recorded categories / original significance',title_fontsize=8,labelspacing=.35,borderaxespad=0)
    key.text(0,.43,'95% conditional case-only bootstrap intervals\nExact original FairFace references fixed\nNot adjusted for selection or multiplicity\nNo independent clinical validation\nCategory labels are not verified ancestry',transform=key.transAxes,fontsize=8,linespacing=1.3,va='top')
    fig.suptitle('Prespecified syndrome–feature effects by recorded category',fontsize=13,y=.985)
    save(fig,output,'alternative_forest_plot')
    (output/'forest_caption.txt').write_text(FOREST_CAPTION.replace('2,000 replications by default',f'{args.bootstrap_reps:,} replications')+'\n')
    return sources+[selection_path],sampling


def regions(primary,cells,selected,denominators,output):
    rows = []
    for entry in selected.itertuples():
        for race in RACES:
            group = primary.loc[primary.disease.eq(entry.disease)&primary.canonical_race.eq(race)]
            for region in REGIONS:
                part = group.loc[group.display_region.eq(region)]
                supported = int(part.supported.sum()) if len(part) else None
                significant = int(part.original_significant.sum()) if len(part) else None
                rows.append({'disease':entry.disease,'canonical_race':race,'display_region':region,
                             'status':'eligible' if len(group) else 'ineligible','case_n':int(group.case_n.iloc[0]) if len(group) else None,
                             'planned_features':denominators[region],'supported_tests':supported,'significant_tests':significant,
                             'percent_significant':100*significant/supported if supported else None})
    result = pd.DataFrame(rows)
    for (d,r),group in result.loc[result.status.eq('eligible')].groupby(['disease','canonical_race']):
        source = primary.loc[primary.disease.eq(d)&primary.canonical_race.eq(r)]
        if group.planned_features.sum()!=120 or group.significant_tests.sum()!=source.original_significant.sum() or group.supported_tests.sum()!=source.supported.sum():
            raise ValueError('BLOCKED: regional partition/count mismatch')
    result.to_csv(output/'anatomical_region_chart_data.csv',index=False)
    fig,axes = plt.subplots(3,1,figsize=(11,9),sharex=True,sharey=True)
    fig.subplots_adjust(left=.075,right=.75,bottom=.15,top=.92,hspace=.34)
    positions = np.arange(7)
    for ax,race in zip(axes,RACES):
        for i,entry in enumerate(selected.itertuples()):
            group = result.loc[result.disease.eq(entry.disease)&result.canonical_race.eq(race)].set_index('display_region').reindex(REGIONS)
            if group.status.eq('ineligible').all():
                continue
            ax.bar(positions+(i-2)*.145,group.percent_significant,width=.135,color=SYNDROME_COLORS[i],zorder=3)
        ax.set_ylim(0,105)
        ax.set_yticks([0,25,50,75,100])
        ax.grid(axis='y',alpha=.25,linewidth=.6,zorder=0)
        ax.spines[['top','right']].set_visible(False)
        ax.set_title(race+' recorded category',loc='left',fontsize=11,pad=6)
        ax.set_ylabel('Significant tests (%)',fontsize=10)
        for i,entry in enumerate(selected.itertuples()):
            part = result.loc[result.disease.eq(entry.disease)&result.canonical_race.eq(race)]
            eligible = part.status.eq('eligible').any()
            name = readable(entry.disease).replace('\n',' ')
            name = {'Hyperphosphatasia–intellectual disability syndrome':'Hyperphosphatasia syndrome'}.get(name,name)
            ax.plot(1.025,.92-i*.21,marker='s',markersize=5,color=SYNDROME_COLORS[i],transform=ax.transAxes,clip_on=False)
            ax.text(1.045,.95-i*.21,name+('\nn = '+str(int(part.case_n.iloc[0])) if eligible else '\nIneligible — absent bars'),
                    transform=ax.transAxes,fontsize=8.4,color='#222222',va='top',linespacing=1.05)
    axes[-1].set_xticks(positions,[r.replace(' / ',' /\n')+f'\n{denominators[r]} features' for r in REGIONS],fontsize=8.5)
    fig.suptitle('Original global-BH discovery rates across seven anatomical groups',fontsize=13,y=.98)
    fig.text(.075,.025,'Five syndromes selected by eligible patient counts and category coverage.\n'
             'Original 9,600-test BH family; correlated measurements are not independent discoveries.',fontsize=9)
    save(fig,output,'alternative_anatomical_region_chart')
    (output/'regions_caption.txt').write_text(REGION_CAPTION+'\n')


def main(args):
    if args.bootstrap_reps<2:
        raise ValueError('At least two bootstrap repetitions required')
    output = common.private_output(args.output_dir)
    protected = common.protected_hashes()
    primary,cells,selected,denominators,manifest,inputs = load()
    selected.to_csv(output/'five_syndrome_selection.csv',index=False)
    font = 'Times New Roman' if 'Times New Roman' in {f.name for f in font_manager.fontManager.ttflist} else 'DejaVu Serif'
    plt.rcParams.update({'font.family':font,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
    sampling = None
    if args.figure in ['forest','both']:
        extra,sampling = forest(primary,cells,selected,manifest,output,args)
        inputs += extra
    if args.figure in ['regions','both']:
        regions(primary,cells,selected,denominators,output)
    if protected!=common.protected_hashes():
        raise ValueError('Protected artifacts changed')
    record = {'study_id':manifest['study_id'],'selection_rule':CRITERION,'selected_syndromes':selected.to_dict('records'),
              'figure_requested':args.figure,'fixed_forest_features':FEATURES,'region_map':REGION_MAP,'region_denominators':denominators,
              'original_primary_tests':9600,'original_discoveries':982,'eligible_cells':80,
              'primary_tests_or_q_recomputed':False,'bootstrap_method':'conditional case-only percentile; exact original reference fixed',
              'seed':args.seed,'bootstrap_reps':args.bootstrap_reps,'sampling_replay':sampling,'font':font,'dpi':300,
              'source_hashes':{str(p):common.sha(p) for p in set(inputs)},'code_sha256':common.sha(Path(__file__)),
              'software':{p:importlib.metadata.version(p) for p in ['numpy','pandas','scipy','matplotlib']},
              'deviations':['Companion script unavailable; implemented specification directly; seven-region mapping explicitly frozen in code'],
              'protected_files_unchanged':True,'patient_rows_exported':False,
              'caveats':[FOREST_CAPTION,REGION_CAPTION],
              'output_sha256':{p.name:common.sha(p) for p in output.iterdir() if p.is_file() and p.name!='alternative_figures_manifest.json'}}
    common.write_json(output/'alternative_figures_manifest.json',record)
    print(selected[['disease','total_eligible_patients','eligible_categories']].to_string(index=False))


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--figure',choices=['forest','regions','both'],default='both')
    parser.add_argument('--bootstrap-reps',type=int,default=2000)
    parser.add_argument('--seed',type=int,default=20261115)
    parser.add_argument('--output-dir',type=Path,default=ROOT/'local_review/main_analysis_alternative_figures')
    args = parser.parse_args()
    main(args)

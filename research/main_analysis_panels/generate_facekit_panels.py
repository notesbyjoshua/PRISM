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
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from revision_results.race_stratified_20261115 import common
from research.race_stratified_sensitivity.core import bh

FEATURES = ['forehead_width_mid', 'eb_inter_distance', 'eb_thickness_mean', 'eye_fissure_length_mean',
            'eye_fissure_height_mean', 'inter_canthal_distance', 'inter_pupillary_distance', 'nose_bridge_width',
            'nose_tip_width', 'nose_length', 'philtrum_length', 'upper_vermilion_height', 'mouth_width', 'chin_width', 'chin_height']
LABELS = ['Mid-forehead width', 'Medial brow spacing', 'Eyebrow polygon area', 'Eye fissure length',
          'Eye fissure height', 'Inner-canthal distance', 'Interpupillary distance', 'Nasal bridge width',
          'Nasal tip width', 'Nose length', 'Philtrum length', 'Upper vermilion height', 'Mouth width', 'Chin width', 'Chin height']
RACES = ['African', 'Asian', 'European']
COLORS = {'African': '#0072B2', 'Asian': '#D55E00', 'European': '#009E73'}
MARKERS = ['o', 's', '^', 'D', 'P']
CRITERION = 'At least two eligible recorded categories; total eligible patients descending, eligible categories descending, exact disease key ascending'
CAPTION_A = (
    'Panel A. All 80 eligible syndrome–recorded-category cells from 67 syndromes in PRISM_RACE_STRATIFIED_20261115. '
    'Each point shows its unique GMDB patient count and the number of significant features among 120 normalized facial measurements. '
    'Significance uses the original single global Benjamini–Hochberg correction of all 9,600 two-sided Mann–Whitney U tests '
    '(q < 0.05; 982 significant associations overall). Colors identify recorded categories. Outlined shapes identify five '
    'syndromes selected using at least two eligible categories and the largest total eligible patient counts, before consulting '
    'discovery counts; the side key reports totals across eligible cells only. All unselected cells remain visible; coincident '
    'points may overlap. Discovery counts depend on sample size, image-source differences and correlated measurements. '
    'They measure neither clinical diagnostic utility nor counts of independent anatomical traits. Dataset categories are '
    'not genetically verified ancestry; no formal cross-category interaction or independent clinical validation is shown.'
)
CAPTION_B = (
    'Panel B. Fifteen prespecified normalized FaceKit measurements for every eligible recorded-category cell of the same '
    'five syndromes selected for Panel A. Fill shows signed Hedges’ g (GMDB syndrome minus same-category FairFace reference), '
    'with a common symmetric −2 to +2 scale; effects outside that range saturate visually without changing stored values. '
    'Dots indicate original global q < 0.05 across all 9,600 primary tests; nonsignificant cells remain visible. Gray indicates '
    'unsupported or undefined effects and is never imputed as zero. Group lines separate syndromes, and n denotes patients. '
    'The 15 displayed measurements are a fixed subset of 120 and were not chosen by q-values. Distances, separations and polygon '
    'areas are normalized 2D image measurements, not physical millimeters, hair density or clinical diagnoses. Correlated metrics, '
    'sample sizes and GMDB/FairFace acquisition differences limit interpretation. Categories are recorded dataset labels, '
    'not genetically verified ancestry; the figure provides no formal interaction test or independent clinical validation.'
)


def readable(key):
    name = key.split('__')[0].split(';')[0]
    if name == 'WILLIAMS-BEUREN SYNDROME':
        return 'Williams–Beuren syndrome'
    if name == 'Hyperphosphatasia with mental retardation syndrome':
        return 'Hyperphosphatasia–intellectual\ndisability syndrome'
    return name


def load(directory):
    path = directory / 'primary_all_tests.csv'
    manifest_path = directory / 'statistics_manifest.json'
    manifest = json.loads(manifest_path.read_text())
    if manifest.get('study_id') != 'PRISM_RACE_STRATIFIED_20261115':
        raise ValueError('Expected the locked November 15 race-stratified study; no historical fallback')
    if common.sha(path) != manifest['output_sha256'][path.name]:
        raise ValueError('Primary input SHA256 does not match its verification manifest')
    frame = pd.read_csv(path, float_precision='round_trip')
    if len(frame) != 9600 or frame[['disease','canonical_race','feature']].duplicated().any():
        raise ValueError('Expected 9,600 unique primary tests')
    schema_path = ROOT / 'publication/feature_definitions/current_120_feature_formulas.csv'
    definitions = pd.read_csv(schema_path).set_index('feature')
    if len(definitions) != 120 or set(definitions.index) != set(frame.feature) or not set(FEATURES).issubset(definitions.index):
        raise ValueError('Current 120-feature definitions do not match the primary table')
    for column in ['formula', 'anatomical_region']:
        if not frame[column].eq(frame.feature.map(definitions[column])).all():
            raise ValueError('Primary feature metadata differ: ' + column)
    if frame[['p','q']].isna().any().any() or not frame.canonical_race.isin(RACES).all():
        raise ValueError('Invalid p/q or recorded categories')
    np.testing.assert_allclose(bh(frame.p), frame.q, rtol=1e-12, atol=1e-15)
    cells = frame.groupby(['disease','canonical_race'], sort=True).agg(
        case_n=('case_n','first'), case_n_variants=('case_n','nunique'), planned_tests=('feature','size'), distinct_features=('feature','nunique')).reset_index()
    if (len(cells) != 80 or cells.disease.nunique() != 67 or not cells.planned_tests.eq(120).all()
            or not cells.distinct_features.eq(120).all() or not cells.case_n_variants.eq(1).all()):
        raise ValueError('Expected 80 cells, 67 syndromes and 120 measurements per cell')
    if not cells.case_n.ge(5).all() or cells.case_n.sum() != 736:
        raise ValueError('Unexpected eligible patient counts')
    selected = select(cells)
    discovery_counts = frame.loc[frame.supported & frame.q.lt(.05)].groupby(['disease','canonical_race']).size()
    cells['significant_features'] = [int(discovery_counts.get((r.disease,r.canonical_race),0)) for r in cells.itertuples()]
    independent = frame.assign(discovery=frame.supported & frame.q.lt(.05)).groupby(['disease','canonical_race']).discovery.sum()
    np.testing.assert_array_equal(cells.significant_features, independent.reindex(pd.MultiIndex.from_frame(cells[['disease','canonical_race']])).to_numpy())
    if cells.significant_features.sum() != 982:
        raise ValueError('Expected 982 original discoveries')
    return frame, cells, selected, [path, manifest_path, schema_path]


def select(cells):
    eligible = cells.groupby('disease', as_index=False).agg(
        total_eligible_patients=('case_n','sum'), eligible_category_count=('canonical_race','size'),
        eligible_categories=('canonical_race', lambda values: '|'.join(r for r in RACES if r in set(values))))
    eligible = eligible.loc[eligible.eligible_category_count.ge(2)].sort_values(
        ['total_eligible_patients','eligible_category_count','disease'], ascending=[False,False,True], kind='stable')
    if len(eligible) < 5:
        raise ValueError('Fewer than five syndromes have at least two eligible categories')
    selected = eligible.head(5).copy()
    selected.insert(0, 'selection_rank', range(1,6))
    selected['criterion'] = CRITERION
    return selected


def save(fig, output, name):
    fig.canvas.draw()
    from matplotlib.text import Text
    renderer = fig.canvas.get_renderer()
    for text in fig.findobj(Text):
        if not text.get_visible() or not text.get_text():
            continue
        box = text.get_window_extent(renderer)
        if box.x0 < -.5 or box.y0 < -.5 or box.x1 > fig.bbox.x1+.5 or box.y1 > fig.bbox.y1+.5:
            raise ValueError('Clipped figure text: ' + text.get_text())
    for extension in ['png','pdf','svg']:
        fig.savefig(output / (name+'.'+extension), dpi=300, facecolor='white')
    plt.close(fig)


def panel_a(cells, selected, output):
    fig, ax = plt.subplots(figsize=(10,5.4))
    fig.subplots_adjust(left=.085, right=.65, bottom=.14, top=.88)
    chosen = set(selected.disease)
    cells = cells.copy()
    cells['selected_syndrome'] = cells.disease.isin(chosen)
    for race in RACES:
        group = cells.loc[cells.canonical_race.eq(race) & ~cells.selected_syndrome]
        ax.scatter(group.case_n,group.significant_features,s=38,c=COLORS[race],alpha=.65,edgecolors='none',zorder=2)
    for i, row in enumerate(selected.itertuples()):
        group = cells.loc[cells.disease.eq(row.disease)]
        ax.scatter(group.case_n,group.significant_features,s=85,marker=MARKERS[i],c=group.canonical_race.map(COLORS),
                   edgecolors='#222222',linewidths=1.15,zorder=3)
    ax.set_xlim(0, max(35,cells.case_n.max()+2))
    ax.set_ylim(-2,70)
    ax.set_xticks(np.arange(0,36,5))
    ax.set_yticks(np.arange(0,71,10))
    ax.set_xlabel('Unique eligible GMDB patients in cell',labelpad=8)
    ax.set_ylabel('Significant facial features (out of 120)',labelpad=8)
    ax.set_title('A   Complete study overview: 80 eligible cells',loc='left',fontsize=12,pad=14)
    ax.spines[['top','right']].set_visible(False)
    ax.grid(alpha=.2,linewidth=.5,zorder=0)
    handles = [Line2D([],[],marker='o',linestyle='',color=COLORS[r],label=r) for r in RACES]
    fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.69,.9),frameon=False,title='Recorded category',fontsize=10)
    fig.text(.70,.61,'Outlined: five selected syndromes',fontsize=10,weight='bold')
    for i,row in enumerate(selected.itertuples()):
        y = .55-i*.102
        fig.add_artist(Line2D([.713],[y+.007],transform=fig.transFigure,marker=MARKERS[i],markersize=8,
                              markerfacecolor='white',markeredgecolor='#222222',linestyle=''))
        fig.text(.735,y,readable(row.disease)+f'\nEligible total n = {row.total_eligible_patients}',fontsize=9,va='center',linespacing=1.1)
    fig.text(.70,.025,'Selection uses patient counts and\ncategory coverage, not discovery counts.',fontsize=8.5)
    cells.to_csv(output/'PRIVATE_panel_A_all_cells.csv',index=False)
    save(fig,output,'panel_A_all_cells_scatter')


def panel_b(primary, cells, selected, output):
    rows = []
    for syndrome in selected.disease:
        for race in RACES:
            cell = cells.loc[cells.disease.eq(syndrome)&cells.canonical_race.eq(race)]
            if len(cell):
                rows.append((syndrome,race,int(cell.case_n.iloc[0])))
    desired = pd.DataFrame([{'disease':d,'canonical_race':r,'feature':f,'display_row':i,'display_column':j}
                            for i,(d,r,n) in enumerate(rows) for j,f in enumerate(FEATURES)])
    values = desired.merge(primary,on=['disease','canonical_race','feature'],how='left',validate='one_to_one',indicator=True)
    if not values['_merge'].eq('both').all():
        raise ValueError('Missing exact Panel B key')
    values = values.drop(columns='_merge')
    values['display_supported'] = values.supported & np.isfinite(values.hedges_g)
    values['significance_dot'] = values.supported & values.q.lt(.05)
    if (values.significance_dot & ~values.display_supported).any():
        raise ValueError('Significant row has undefined effect; requires explicit display review')
    values['display_g'] = values.hedges_g.where(values.display_supported)
    values['color_scale_saturated'] = values.display_supported & values.hedges_g.abs().gt(2)
    matrix = values.pivot(index='display_row',columns='display_column',values='display_g').reindex(index=range(len(rows)),columns=range(15)).to_numpy()
    fig,ax = plt.subplots(figsize=(11.6,7.5))
    fig.subplots_adjust(left=.325,right=.86,bottom=.28,top=.90)
    cmap = plt.get_cmap('RdBu_r').copy()
    cmap.set_bad('#D3D3D3')
    im = ax.pcolormesh(np.arange(16)-.5,np.arange(len(rows)+1)-.5,np.ma.masked_invalid(matrix),
                       cmap=cmap,norm=Normalize(-2,2),shading='flat',edgecolors='white',linewidth=.35)
    ax.set_ylim(len(rows)-.5,-.5)
    dots = values.loc[values.significance_dot]
    ax.scatter(dots.display_column,dots.display_row,s=15,c='black',edgecolors='white',linewidths=.45)
    ax.set_xticks(range(15),LABELS,rotation=58,ha='right',fontsize=9)
    ax.set_yticks(range(len(rows)),[f'{r} · n = {n}' for d,r,n in rows],fontsize=9)
    ax.tick_params(length=0,pad=6)
    ax.set_title('B   Five syndromes × 15 prespecified measurements',fontsize=12,pad=15,loc='left')
    start = 0
    for syndrome in selected.disease:
        count = sum(d==syndrome for d,r,n in rows)
        ax.text(-.19, (start+(count-1)/2), readable(syndrome), transform=ax.get_yaxis_transform(),
                ha='right',va='center',fontsize=10)
        if start:
            ax.axhline(start-.5,color='#222222',linewidth=1.3)
        start += count
    ax.set_xticks(np.arange(-.5,15,1),minor=True)
    ax.set_yticks(np.arange(-.5,len(rows),1),minor=True)
    ax.grid(which='minor',color='white',linewidth=.35)
    ax.tick_params(which='minor',length=0)
    bar_ax = fig.add_axes([.89,.37,.018,.40])
    bar = fig.colorbar(im,cax=bar_ax,extend='both',ticks=[-2,-1,0,1,2])
    bar.solids.set_rasterized(False)
    bar.set_label('Hedges’ g (GMDB − FairFace)',fontsize=10,labelpad=9)
    fig.text(.89,.28,'● Original\nglobal q < 0.05',fontsize=9)
    fig.text(.89,.19,'Gray: unsupported\nor undefined',fontsize=8.5)
    fig.text(.325,.025,'Normalized 2D metrics; category labels are recorded, not verified ancestry.\n'
             'Color saturation beyond ±2 preserves the original numerical effects.',fontsize=9)
    values.to_csv(output/'PRIVATE_panel_B_heatmap_values.csv',index=False)
    save(fig,output,'panel_B_five_syndromes_heatmap')
    return values


def main(args):
    output = common.private_output(args.output_dir)
    old_path = output/'run_manifest.json'
    old = json.loads(old_path.read_text()) if old_path.exists() else {}
    protected = common.protected_hashes()
    primary,cells,selected,inputs = load(args.input_dir)
    input_hashes = {str(path.relative_to(ROOT)):common.sha(path) for path in inputs}
    if old and old['input_sha256'] != input_hashes:
        raise ValueError('Output directory already belongs to different inputs')
    fonts = {f.name for f in font_manager.fontManager.ttflist}
    font = 'Times New Roman' if 'Times New Roman' in fonts else 'DejaVu Serif'
    plt.rcParams.update({'font.family':font,'font.size':10,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
    selected.to_csv(output/'PRIVATE_five_syndromes_selected.csv',index=False)
    if args.panel in ['a','both']:
        panel_a(cells,selected,output)
    if args.panel in ['b','both']:
        panel_b(primary,cells,selected,output)
    (output/'FIGURE_CAPTIONS.txt').write_text(CAPTION_A+'\n\n'+CAPTION_B+'\n')
    if protected != common.protected_hashes() or any(common.sha(ROOT/p)!=v for p,v in input_hashes.items()):
        raise ValueError('Protected inputs changed during plotting')
    manifest = {'study_id':'PRISM_RACE_STRATIFIED_20261115','input_sha256':input_hashes,
                'rows':len(primary),'eligible_cells':len(cells),'syndromes':cells.disease.nunique(),
                'features_per_cell':120,'original_discoveries':int(cells.significant_features.sum()),
                'selected_syndromes':selected.to_dict('records'),'criterion':CRITERION,'panel_B_features':FEATURES,
                'selection_uses_significance':False,'panel_B_scale':[-2,2],'primary_q_recalculated_for_display':False,
                'original_global_BH_independently_verified':True,'font':font,'png_dpi':300,
                'software':{x:importlib.metadata.version(x) for x in ['matplotlib','numpy','pandas','scipy']},
                'python':sys.version,'git_revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                'script_sha256':common.sha(Path(__file__)),'last_command_panel':args.panel,
                'deviations':['Referenced attached starter script was absent; implemented specification directly'],
                'protected_files_unchanged':True,'real_outputs_private':True,
                'output_sha256':{p.name:common.sha(p) for p in output.iterdir() if p.is_file() and p.name!='run_manifest.json'}}
    common.write_json(old_path,manifest)
    print(selected[['selection_rank','disease','total_eligible_patients','eligible_categories']].to_string(index=False))


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--panel',choices=['a','b','both'],default='both')
    parser.add_argument('--input-dir',type=Path,default=common.PRIVATE)
    parser.add_argument('--output-dir',type=Path,default=ROOT/'local_review/main_analysis_panels')
    args = parser.parse_args()
    main(args)

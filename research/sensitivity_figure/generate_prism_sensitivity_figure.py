import argparse
import hashlib
import importlib.metadata
import json
import subprocess
import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.text import Text
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from revision_results.race_stratified_20261115.common import private_output

SOURCE = ROOT/'research/race_stratified_sensitivity_results_20261115/summary.json'
EXPECTED_SOURCE_SHA256 = 'bdf73bc8f210e913a9676d6521bd46d459ae7cc37e45321833f56b05dc7d84df'
RELEASE_INPUT = Path(__file__).resolve().parent/'results_20261117/reference_comparison_inputs.json'
EXPECTED_RELEASE_INPUT_SHA256 = '69135fb782ab75bcffa2a12ded89f66de7883726b289410336a224f2b2cd662d'
BLUE, ORANGE, GRAY = '#315C7A','#D55E00','#BAC5CF'
CAPTION = (
    'Sensitivity of race-stratified FaceKit syndrome–feature associations. '
    '(A) Among {B:,} FairFace-only placebo iterations, {zero:,} produced no globally BH-significant discoveries and '
    '{nonzero} produced at least one. Only {supported:,}/{planned:,} planned hypotheses were supported; European '
    'placebo contrasts were unsupported because the reference sample exhausted the available pool. The mean was '
    '{mean:.3f} discoveries per run (maximum {maximum}). This same-source placebo does not test GMDB–FairFace '
    'source effects or establish a universally calibrated false-positive rate. '
    '(B) {permutations:,} permutations per null yielded {unblocked} empirical point-estimate BH discoveries under '
    'category-only shuffling and {blocked} under category/sex/age-blocked shuffling, with {overlap_unblocked} and '
    '{overlap_blocked} respective overlaps with primary findings in the same direction. '
    '{cross} associations met cross-design point-estimate criteria; {certain} met the conservative simultaneous '
    'Monte Carlo upper-bound criterion. Each null uses the planned {planned:,}-comparison correction family. '
    'The null designs are not independent replications or nested sets; their differing counts do not establish '
    'a causal effect of blocking. '
    '(C) Reference resampling retained {african_stable}/{african_total} African and {asian_stable}/{asian_total} '
    'primary associations. Stable means original-effect sign agreement in ≥95% and global BH q<0.05 in ≥80% '
    'of {reference_B:,} draws. Among groups with varied reference photographs, {varied_stable}/{varied_total} '
    '({varied_percent:.1f}%) met the rule. European reference photographs could not vary; '
    '{european_stable}/{european_total} ({european_percent:.1f}%) recurrence therefore does not measure '
    'European reference-selection sensitivity. The overall {overall_stable}/{original} ({overall_percent:.1f}%) '
    'includes that fixed category. Overlapping reference draws are not independent patient-level replication. '
    'These analyses evaluate internal sensitivity and establish neither independent clinical validation nor '
    'elimination of image-source confounding. Demographic labels are recorded annotations, not verified genetic ancestry.'
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate(data):
    if (data['status']!='COMPLETED' or data['study_id']!='PRISM_THREE_RACE_SENSITIVITY_20261115'
            or data['baseline_id']!='PRISM_RACE_STRATIFIED_20261115'):
        raise ValueError('Completed November 15 race-stratified sensitivity summary required')
    if data['planned_hypotheses']!=9600 or data['eligible_cells']!=80 or data['observed_main_discoveries']!=982:
        raise ValueError('Unexpected primary study dimensions')
    placebo = data['placebo']
    if placebo['B']!=1000 or placebo['supported_hypotheses']+placebo['unsupported_hypotheses']!=data['planned_hypotheses']:
        raise ValueError('Inconsistent placebo family/iterations')
    nonzero = placebo['fraction_at_least_one']*placebo['B']
    if not np.isclose(nonzero,round(nonzero)) or not 0<=nonzero<=placebo['B']:
        raise ValueError('Invalid placebo nonzero-run count')
    if not nonzero<=placebo['mean']*placebo['B']<=nonzero*placebo['maximum']:
        raise ValueError('Placebo mean/max cannot agree with discovery-run count')
    group_rows = data['reference_stability_by_race']
    groups = {r['canonical_race']:r for r in group_rows}
    if len(group_rows)!=3 or set(groups)!={'African','Asian','European'}:
        raise ValueError('Three unique recorded reference categories required')
    if not (groups['African']['reference_images_vary'] is True and groups['Asian']['reference_images_vary'] is True
            and groups['European']['reference_images_vary'] is False):
        raise ValueError('Unexpected fixed/varied reference status')
    for row in groups.values():
        if not 0<=row['descriptive_stable']<=row['original_discoveries']:
            raise ValueError('Invalid subgroup stability totals')
    ref = data['reference_resampling']
    if (ref['B']!=1000 or sum(r['original_discoveries'] for r in groups.values())!=data['observed_main_discoveries']
            or sum(r['descriptive_stable'] for r in groups.values())!=ref['stable_main_associations']
            or ref['stable_main_associations']+ref['unstable_main_associations']!=data['observed_main_discoveries']):
        raise ValueError('Reference subgroup/global totals disagree')
    if sum(data['reference_quotas'].values())!=data['reference_X'] or data['reference_quotas']['European']!=647:
        raise ValueError('Frozen European reference allocation differs')
    for spec in ['within_race','within_race_sex_age']:
        row = data['gmdb_permutation'][spec]
        if row['B']!=100000 or not np.isclose(row['p_floor'],1/(row['B']+1)):
            raise ValueError('Unexpected permutation repetitions/p resolution')
        if not 0<=row['overlap_main_same_direction']<=row['overlap_main']<=row['discoveries']:
            raise ValueError('Invalid empirical discovery/overlap counts')
        if row['MC_upper_q_below_005']!=0:
            raise ValueError('Expected zero simultaneous Monte Carlo upper-bound discoveries')
    if data['simultaneous_MC_certain_cross_design_candidates']!=0:
        raise ValueError('Cannot present Monte Carlo-certain cross-design discoveries')
    if not 0<=data['point_estimate_cross_design_candidates']<=min(r['overlap_main_same_direction'] for r in data['gmdb_permutation'].values()):
        raise ValueError('Invalid cross-design point candidate count')
    return groups


def values(data):
    groups = validate(data)
    placebo = data['placebo']
    unblocked,blocked = [data['gmdb_permutation'][x] for x in ['within_race','within_race_sex_age']]
    varied = [r for r in groups.values() if r['reference_images_vary']]
    stable,total = sum(r['descriptive_stable'] for r in varied),sum(r['original_discoveries'] for r in varied)
    return dict(B=placebo['B'],nonzero=round(placebo['B']*placebo['fraction_at_least_one']),
                zero=placebo['B']-round(placebo['B']*placebo['fraction_at_least_one']),
                supported=placebo['supported_hypotheses'],planned=data['planned_hypotheses'],mean=placebo['mean'],maximum=placebo['maximum'],
                permutations=unblocked['B'],unblocked=unblocked['discoveries'],blocked=blocked['discoveries'],
                overlap_unblocked=unblocked['overlap_main_same_direction'],overlap_blocked=blocked['overlap_main_same_direction'],
                other_unblocked=unblocked['discoveries']-unblocked['overlap_main_same_direction'],
                other_blocked=blocked['discoveries']-blocked['overlap_main_same_direction'],
                cross=data['point_estimate_cross_design_candidates'],certain=data['simultaneous_MC_certain_cross_design_candidates'],
                reference_B=data['reference_resampling']['B'],varied_stable=stable,varied_total=total,varied_percent=100*stable/total,
                overall_stable=data['reference_resampling']['stable_main_associations'],original=data['observed_main_discoveries'],
                overall_percent=100*data['reference_resampling']['stable_main_associations']/data['observed_main_discoveries'],
                **{f'{race.lower()}_{key}':r[field] for race,r in groups.items()
                   for key,field in [('stable','descriptive_stable'),('total','original_discoveries')]},
                european_percent=100*groups['European']['descriptive_stable']/groups['European']['original_discoveries'])


def ticks(ax, y, maximum, positions):
    ax.plot([0,maximum],[y,y],color='#888888',linewidth=.6)
    for x in positions:
        ax.plot([x,x],[y,y-.018],color='#888888',linewidth=.6)
        ax.text(x,y-.035,f'{x:g}',ha='center',va='top',fontsize=8)


def setup(ax,title,maximum):
    ax.set_xlim(0,maximum)
    ax.set_ylim(0,1)
    ax.axis('off')
    ax.set_title(title,loc='left',fontsize=11,pad=14,weight='bold')


def panel_a(ax,v):
    setup(ax,'A   FairFace-only same-source placebo',100)
    ax.text(0,.97,f"{v['B']:,} iterations",transform=ax.transAxes,fontsize=9)
    zero_percent=100*v['zero']/v['B']
    nonzero_percent=100*v['nonzero']/v['B']
    ax.barh(.76,zero_percent,height=.15,color=BLUE)
    ax.barh(.76,nonzero_percent,left=zero_percent,height=.15,color=ORANGE)
    ax.text(zero_percent/2,.76,f"{v['zero']:,} runs · zero discoveries ({zero_percent:.1f}%)",ha='center',va='center',color='white',fontsize=10)
    ticks(ax,.64,100,[0,25,50,75,100])
    ax.text(0,.49,f"≥1 discovery: {v['nonzero']} runs ({nonzero_percent:.1f}%)",transform=ax.transAxes,fontsize=10,color=ORANGE,weight='bold')
    ax.text(0,.37,f"Mean {v['mean']:.3f} discoveries per run; maximum {v['maximum']}",transform=ax.transAxes,fontsize=9)
    ax.text(0,.25,f"Supported: {v['supported']:,} / {v['planned']:,} planned hypotheses",transform=ax.transAxes,fontsize=10,weight='bold')
    ax.text(0,.13,'European pseudo-cases unavailable: all 647 reference images used.',transform=ax.transAxes,fontsize=9)
    ax.text(0,.02,'Does not test the GMDB–FairFace source effect.',transform=ax.transAxes,fontsize=9)


def panel_b(ax,v):
    setup(ax,'B   Within-GMDB label permutations',120)
    ax.text(0,.98,f"{v['permutations']:,} permutations per null · {v['planned']:,}-test correction families",transform=ax.transAxes,fontsize=9)
    for y,label,overlap,other,total in [(.81,'Within recorded category',v['overlap_unblocked'],v['other_unblocked'],v['unblocked']),
                                      (.57,'Within category + sex/age blocks',v['overlap_blocked'],v['other_blocked'],v['blocked'])]:
        ax.text(0,y+.075,label,fontsize=9)
        ax.barh(y,overlap,height=.11,color=BLUE)
        ax.barh(y,other,left=overlap,height=.11,color=GRAY)
        ax.text(overlap/2,y,str(overlap),ha='center',va='center',color='white',fontsize=10)
        ax.text(overlap+other/2,y,str(other),ha='center',va='center',color='#222222',fontsize=10)
        ax.text(total+2,y,str(total),va='center',fontsize=10,weight='bold')
    ticks(ax,.46,120,[0,25,50,75,100,120])
    ax.text(.50,.28,'Other empirical discoveries',transform=ax.transAxes,fontsize=9)
    ax.plot(.475,.295,marker='s',markersize=7,color=GRAY,transform=ax.transAxes)
    ax.text(.025,.28,'Overlap with primary, same direction',transform=ax.transAxes,fontsize=9)
    ax.plot(0,.295,marker='s',markersize=7,color=BLUE,transform=ax.transAxes,clip_on=False)
    ax.text(0,.14,f"{v['cross']} cross-design point candidates; {v['certain']} pass simultaneous MC upper bounds.",transform=ax.transAxes,fontsize=10,weight='bold',color='#8B2F3F')
    ax.text(0,.025,'Empirical point-estimate BH discoveries; null designs are not replications.',transform=ax.transAxes,fontsize=9)


def comparison_values(data, summary, plan):
    if plan['baseline_id'] != data['baseline_id'] or plan['planned_family'] != data['planned_hypotheses']:
        raise ValueError('Reference diagnostic baseline/family mismatch')
    if plan['requested_repetitions'] != 1000 or plan['original_quotas'] != data['reference_quotas']:
        raise ValueError('Reference diagnostic repetitions/original allocation mismatch')
    if plan['stability_rule'] != 'q<.05 in >=80% and same finite nonzero original-g sign in >=95% of draws':
        raise ValueError('Reference stability rule differs')
    if summary[['design', 'canonical_race']].duplicated().any():
        raise ValueError('Duplicate reference summary keys')
    if plan['ladder'][0] != data['reference_quotas']:
        raise ValueError('Reference ladder does not start at the original allocation')
    candidates = [(i, sizes) for i, sizes in enumerate(plan['ladder']) if i > 0
                  and all(data['reference_quotas'][r] < sizes[r] < plan['pools'][r]
                          for r in ['African', 'Asian'])
                  and sizes['European'] == data['reference_quotas']['European']]
    if not candidates:
        raise ValueError('No larger resampleable paired configuration available')
    index, sizes = max(candidates, key=lambda item: sum(item[1][r] for r in ['African', 'Asian']))
    result = {'design': f'C_stage_{index}', 'reference_n': sizes, 'pool_n': plan['pools'], 'repetitions': 1000,
              'selection_rule': 'Largest paired resampleable configuration by reference size, not by outcomes', 'groups': {}}
    original = {r['canonical_race']: r for r in data['reference_stability_by_race']}
    for race in ['African', 'Asian']:
        first = summary.loc[summary.design.eq('A_original') & summary.canonical_race.eq(race)]
        last = summary.loc[summary.design.eq(result['design']) & summary.canonical_race.eq(race)]
        if len(first) != 1 or len(last) != 1:
            raise ValueError('Missing original/larger reference rows')
        first, last = first.iloc[0], last.iloc[0]
        if (first.reference_n != data['reference_quotas'][race]
                or first.stable_original != original[race]['descriptive_stable']
                or first.original_discoveries != original[race]['original_discoveries']
                or first.repetitions != 1000 or not first.reference_selection_varies
                or first.sampling_scheme != 'original_without_replacement_subset'):
            raise ValueError('Original reference results do not match public summary')
        if (last.reference_n != sizes[race] or last.repetitions != 1000 or not last.reference_selection_varies
                or last.sampling_scheme != 'without_replacement_subset'
                or last.original_discoveries != first.original_discoveries
                or not 0 <= last.stable_original <= last.original_discoveries
                or not float(last.stable_original).is_integer()):
            raise ValueError('Larger comparison is not a valid same-denominator subset experiment')
        if any(not np.isclose(row.stable_original_percent, 100 * row.stable_original / row.original_discoveries)
               for row in [first, last]):
            raise ValueError('Larger-reference percentage/count mismatch')
        result['groups'][race] = {'stable': int(last.stable_original), 'total': int(last.original_discoveries)}
    result['stable'] = sum(r['stable'] for r in result['groups'].values())
    result['total'] = sum(r['total'] for r in result['groups'].values())
    result['percent'] = 100 * result['stable'] / result['total']
    return result


def load_comparison(data, source='private'):
    if source == 'released':
        if sha(RELEASE_INPUT) != EXPECTED_RELEASE_INPUT_SHA256:
            raise ValueError('Released reference plotting input changed; explicit revalidation required')
        bundle = json.loads(RELEASE_INPUT.read_text())
        if (bundle['schema_version'] != 1 or bundle['status'] != 'COMPLETED'
                or bundle['original_summary_sha256'] != EXPECTED_SOURCE_SHA256):
            raise ValueError('Released reference plotting input provenance mismatch')
        result = comparison_values(data, pd.DataFrame(bundle['summary_rows']), bundle['plan'])
        result['input_sha256'] = {str(RELEASE_INPUT.relative_to(ROOT)): sha(RELEASE_INPUT)}
        return result
    if source != 'private':
        raise ValueError('Unknown reference comparison source')
    directory = ROOT / 'local_review/reference_stability_diagnostic'
    paths = [directory / 'source_and_execution_manifest.json', directory / 'reference_size_ladder_summary.csv',
             directory / 'verification_report.json']
    if any(not path.is_file() for path in paths):
        raise FileNotFoundError('Verified reference diagnostic aggregates are required for updated Panel C')
    manifest = json.loads(paths[0].read_text())
    verification = json.loads(paths[2].read_text())
    if manifest['status'] != 'COMPLETED' or verification['final_status'] != 'PASS':
        raise ValueError('Reference diagnostic is incomplete or unverified')
    for path in paths[1:]:
        if sha(path) != manifest['output_sha256'][path.name]:
            raise ValueError('Reference diagnostic aggregate hash mismatch: ' + path.name)
    result = comparison_values(data, pd.read_csv(paths[1], float_precision='round_trip'), manifest['plan'])
    result['input_sha256'] = {str(path.relative_to(ROOT)): sha(path) for path in paths}
    return result


def panel_c(ax, v, comparison):
    setup(ax, 'C   Reference stability: original versus larger subsets', 100)
    ax.text(0, .98, f"{comparison['repetitions']:,} without-replacement draws per design; same original discoveries", transform=ax.transAxes, fontsize=9)
    expanded = '#117733'
    for heading, y1, y2, race in [(.91, .81, .67, 'African'), (.55, .45, .31, 'Asian')]:
        total = v[race.lower() + '_total']
        initial = v[race.lower() + '_stable']
        later = comparison['groups'][race]['stable']
        ax.text(0, heading, race + ' recorded category', fontsize=10, weight='bold')
        for y, count, reference_n, color, label in [
                (y1, initial, v[race.lower() + '_reference_n'], BLUE, 'Original'),
                (y2, later, comparison['reference_n'][race], expanded, 'Larger subset')]:
            pct = 100 * count / total
            ax.barh(y, pct, height=.105, color=color)
            ax.text(pct / 2, y, f'{label}: {count}/{total} ({pct:.1f}%) · ref n = {reference_n}',
                    ha='center', va='center', color='white', fontsize=9)
    ticks(ax, .22, 100, [0, 25, 50, 75, 100])
    ax.text(0, .115, f"Combined: {v['varied_stable']}/{v['varied_total']} ({v['varied_percent']:.1f}%) → "
            f"{comparison['stable']}/{comparison['total']} ({comparison['percent']:.1f}%).", transform=ax.transAxes, fontsize=10, weight='bold')
    ax.text(0, .005, 'Stable: global q < 0.05 in ≥80% of draws; original g sign in ≥95%.\n'
            'European references fixed. Full-pool and bootstrap results are separate.',
            transform=ax.transAxes, fontsize=8.5, linespacing=1.2)


PANELS={'a':panel_a,'b':panel_b,'c':panel_c}


def save(fig,path):
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    for item in fig.findobj(Text):
        if item.get_visible() and item.get_text():
            box=item.get_window_extent(renderer)
            if box.x0<-.5 or box.y0<-.5 or box.x1>fig.bbox.x1+.5 or box.y1>fig.bbox.y1+.5:
                raise ValueError('Clipped text: '+item.get_text())
    for ext in ['png','pdf','svg']:
        fig.savefig(path.with_suffix('.'+ext),dpi=300,facecolor='white')
    plt.close(fig)


def main(args):
    if not SOURCE.is_file():
        raise FileNotFoundError('BLOCKED: authoritative public sensitivity summary missing')
    source_hash=sha(SOURCE)
    if source_hash!=EXPECTED_SOURCE_SHA256:
        raise ValueError('BLOCKED: source summary changed; explicit revalidation required')
    data=json.loads(SOURCE.read_text())
    v=values(data)
    comparison = load_comparison(data,args.comparison_source) if args.panel in ['c', 'composite', 'all'] else None
    for race in ['African', 'Asian']:
        v[race.lower() + '_reference_n'] = data['reference_quotas'][race]
    output=private_output(args.output_dir)
    archive = output / 'archive_original_panel_C'
    if comparison is not None and not archive.exists() and (output / 'sensitivity_panel_C.png').exists():
        archive.mkdir()
        for name in ['sensitivity_panel_C', 'sensitivity_ABC_composite']:
            for extension in ['png', 'pdf', 'svg']:
                path = output / (name + '.' + extension)
                if path.exists():
                    shutil.copy2(path, archive / path.name)
        for name in ['FIGURE_CAPTION.txt', 'figure_manifest.json']:
            if (output / name).exists():
                shutil.copy2(output / name, archive / name)
    git_before=subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True)
    font='Times New Roman' if 'Times New Roman' in {f.name for f in font_manager.fontManager.ttflist} else 'DejaVu Serif'
    plt.rcParams.update({'font.family':font,'font.size':9,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
    generated=[]
    if args.panel in ['composite','all']:
        fig,axes=plt.subplots(3,1,figsize=(7.5,10.7),gridspec_kw={'height_ratios':[1,1,1.35]})
        fig.subplots_adjust(left=.085,right=.96,bottom=.065,top=.925,hspace=.34)
        for key,ax in zip(PANELS,axes):
            panel_c(ax,v,comparison) if key == 'c' else PANELS[key](ax,v)
        fig.suptitle('Sensitivity of race-stratified FaceKit associations',fontsize=13,y=.985)
        fig.text(.085,.012,'Internal sensitivity only; no independent clinical validation\nor elimination of image-source confounding.',fontsize=9,linespacing=1.2)
        save(fig,output/'sensitivity_ABC_composite')
        generated.append('sensitivity_ABC_composite')
    for key,draw in PANELS.items():
        if args.panel in [key,'all']:
            fig,ax=plt.subplots(figsize=(7.5,4.5) if key == 'c' else (7.5,3.4))
            fig.subplots_adjust(left=.085,right=.96,bottom=.10,top=.83)
            panel_c(ax,v,comparison) if key == 'c' else draw(ax,v)
            save(fig,output/('sensitivity_panel_'+key.upper()))
            generated.append('sensitivity_panel_'+key.upper())
    caption = CAPTION.format(**v)
    if comparison is not None:
        caption = caption.split('(C)')[0] + (
                    f"(C) Original versus larger-reference subset stability among the same original discoveries. "
                    f"African: original ref n={v['african_reference_n']}, {v['african_stable']}/{v['african_total']} "
                    f"({100*v['african_stable']/v['african_total']:.1f}%) stable; larger ref n={comparison['reference_n']['African']}, "
                    f"{comparison['groups']['African']['stable']}/{comparison['groups']['African']['total']} "
                    f"({100*comparison['groups']['African']['stable']/comparison['groups']['African']['total']:.1f}%). "
                    f"Asian: original ref n={v['asian_reference_n']}, {v['asian_stable']}/{v['asian_total']} "
                    f"({100*v['asian_stable']/v['asian_total']:.1f}%) stable; larger ref n={comparison['reference_n']['Asian']}, "
                    f"{comparison['groups']['Asian']['stable']}/{comparison['groups']['Asian']['total']} "
                    f"({100*comparison['groups']['Asian']['stable']/comparison['groups']['Asian']['total']:.1f}%). "
                    f"Combined: {v['varied_stable']}/{v['varied_total']} ({v['varied_percent']:.1f}%) versus "
                    f"{comparison['stable']}/{comparison['total']} ({comparison['percent']:.1f}%). "
                    f"Each displayed design uses {comparison['repetitions']:,} draws. Stable means global BH q<0.05 "
                    'in ≥80% of draws and original-effect g-sign agreement in ≥95%. The larger allocation is the largest '
                    'paired ladder configuration that still varies reference selection in both categories. '
                    'This configuration was selected by reference size, not stability outcomes. Both displayed designs use without-replacement subsets and fixed cases; global q is recalculated across all 9,600 tests in the source experiments. '
                    'African/Asian sizes increase jointly, so global BH competition can contribute to q changes. Near-full subsets reduce finite-pool selection variation; this is not independent replication or a patient confidence interval. '
                    'European reference photographs remain fixed at n=647 and are omitted from the stability bars. '
                    'The deterministic full-pool diagnostic and its distinct with-replacement bootstrap are not plotted. '
                    'These analyses evaluate internal sensitivity and establish neither independent clinical validation nor '
                    'elimination of image-source confounding. Demographic labels are recorded annotations, not verified genetic ancestry.')
        plot_rows = []
        for race in ['African', 'Asian']:
            for label, reference_n, stable in [('Original', v[race.lower()+'_reference_n'], v[race.lower()+'_stable']),
                                               ('Larger subset', comparison['reference_n'][race], comparison['groups'][race]['stable'])]:
                total = v[race.lower()+'_total']
                plot_rows.append({'canonical_race':race, 'design':label, 'reference_n':reference_n, 'original_discoveries':total,
                                  'stable_original':stable, 'stable_percent':100*stable/total, 'repetitions':comparison['repetitions'],
                                  'sampling_scheme':'without_replacement_subset'})
        pd.DataFrame(plot_rows).to_csv(output/'panel_C_reference_comparison_values.csv',index=False)
    (output/'FIGURE_CAPTION.txt').write_text(caption+'\n')
    if sha(SOURCE)!=source_hash or git_before!=subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True):
        raise ValueError('Source or checkout changed during rendering')
    if comparison is not None and any(sha(ROOT / name) != expected for name, expected in comparison['input_sha256'].items()):
        raise ValueError('Reference diagnostic source changed during rendering')
    manifest={'reference_comparison':comparison,'original_panel_C_archive':str(archive.relative_to(output)) if archive.exists() else None,'source':str(SOURCE.relative_to(ROOT)),'source_sha256':source_hash,'study_id':data['study_id'],
              'baseline_id':data['baseline_id'],'status':data['status'],'script_sha256':sha(Path(__file__)),
              'panel_requested':args.panel,'figures_generated_this_run':generated,'data_values':v,
              'software':{k:importlib.metadata.version(k) for k in ['numpy','matplotlib']},'python':sys.version,
              'font':font,'PNG_dpi':300,'protected_checkout_status_unchanged':True,'source_unchanged':True,
              'raw_patient_data_read':False,'inference_rerun':False,'outputs_private':True,
              'source_limitations':['Placebo supports 2,760 hypotheses, not all 9,600','European reference photographs fixed',
                                    'No simultaneous Monte Carlo upper-bound discoveries','Overlapping resamples are not replication'],
              'deviations':['Companion script unavailable; implemented the supplied specification directly'],
              'output_sha256':{p.name:sha(p) for p in output.iterdir() if p.is_file() and p.name!='figure_manifest.json'}}
    (output/'figure_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'generated':generated,'source_sha256':source_hash,'plotted_values':v},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--panel',choices=['all','composite','a','b','c'],default='all')
    parser.add_argument('--output-dir',type=Path,default=ROOT/'local_review/sensitivity_summary_figure')
    parser.add_argument('--comparison-source',choices=['released','private'],default='released')
    args=parser.parse_args()
    main(args)

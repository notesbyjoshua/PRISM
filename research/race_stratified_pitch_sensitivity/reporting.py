import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.text import Text
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from research.race_stratified_pitch_sensitivity.design import ARMS, METRICS, RACES, full_family, locked

COLORS = ['#0072B2', '#D55E00']
LABELS = ['Random FairFace', 'Pitch-imbalanced FairFace']


def describe(values):
    values = np.asarray(values, float)
    return {'runs': len(values), 'runs_with_any': int(np.sum(values > 0)), 'fraction_with_any': float(np.mean(values > 0)),
            'mean': float(np.mean(values)), 'median': float(np.median(values)), 'q25': float(np.quantile(values, .25)),
            'q75': float(np.quantile(values, .75)), 'q95': float(np.quantile(values, .95)),
            'minimum': int(np.min(values)), 'maximum': int(np.max(values))}


def collect(study, plan, checkpoint):
    support, cells = study['support'], study['cells']
    supported = study['supported']
    supported_rows = support.loc[support.supported]
    support_races = supported_rows.canonical_race.to_numpy()
    feature_regions = study['definitions'].anatomical_region.to_numpy()
    regions = sorted(set(feature_regions))
    repetitions, category_gaps, pair_gaps, counts, feature_any, region_counts = [], [], [], [], [], []
    hits = np.zeros((2, 9600), dtype=int)
    for chunk in checkpoint.state['chunks']:
        with np.load(checkpoint.directory/chunk['file']) as arrays:
            for offset, rep in enumerate(range(chunk['start'], chunk['stop'])):
                for arm_index, arm in enumerate(ARMS):
                    p, g, q, found = full_family(arrays['p'][offset, arm_index], arrays['g'][offset, arm_index], supported, study['features'])
                    pose = arrays['pose'][offset, arm_index]
                    count = found.reshape(80, 120).sum(1)
                    unique = found.reshape(80, 120).any(0)
                    by_region = np.asarray([found.reshape(80, 120)[:, feature_regions == region].sum() for region in regions])
                    hits[arm_index] += found
                    counts.append(count)
                    feature_any.append(unique)
                    region_counts.append(by_region)
                    row = {'repetition': rep, 'arm': arm, 'planned_family': 9600,
                           'supported_cells': int(support.supported.sum()), 'supported_tests': int(supported.sum()),
                           'discoveries': int(found.sum()), 'unique_feature_names_discovered': int(unique.sum()),
                           'minimum_q': float(q.min()), 'minimum_supported_q': float(q[supported].min()),
                           'discoveries_abs_g_ge_0_8': int(np.sum(found & (np.abs(g) >= .8))),
                           'undefined_g_supported': int(np.sum(supported & ~np.isfinite(g))),
                           'mean_cell_abs_pitch_gap': float(pose[:, 12].mean()),
                           'median_cell_abs_pitch_gap': float(np.median(pose[:, 13])),
                           'mean_cell_signed_pitch_gap': float(pose[:, 17].mean()),
                           'mean_cell_signed_yaw_gap': float(pose[:, 19].mean()),
                           'mean_cell_abs_yaw_gap': float(pose[:, 20].mean()),
                           'mean_cell_signed_roll_gap': float(pose[:, 21].mean()),
                           'mean_cell_abs_roll_gap': float(pose[:, 22].mean())}
                    for race in RACES:
                        take = support.canonical_race.eq(race).to_numpy()
                        row[race+'_discoveries'] = int(count[take].sum())
                        take = support_races == race
                        category_gaps.append({'repetition': rep, 'arm': arm, 'canonical_race': race,
                                              'supported_cells': int(np.sum(take)),
                                              **{metric: float(pose[take, i].mean()) if take.any() else np.nan
                                                 for i, metric in enumerate(METRICS)}})
                    repetitions.append(row)
                    for i, cell in enumerate(supported_rows.itertuples()):
                        pair_gaps.append({'repetition': rep, 'arm': arm, 'cell_index': cell.cell_index,
                                          'disease': cell.disease, 'canonical_race': cell.canonical_race,
                                          **dict(zip(METRICS, pose[i]))})
    return (pd.DataFrame(repetitions), pd.DataFrame(category_gaps), pd.DataFrame(pair_gaps),
            np.asarray(counts), np.asarray(feature_any), np.asarray(region_counts), hits, regions)


def save_figure(fig, output, name):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    outside_ticks = set()
    for ax in fig.axes:
        for axis in [ax.xaxis, ax.yaxis]:
            lower, upper = sorted(axis.get_view_interval())
            for tick in axis.get_major_ticks()+axis.get_minor_ticks():
                if not lower <= tick.get_loc() <= upper:
                    outside_ticks.update([tick.label1, tick.label2])
    for item in fig.findobj(Text):
        if item.get_visible() and item.get_text() and item not in outside_ticks:
            box = item.get_window_extent(renderer)
            if box.x0 < -.5 or box.y0 < -.5 or box.x1 > fig.bbox.x1+.5 or box.y1 > fig.bbox.y1+.5:
                raise ValueError('Clipped figure text: '+item.get_text())
    for extension in ['png', 'pdf', 'svg']:
        fig.savefig(output/(name+'.'+extension), dpi=300, facecolor='white')
    plt.close(fig)


def figures(output, repetitions, categories, pose_summary, supported_cells, supported_tests, european_reference):
    font = 'Times New Roman' if 'Times New Roman' in {f.name for f in font_manager.fontManager.ttflist} else 'DejaVu Serif'
    plt.rcParams.update({'font.family': font, 'font.size': 10, 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    distributions = [repetitions.loc[repetitions.arm.eq(arm), 'discoveries'].to_numpy() for arm in ARMS]
    maximum = max(values.max() for values in distributions)
    fig, ax = plt.subplots(figsize=(7.4, 4.9))
    jitter = np.random.default_rng(88341)
    boxes = ax.boxplot(distributions, positions=[0, 1], widths=.42, patch_artist=True, showfliers=False,
                       medianprops={'color': 'black', 'linewidth': 1.4})
    for i, (values, box) in enumerate(zip(distributions, boxes['boxes'])):
        box.set_facecolor(COLORS[i])
        box.set_alpha(.4)
        ax.scatter(i+jitter.uniform(-.16, .16, len(values)), values, s=7, alpha=.12, color=COLORS[i], edgecolors='none', rasterized=False)
        ax.scatter(i, values.mean(), marker='D', s=24, color='black', zorder=4)
        ax.text(i, max(1, maximum)*1.09,
                f'n = {len(values):,}\nMean {values.mean():.3f} · median {np.median(values):g}\n95th percentile {np.quantile(values, .95):g}\nAny discovery: {np.sum(values > 0)}/{len(values)}',
                ha='center', va='bottom', fontsize=9)
    ax.set_xticks([0, 1], LABELS)
    ax.set(ylabel='Global BH-significant pseudo-cohort tests', ylim=(-max(1, maximum)*.025, max(1, maximum)*1.34))
    ax.set_title('Paired same-source placebo discovery distributions', fontsize=12, pad=12)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.spines[['top', 'right']].set_visible(False)
    fig.text(.5, .025, f'Both arms: {supported_cells} supported cells · {supported_tests:,} / 9,600 planned tests\nEuropean contrasts unsupported at original reference n = {european_reference}', ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .10, 1, 1))
    save_figure(fig, output, 'pitch_placebo_distribution')
    fig, axes = plt.subplots(1, 2, figsize=(9.1, 4.5))
    locations = np.arange(3)
    for arm_index, arm in enumerate(ARMS):
        for race_index, race in enumerate(RACES):
            row = categories.loc[categories.arm.eq(arm) & categories.canonical_race.eq(race)].iloc[0]
            if not row.supported_tests:
                continue
            gap = pose_summary.loc[pose_summary.arm.eq(arm) & pose_summary.canonical_race.eq(race) & pose_summary.metric.eq('mean_abs_pitch_gap'), 'mean'].iloc[0]
            x = race_index+(arm_index-.5)*.3
            axes[0].bar(x, gap, width=.28, color=COLORS[arm_index], label=LABELS[arm_index] if race_index == 0 else None)
            axes[0].annotate(f'{gap:.2f}', (x, gap), xytext=(0, 4 if gap >= 0 else -12), textcoords='offset points', ha='center', fontsize=8)
            value = 100*row.fraction_with_any
            axes[1].bar(x, value, width=.28, color=COLORS[arm_index])
            axes[1].annotate(f'{int(row.runs_with_any)}/{int(row.runs)}', (x, value), xytext=(0, 5), textcoords='offset points', ha='center', fontsize=8)
    axes[0].set_title('(a) Achieved absolute-pitch separation', fontsize=11)
    axes[0].set_ylabel('Mean cell |pitch| gap (degrees)\n(pseudo-case − reference)')
    axes[0].axhline(0, color='#777777', linewidth=.6)
    axes[0].margins(y=.2)
    axes[1].set_title('(b) Repetitions with any discovery', fontsize=11)
    axes[1].set(ylabel='Repetitions with ≥1 discovery (%)', ylim=(0, 118), yticks=[0, 25, 50, 75, 100])
    for ax in axes:
        ax.set_xticks(locations, ['African', 'Asian', 'European'])
        ax.set_xlim(-.55, 2.55)
        ax.spines[['top', 'right']].set_visible(False)
        ax.axvspan(1.65, 2.45, facecolor='#F2F2F2', hatch='///', edgecolor='#CCCCCC', alpha=.5)
        ax.text(2, .48, f'N/A\nreference n = {european_reference}', transform=ax.get_xaxis_transform(), ha='center', va='center', fontsize=9)
    axes[0].legend(loc='upper center', bbox_to_anchor=(.5, -.13), frameon=False, fontsize=9)
    fig.text(.5, .02, 'Common support mask and original 9,600-test global BH family in both arms; N/A means unsupported.', ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .12, 1, 1))
    save_figure(fig, output, 'pitch_imbalance_diagnostics')
    return font


def write_reports(study, plan, output, checkpoint):
    reps, gaps, pairs, counts, any_feature, region_counts, hits, regions = collect(study, plan, checkpoint)
    reps.to_csv(output/'repetition_summary.csv', index=False)
    gaps.to_csv(output/'pose_gap_by_repetition_category.csv', index=False)
    pairs.to_csv(output/'repetition_pose_gaps.csv', index=False)
    summaries = {arm: describe(reps.loc[reps.arm.eq(arm), 'discoveries']) for arm in ARMS}
    summary = pd.DataFrame([{'arm': arm, **summaries[arm]} for arm in ARMS])
    summary.to_csv(output/'discovery_summary_overall.csv', index=False)
    pivot = reps.pivot(index='repetition', columns='arm', values='discoveries')
    paired = pd.DataFrame({'repetition': pivot.index, 'random': pivot.random, 'pitch_imbalanced': pivot.pitch_imbalanced})
    paired['difference_B_minus_A'] = paired.pitch_imbalanced-paired.random
    paired.to_csv(output/'paired_discovery_differences.csv', index=False)
    category_rows, region_rows, cell_rows, feature_rows, correlation_rows, pose_correlation_rows = [], [], [], [], [], []
    pose_rows = []
    for (arm, race), frame in gaps.groupby(['arm', 'canonical_race']):
        for metric in METRICS:
            values = frame[metric].dropna().to_numpy()
            pose_rows.append({'arm': arm, 'canonical_race': race, 'metric': metric, 'units': 'fraction' if 'fraction' in metric else 'degrees',
                              'runs': len(values), **{stat: float(function(values)) if len(values) else np.nan
                                                     for stat, function in [('mean', np.mean), ('median', np.median),
                                                                            ('q25', lambda x: np.quantile(x, .25)),
                                                                            ('q75', lambda x: np.quantile(x, .75)),
                                                                            ('minimum', np.min), ('maximum', np.max)]}})
    pose_summary = pd.DataFrame(pose_rows)
    pose_summary.to_csv(output/'pose_gap_summary.csv', index=False)
    for arm_index, arm in enumerate(ARMS):
        take = reps.arm.eq(arm).to_numpy()
        for race in RACES:
            supported = study['support'].canonical_race.eq(race) & study['support'].supported
            n = int(supported.sum())
            row = {'arm': arm, 'canonical_race': race, 'eligible_cells': int(study['support'].canonical_race.eq(race).sum()),
                   'supported_cells': n, 'supported_tests': n*120, 'status': 'supported' if n else 'UNSUPPORTED',
                   'n_reference': study['quotas'][race]}
            row.update(describe(reps.loc[take, race+'_discoveries']) if n else {k: np.nan for k in describe([0])})
            category_rows.append(row)
            if n:
                own_gaps = gaps.loc[gaps.arm.eq(arm) & gaps.canonical_race.eq(race)].sort_values('repetition')
                discovery_counts = reps.loc[take].sort_values('repetition')[race+'_discoveries'].to_numpy()
                for metric in ['mean_abs_pitch_gap', 'median_abs_pitch_gap', 'abs_yaw_mean_gap', 'abs_roll_mean_gap']:
                    values = own_gaps[metric].to_numpy()
                    varied = np.ptp(values) > 0 and np.ptp(discovery_counts) > 0
                    pose_correlation_rows.append({'arm': arm, 'canonical_race': race, 'metric': metric,
                                                  'spearman_rho': float(spearmanr(values, discovery_counts).statistic) if varied else np.nan,
                                                  'status': 'descriptive_only' if varied else 'undefined_constant_vector'})
        for i, region in enumerate(regions):
            n_features = int(study['definitions'].anatomical_region.eq(region).sum())
            region_rows.append({'arm': arm, 'anatomical_region': region, 'unique_features_in_region': n_features,
                                'supported_tests': int(study['support'].supported.sum())*n_features,
                                'total_discoveries_across_repetitions': int(region_counts[take, i].sum()),
                                **describe(region_counts[take, i])})
        for i, row in enumerate(study['support'].itertuples()):
            cell_rows.append({'arm': arm, 'disease': row.disease, 'canonical_race': row.canonical_race,
                              'n_case': row.n_case, 'n_reference': row.n_reference, 'supported': row.supported,
                              **(describe(counts[take, i]) if row.supported else {key: np.nan for key in describe([0])})})
        for i, feature in enumerate(study['features']):
            feature_rows.append({'arm': arm, 'feature': feature, 'anatomical_region': study['definitions'].anatomical_region.iloc[i],
                                 'total_discovery_tests_across_repetitions': int(hits[arm_index].reshape(80, 120)[:, i].sum()),
                                 'runs_with_feature_discovered': int(any_feature[take, i].sum()),
                                 'fraction_runs_with_feature_discovered': float(any_feature[take, i].mean())})
        with np.errstate(invalid='ignore', divide='ignore'):
            correlations = np.corrcoef(any_feature[take].astype(float), rowvar=False)
        for i in range(120):
            for j in range(i+1, 120):
                correlation_rows.append({'arm': arm, 'feature_1': study['features'][i], 'feature_2': study['features'][j],
                                         'discovery_indicator_correlation': correlations[i, j],
                                         'status': 'descriptive_only' if np.isfinite(correlations[i, j]) else 'undefined_constant_indicator'})
    categories = pd.DataFrame(category_rows)
    categories.to_csv(output/'discovery_summary_by_category.csv', index=False)
    pd.DataFrame(region_rows).to_csv(output/'discovery_summary_by_region.csv', index=False)
    pd.DataFrame(cell_rows).to_csv(output/'discovery_summary_by_cell.csv', index=False)
    pd.DataFrame(feature_rows).to_csv(output/'feature_discovery_frequency.csv', index=False)
    pd.DataFrame(correlation_rows).to_csv(output/'feature_discovery_indicator_correlations.csv', index=False)
    pd.DataFrame(pose_correlation_rows).to_csv(output/'pose_discovery_correlations.csv', index=False)
    recurrence = pd.concat([study['primary'][locked.KEYS].assign(arm=arm, supported=study['supported'],
                          discovery_count=hits[i], discovery_frequency=np.where(study['supported'], hits[i]/plan['repetitions_per_arm'], np.nan))
                          for i, arm in enumerate(ARMS)], ignore_index=True)
    recurrence.to_csv(output/'test_discovery_frequencies.csv', index=False)
    reps[['repetition', 'arm', 'discoveries']].to_csv(output/'pitch_placebo_distribution_plot_values.csv', index=False)
    categories.to_csv(output/'pitch_imbalance_diagnostics_category_plot_values.csv', index=False)
    pose_summary.loc[pose_summary.metric.eq('mean_abs_pitch_gap')].to_csv(output/'pitch_imbalance_diagnostics_gap_plot_values.csv', index=False)
    font = figures(output, reps, categories, pose_summary, plan['supported_cells'], plan['supported_tests'], study['quotas']['European'])
    a, b = summaries['random'], summaries['pitch_imbalanced']
    difference = paired.difference_B_minus_A.to_numpy()
    history = study['benchmark']
    support_counts = study['support'].groupby('canonical_race').supported.sum().astype(int).to_dict()
    methods = (f'Using the locked November 15 race-stratified PRISM baseline (67 syndromes, 80 eligible cells, 736 eligible patients and 120 ordered measurements), '
               f'we ran {plan["repetitions_per_arm"]:,} repetitions per same-source FairFace arm with seed {plan["seed"]}. Pseudo-case sizes matched eligible syndrome/category counts; '
               'no GMDB feature vectors entered either arm. The random arm used disjoint random samples; the imbalanced arm sampled pseudo-cases from the upper half and references '
               'from the lower half of absolute estimated pitch, with deterministic lexical file ties and the odd observation assigned to the upper half. Original reference quotas '
               f'were 52 African, 190 Asian and 647 European. A common mask supported {plan["supported_cells"]} cells and {plan["supported_tests"]:,} tests in both arms. '
               'Each category used a shared reference disjoint from all its pseudo-cases where capacity permitted; any per-pair fallback is recorded in the support audit. '
               'Each arm/repetition used the original per-feature two-sided exact/asymptotic Mann–Whitney rules, exact-gamma-corrected Hedges g and one global BH family of '
               '9,600 p-values, with unsupported comparisons fixed at p=1. Significance was supported q<0.05.')
    results = (f'The random arm had a median of {a["median"]:g}, mean {a["mean"]:.3f} and maximum {a["maximum"]} significant pseudo-cohort tests per repetition; '
               f'{a["runs_with_any"]}/{a["runs"]} ({100*a["fraction_with_any"]:.1f}%) repetitions had at least one. The pitch-imbalanced arm had median {b["median"]:g}, '
               f'mean {b["mean"]:.3f} and maximum {b["maximum"]}; {b["runs_with_any"]}/{b["runs"]} ({100*b["fraction_with_any"]:.1f}%) had at least one. '
               f'The paired B−A discovery-count difference had median {np.median(difference):g}, mean {difference.mean():.3f}, IQR '
               f'{np.quantile(difference,.25):g}–{np.quantile(difference,.75):g} and range {difference.min()}–{difference.max()}. '
               'Category-specific discovery distributions and achieved pitch gaps are reported below. European contrasts were unsupported, not negative evidence.')
    limitations = ('These are same-source pseudo-cohort associations under a deliberate image-selection procedure, not syndrome associations or clinical discoveries. '
                   'The signed pose estimates are small-angle XYZ Tait–Bryan angles in degrees from the pinned MediaPipe facial transformation matrix. Absolute pitch is meaningful '
                   'as deviation from the extractor zero, but physical up/down orientation and pose-estimate accuracy are not independently validated. All photographs already passed '
                   'the original pose/QC gate; no new extraction or alternative normalization was performed. Pitch selection can also change yaw, roll, recorded sex, gaze, expression, '
                   'lighting and other capture attributes; it does not isolate a causal pitch effect. All FairFace ages lie in the same recorded 3–9 bin; comparable exact ages and '
                   'several acquisition covariates are unavailable. Dataset demographic categories do not establish genetic ancestry; deduplicated photographs do not establish unique '
                   'people or clinically healthy controls. Shared references, correlated features and overlapping draws induce dependence. Simulation quantiles describe this sampling '
                   'design, not patient-population confidence intervals. This experiment supplies neither a p-value for the original GMDB–FairFace study nor independent clinical validation.')
    lines = ['# Race-stratified head-pitch imbalance sensitivity', '',
             'Status: AWAITING_INDEPENDENT_VALIDATION' if plan['repetitions_per_arm'] == 1000 else 'Status: PILOT_AWAITING_VALIDATION', '',
             '## Manuscript-ready Methods', '', methods, '', '## Manuscript-ready Results', '', results, '',
             '## Manuscript-ready Limitations', '', limitations, '', '## Pose and common support', '',
             f'Finite stored pitch is present in every accepted category pool; no additional pose exclusions were required. Supported cells: {support_counts}. '
             'All 57 European cells (6,840 tests) retain p=1 because reference n=647 exhausts its pool; the low half additionally cannot supply that quota. '
             f'Newly unsupported relative to nominal random placebo: {int(study["support"].newly_unsupported_after_split.sum())}. '
             'Explicit cell identifiers and reasons remain private in support_audit.csv.', '',
             'Both supported categories use shared category references and globally disjoint pseudo-cases in each arm. Reuse across repetitions is expected. '
             'The two arms have independent deterministic category substreams; pairing refers to repetition indices and the locked design, not reused photo selections.', '',
             '## Overall discovery distributions', '', summary.to_markdown(index=False), '',
             '## Category-specific discovery distributions', '', categories.to_markdown(index=False), '',
             '## Achieved pose gaps', '',
             'Category metrics first average the cell-specific gap within each repetition, then summarize across repetitions. Gaps are pseudo-case minus reference. '
             'Signed and absolute yaw/roll gaps, sex-fraction differences, pitch quartiles, minimum separation and range overlap are retained in the CSVs. '
             'Unsupported categories have unavailable statistics, not zero effects.', '']
    for metric in ['mean_abs_pitch_gap', 'median_abs_pitch_gap', 'signed_yaw_mean_gap', 'abs_yaw_mean_gap', 'signed_roll_mean_gap', 'abs_roll_mean_gap']:
        lines.extend([f'### {metric} (degrees)', '', pose_summary.loc[pose_summary.metric.eq(metric), ['arm','canonical_race','mean','median','q25','q75','minimum','maximum']].to_markdown(index=False), ''])
    lines.extend(['## Historical benchmark and interpretation', '',
                  f'The separately completed unbiased race-stratified placebo had {round(history["fraction_at_least_one"]*history["B"])}/{history["B"]} runs with any discovery, '
                  f'mean {history["mean"]:.3f}, median {history["median"]:g}, maximum {history["maximum"]}, and {history["supported_hypotheses"]:,} supported tests. '
                  'The new random arm uses the same support and quotas but a new independent seed/substream schedule; its finite-run counts need not match exactly. '
                  'The benchmark is a reproducibility comparison, not an empirical p-value for the original clinical contrast.', '',
                  f'The observed mean-count difference was {b["mean"]-a["mean"]:.3f} tests per repetition; the difference in the fraction with any finding was '
                  f'{100*(b["fraction_with_any"]-a["fraction_with_any"]):.1f} percentage points. '
                  +(f'The mean-count ratio was {b["mean"]/a["mean"]:.2f} (imbalanced/random). ' if a['mean'] > 0 else 'A finite mean-count ratio is unavailable because the random mean was zero. ')
                  +'These are descriptive contrasts conditional on the frozen sampling design. Absolute counts and zero-discovery repetitions must be considered alongside relative increases; '
                  'the experiment does not estimate how many original GMDB/FairFace associations are explained by pose.', '',
                  'Earlier pooled pose-selection work also used same-source image-selection controls, but its cohort definitions, reference allocations and correction family differed. '
                  'Its numerical discoveries are not reused or treated as results of this study.', '',
                  'Feature/region frequencies and feature discovery-indicator correlations are descriptive. They do not identify a biological phenotype, syndrome specificity or independent anatomical signals. '
                  'Pose/discovery Spearman correlations are exploratory, without causal interpretation or inferential p-values.', '',
                  '## Figure captions', '',
                  '**Figure 1. Paired placebo discovery distributions.** Each point is one randomized repetition, not one patient. Boxes show median and interquartile range; '
                  'whiskers follow the 1.5-IQR rule, all repetitions are shown as jittered points, and diamonds mark means. Counts use supported global BH q<0.05, the same common mask '
                  f'({plan["supported_tests"]:,}/9,600 tests) and the original correction family in both arms. Mean, median, 95th percentile and repetition budgets are annotated. '
                  'The vertical scale is linear and includes the complete observed range.', '',
                  '**Figure 2. Diagnostic decomposition.** (a) Mean absolute-pitch gap, averaging cell-specific pseudo-case/reference mean differences within each repetition and then across '
                  'repetitions. (b) Fraction of repetitions with at least one supported BH discovery in each category, with exact numerator/denominator labels. '
                  'European contrasts are shaded N/A at their unchanged reference quota of 647. Both arms use the same support and 9,600-test BH family. '
                  'The contrast reflects deliberate pose-correlated selection, not an isolated causal manipulation of head pitch.', '',
                  '## Files and reproducibility', '',
                  'analysis_plan.json freezes the seed, numerical code/source hashes, features, support and allocation rules before outcomes. support_audit.csv contains every eligible cell. '
                  'Checkpoints store supported p/g vectors, achieved pose metrics and selection digests; unsupported p-values reconstruct as 1. All detailed intermediates remain private. '
                  'No selector, prior result, manuscript, publication figure or source dataset is modified. No result upload or result commit is performed.', '',
                  f'Rendering font: {font}. PNG 300 DPI; PDF/SVG vector. See validation_report.json and source_and_execution_manifest.json for independent verification and hashes.', ''])
    (output/'RESULTS_AND_LIMITATIONS.md').write_text('\n'.join(lines))
    locked.write_json(output/'aggregate_summary.json', {'arms': summaries,
                      'paired_difference': {'mean': float(difference.mean()), 'median': float(np.median(difference)),
                                            'minimum': int(difference.min()), 'maximum': int(difference.max())},
                      'supported_cells': plan['supported_cells'], 'supported_tests': plan['supported_tests'],
                      'planned_tests': 9600, 'repetitions_per_arm': plan['repetitions_per_arm'], 'font': font})

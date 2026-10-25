#!/usr/bin/env python3
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
COLORS = {'neither': '#a5a8ac', 'FairFace_only': '#b65a18', 'within_GMDB_only': '#6b4792', 'both': '#2166ac'}


def read(name):
    return pd.read_csv(OUT/name)


def label(syndrome):
    return re.sub(r'__[0-9a-f]{8}$', '', syndrome)


def save(fig, name):
    folder = OUT/'figures'
    folder.mkdir(exist_ok=True)
    fig.savefig(folder/f'{name}.pdf', bbox_inches='tight', metadata={'CreationDate': None, 'ModDate': None})
    fig.savefig(folder/f'{name}.png', dpi=350, bbox_inches='tight')
    plt.close(fig)


def figures():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8, 'axes.labelsize': 8, 'axes.titlesize': 9,
                         'xtick.labelsize': 7, 'ytick.labelsize': 7, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    joined = read('within_gmdb_vs_fairface_comparison.csv')
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for category, color in COLORS.items():
        frame = joined.loc[joined.discovery_category.eq(category)]
        ax.scatter(frame.fairface_g, frame.hedges_g, s=9 if category == 'neither' else 15,
                   c=color, alpha=.45 if category == 'neither' else .85, linewidths=0,
                   label=f'{category.replace("_", " ")} ({len(frame):,})')
    limit = max(joined.fairface_g.abs().max(), joined.hedges_g.abs().max())*1.05
    ax.plot([-limit, limit], [-limit, limit], color='#777777', ls='--', lw=.8)
    ax.axhline(0, lw=.6, color='#bbbbbb')
    ax.axvline(0, lw=.6, color='#bbbbbb')
    ax.set(xlabel="Syndrome minus FairFace: Hedges' g", ylabel="Syndrome minus other GMDB: Hedges' g",
           xlim=(-limit, limit), ylim=(-limit, limit))
    ax.set_title('Comparison of reference definitions across 6,125 hypotheses', loc='left')
    ax.legend(title='Global BH q < 0.05', frameon=False, loc='upper left', title_fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout(rect=(0, .07, 1, 1))
    fig.text(.09, .025, 'Each comparison uses its own 49 x 125 BH family. The comparator populations answer different questions.', fontsize=7)
    save(fig, 'within_gmdb_vs_fairface')
    highlights = read('highlighted_examples.csv')
    iterations = int(highlights.bootstrap_iterations.iloc[0])
    frozen = pd.read_csv(ROOT/'revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz')
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    for i, row in enumerate(highlights.itertuples()):
        old = frozen.loc[frozen.disease.eq(row.syndrome) & frozen.feature.eq(row.feature)].iloc[0]
        for y, g, lower, upper, color in [
            (i-.12, old.hedges_g, old.hedges_g_ci_95_low, old.hedges_g_ci_95_high, '#b65a18'),
            (i+.12, row.hedges_g, row.ci_95_lower, row.ci_95_upper, '#2166ac')]:
            ax.plot([lower, upper], [y, y], color=color, lw=1.5)
            ax.scatter(g, y, c=color, s=24, zorder=3)
        ax.text(2.35, i, f'{row.n_target}/{row.n_comparator}', va='center', fontsize=7)
    ax.set_yticks(range(4), [label(x.syndrome).replace('WILLIAMS-BEUREN SYNDROME; WBS', 'Williams-Beuren syndrome')+
                            '\n'+x.feature.replace('_', ' ') for x in highlights.itertuples()])
    ax.invert_yaxis()
    ax.axvline(0, color='#777777', lw=.8, ls='--')
    ax.set(xlabel="Hedges' g and marginal bootstrap 95% interval", xlim=(-2.4, 2.9), ylim=(3.5, -.65))
    ax.text(2.35, -.5, 'GMDB n\ntarget/other', fontsize=7)
    ax.set_title('Four prespecified measurement examples', loc='left')
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], color='#2166ac', marker='o', label='Other GMDB comparator'),
                       Line2D([], [], color='#b65a18', marker='o', label='FairFace comparator')],
              loc='lower right', frameon=False, fontsize=7)
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout(rect=(0, .08, 1, 1))
    fig.text(.04, .02, f'GMDB: {iterations:,} patient bootstraps. FairFace: frozen intervals (500 resamples). Normalized measures require clinical review.', fontsize=6.5)
    save(fig, 'highlighted_comparators')
    photo = read('sensitivity_effect_sizes.csv')
    random = read('random_photo_selection_summary.csv')
    thresholds = read('subgroup_threshold_sensitivity.csv')
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.5))
    ax = axes[0]
    for discordant, color in [(False, '#a5a8ac'), (True, '#b65a18')]:
        frame = photo.loc[photo.discordant_call.eq(discordant)]
        ax.scatter(frame.hedges_g_primary, frame.hedges_g_alternative, s=5, color=color, alpha=.5,
                   label=f'{"Discordant" if discordant else "Same"} call ({len(frame):,})')
    limits = [photo.hedges_g_primary.min()-.1, photo.hedges_g_primary.max()+.1]
    ax.plot(limits, limits, color='#777777', lw=.6, ls='--')
    ax.set(xlabel='Primary g', ylabel='Alternative g', title='A  Photograph rule')
    ax.legend(frameon=False, fontsize=5.5, loc='upper left')
    ax = axes[1]
    ax.hist(random.significant_frequency, bins=np.linspace(0, 1, 11), color='#2166ac', edgecolor='white')
    ax.set(xlabel='Significance frequency', ylabel='Hypotheses', title=f'B  {int(random.replicates.iloc[0])} random selections', xlim=(0, 1))
    ax = axes[2]
    for (analysis, layer), frame in thresholds.groupby(['analysis', 'layer']):
        frame = frame.sort_values('minimum_disease_n')
        ax.plot(frame.minimum_disease_n, 100*frame.significant_recorrected/frame.hypotheses,
                marker='o', ms=3, lw=1, ls='-' if layer == 'broad' else '--',
                color='#2166ac' if analysis == 'B' else '#b65a18', label=f'{analysis} {layer}')
    ax.set(xlabel='Minimum disease n/group', ylabel='Significant (%)', title='C  Restricted BH families')
    ax.set_xticks([2, 5, 10, 20])
    ax.legend(frameon=False, fontsize=5.5)
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
        ax.title.set_fontsize(8)
    fig.tight_layout(rect=(0, .14, 1, 1), w_pad=1.0)
    fig.text(.04, .035, 'A/B: same patients and eligibility; selection sensitivity. C: saved fits, restricted-family BH re-correction; no refitting.\n'
             'Counts and test denominators are provided in the accompanying source CSVs.', fontsize=6.5)
    save(fig, 'robustness_summary')


def reports():
    r, j = read('within_gmdb_results.csv'), read('within_gmdb_vs_fairface_comparison.csv')
    restricted, adjusted = read('within_gmdb_eligible49_sensitivity.csv'), read('within_gmdb_age_sex_adjusted.csv')
    h = read('highlighted_examples.csv')
    photos = read('sensitivity_effect_size_summary.csv')
    changes = read('photo_change_by_syndrome.csv')
    thresholds = read('subgroup_threshold_sensitivity.csv')
    random = read('random_photo_replicate_counts.csv')
    random_hyp = read('random_photo_selection_summary.csv')
    metadata = read('metadata_availability.csv')
    direction = read('feature_direction_summary.csv')
    categories = j.discovery_category.value_counts().to_dict()
    comparison = r[['syndrome', 'feature', 'hedges_g', 'q_global_bh']].merge(
        restricted[['syndrome', 'feature', 'hedges_g', 'q_global_bh']], on=['syndrome','feature'], suffixes=('_main', '_eligible49'),validate='one_to_one')
    summary = {'bootstrap_iterations': int(r.bootstrap_iterations.iloc[0]), 'photo_replicates': len(random),
               'target_missing_fraction_max': float(r.target_missing_fraction.max()),
               'comparator_missing_fraction_max': float(r.comparator_missing_fraction.max()), 'within_significant': int(r.significant.sum()), 'hypotheses': len(r), 'estimable': int(r.estimable.sum()),
               'positive_significant': int((r.significant&r.hedges_g.gt(0)).sum()),
               'negative_significant': int((r.significant&r.hedges_g.lt(0)).sum()),
               'pearson_g_between_comparators': j.hedges_g.corr(j.fairface_g),
               'spearman_g_between_comparators': j.hedges_g.corr(j.fairface_g, method='spearman'),
               'direction_agreement': j.same_direction.mean(), 'discovery_categories': categories,
               'bootstrap_ci_width_median': (r.ci_95_upper-r.ci_95_lower).median(),
               'bootstrap_finite_min': int(r.bootstrap_finite.min()), 'bootstrap_finite_max': int(r.bootstrap_finite.max()),
               'target_n_min': int(r.n_target.min()), 'target_n_max': int(r.n_target.max()),
               'other_n_min': int(r.n_comparator.min()), 'other_n_max': int(r.n_comparator.max()),
               'eligible49_significant': int(restricted.significant.sum()),
               'eligible49_vs_main_g_correlation': comparison.hedges_g_main.corr(comparison.hedges_g_eligible49),
               'eligible49_vs_main_discordant': int(comparison.q_global_bh_main.lt(.05).ne(comparison.q_global_bh_eligible49.lt(.05)).sum()),
               'adjusted_significant_four_model_family': int(adjusted.q_four_models.lt(.05).sum()),
               'syndromes_any_changed_photo': int(changes.patients_photo_changed.gt(0).sum()),
               'main_eligible_syndromes_any_changed_photo': int(changes.loc[changes.syndrome.isin(r.syndrome), 'patients_photo_changed'].gt(0).sum()),
               'random_photo_significant_min': int(random.significant.min()), 'random_photo_significant_median': float(random.significant.median()),
               'random_photo_significant_max': int(random.significant.max()),
               'random_photo_always_significant': int(random_hyp.significant_frequency.eq(1).sum()),
               'random_photo_ever_significant': int(random_hyp.significant_frequency.gt(0).sum())}
    (OUT/'RUN_SUMMARY.json').write_text(json.dumps(summary, indent=2))
    lines = ['# New methods and results for author review', '', 'Local-only aggregates; no public release or clinical interpretation has been approved.', '',
             '## Verified existing behavior', '',
             'The reconstructed primary cohort has 1,850 distinct recorded patient IDs and 379 syndrome labels. All 11 cohort-flow targets pass. '
             'The original disease table has zero missing recorded patient IDs. The reference has 942 unique normalized image paths after 58 duplicate rows '
             'are removed from 1,000. File uniqueness is not participant uniqueness, and recorded patient-ID uniqueness cannot detect duplicate people with different IDs.', '',
             'Frozen significance: 811/6,125 pooled comparisons, including 379 with |g| >= 0.8; alternative selection has 834 with 770 shared. '
             'Independent numerical checks are described in the unit-test and validation files; frozen bootstrap intervals were read, not re-created from the original RNG stream.', '',
             '## New implemented within-GMDB results', '',
             f'{summary["estimable"]:,} of {len(r):,} planned comparisons were estimable, with {summary["within_significant"]} meeting BH q < 0.05 '
             f'({summary["positive_significant"]} positive, {summary["negative_significant"]} negative). '
             f'Target n ranges {summary["target_n_min"]}-{summary["target_n_max"]}; the other-GMDB comparator ranges '
             f'{summary["other_n_min"]}-{summary["other_n_max"]}. Maximum feature missingness: '
             f'target {summary["target_missing_fraction_max"]:.2%}, other {summary["comparator_missing_fraction_max"]:.2%}. '
             f'The {summary["bootstrap_iterations"]:,}-resample marginal percentile interval width has median {summary["bootstrap_ci_width_median"]:.6g}; '
             f'finite bootstrap estimates range {summary["bootstrap_finite_min"]}-{summary["bootstrap_finite_max"]} per comparison.', '',
             'The main comparator contains every other age-eligible selected GMDB patient, including patients in smaller syndrome groups. '
             'Independent patient resampling within target/comparator groups keeps all 125 measurements together within each patient. '
             'BH is applied once to the 49 x 125 primary hypotheses. Statistical dependence among features and shared comparators remains; '
             'the test count is not a count of independent phenotypic discoveries. Percentile intervals are marginal and do not receive simultaneous coverage adjustment.', '',
             f'Across the same 6,125 frozen FairFace hypotheses, g Pearson r={summary["pearson_g_between_comparators"]:.6f}, '
             f'Spearman rho={summary["spearman_g_between_comparators"]:.6f}, and sign agreement={summary["direction_agreement"]:.2%}. '
             f'Significance categories: {categories}. These correlations are descriptive; no independence-based p-value is assigned.', '',
             f'Restricting the comparator to other patients in the 49 eligible syndromes gives {summary["eligible49_significant"]}/6,125 significant tests '
             f'under its own BH family, with g correlation {summary["eligible49_vs_main_g_correlation"]:.6f} and '
             f'{summary["eligible49_vs_main_discordant"]} different significance calls relative to the main comparator. '
             'This sensitivity has point estimates and q-values; no duplicate bootstrap interval analysis is claimed.', '',
             '### Four prespecified examples', '']
    for row in h.itertuples():
        lines.append(f'- {label(row.syndrome)} / {row.feature}: GMDB n={row.n_target}/{row.n_comparator}, '
                     f'g={row.hedges_g:.6g}, 95% bootstrap interval [{row.ci_95_lower:.6g}, {row.ci_95_upper:.6g}], '
                     f'global q={row.q_global_bh:.6g}; frozen FairFace g={row.fairface_g:.6g}, q={row.fairface_q:.6g}.')
    lines += ['', '### Recorded-age and sex sensitivity', '',
              'Four complete-case OLS models use target indicator, recorded age centered at 6, and Male/Female source label. '
              'Each target is compared with every other GMDB patient with complete included fields. '
              'HC3 standard errors and normal-reference confidence intervals/p-values are used; BH spans these four planned models. '
              'The coefficients are differences in normalized feature units and are not standardized g values. Unknown sex labels are excluded. '
              'This is a limited linear age specification, with no claim that inclusion removes residual age effects, sex-related imbalance, syndrome acquisition effects, or unmeasured confounding.', '']
    for row in adjusted.itertuples():
        lines.append(f'- {label(row.syndrome)} / {row.feature}: complete n={row.n_complete} (target {row.n_target}, other {row.n_other}), '
                     f'excluded {row.excluded}; beta={row.target_beta:.6g}, HC3 95% interval '
                     f'[{row.ci_95_lower:.6g}, {row.ci_95_upper:.6g}], q(4)={row.q_four_models:.6g}; '
                     f'age range overlap [{row.age_overlap_low:.3g}, {row.age_overlap_high:.3g}], max leverage={row.maximum_leverage:.4g}.')
    lines += ['', '## Photograph-selection sensitivity', '',
              f'Of 379 retained syndromes, {summary["syndromes_any_changed_photo"]} have at least one changed selected photograph; '
              f'{summary["main_eligible_syndromes_any_changed_photo"]} of the 49 main-eligible syndromes do. '
              'The exact two-rule comparison retains the original patients and has 105 discordant significance calls.', '']
    for row in photos.itertuples():
        lines.append(f'- {row.subset}: n={row.n_hypotheses}; median |delta g|={row.median_absolute_difference:.6g}, '
                     f'90th/95th/99th percentiles={row.q90_absolute_difference:.6g}/{row.q95_absolute_difference:.6g}/{row.q99_absolute_difference:.6g}; '
                     f'max={row.max_absolute_difference:.6g}; mean signed delta={row.mean_signed_difference:.6g}; '
                     f'Pearson r={row.pearson_r:.6f}, Spearman rho={row.spearman_rho:.6f}, sign agreement={row.sign_agreement_fraction:.2%}.')
    lines += ['', 'Near-threshold is prespecified as either frozen global q in [0.025,0.10]. '
              'Feature and syndrome breakdowns are provided in separate CSVs; these describe shifts rather than test a difference between correlated estimates.', '',
              f'Across {summary["photo_replicates"]} uniform one-photo-per-recorded-patient random selections, global-BH-significant counts range '
              f'{summary["random_photo_significant_min"]}-{summary["random_photo_significant_max"]}, median '
              f'{summary["random_photo_significant_median"]:g}. {summary["random_photo_always_significant"]} hypotheses are significant in every sampled selection; '
              f'{summary["random_photo_ever_significant"]} are significant at least once. '
              'The empirical g 2.5th/97.5th quantiles describe variation over available photograph choices, not population-bootstrap confidence intervals. '
              'Each replicate keeps 1,850 patients, the 49 targets, 942 reference paths and the 6,125-test BH family. '
              f'With {summary["photo_replicates"]} draws, estimated selection-frequency Monte Carlo SE can reach '
              f'{0.5 / np.sqrt(summary["photo_replicates"]):.4f}; frequency 0 or 1 is not proof of probability 0 or 1. '
              'The same observed photographs and patients are reused, so this is not independent replication.', '',
              '## Subgroup size restriction', '',
              'The restricted-family table uses the existing finite-observation sample sizes and saved raw p-values. '
              'Each threshold/analysis/layer has a separate recomputed BH family. This is a restriction and re-correction of existing fits, '
              'not a model refit; no retained effect or HC3 interval changes. Median widths can change because the retained hypotheses differ. '
              'Original p-values, old global q, and new restricted q are all retained. Sparse disease cells and high leverage can make HC3 unstable even when estimation succeeds. '
              'Analysis C interval-width medians pool feature-specific scales, including angles and normalized ratios; they are not a comparable cross-feature precision metric. Prefer within-feature intervals.', '',
              '| Analysis | Layer | Minimum disease n | Tests | Significant/new BH | Significant/old q filtered | Median interval width |',
              '|---|---|---:|---:|---:|---:|---:|']
    for row in thresholds.itertuples():
        lines.append(f'| {row.analysis} | {row.layer} | {row.minimum_disease_n} | {row.hypotheses} | '
                     f'{row.significant_recorrected} | {row.significant_original_q_filtered} | {row.median_interval_width:.6g} |')
    lines += ['', '## Changes to scientific conclusions', '',
              'Within-GMDB comparisons reduce the direct GMDB-versus-FairFace source contrast, and the much smaller discovery set shows that the '
              'chosen reference population materially affects which associations are detected. This does not identify which FairFace associations are technical artifacts: '
              'both syndrome and comparator composition change. GMDB controls have other diagnoses, are not a healthy normative population, and can share '
              'genetic/clinical features or acquisition practices. These results support differences relative to an assembled syndrome comparator only. '
              'No diagnostic performance, clinical phenotype equivalence, novelty, or external validation is established.', '',
              '## Unavailable fields and proposed future research', '',
              'No verified image-width/height, expression, occlusion, clinical gaze labels, raw historical landmarks, or explicit diagnosis-confirmation field '
              'was identified in the master tables. Iris-offset features are geometric proxies, not validated gaze classifications. '
              'Genetic and HPO fields are inventoried by nonempty count; no genetic confirmation is inferred from syndrome names, gene names, or nonempty variants. '
              'See component/imaging reports for blockers and a future authorized validation design. HPO mappings are review candidates only.', '']
    (OUT/'METHODS_AND_RESULTS_FOR_AUTHOR_REVIEW.md').write_text('\n'.join(lines))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    figures()
    reports()

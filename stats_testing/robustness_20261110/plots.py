import argparse
import hashlib
import json
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .common import ROOT, DEFAULT_OUT, case_rows, features, private_output, protect_baseline, require, sha, write_aggregate, write_json

plt.rcParams.update({'font.size': 9, 'font.family': 'DejaVu Sans', 'svg.fonttype': 'none',
                     'svg.hashsalt': 'PRISM_ROBUSTNESS_20261110', 'axes.spines.top': False, 'axes.spines.right': False})
COLORS = {'shared': '#1B766E', 'FairFace-only': '#C87832', 'within-GMDB-only': '#775B97', 'neither': '#B7C0C8'}


def save(fig, out, name, sources, caption):
    folder = out / 'figures'
    folder.mkdir(exist_ok=True, mode=0o700)
    outputs = {}
    for extension in ['svg', 'png']:
        path = folder / f'{name}.{extension}'
        fig.savefig(path, bbox_inches='tight', dpi=300, metadata={'Date': None} if extension == 'svg' else None)
        if extension == 'svg':
            path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines()) + '\n')
        outputs[path.name] = sha(path)
    plt.close(fig)
    write_json({'private_review_draft': True, 'source_sha256': {str(path.relative_to(out)): sha(path) for path in sources},
                'code_sha256': sha(Path(__file__)), 'outputs': outputs, 'caption_and_alt_text': caption}, folder / (name + '.json'))


def comparator_plots(out):
    path = out / 'evidence/P0_exhaustive_evidence_matrix.csv'
    matrix = pd.read_csv(path, low_memory=False)
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.6), sharex=True, sharey=True)
    bands = [('10–19', matrix.ff_n_disease.between(10, 19)), ('20–39', matrix.ff_n_disease.between(20, 39)), ('40+', matrix.ff_n_disease.ge(40))]
    limit = np.ceil(np.abs(matrix[['ff_hedges_g', 'gmdb_hedges_g']]).max().max())
    for ax, (label, mask) in zip(axes, bands):
        part = matrix.loc[mask]
        for category in ['neither', 'FairFace-only', 'within-GMDB-only', 'shared']:
            group = part.loc[part.category.eq(category)]
            ax.scatter(group.ff_hedges_g, group.gmdb_hedges_g, s=8, alpha=.45, color=COLORS[category], linewidths=0, label=category)
        ax.plot([-limit, limit], [-limit, limit], '--', color='#697781', lw=.8)
        ax.axhline(0, color='#AAB4BC', lw=.5)
        ax.axvline(0, color='#AAB4BC', lw=.5)
        ax.set(xlim=(-limit, limit), ylim=(-limit, limit), xlabel="Syndrome − FairFace: Hedges’ g", title=f'Real syndrome n={label}\n{len(part):,} hypotheses')
    axes[0].set_ylabel("Syndrome − other GMDB: Hedges’ g")
    axes[-1].legend(fontsize=7, loc='upper left')
    fig.tight_layout()
    save(fig, out, 'P0_comparator_scatter', [path], 'All 5,760 syndrome-feature aggregates by patient-count band. '
         'Colors show separate global-BH significance categories, dashed line equality. Cases overlap and comparator estimands differ; agreement is not independent replication. No concordance uncertainty intervals are shown.')
    categories = out / 'aggregate_candidates/P0_overlap_categories.csv'
    frame = pd.read_csv(categories)
    fig, ax = plt.subplots(figsize=(6.8, 3.6))
    ax.barh(frame.category, frame.hypotheses, color=[COLORS[label] for label in frame.category])
    for y, n in enumerate(frame.hypotheses):
        ax.text(n + 35, y, str(n), va='center')
    ax.set(xlim=(0, 5300), xlabel='Primary syndrome-feature hypotheses (total 5,760)')
    fig.tight_layout()
    save(fig, out, 'P0_comparator_categories', [categories], 'Four exhaustive categories: shared 173, FairFace-only 782, within-GMDB-only 49, neither 4,756. '
         'Primary and all-other-GMDB each use their own 5,760-hypothesis global-BH family. Counts are not false-positive estimates.')
    cases = out / 'aggregate_candidates/P0_discordant_candidates.csv'
    frame = pd.read_csv(cases).head(8).iloc[::-1]
    fig, ax = plt.subplots(figsize=(9.3, 6.4))
    y = np.arange(len(frame))
    for offset, prefix, color, lower, upper in [(-.12, 'ff_', COLORS['FairFace-only'], 'hedges_g_ci_95_low', 'hedges_g_ci_95_high'),
                                                (.12, 'gmdb_', COLORS['shared'], 'ci_95_lower', 'ci_95_upper')]:
        point = frame[prefix + 'hedges_g']
        ax.errorbar(point, y + offset, xerr=np.vstack([point - frame[prefix + lower], frame[prefix + upper] - point]),
                    fmt='o', markersize=4, color=color, capsize=2, label='FairFace' if prefix == 'ff_' else 'Other GMDB')
    ax.set_yticks(y, [textwrap.shorten(s.split('__')[0], width=42, placeholder='…') + f'\nID {identifier} · {feature}'
                      for s, identifier, feature in zip(frame.syndrome, frame.canonical_syndrome_id, frame.feature)], fontsize=8)
    ax.axvline(0, color='#AAB4BC', lw=.8)
    ax.set_xlabel("Syndrome-minus-comparator Hedges’ g")
    ax.legend(loc='best', fontsize=8)
    fig.tight_layout()
    save(fig, out, 'P0_discordant_cases', [cases], 'Eight largest absolute comparator discrepancies among primary q<.05, |g|>=.8 pairs with sign reversal or at least halved magnitude. '
         'Source-specific marginal intervals are shown for context; they are dependent and do not constitute a confidence interval for the difference. Selection is explicitly exploratory and emphasizes weakening cases.')


def source_plots(out):
    folder = out / 'aggregate_candidates'
    path = folder / 'P0_negative_control_repetitions.csv'
    frame = pd.read_csv(path)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1), sharey=True)
    for ax, (design, part) in zip(axes, frame.groupby('design', sort=True)):
        ax.hist(part.discoveries, bins=np.arange(-.5, max(frame.discoveries) + 5, 5), color='#48778C', edgecolor='white')
        ax.set(title=design.replace('_', ' '), xlabel='Discoveries per 5,760-test repetition')
    axes[0].set_ylabel('Repetitions (100 per design)')
    fig.tight_layout()
    save(fig, out, 'P0_negative_control_distribution', [path], 'Within-FairFace random label-stratified and deliberate absolute-pitch-imbalance designs. '
         'Each uses 48 disjoint pseudo-case groups, 872 paths total, a disjoint shared reference half, 120 features and global BH across 5,760 tests. '
         'The observed discovery-count distribution is workflow sensitivity, not a clinical false-positive rate; identities are not verified people.')
    path = folder / 'P0_feature_quality_correlations.csv'
    frame = pd.read_csv(path)
    frame = frame.loc[frame.recorded_group.eq('all')]
    names = ['eb_lateral_thickness_mean', 'eb_lateral_thickness_l', 'eb_lateral_thickness_r', 'forehead_height', 'chin_height']
    proxies = ['abs_yaw', 'abs_pitch', 'abs_roll', 'pose_score']
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), sharey=True)
    for ax, source in zip(axes, ['GMDB', 'FairFace']):
        values = frame.loc[frame.source.eq(source)].pivot(index='feature', columns='proxy', values='rho').loc[names, proxies]
        image = ax.imshow(values, cmap='RdBu_r', vmin=-.5, vmax=.5, aspect='auto')
        ax.set_xticks(range(4), ['|Yaw|', '|Pitch|', '|Roll|', 'Pose score'], rotation=35, ha='right')
        ax.set_yticks(range(5), [n.replace('_', ' ') for n in names])
        ax.set_title(source)
        for i in range(5):
            for j in range(4):
                ax.text(j, i, f'{values.iloc[i,j]:.2f}', ha='center', va='center', fontsize=8)
    fig.subplots_adjust(left=.32, bottom=.28, right=.88, wspace=.15)
    fig.colorbar(image, cax=fig.add_axes([.91, .28, .02, .6]), label='Spearman rho')
    save(fig, out, 'P0_recurrent_feature_quality', [path], 'All-source rank correlations between predeclared recurrent features and observed pose proxies. '
         'These are associations, not identified artifacts or biological effects. Source/group correlation families are separate; full 120-feature tables are supplied. '
         'Processed dimensions are constant and cannot identify original acquisition resolution.')
    path = folder / 'P0_acquisition_proxy_balance.csv'
    frame = pd.read_csv(path)
    groups = ['all', 'African', 'Asian', 'European', 'Others']
    fig, axes = plt.subplots(1, 3, figsize=(8.8, 3.5), sharey=True)
    for ax, proxy in zip(axes, ['pose_yaw', 'pose_pitch', 'pose_roll']):
        for source, offset, color in [('GMDB_selected', -.12, '#1B766E'), ('FairFace_reference', .12, '#C87832')]:
            part = frame.loc[frame.source.eq(source) & frame.proxy.eq(proxy)].set_index('group').reindex(groups)
            require(part[['q1', 'median', 'q3']].notna().all().all(), 'Acquisition distribution source/group is missing')
            y = np.arange(len(groups)) + offset
            ax.hlines(y, part.q1, part.q3, color=color, linewidth=2)
            ax.scatter(part['median'], y, color=color, s=12, label=source.split('_')[0])
        ax.axvline(0, color='#BDC7CC', lw=.7)
        ax.set(xlabel=proxy.replace('_', ' '), yticks=range(len(groups)), yticklabels=groups)
    axes[-1].legend(fontsize=8)
    axes[0].invert_yaxis()
    fig.tight_layout()
    save(fig, out, 'P0_acquisition_proxy_distributions', [path], 'Recorded signed-pose medians and interquartile ranges by source and recorded broad group. '
         'Bars describe distributions, not confidence intervals. Exact group counts and missing-label summaries are in the source table. '
         'Recorded groups are dataset annotations, not biological ancestry. No clinical source/disease effect is identified.')


def dependence_plots(out):
    folder = out / 'aggregate_candidates'
    path = folder / 'P1_source_correlations.csv'
    frame = pd.read_csv(path)
    names = pd.read_csv(ROOT / 'publication/feature_definitions/current_120_feature_formulas.csv').feature.tolist()
    fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.5))
    matrices = []
    for source in ['FairFace', 'GMDB_selected']:
        matrix = np.eye(120)
        for row in frame.loc[frame.source.eq(source)].itertuples():
            i, j = names.index(row.feature_a), names.index(row.feature_b)
            matrix[i, j] = matrix[j, i] = row.spearman
        matrices.append(matrix)
    for ax, values, title in zip(axes, [*matrices, matrices[1] - matrices[0]], ['FairFace (n=3,409)', 'Selected GMDB (n=1,845)', 'GMDB − FairFace']):
        image = ax.imshow(values, cmap='RdBu_r', vmin=-1, vmax=1)
        ax.set_title(title, fontsize=9)
        ax.set_xticks([0, 29, 59, 89, 119], [1, 30, 60, 90, 120])
        ax.set_yticks([0, 29, 59, 89, 119], [1, 30, 60, 90, 120])
        ax.set_xlabel('Ordered S6 feature index')
    fig.subplots_adjust(left=.05, right=.87, bottom=.18, wspace=.3)
    fig.colorbar(image, cax=fig.add_axes([.9, .18, .015, .67]), label='Spearman rho / difference')
    save(fig, out, 'P1_source_correlation_heatmap', [path], 'Separate 120-by-120 Spearman matrices in canonical S6 order and their difference. '
         'Mixture, population, sample size and acquisition can all change correlations. Clusters do not estimate independent anatomical dimensions; feature IDs and assignments are in the source tables.')
    path = folder / 'P1_feature_dependency_edges.csv'
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.5))
    for ax in axes:
        ax.axis('off')
        ax.set(xlim=(0, 1), ylim=(0, 1))
    ax = axes[0]
    positions = {'right area': (.15, .8), 'left area': (.85, .8), 'bilateral mean': (.5, .47), 'relative asymmetry': (.5, .12)}
    for label, (x, y) in positions.items():
        ax.text(x, y, label, ha='center', va='center', bbox={'boxstyle': 'round,pad=.4', 'facecolor': '#E5EDF0', 'edgecolor': '#657984'}, fontsize=8)
    for a, b in [('right area', 'bilateral mean'), ('left area', 'bilateral mean'), ('right area', 'relative asymmetry'), ('left area', 'relative asymmetry'), ('bilateral mean', 'relative asymmetry')]:
        ax.annotate('', xy=positions[b], xytext=positions[a], arrowprops={'arrowstyle': '->', 'color': '#48778C', 'shrinkA': 22, 'shrinkB': 22})
    ax.set_title('A  Explicit algebraic dependence', fontsize=10)
    ax = axes[1]
    ax.text(.5, .82, 'Bizygomatic width', ha='center', bbox={'boxstyle': 'round,pad=.4', 'facecolor': '#FAE6CE', 'edgecolor': '#A37E54'})
    for x, label in zip([.12, .5, .88], ['Eye span', 'Chin height', 'Forehead height']):
        ax.text(x, .25, label, ha='center', va='center', rotation=25, fontsize=8)
        ax.annotate('', xy=(x, .35), xytext=(.5, .73), arrowprops={'arrowstyle': '->', 'linestyle': '--', 'color': '#A37E54'})
    ax.set_title('B  Shared denominator', fontsize=10)
    fig.text(.5, .02, 'Shared normalization can couple indices; it is not an algebraic identity between outputs.', ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .08, 1, 1))
    save(fig, out, 'P1_feature_dependency_network', [path], 'Original schematic of exact eyebrow right/left/mean/asymmetry dependence and shared-normalizer connections. '
         'The complete edge CSV is parsed from exact S6 formulas. Shared denominators are labeled separately from explicit output-to-output algebraic references; no independence count is inferred.')
    path = folder / 'P1_signature_source.csv'
    frame = pd.read_csv(path, low_memory=False)
    fig, axes = plt.subplots(2, 2, figsize=(8, 6.5), sharex=True)
    for ax, (syndrome, part) in zip(axes.flat, frame.groupby('syndrome', sort=True)):
        part = part.assign(magnitude=part.ff_hedges_g.abs()).nlargest(5, 'magnitude').iloc[::-1]
        y = np.arange(len(part))
        ax.hlines(y, part.ff_hedges_g_ci_95_low, part.ff_hedges_g_ci_95_high, color='#48778C')
        ax.plot(part.ff_hedges_g, y, 'o', color='#48778C', markersize=4)
        ax.set_yticks(y, [textwrap.fill(f.replace('_', ' '), 22) for f in part.feature], fontsize=7)
        ax.axvline(0, color='#C0C8CD', lw=.7)
        ax.set_title(syndrome.split('__')[0].split(';')[0], fontsize=9)
        ax.set_xlabel("Syndrome − FairFace Hedges’ g")
    fig.tight_layout(w_pad=2, h_pad=1.5)
    save(fig, out, 'P1_reference_cluster_signatures', [path], 'Four predeclared n>=20 syndrome cohorts. Reference-only representatives are selected by formula simplicity without outcome selection. '
         'The five largest |g| representatives per panel are displayed for readability, an outcome-selected descriptive illustration. Intervals are frozen marginal primary intervals; all representatives and nonsignificant cells remain in the source.')


def stability_and_geometry_plots(out):
    folder = out / 'aggregate_candidates'
    path = out / 'evidence/P0_exhaustive_evidence_matrix.csv'
    frame = case_rows(pd.read_csv(path, low_memory=False))
    fig, axes = plt.subplots(3, 2, figsize=(8.4, 9))
    specs = [('ff_', 'Primary FF', 'hedges_g_ci_95_low', 'hedges_g_ci_95_high'),
             ('gmdb_', 'Other GMDB', 'ci_95_lower', 'ci_95_upper'),
             ('alternate_', 'Alternate photo', 'hedges_g_ci_95_low', 'hedges_g_ci_95_high'),
             ('all_reference_', 'All eligible FF', 'hedges_g_ci_95_low', 'hedges_g_ci_95_high'),
             ('sex_reference_', 'Category × sex FF', 'hedges_g_ci_95_low', 'hedges_g_ci_95_high')]
    for ax, row in zip(axes.flat, frame.itertuples()):
        for y, (prefix, label, lower, upper) in enumerate(specs):
            ax.hlines(y, getattr(row, prefix + lower), getattr(row, prefix + upper), color='#5A8192')
            ax.plot(getattr(row, prefix + 'hedges_g'), y, 'o', markersize=4, color='#C47836' if prefix == 'gmdb_' else '#48778C')
        ax.set_yticks(range(5), [x[1] for x in specs], fontsize=7)
        ax.axvline(0, color='#B5C0C7', lw=.8)
        ax.set_title(row.syndrome.split('__')[0].split(';')[0] + '\n' + row.feature.replace('_', ' '), fontsize=9)
        ax.set_xlabel("Hedges’ g, comparator-specific estimand")
        ax.invert_yaxis()
    fig.tight_layout(h_pad=1.8, w_pad=2)
    save(fig, out, 'P1_predeclared_case_specifications', [path], 'Six predeclared handout cases across frozen primary/reference/photo and within-GMDB contrasts. '
         'Orange marks other-GMDB effects. Each specification keeps its own original 5,760-test BH family; intervals are source-specific marginal intervals, not joint difference intervals. '
         'The Rubinstein–Taybi malar and Cornelia de Lange lip-peak-width contrasts weaken against other GMDB conditions.')
    path = folder / 'P1_normalization_cases.csv'
    frame = pd.read_csv(path)
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.7), sharex=True)
    for ax, (feature, part) in zip(axes, frame.groupby('feature', sort=True)):
        definitions = part.definition.drop_duplicates().tolist()
        for offset, comparator, color in [(-.12, 'FairFace', '#48778C'), (.12, 'other_GMDB', '#C47836')]:
            group = part.loc[part.comparator.eq(comparator)].set_index('definition').loc[definitions]
            y = np.arange(3) + offset
            ax.hlines(y, group.case_resampling_q025, group.case_resampling_q975, color=color)
            ax.plot(group.hedges_g, y, 'o', color=color, markersize=4, label=comparator)
        ax.set_yticks(range(3), [textwrap.fill(s.replace('_', ' '), 18) for s in definitions], fontsize=8)
        ax.axvline(0, color='#B5C0C7', lw=.8)
        ax.set(title=feature.replace('_', ' '), xlabel="Hedges’ g for the stated ratio")
    axes[-1].legend(fontsize=8)
    fig.tight_layout(w_pad=1.5)
    save(fig, out, 'P1_normalization_sensitivity', [path], 'Original and two legitimate alternate ratios for each predeclared case. '
         'Bars are conditional 200-case-patient 2.5–97.5 percentile ranges with fixed comparators, not calibrated clinical confidence intervals. '
         'Alternate ratios change the estimand; a changed point sign is not proof of contradictory anatomy. New p/q values belong only to the declared 12-cell exploratory family.')
    path = folder / 'P1_ratio_geometry_simulation.csv'
    frame = pd.read_csv(path)
    fig, ax = plt.subplots(figsize=(6.3, 3.7))
    for change, part in frame.groupby('simulated_numerator_change', sort=True):
        ax.plot(100 * part.simulated_denominator_change, 100 * part.relative_ratio_change, label=f'Numerator change {change:+.0%}')
    ax.axhline(0, color='#B5C0C7', lw=.8)
    ax.axvline(0, color='#B5C0C7', lw=.8)
    ax.set(xlabel='Simulated denominator change (%)', ylabel='Relative ratio change (%)', title='Numerical normalization mechanism; no patient data')
    ax.legend(fontsize=8)
    fig.tight_layout()
    save(fig, out, 'P1_ratio_simulation', [path], 'Numerical geometry simulation of numerator/denominator changes. '
         'Holding the numerator fixed while changing the denominator shifts the ratio. Curves illustrate algebra only and do not attribute an observed disease effect to a particular component.')
    path = folder / 'P1_bootstrap_coverage_simulation.csv'
    frame = pd.read_csv(path)
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    y = np.arange(len(frame))
    ax.hlines(y, frame.coverage_mc_wilson_low, frame.coverage_mc_wilson_high, color='#48778C')
    ax.plot(frame.empirical_coverage, y, 'o', color='#48778C')
    ax.axvline(.95, ls='--', color='#C47836', label='Nominal .95')
    ax.set_yticks(y, [f'{r.distribution}, case n={r.case_n}' for r in frame.itertuples()])
    ax.set(xlim=(.82, .98), xlabel='Coverage across 300 independent numerical datasets per regime')
    ax.legend(fontsize=8)
    fig.tight_layout()
    save(fig, out, 'P1_conditional_bootstrap_coverage', [path], 'Coverage diagnostics for the new case-only, fixed-reference percentile procedure under normal/scaled-t5 simulations. '
         'Intervals describe Monte Carlo uncertainty across independent datasets. Undercoverage at small n motivates descriptive labeling of empirical ranges; '
         'this is not a replay or automatic invalidation of the frozen two-sample primary bootstrap.')


def candidate_plot(out):
    path = out / 'aggregate_candidates/P1_ranking_rule_overlap.csv'
    frame = pd.read_csv(path)
    labels = ['high_magnitude', 'high_stability', 'conservative_cross_comparator', 'unverified_geometry_hypotheses', 'primary_q_baseline', 'adjusted_q_baseline']
    matrix = frame.pivot(index='first_rule', columns='second_rule', values='shared_top20').loc[labels, labels]
    fig, ax = plt.subplots(figsize=(6.8, 5.2))
    image = ax.imshow(matrix, vmin=0, vmax=20, cmap='Blues')
    names = ['Magnitude', 'Stability', 'Cross-comparator', 'Unmapped geometry', 'Primary q', 'Adjusted q']
    ax.set_xticks(range(6), names, rotation=35, ha='right')
    ax.set_yticks(range(6), names)
    for i in range(6):
        for j in range(6):
            ax.text(j, i, str(matrix.iloc[i, j]), ha='center', va='center', color='white' if matrix.iloc[i,j] > 12 else '#1F3541')
    fig.colorbar(image, ax=ax, label='Shared top-20 pairs')
    fig.tight_layout()
    save(fig, out, 'P1_candidate_rule_sensitivity', [path], 'Overlap of transparent exploratory top-20 ranking rules. '
         'Rules use the same patients and data; overlap does not validate a method or estimate held-out performance. '
         'Unmapped clinical correspondence is not evidence of novelty, and unsupported adjustments are not proven nulls.')


def geometry_diagram(out):
    path = out / 'aggregate_candidates/P1_geometry_diagram_definitions.csv'
    definitions = features().set_index('feature')
    frame = definitions.loc[['outer_canthal_distance', 'philtrum_width', 'face_aspect_ratio']].reset_index()
    frame['illustration_scope'] = 'original invented coordinates; no real photograph, landmark record or patient-derived shape'
    write_aggregate(frame, path)
    fig, axes = plt.subplots(1, 2, figsize=(8.3, 4.3))
    for ax, feature, title in zip(axes, ['outer_canthal_distance', 'philtrum_width'], ['A  Outer-eye-corner span', 'B  Upper-lip-peak span']):
        t = np.linspace(0, 2 * np.pi, 160)
        ax.plot(.5 + .33 * np.cos(t), .55 + .4 * np.sin(t), color='#BDC7CC', lw=1.2)
        ax.plot([.23, .3, .37], [.66, .69, .66], color='#81939C')
        ax.plot([.63, .7, .77], [.66, .69, .66], color='#81939C')
        ax.plot([.5, .46, .54], [.6, .45, .45], color='#81939C')
        ax.plot([.38, .45, .5, .55, .62], [.34, .38, .35, .38, .34], color='#81939C')
        ax.plot([.17, .83], [.54, .54], '--', color='#A67A4B', lw=1)
        ax.scatter([.17, .83], [.54, .54], s=12, color='#A67A4B')
        ax.text(.5, .51, 'B: bizygomatic width', ha='center', va='top', fontsize=8, color='#855C31')
        a, b, y = (.23, .77, .66) if feature == 'outer_canthal_distance' else (.45, .55, .38)
        ax.scatter([a, b], [y, y], s=20, color='#1B766E', zorder=4)
        ax.annotate('', xy=(a, y + .035), xytext=(b, y + .035), arrowprops={'arrowstyle': '<->', 'color': '#1B766E'})
        ax.text(.5, y + .075 if feature == 'outer_canthal_distance' else .25, 'N: numerator distance', ha='center', color='#1B766E', fontsize=8)
        ax.text(.5, .06, 'Original index = N / (B + ε)', ha='center', fontsize=10)
        ax.set(title=title, xlim=(0, 1), ylim=(0, 1.03), aspect='equal')
        ax.axis('off')
    fig.text(.5, .075, 'Height alternative ≈ (N / B) × (B / H) = N / H; each denominator changes the estimand.', ha='center', fontsize=9)
    fig.text(.5, .025, 'Synthetic schematic, not anatomical calibration. Lip peaks do not establish philtral-column width.', ha='center', fontsize=8)
    fig.subplots_adjust(bottom=.14, top=.94, left=.02, right=.98, wspace=.1)
    save(fig, out, 'P1_geometry_landmark_schematic', [path], 'Original synthetic schematic with invented coordinates for outer-eye-corner and upper-lip-peak numerators '
         'and bizygomatic denominator. Positions illustrate the construct only and are not exported MediaPipe coordinates. '
         'Multiplying by the pinned face-aspect ratio gives the height alternative up to the extractor epsilon convention. '
         'No patient photograph, saved landmark or real shape is represented; no millimeter inference is possible.')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    out = private_output(args.output_dir)
    protect_baseline(out)
    comparator_plots(out)
    if (out / 'aggregate_candidates/P0_negative_control_summary.csv').exists():
        source_plots(out)
        dependence_plots(out)
        stability_and_geometry_plots(out)
        candidate_plot(out)
        geometry_diagram(out)
    print('Created private aggregate research figures; no public figure paths changed')


if __name__ == '__main__':
    main()

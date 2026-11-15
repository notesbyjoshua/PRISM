import json

import numpy as np
import pandas as pd

from core import KEYS, RACES, bh, empirical_p, write_json


def table(frame):
    lines = ['| ' + ' | '.join(frame.columns) + ' |', '| ' + ' | '.join(['---'] * len(frame.columns)) + ' |']
    for row in frame.itertuples(index=False, name=None):
        cells = [f'{v:.4g}' if isinstance(v, float) else str(v) for v in row]
        lines.append('| ' + ' | '.join(v.replace('|', '\\|') for v in cells) + ' |')
    return '\n'.join(lines)


def count_distribution(values):
    values = np.asarray(values)
    return {'mean': float(values.mean()), 'median': float(np.median(values)), 'q95': float(np.quantile(values, .95)),
            'minimum': int(values.min()), 'maximum': int(values.max()), 'fraction_at_least_one': float((values > 0).mean())}


def figures(study, plan, output, placebo, reference, repetitions, permutation):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    observed = int(study.primary.q.lt(.05).sum())
    caption = f'T={study.total:,}; real GMDB eligible n={int(study.cells.n_current_gmdb.sum())}; locked FF n={study.sampling["X"]}'
    if placebo is not None:
        fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')
        values = placebo.discoveries.to_numpy()
        axes[0].hist(values, bins=np.arange(values.max() + 2) - .5, color='#28748c')
        axes[0].set(xlabel='BH discoveries per placebo repetition', ylabel='Repetitions')
        ordered = np.sort(values)
        axes[1].step(ordered, np.arange(1, len(values) + 1) / len(values), where='post', color='#28748c')
        axes[1].set(xlabel='BH discoveries', ylabel='ECDF')
        fig.suptitle(f'FairFace-only placebo: B={len(values)}, seed={plan["analysis_seeds"]["placebo"]}\n'
                     f'Main observed={observed} (different-source context, not a permutation p-value)')
        fig.supxlabel(caption + '; European placebo unsupported; full BH denominator retained', fontsize=8)
        fig.savefig(output / 'fairface_placebo_discovery_distribution.png', dpi=180)
        plt.close(fig)
    if reference is not None:
        fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')
        for ax, column, title in zip(axes, ['discoveries', 'overlap_with_main'], ['Discoveries', 'Overlap with main discoveries']):
            ax.hist(repetitions[column], bins=25, color='#28748c')
            ax.set(xlabel=title, ylabel='Reference draws')
        fig.suptitle(f'Reference selection stability: B={len(repetitions)}, seed={plan["analysis_seeds"]["reference-resampling"]}')
        fig.supxlabel(caption + '; shared overlapping draws; not a null or independent replication', fontsize=8)
        fig.savefig(output / 'reference_resampling_discovery_overlap.png', dpi=180)
        plt.close(fig)
        fig, ax = plt.subplots(figsize=(7, 5), layout='constrained')
        for race in RACES:
            part = reference.loc[reference.canonical_race.eq(race)]
            ax.scatter(part.sign_agreement_frequency, part.significant_frequency, s=7, alpha=.35, label=race)
        ax.axhline(.8, color='grey', linestyle='--')
        ax.axvline(.95, color='grey', linestyle='--')
        ax.set(xlabel='Agreement with main effect sign', ylabel='Fraction with global BH q<.05', xlim=(-.02, 1.02), ylim=(-.02, 1.02))
        ax.legend(fontsize=8)
        ax.set_title('Descriptive reference stability; thresholds prespecified')
        fig.supxlabel(caption + f'; B={len(repetitions)}; reference seed={plan["analysis_seeds"]["reference-resampling"]}', fontsize=8)
        fig.savefig(output / 'reference_effect_sign_stability.png', dpi=180)
        plt.close(fig)
        fig, ax = plt.subplots(figsize=(6, 5), layout='constrained')
        for race in RACES:
            part = reference.loc[reference.canonical_race.eq(race)]
            ax.scatter(part.observed_main_g, part.median_g, s=6, alpha=.35, label=race)
        low, high = np.nanmin(reference[['observed_main_g', 'median_g']]), np.nanmax(reference[['observed_main_g', 'median_g']])
        ax.plot([low, high], [low, high], '--', color='grey')
        ax.set(xlabel='Original main Hedges g', ylabel='Median g across reference draws', title='Selection quantiles are not patient confidence intervals')
        ax.legend(fontsize=8)
        fig.supxlabel(caption + f'; B={len(repetitions)}; reference seed={plan["analysis_seeds"]["reference-resampling"]}', fontsize=8)
        fig.savefig(output / 'reference_observed_vs_median_effect.png', dpi=180)
        plt.close(fig)
    if permutation is not None:
        names = permutation.null_specification.unique()
        fig, axes = plt.subplots(1, len(names), figsize=(6 * len(names), 5), squeeze=False, layout='constrained')
        for ax, name in zip(axes.flat, names):
            part = permutation.loc[permutation.null_specification.eq(name)]
            ax.scatter(part.main_g, part.observed_hedges_g, c=np.where(part.q_global_BH_empirical.lt(.05), '#c05a43', '#28748c'),
                       s=5, alpha=.4)
            ax.axhline(0, color='grey', linewidth=.7)
            ax.axvline(0, color='grey', linewidth=.7)
            ax.set(xlabel='Main syndrome minus FairFace g', ylabel='Syndrome minus other GMDB g', title=name)
        fig.suptitle(f'Different comparators, shared cases; orange: empirical BH q<.05\n'
                     f'B={int(permutation.B.iloc[0]):,}; permutation seed={plan["analysis_seeds"]["permutation"]}')
        fig.supxlabel(caption + '; other-GMDB pool includes eligible and ineligible syndromes; not external validation', fontsize=8)
        fig.savefig(output / 'gmdb_permutation_vs_main_effect.png', dpi=180)
        plt.close(fig)


def report(study, plan, output):
    def load(name):
        path = output / name
        return pd.read_csv(path, float_precision='round_trip') if path.exists() else None

    placebo = load('fairface_placebo_summary.csv')
    reference = load('reference_resampling_aggregate.csv')
    repetitions = load('reference_resampling_repetitions.csv')
    permutation = load('gmdb_permutation_aggregate.csv')
    observed = int(study.primary.q.lt(.05).sum())
    summary = {'status': 'PARTIAL', 'planned_family': study.total, 'observed_main_discoveries': observed,
               'baseline': 'PRISM_RACE_STRATIFIED_20261115', 'public_export': 'NOT RUN; private review candidates only'}
    if placebo is not None:
        summary['placebo'] = {'B': len(placebo), **count_distribution(placebo.discoveries),
                               'supported_hypotheses': int(placebo.supported_hypotheses.iloc[0]),
                               'unsupported_hypotheses': study.total - int(placebo.supported_hypotheses.iloc[0]),
                               'primary_exceedance_fraction_descriptive_only': float(placebo.discoveries.ge(observed).mean())}
    if reference is not None:
        summary['reference_resampling'] = {'B': len(repetitions), **count_distribution(repetitions.discoveries),
                                            'overlap_with_main': count_distribution(repetitions.overlap_with_main),
                                            'stable_main_associations': int((reference.observed_main_q.lt(.05) & reference.descriptive_stable).sum()),
                                            'unstable_main_associations': int((reference.observed_main_q.lt(.05) & ~reference.descriptive_stable).sum()),
                                            'any_sign_changes': int(reference.sign_agreement_frequency.lt(1).sum())}
    if permutation is not None:
        summary['gmdb_permutation'] = {}
        for name, part in permutation.groupby('null_specification'):
            np.testing.assert_allclose(part.p_empirical_permutation, empirical_p(part.exceedances, int(part.B.iloc[0])))
            np.testing.assert_allclose(part.q_global_BH_empirical, bh(part.p_empirical_permutation), rtol=1e-12)
            sig = part.q_global_BH_empirical.lt(.05)
            main = part.main_q.lt(.05)
            summary['gmdb_permutation'][name] = {'B': int(part.B.iloc[0]), 'p_floor': float(part.p_value_floor.iloc[0]),
                'discoveries': int(sig.sum()), 'overlap_main': int((sig & main).sum()),
                'overlap_main_same_direction': int((sig & main & part.same_effect_direction_as_main).sum()),
                'MC_decision_ambiguous': int(part.MC_FDR_decision_ambiguous.sum()),
                'MC_upper_q_below_005': int(part.q_MC_simultaneous_upper.lt(.05).sum())}
    if placebo is not None and reference is not None and permutation is not None:
        complete = len(placebo) >= 1000 and len(repetitions) >= 1000 and permutation.B.min() >= 100000
        summary['status'] = 'COMPLETED' if complete else 'EXPLORATION_COMPLETE'
    figures(study, plan, output, placebo, reference, repetitions, permutation)
    robust = None
    if reference is not None and permutation is not None:
        robust = reference.loc[reference.observed_main_q.lt(.05) & reference.descriptive_stable].copy()
        for name, part in permutation.groupby('null_specification'):
            columns = KEYS + ['q_global_BH_empirical', 'same_effect_direction_as_main', 'MC_FDR_decision_ambiguous']
            selected = part[columns].rename(columns={key: f'{name}_{key}' for key in columns if key not in KEYS})
            robust = robust.merge(selected, on=KEYS, validate='one_to_one')
            robust = robust.loc[robust[f'{name}_q_global_BH_empirical'].lt(.05) & robust[f'{name}_same_effect_direction_as_main']]
        q_columns = [column for column in robust if column.endswith('_q_global_BH_empirical')]
        robust['worst_empirical_q'] = robust[q_columns].max(axis=1)
        robust = robust.sort_values(['worst_empirical_q', 'observed_main_q'] + KEYS)
        robust.to_csv(output / 'cross_design_robust_candidates.csv', index=False)
        summary['descriptive_cross_design_robust_associations'] = len(robust)
    write_json(output / 'overall_summary.json', summary)
    lines = ['# Three-category PRISM sensitivities: results and limitations', '',
             f'Status: **{summary["status"]}**. The locked baseline has {study.total:,} hypotheses, '
             f'{len(study.cells)} eligible cells, {study.cells.disease.nunique()} syndromes, '
             f'{study.cells.n_current_gmdb.sum()} patients in eligible cells and {observed} discoveries. '
             'No ML, image extraction, clinical/HPO interpretation, manuscript edits or synthetic PDIDB analysis was performed.', '',
             '## Locked design and cohort audit', '',
             'The baseline is the completed November 15 race-stratified study, not the historical pooled 48-syndrome family. '
             'Exactly 120 ordered finite features, unique recorded GMDB patients and selected images, pinned extraction hashes, '
             'same-category matching and the fixed eligible-cell list are verified against current private inputs. '
             'All 1,429 mapped GMDB patients are available to the within-GMDB comparator, including patients in ineligible cells. '
             'Unknown/Others GMDB annotations and unsupported FairFace annotations are excluded. The original conditional ages '
             'and QC are retained. The broad Asian mapping still includes Middle Eastern/West Asian; these are recorded labels, '
             'not verified genetic ancestry.', '',
             table(pd.DataFrame([{'race': race, **values} for race, values in plan['races'].items()])), '',
             f'The reference remains X={study.sampling["X"]} with quotas {json.dumps(study.sampling["quotas"])}. '
             'The accepted FairFace content-QC exclusions and source crosswalk match the main study exactly. '
             'Every A/C repetition corrects the full planned family. Each GMDB empirical null has a separate global BH family. '
             'Unsupported comparisons retain p=1; no family is reduced after inspecting results.', '',
             '## A. FairFace-only negative controls', '']
    if placebo is None:
        lines += ['NOT RUN in this invocation.']
    else:
        a = summary['placebo']
        lines += [f'Executed B={a["B"]:,}; seed={plan["analysis_seeds"]["placebo"]}. Mean discoveries={a["mean"]:.4g}, '
                  f'median={a["median"]:.4g}, 95th percentile={a["q95"]:.4g}, maximum={a["maximum"]}; '
                  f'fraction with at least one={a["fraction_at_least_one"]:.4g}. Main observed count={observed}.', '',
                  f'{a["supported_hypotheses"]:,} hypotheses have disjoint same-source case/reference support; '
                  f'{a["unsupported_hypotheses"]:,} retain p=1. African and Asian groups use one shared reference and '
                  'a disjoint multiclass pseudo-case pool per repetition. European capacity equals the locked reference '
                  'size, so all European placebo contrasts are impossible even under the per-pair fallback. '
                  'No replacement, overlap or shrinking is used. This substantially limits whole-study calibration.', '',
                  'This is a same-source placebo experiment. The descriptive fraction of placebo counts exceeding the '
                  'real main count is not a valid p-value for syndrome effects across GMDB/FairFace sources. '
                  'Maximum effect magnitudes can be noisy with five artificial cases.', '',
                  f'![Placebo discoveries]({output / "fairface_placebo_discovery_distribution.png"})']
    lines += ['', '## B. Within-GMDB syndrome-label permutations', '']
    if permutation is None:
        lines += ['NOT RUN in this invocation.']
    else:
        lines += [table(pd.DataFrame([{'null': name, **values} for name, values in summary['gmdb_permutation'].items()])), '',
                  'Each repetition shuffles one multiclass syndrome assignment across whole patients within a race, '
                  'shared across every syndrome/feature, preserving all syndrome counts. Cases are compared with all other '
                  'GMDB syndrome patients in that race. The unblocked statistic is |U−n_case*n_other/2|. '
                  'Empirical p=(1+exceedances)/(B+1), including ties; no zero p-values are reported. '
                  'The fixed 1,000→10,000→100,000 ladder is planned before outcomes, without outcome-dependent stopping. '
                  'The exact achieved B and p-value floor are in the table and each output row.', '',
                  f'Age/sex imbalance triggered the blocked alternative: {sum(item["substantial_imbalance"] for item in study.covariates)} '
                  'eligible cells met prespecified age |SMD|≥.5 or sex fraction gap≥.20. Blocks are recorded sex, '
                  'including Unknown, crossed with ages [3,5), [5,7), [7,10). Only multiclass blocks shuffle. '
                  'Blocked statistics are centered at the exact conditional expected U from fixed per-block syndrome counts '
                  'and mean feature ranks. Immobile counts are in cohort_audit.json. These are different null assumptions.', '',
                  'q_MC_simultaneous_lower/upper use exact binomial tail-probability bounds with Bonferroni coverage across '
                  'the planned Monte Carlo family, then BH; they assess Monte Carlo decision sensitivity, not clinical '
                  'effect uncertainty. Point-estimate FDR decisions flagged ambiguous require more Monte Carlo precision. '
                  'Blocked labels may still be nonexchangeable because collection conditions and unmeasured covariates '
                  'can be tied to diagnosis. Neither comparator establishes external validation.', '',
                  f'![GMDB and main effects]({output / "gmdb_permutation_vs_main_effect.png"})']
    lines += ['', '## C. Reference-resampling stability', '']
    if reference is None:
        lines += ['NOT RUN in this invocation.']
    else:
        c = summary['reference_resampling']
        lines += [f'Executed B={c["B"]:,}; seed={plan["analysis_seeds"]["reference-resampling"]}. Mean discoveries={c["mean"]:.4g}, '
                  f'median={c["median"]:.4g}, 95th percentile={c["q95"]:.4g}, range={c["minimum"]}–{c["maximum"]}. '
                  f'Mean overlap with original discoveries={c["overlap_with_main"]["mean"]:.4g}. '
                  f'{c["stable_main_associations"]} original discoveries meet the prespecified descriptive stability rule; '
                  f'{c["unstable_main_associations"]} do not. {c["any_sign_changes"]} planned associations change effect sign in at least one draw.', '',
                  'Stability requires global BH q<.05 in at least 80% of draws and at least 95% agreement with the original '
                  'g sign. The European reference includes its entire pool, so its sampled photographs and g values cannot '
                  'vary. European q-values can still change through the global BH ranks of other groups. '
                  'Consequently, apparent European stability does not test sensitivity to unobserved European controls.', '',
                  'The 2.5th/97.5th g quantiles describe selection variation conditional on fixed GMDB cases and the '
                  'available FairFace pool. They are not patient-sampling confidence intervals. Shared overlapping reference '
                  'draws are neither independent biological replications nor a permutation p-value.', '',
                  f'![Reference discoveries and overlap]({output / "reference_resampling_discovery_overlap.png"})', '',
                  f'![Effect sign stability]({output / "reference_effect_sign_stability.png"})', '',
                  f'![Observed and median effects]({output / "reference_observed_vs_median_effect.png"})']
    lines += ['', '## Descriptive cross-design candidates', '']
    if robust is None:
        lines += ['NOT APPLICABLE until both reference resampling and GMDB nulls are executed.']
    elif robust.empty:
        lines += ['No association meets all prespecified cross-design criteria. No ten-association list is manufactured.']
    else:
        columns = KEYS + ['n_case', 'observed_main_g', 'observed_main_q', 'significant_frequency', 'sign_agreement_frequency', 'worst_empirical_q']
        lines += [f'{len(robust)} associations have main q<.05, descriptive reference stability, empirical q<.05 in '
                  'each enabled GMDB null and concordant g direction. The ten strongest by worst empirical q, then main q, '
                  'are shown below. These are descriptive cross-design candidates, not independent replications.', '',
                  table(robust[columns].head(10)), '',
                  'Monte Carlo ambiguity flags for each null are retained in cross_design_robust_candidates.csv.']
    lines += ['', '## Reproducibility, replacement and limits', '',
              'SENSITIVITY_ANALYSIS.json and the former stats_testing/reference_sensitivity.py entrypoint now select this '
              'three-experiment workflow. Historical pooled primary/sensitivity numerical artifacts and publication/manuscript '
              'files remain preserved. ACTIVE_ANALYSIS.json remains unchanged. Existing untracked November 10 work is untouched.', '',
              f'Input/source/code hashes, 120-feature order, Git revision, versions, seeds, nulls, capacity checks and thresholds '
              f'are in analysis_plan.json. Execution requests and elapsed times are in execution_commands.json; '
              'checkpoints record contiguous completed repetitions, signatures and checksums. SeedSequence derivation permits '
              'resume and extension without duplicate iterations or scheduling dependence. No new dependencies were installed.', '',
              'All candidate tables, figures, source identifiers, intermediate arrays and reports remain in verified Git-ignored '
              'private storage. No raw records or small cells were published. Cohort identity is based on recorded patient keys '
              'and image-content checks, without comprehensive person linkage. Source acquisition, hidden covariates, broad '
              'category correspondence, conditional age semantics and small cells remain substantive limitations.', '',
              'The data dictionary is DATA_DICTIONARY.md beside this report. Reproduction uses the commands in the code README '
              'and execution_commands.json. Synthetic fixtures verify null/shift behavior and optimized numerical utilities; '
              'test and independent validation outputs are saved privately. No clinical, HPO or diagnostic utility claims are made.']
    (output / 'RESULTS_AND_LIMITATIONS.md').write_text('\n'.join(lines) + '\n')
    return summary

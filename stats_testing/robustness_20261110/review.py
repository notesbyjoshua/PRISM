import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

from stats_testing.cohort_audit import age_diagnostics, clean_id
from .common import (ROOT, DEFAULT_OUT, KEY, SEED, bh, canonicalize, features, load_private, offline,
                     private_output, protect_baseline, require, write_aggregate, write_json)


def photo_change_details(real, out):
    config = json.loads((ROOT / 'ACTIVE_ANALYSIS.json').read_text())
    raw = pd.read_csv(ROOT / config['private_run'] / 'COHORTS/all_gmdb_candidates_master.csv',
                      dtype={'patient_id': 'string', 'image_id': 'string'})
    age = age_diagnostics(raw).conditional_age_years
    eligible = raw.loc[raw.frontal_ok.eq(True) & age.ge(3) & age.lt(10)].copy()
    eligible['patient_id'] = clean_id(eligible.patient_id)
    groups = list(eligible.groupby('patient_id', sort=True).indices.values())
    targets = set(real.disease.value_counts().loc[lambda n: n >= 10].index)
    baseline = dict(zip(real.patient_id.astype(str), real.image_id.astype(str)))
    keys = eligible.patient_id.astype(str).to_numpy()
    images, labels = eligible.image_id.astype(str).to_numpy(), eligible.disease.to_numpy()
    rng = np.random.default_rng(20261025)
    rows, by_syndrome = [], []
    for rep in range(100):
        chosen = np.array([indices[rng.integers(len(indices))] for indices in groups])
        changed = np.array([images[i] != baseline[keys[i]] for i in chosen])
        primary = np.isin(labels[chosen], list(targets))
        rows.append(dict(repetition=rep, changed_selected_photos=int(changed.sum()), changed_primary_photos=int((changed & primary).sum()),
                         primary_patients=872, selected_patients=1845))
        for syndrome in sorted(targets):
            mask = labels[chosen] == syndrome
            by_syndrome.append(dict(syndrome=syndrome, repetition=rep, patients=int(mask.sum()), changed=int(changed[mask].sum())))
    write_aggregate(pd.DataFrame(rows), out / 'aggregate_candidates/P1_random_photo_change_counts.csv')
    detail = canonicalize(pd.DataFrame(by_syndrome))
    summary = detail.groupby(['canonical_syndrome_id', 'syndrome'], sort=True).agg(patients=('patients', 'first'),
        mean_changed=('changed', 'mean'), median_changed=('changed', 'median'), q025_changed=('changed', lambda x: x.quantile(.025)),
        q975_changed=('changed', lambda x: x.quantile(.975))).reset_index()
    write_aggregate(summary, out / 'aggregate_candidates/P1_random_photo_change_by_syndrome.csv')


def candidate_rankings(out):
    folder = out / 'aggregate_candidates'
    matrix = pd.read_csv(out / 'evidence/P0_exhaustive_evidence_matrix.csv', dtype={'canonical_syndrome_id': 'string'}, low_memory=False)
    photo = pd.read_csv(folder / 'P1_exact_photo_stability.csv', dtype={'canonical_syndrome_id': 'string'})
    matrix = matrix.merge(photo[KEY + ['rank_q025', 'rank_median', 'rank_q975', 'top10_frequency', 'raw_p_q025',
                                      'raw_p_median', 'raw_p_q975', 'global_q_q025', 'global_q_median', 'global_q_q975']],
                          on=KEY, validate='one_to_one')
    for comparator, prefix in [('FairFace', 'patient_ff_'), ('other_GMDB', 'patient_gmdb_')]:
        frame = pd.read_csv(folder / 'P1_patient_resampling.csv', dtype={'canonical_syndrome_id': 'string'})
        frame = frame.loc[frame.comparator.eq(comparator)]
        fields = ['same_direction_frequency', 'bootstrap_q025', 'bootstrap_q975', 'subsample_supported',
                  'subsample_direction_frequency', 'subsample_q025', 'subsample_q975']
        matrix = matrix.merge(frame[KEY + fields].rename(columns={c: prefix + c for c in fields}), on=KEY, validate='one_to_one')
    clusters = pd.read_csv(folder / 'P1_cluster_assignments.csv')
    clusters = clusters.loc[clusters.source.eq('FairFace'), ['feature', 'cluster_id', 'representative']]
    matrix = matrix.merge(clusters, on='feature', validate='many_to_one')
    matrix['magnitude_primary'] = matrix.ff_hedges_g.abs()
    matrix['minimum_direction_frequency'] = matrix[['same_direction_photo_frequency', 'patient_ff_same_direction_frequency',
                                                    'patient_gmdb_same_direction_frequency']].min(axis=1)
    matrix['minimum_absolute_effect'] = matrix[['ff_hedges_g', 'gmdb_hedges_g', 'all_reference_hedges_g', 'sex_reference_hedges_g']].abs().min(axis=1)
    matrix['screen_supported'] = matrix.screen_status.eq('ok')
    matrix['classification'] = np.where(matrix.clinical_clinical_support.notna(), 'related geometry; assess exact source and denominator',
                                        'unmapped clinical correspondence; NOT evidence of novelty')
    supported = matrix.ff_n_disease.ge(20)
    selected_primary = supported & matrix.ff_mannwhitney_q_fdr_global.lt(.05)
    rules = {
        'high_magnitude': (selected_primary, ['magnitude_primary'], [False]),
        'high_stability': (selected_primary & matrix.reference_sign_consistent & matrix.minimum_direction_frequency.ge(.95),
                           ['minimum_direction_frequency', 'minimum_absolute_effect'], [False, False]),
        'conservative_cross_comparator': (selected_primary & matrix.gmdb_q_global_bh.lt(.05) & matrix.same_direction
            & matrix.screen_supported & matrix.screen_q_global_bh.lt(.05) & matrix.reference_sign_consistent
            & matrix.same_direction_photo_frequency.ge(.95), ['minimum_absolute_effect'], [False]),
        'unverified_geometry_hypotheses': (selected_primary & matrix.clinical_clinical_support.isna(), ['magnitude_primary'], [False]),
        'primary_q_baseline': (selected_primary, ['ff_mannwhitney_q_fdr_global'], [True]),
        'adjusted_q_baseline': (selected_primary & matrix.screen_supported, ['screen_q_global_bh'], [True])}
    rankings, top_keys = [], {}
    for name, (mask, columns, order) in rules.items():
        part = matrix.loc[mask].sort_values(columns + KEY, ascending=order + [True, True]).copy()
        part['ranking_rule'], part['rank'] = name, np.arange(1, len(part) + 1)
        rankings.append(part)
        top_keys[name] = set(map(tuple, part.head(20)[KEY].to_numpy()))
    ranked = pd.concat(rankings, ignore_index=True)
    write_aggregate(matrix, folder / 'P1_all_candidate_evidence_private.csv')
    write_aggregate(ranked.loc[ranked['rank'].le(20)], out / 'P1_ranked_candidate_evidence.csv')
    overlaps = []
    for a in rules:
        for b in rules:
            overlaps.append(dict(first_rule=a, second_rule=b, first_top=len(top_keys[a]), second_top=len(top_keys[b]),
                                  shared_top20=len(top_keys[a] & top_keys[b]),
                                  jaccard=len(top_keys[a] & top_keys[b]) / max(1, len(top_keys[a] | top_keys[b]))))
    write_aggregate(pd.DataFrame(overlaps), folder / 'P1_ranking_rule_overlap.csv')
    dependence = []
    for rule, group in ranked.groupby('ranking_rule', sort=True):
        for reduction in ['raw', 'strict_reference_representative', 'one_per_family', 'one_per_reference_cluster']:
            part = group.sort_values('rank')
            if reduction == 'strict_reference_representative':
                part = part.loc[part.representative]
            elif reduction == 'one_per_family':
                part = part.drop_duplicates(['canonical_syndrome_id', 'feature_family'])
            elif reduction == 'one_per_reference_cluster':
                part = part.drop_duplicates(['canonical_syndrome_id', 'cluster_id'])
            selected = set(map(tuple, part.head(20)[KEY].to_numpy()))
            dependence.append(dict(ranking_rule=rule, reduction=reduction, selected=len(selected),
                                   shared_with_raw_top20=len(selected & top_keys[rule]),
                                   role='descriptive selection sensitivity; no new hypothesis family or independent dimension claim'))
    write_aggregate(pd.DataFrame(dependence), folder / 'P1_ranking_dependence_sensitivity.csv')
    photo_data = np.load(out / 'checkpoints/exact_photo_recovery.npz')
    config = json.loads((ROOT / 'ACTIVE_ANALYSIS.json').read_text())
    targets = json.loads((ROOT / config['private_run'] / 'ANALYSIS_PLAN.json').read_text())['targets']
    names = features().feature.tolist()
    ordered = pd.MultiIndex.from_product([targets, names], names=['syndrome', 'feature'])
    lookup = matrix.set_index(['syndrome', 'feature']).loc[ordered].reset_index()
    eligible = lookup.ff_n_disease.ge(20).to_numpy()
    baseline = np.abs(lookup.ff_hedges_g.to_numpy())
    original = set(np.flatnonzero(eligible)[np.argsort(-baseline[eligible], kind='stable')[:20]])
    changes = []
    for rep, effects in enumerate(photo_data['g']):
        selected = set(np.flatnonzero(eligible)[np.argsort(-np.abs(effects[eligible]), kind='stable')[:20]])
        changes.append(dict(perturbation='existing_photo_draw', repetition=rep, top20_shared_with_primary=len(original & selected),
                            evaluation_status='exploratory same-patient perturbation; no unbiased predictive performance'))
    write_aggregate(pd.DataFrame(changes), folder / 'P1_top20_photo_variation.csv')
    ecological = pd.read_csv(folder / 'P0_ecological_quality_overlap.csv')
    correlations = []
    for proxy, part in ecological.groupby('proxy', sort=True):
        rho, p = spearmanr(part.cohort_mean, part.shared / 120)
        correlations.append(dict(proxy=proxy, syndrome_units=len(part), rho=rho, p_value=p,
                                  role='ecological association with shared-cell fraction; not independent feature observations or causal source evidence'))
    correlations = pd.DataFrame(correlations)
    correlations['exploratory_q_six_proxies'] = bh(correlations.p_value)
    write_aggregate(correlations, folder / 'P0_ecological_quality_associations.csv')
    text = ('# P1 transparent candidate selection rules\n\n'
        'All rankings are exploratory descriptions of the same patients and references, not an unbiased prediction benchmark. '
        'No weighted composite score is created or validated. Headline tables require original case n>=20 to avoid optimistic ranking of n≈10 cohorts; '
        'the full 5,760-row evidence matrix retains all smaller cohorts with uncertainty/support labels. Unsupported models are missing, not a proved null.\n\n'
        'High magnitude orders original globally significant pairs by |primary g|. High stability requires primary q<.05, reference sign consistency and '
        'minimum case-bootstrap/photo direction frequency >=.95, then orders that minimum frequency and worst absolute effect across comparators/references. '
        'Conservative cross-comparator candidates require separate primary and GMDB q<.05, same direction, supported full-screen q<.05, '
        'reference sign consistency and >=.95 photo direction frequency; ordering uses the smaller effect magnitude. '
        'Unverified geometry lists unmapped primary associations without claiming novelty. The q baseline and adjusted-model-q baseline are also shown; '
        'the latter avoids ranking raw beta across incompatible feature scales. All ties use canonical ID and feature.\n\n'
        'These specific filtering rules were specified transparently after reviewing P0/P1 results and are post hoc, NOT preregistered. '
        'No advantage over the baselines is claimed. Overlap of top-20 lists, 100 existing photo-draw rankings and family/cluster collapsing quantify selection sensitivity. '
        'Representative choice is reference-only; choosing the highest ranked feature within each family/cluster is outcome-dependent and labeled descriptive. '
        'Equal family representation is a selection sensitivity, not a recalculated primary FDR family.\n\n'
        'P1_ranked_candidate_evidence.csv stays PRIVATE pending disclosure/PI review. It contains aggregate hypothesis evidence and full canonical syndrome IDs, '
        'not patient/image identifiers. Do not publish it merely because names are absent.\n')
    (out / 'P1_selection_rules.md').write_text(text)
    return matrix, ranked


def method_decision(out):
    text = ('# P2 methodological decision: NOT_JUSTIFIED\n\n'
        'P0/P1 expose real comparator/normalization sensitivity, correlated geometric features and limited photograph perturbation opportunities. '
        'These support an exploratory robustness study, not a new extractor or a validated composite ranking method. '
        'In particular, perfect source/case alignment cannot be repaired or identified by reweighting the same observational effects.\n\n'
        'Objective target for a possible future method would be held-out patient effect direction/magnitude, or known-signal recovery under independent numerical simulation regimes. '
        'It must not be agreement with the proposed score itself. Simple baselines are primary q, minimum-n |g|, scale-aware adjusted evidence, '
        'same-direction significance across comparator families and fixed reference-cluster representatives. Current tables describe these baselines without claiming one wins.\n\n'
        'A future leakage-free benchmark would freeze all rules on discovery patients/perturbations, retain patient grouping, and evaluate only untouched partitions/settings; '
        'n≈10 cohorts cannot support credible split-sample rankings. Internal partitions would still not be independent clinical replication. '
        'An improvement claim would be falsified by no advantage in held-out effect agreement or simulated known-signal precision/recall/calibration, '
        'or loss of any gain after removing individual components in ablations. Confidence intervals would need independent benchmark repetitions/partitions.\n\n'
        'No novel score, optimization, classifier, generator training, held-out performance claim or ablation study was implemented. '
        'Those conditional Milestone D method-development deliverables are NOT_JUSTIFIED rather than falsely marked complete. '
        'The negative-control and coverage simulations test declared workflow behavior under assumptions; they do not establish superiority or clinical validity.\n')
    (out / 'P2_method_decision.md').write_text(text)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    offline()
    out = private_output(args.output_dir)
    protect_baseline(out)
    real, reference, raw, alternate, paths = load_private()
    photo_change_details(real, out)
    candidate_rankings(out)
    method_decision(out)
    protect_baseline(out)
    print('Private exploratory rankings, photo-change context and NOT_JUSTIFIED method decision complete')


if __name__ == '__main__':
    main()

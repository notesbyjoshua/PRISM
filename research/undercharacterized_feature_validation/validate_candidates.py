import argparse
import importlib.metadata
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from research.race_stratified_sensitivity import core

HERE = Path(__file__).resolve().parent
KEYS = core.KEYS


def locked_candidates(document):
    rows = []
    for line in document.read_text().splitlines():
        if re.match(r'\| \d+ \| `', line):
            parts = [p.strip().strip('`') for p in line.strip('|').split('|')]
            rows.append(dict(candidate_rank=int(parts[0]), disease=parts[1], canonical_race=parts[2],
                             feature=parts[3], target_n=int(parts[4]), target_g=float(parts[5]), target_q=float(parts[6])))
    if len(rows) != 10 or [r['candidate_rank'] for r in rows] != list(range(1, 11)):
        raise ValueError('Exactly ten locked candidates in supplied order required')
    return pd.DataFrame(rows)


def percentiles(x, y):
    y = np.sort(y)
    return 50 * (np.searchsorted(y, x, side='left') + np.searchsorted(y, x, side='right')) / len(y)


def check_patients(frame):
    if frame.patient_id.isna().any() or not frame.patient_id.is_unique:
        raise ValueError('Distinct recorded patient required for each case row')


def loo(x, y):
    return pd.DataFrame([{'g': float(core.effect(np.delete(x, i), y)),
                          'mean': float(np.delete(x, i).mean()),
                          'p_raw': float(mannwhitneyu(np.delete(x, i), y, alternative='two-sided', method='auto').pvalue)}
                         for i in range(len(x))])


def decomposition(numerator, denominator, frozen, epsilon=0):
    if numerator is None or denominator is None:
        return {'status': 'not_evaluable', 'reason': 'raw aligned numerator/denominator not stored'}
    numerator, denominator = np.asarray(numerator), np.asarray(denominator)
    if not np.isfinite(numerator).all() or not np.isfinite(denominator).all() or (denominator <= 1e-10).any():
        raise ValueError('Missing/nonfinite components or near-zero denominator')
    np.testing.assert_allclose(numerator / (denominator + epsilon), frozen, rtol=1e-10, atol=1e-12)
    return {'status': 'formula_verified', 'absolute_anatomy_established': False}


def json_value(value):
    return json.dumps(value, default=lambda x: x.item() if isinstance(x, np.generic) else str(x), sort_keys=True)


def run(document, audit_directory, output):
    output = core.baseline.private_output(output)
    plan_path = HERE / 'analysis_plan.json'
    plan = json.loads(plan_path.read_text())
    core.write_json(output / 'analysis_plan_frozen.json', plan)
    protected = core.baseline.protected_hashes()
    loaded = core.load_baseline()
    disease, reference, features, definitions, source_audit, counts, primary, sampling, sources = loaded
    candidates = locked_candidates(document)
    baseline_dir = ROOT / json.loads((ROOT / 'SENSITIVITY_ANALYSIS.json').read_text())['baseline_private_output']
    reference = core.baseline.sample_reference(reference, sampling['quotas'], sampling['seed'])
    audit_path = audit_directory / 'all_982_association_literature_audit.csv'
    report_path = audit_directory / 'LITERATURE_COMPARISON_REPORT.md'
    if not audit_path.is_file() or not report_path.is_file():
        raise FileNotFoundError('Exact supplied audit CSV and report required; use --audit-directory')
    literature = pd.read_csv(audit_path, float_precision='round_trip')
    if len(literature) != 982 or literature[KEYS].duplicated().any():
        raise ValueError('Unexpected literature inventory')
    selected = candidates.merge(primary, on=KEYS, validate='one_to_one').merge(
        literature[KEYS + ['literature_status', 'targeted_narrative_review_performed', 'targeted_narrative_sources', 'comparability_note']],
        on=KEYS, validate='one_to_one')
    if len(selected) != 10 or not selected.case_n.eq(selected.target_n).all():
        raise ValueError('Locked candidates missing or sample counts differ')
    np.testing.assert_allclose(selected.hedges_g, selected.target_g, rtol=0, atol=5.1e-7)
    np.testing.assert_allclose(selected.q, selected.target_q, rtol=5e-8, atol=1e-14)
    if not selected.literature_status.eq('not_identified_in_reviewed_annotation_catalogs').all():
        raise ValueError('Candidates differ from the supplied undercharacterized audit class')
    recalculated = []
    for (syndrome, race), group in primary.groupby(['disease', 'canonical_race'], sort=False):
        case = disease.loc[disease.disease.eq(syndrome) & disease.canonical_race.eq(race)]
        check_patients(case)
        ref = reference.loc[reference.canonical_race.eq(race)]
        u, p, g = core.tests_direct(case[features].to_numpy(float), ref[features].to_numpy(float))
        replay = pd.DataFrame({'disease': syndrome, 'canonical_race': race, 'feature': features, 'replayed_p': p, 'replayed_g': g, 'replayed_u': u})
        recalculated.append(replay)
    replay = primary.merge(pd.concat(recalculated), on=KEYS, validate='one_to_one')
    for old, new in [('p', 'replayed_p'), ('hedges_g', 'replayed_g'), ('mannwhitney_u', 'replayed_u')]:
        np.testing.assert_allclose(replay[old], replay[new], rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(core.bh(replay.replayed_p), replay.q, rtol=1e-10, atol=1e-12)
    sensitivity_dir = core.DEFAULT_OUTPUT
    sensitivity_hashes = json.loads((sensitivity_dir / 'artifact_sha256.json').read_text())
    sensitivity = {}
    source_hashes_at_start = {str(p):core.baseline.sha(p) for p in sources}
    for name in ['reference_resampling_aggregate.csv', 'gmdb_permutation_aggregate.csv']:
        path = sensitivity_dir / name
        if core.baseline.sha(path) != sensitivity_hashes[name]:
            raise ValueError('Sensitivity source hash mismatch')
        sensitivity[name] = pd.read_csv(path, float_precision='round_trip')
        sources.append(path)
    raw_path = core.baseline.EXTRACTION / 'gmdb/gmdb_120_all_records_measurements.parquet'
    raw_columns = pd.read_parquet(raw_path).columns.tolist()
    extractor = ROOT / 'FaceKit/src/facekit/core/geometric/extractor.py'
    if core.baseline.sha(extractor) != definitions.extractor_sha256.iloc[0]:
        raise ValueError('Pinned extraction implementation differs')
    definitions.to_csv(output / 'measurement_definitions_reviewed.csv', index=False)
    sources += [document, audit_path, report_path, plan_path, extractor, raw_path,
                audit_directory / 'review_manifest.json', sensitivity_dir / 'artifact_sha256.json']
    rows, cross_rows, loo_rows, context_rows = [], [], [], []
    for entry in selected.sort_values('candidate_rank').itertuples():
        case = disease.loc[disease.disease.eq(entry.disease) & disease.canonical_race.eq(entry.canonical_race)].sort_values('patient_id')
        ref = reference.loc[reference.canonical_race.eq(entry.canonical_race)]
        check_patients(case)
        x, y = case[entry.feature].to_numpy(float), ref[entry.feature].to_numpy(float)
        np.testing.assert_allclose([x.mean(), y.mean()], [entry.case_mean, entry.reference_mean], rtol=1e-12)
        if not np.isfinite(x).all() or not np.isfinite(y).all():
            raise ValueError('Unexpected missing/nonfinite measurements')
        direct = mannwhitneyu(x, y, alternative='two-sided', method='auto')
        np.testing.assert_allclose(direct.pvalue, entry.p, rtol=1e-10, atol=1e-15)
        sign = np.sign(entry.hedges_g)
        pct = percentiles(x, y)
        beyond = (x > np.median(y)) if sign > 0 else (x < np.median(y))
        perturb = loo(x, y)
        perturb['candidate_rank'] = entry.candidate_rank
        perturb['omission_index_nonidentifying'] = range(len(x))
        loo_rows.append(perturb)
        agreement = float(np.mean(np.sign(perturb.g) == sign))
        delta = float(np.max(np.abs(perturb.g - entry.hedges_g)))
        rng = np.random.default_rng(np.random.SeedSequence([plan['seed'], entry.candidate_rank]))
        draws = x[rng.integers(0, len(x), size=(plan['case_bootstrap_B'], len(x)))]
        gs = core.effect_moments(len(x), len(y), draws.mean(1), draws.var(1, ddof=1), y.mean(), y.var(ddof=1))
        lo, hi = np.quantile(gs, [.025, .975])
        row = {key: getattr(entry, key) for key in KEYS}
        row.update(candidate_rank=entry.candidate_rank, case_n=len(x), reference_n=len(y), case_mean=x.mean(), reference_mean=y.mean(),
                   case_median=np.median(x), reference_median=np.median(y), case_sd=x.std(ddof=1), reference_sd=y.std(ddof=1),
                   case_iqr=np.ptp(np.quantile(x, [.25, .75])), reference_iqr=np.ptp(np.quantile(y, [.25, .75])),
                   raw_mean_difference=x.mean()-y.mean(), hedges_g=entry.hedges_g, mannwhitney_u=direct.statistic,
                   p_primary=entry.p, q_primary_global_9600=entry.q, primary_MW_method=entry.mw_method,
                   case_fraction_expected_side_of_reference_median=beyond.mean(), cases_expected_side=int(beyond.sum()),
                   cases_below_reference_median=int((x < np.median(y)).sum()), cases_equal_reference_median=int((x == np.median(y)).sum()),
                   cases_above_reference_median=int((x > np.median(y)).sum()), median_case_reference_percentile=np.median(pct),
                   case_percentile_iqr=np.ptp(np.quantile(pct, [.25, .75])), case_percentile_q25=np.quantile(pct,.25), case_percentile_q75=np.quantile(pct,.75),
                   cases_percentile_below_10=int((pct < 10).sum()), cases_percentile_below_25=int((pct < 25).sum()),
                   cases_percentile_above_75=int((pct > 75).sum()), cases_percentile_above_90=int((pct > 90).sum()),
                   probability_case_exceeds_reference=direct.statistic/(len(x)*len(y)),
                   rank_biserial=2*direct.statistic/(len(x)*len(y))-1,
                   loo_g_min=perturb.g.min(), loo_g_median=perturb.g.median(), loo_g_max=perturb.g.max(), loo_sign_agreement=agreement,
                   loo_max_abs_delta_g=delta, loo_raw_p_min=perturb.p_raw.min(), loo_raw_p_max=perturb.p_raw.max(),
                   loo_sign_reversals=int((np.sign(perturb.g)!=sign).sum()), case_bootstrap_g_lower=lo, case_bootstrap_g_upper=hi,
                   case_bootstrap_B=plan['case_bootstrap_B'], bootstrap_scope='case-only conditional; fixed reference; not adjusted for selection',
                   missing_case_measurements=0, missing_reference_measurements=0, nonpositive_case_values=int((x<=0).sum()),
                   nonpositive_reference_values=int((y<=0).sum()), extreme_value_dependence=delta>.25*abs(entry.hedges_g),
                   numerator_status='not_evaluable', numerator_case_vs_reference='not_evaluable', denominator_case_vs_reference='not_evaluable',
                   normalization_assessment='raw aligned landmark numerators, bizygomatic_width and midface_height not stored; ratios cannot identify contributions',
                   patient_consistency_score=2 if beyond.mean()>=.8 else 1 if beyond.mean()>=.6 else 0,
                   loo_robustness_score=2 if agreement==1 and delta<=.25*abs(entry.hedges_g) else 1 if agreement==1 else 0,
                   normalization_interpretability_score=np.nan,
                   clinical_interpretability_score=int(entry.feature in ['chin_width','chin_height','forehead_width_mid']),
                   clinical_evidence_status='targeted_sources_reviewed; candidate_undercharacterized; no_independent_validation' if entry.targeted_narrative_review_performed else 'catalog_only; candidate_undercharacterized; no_independent_validation',
                   additional_targeted_review_required=not entry.targeted_narrative_review_performed,
                   missing_analyses='numerator/denominator components unavailable; independent clinical/landmark corroboration unavailable; optional age/sex matching not run',
                   key_cautions=entry.comparability_note + '; GMDB/FairFace source confounding; conditional age; recorded labels; correlated anatomical families')
        rs = sensitivity['reference_resampling_aggregate.csv']
        rs = rs.loc[(rs[KEYS] == pd.Series({k:getattr(entry,k) for k in KEYS})).all(axis=1)]
        if len(rs)!=1:
            raise ValueError('Missing/duplicate reference sensitivity candidate')
        rs = rs.iloc[0]
        np.testing.assert_allclose([rs.observed_main_g,rs.observed_main_q],[entry.hedges_g,entry.q],rtol=1e-12)
        row.update(reference_resampling_B=int(rs.repetitions), reference_varied=bool(rs.reference_varies),
                   reference_g_sign_agreement=rs.sign_agreement_frequency, reference_q_lt_0_05_frequency=rs.significant_frequency,
                   reference_median_g=rs.median_g, reference_g_resampling_q025=rs.g_resampling_q025, reference_g_resampling_q975=rs.g_resampling_q975,
                   reference_stability_assessment='varying references; descriptive only' if rs.reference_varies else 'no reference variation; European eligible pool exhausted; q may change via other groups',
                   comparator_stability_score=(2 if rs.sign_agreement_frequency>=.95 and rs.significant_frequency>=.8 else 1 if rs.sign_agreement_frequency>=.95 else 0) if rs.reference_varies else np.nan)
        perms = sensitivity['gmdb_permutation_aggregate.csv']
        perms = perms.loc[(perms[KEYS] == pd.Series({k:getattr(entry,k) for k in KEYS})).all(axis=1)]
        statuses = []
        for label,spec in [('unblocked','within_race'),('blocked','within_race_sex_age')]:
            pr = perms.loc[perms.null_specification.eq(spec)]
            if len(pr)!=1:
                raise ValueError('Missing/duplicate permutation candidate')
            pr=pr.iloc[0]
            for target,field in [('p_empirical','p_empirical_permutation'),('q','q_global_BH_empirical'),('g_vs_other_GMDB','observed_hedges_g'),
                                 ('B','B'),('p_resolution','p_value_floor'),('q_mc_simultaneous_upper','q_MC_simultaneous_upper'),('q_mc_simultaneous_lower','q_MC_simultaneous_lower')]:
                row[f'permutation_{label}_{target}']=pr[field]
            statuses.append(f'{label}: '+('upper-bound significant' if pr.q_MC_simultaneous_upper<.05 else 'not significant under simultaneous MC upper bound'))
        row['permutation_mc_conservative_status']='; '.join(statuses)
        same = primary.loc[primary.disease.eq(entry.disease) & primary.feature.eq(entry.feature)]
        others = same.loc[~same.canonical_race.eq(entry.canonical_race) & same.supported]
        row['other_category_effects']=json_value(same[['canonical_race','case_n','hedges_g','q','supported']].to_dict('records'))
        row['cross_category_score']=np.nan if len(others)==0 else 0 if (np.sign(others.hedges_g)!=sign).any() else 2 if others.q.lt(.05).any() else 1
        for race in core.RACES:
            cell=same.loc[same.canonical_race.eq(race)]
            cross_rows.append({'candidate_rank':entry.candidate_rank,'canonical_race':race,'status':'eligible' if len(cell) else 'not_eligible',
                               'n':int(cell.case_n.iloc[0]) if len(cell) else None,'g':float(cell.hedges_g.iloc[0]) if len(cell) else None,
                               'q_global_9600':float(cell.q.iloc[0]) if len(cell) else None,
                               'direction_vs_candidate':('same' if np.sign(cell.hedges_g.iloc[0])==sign else 'opposite') if len(cell) else 'not_evaluable'})
        other=primary.loc[primary.canonical_race.eq(entry.canonical_race)&primary.feature.eq(entry.feature)&~primary.disease.eq(entry.disease)]
        distinct={'eligible_other_syndromes':len(other),'significant_same_direction':int((other.q.lt(.05)&(np.sign(other.hedges_g)==sign)).sum()),
                  'significant_opposite_direction':int((other.q.lt(.05)&(np.sign(other.hedges_g)!=sign)).sum()),'nonsignificant':int(other.q.ge(.05).sum()),
                  'candidate_directional_g_rank':1+int((sign*other.hedges_g>abs(entry.hedges_g)).sum()),'rank_denominator':len(other)+1}
        row['other_syndrome_same_feature_summary']=json_value(distinct)
        row.update(distinct)
        components=[]
        if entry.feature.endswith('_mean'):
            for side in ['r','l']:
                name=entry.feature[:-4]+side
                result=core.effect(case[name].to_numpy(float),ref[name].to_numpy(float))
                components.append({'feature':name,'g':float(result),'same_direction':bool(np.sign(result)==sign)})
        row['bilateral_component_context']=json_value(components)
        for feature in ['forehead_width_upper','face_triangularity']:
            values=primary.loc[primary.disease.eq(entry.disease)&primary.canonical_race.eq(entry.canonical_race)&primary.feature.eq(feature)]
            for value in values.itertuples():
                context_rows.append({'candidate_rank':entry.candidate_rank,'context_feature':feature,'g':value.hedges_g,'q_original_global_9600':value.q,'status':'context only; not independent discovery'})
        score_fields=['patient_consistency_score','loo_robustness_score','normalization_interpretability_score','comparator_stability_score','cross_category_score','clinical_interpretability_score']
        available=[row[k] for k in score_fields if pd.notna(row[k])]
        row.update(score_earned=sum(available),score_possible=2*len(available),completeness=len(available)/6,score_fraction=sum(available)/(2*len(available)))
        rows.append(row)
    summary=pd.DataFrame(rows).sort_values('candidate_rank')
    summary.to_csv(output/'candidate_validation_summary.csv',index=False)
    pd.DataFrame(cross_rows).to_csv(output/'cross_category_evidence.csv',index=False)
    pd.DataFrame(context_rows).to_csv(output/'related_feature_context.csv',index=False)
    pd.concat(loo_rows).to_csv(output/'PRIVATE_loo_perturbations.csv',index=False)
    evidence=summary.melt(id_vars=['candidate_rank']+KEYS,var_name='metric',value_name='value')
    evidence['status']=np.where(evidence.value.isna(),'not_evaluable','reported')
    evidence['provenance']='locked primary replay; frozen analysis plan; existing sensitivity tables; see run_manifest.json'
    evidence.to_csv(output/'candidate_evidence_long.csv',index=False)
    ranked=summary.sort_values(['score_fraction','score_possible','candidate_rank'],ascending=[False,False,True]).copy()
    ranked.insert(0,'followup_priority_rank',range(1,11))
    ranked.to_csv(output/'ranked_followup_priorities.csv',index=False)
    report(ranked,output,plan)
    figure(summary,output)
    for name, expected in source_hashes_at_start.items():
        if core.baseline.sha(Path(name)) != expected:
            raise ValueError('Input changed during run')
    if protected!=core.baseline.protected_hashes():
        raise ValueError('Protected artifact changed')
    manifest={'study_id':plan['study_id'],'active_manuscript_study_unchanged':'PRISM_120_20261028',
              'baseline_verified_tests':len(replay),'baseline_discoveries':int(primary.q.lt(.05).sum()),
              'eligible_patients':int(counts.loc[counts.eligible,'n_current_gmdb'].sum()),'sampling':sampling,
              'source_audit':source_audit,'seed':plan['seed'],'bootstrap_B':plan['case_bootstrap_B'],
              'python_version':sys.version, 'input_schema':{'case_columns':disease.columns.tolist(), 'reference_columns':reference.columns.tolist(), 'features':features},
              'software':{k:importlib.metadata.version(k) for k in ['numpy','pandas','scipy','matplotlib']},
              'git_revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
              'input_sha256':{str(p):core.baseline.sha(p) for p in set(sources)},
              'code_sha256':{str(p.relative_to(ROOT)):core.baseline.sha(p) for p in HERE.rglob('*') if p.suffix in ['.py','.json']},
              'raw_extraction_columns':raw_columns,'raw_components_stored':False,'protected_files_unchanged':True,
              'decomposition':'not_evaluable; extraction stores final metrics only; no image access or reconstruction attempted',
              'optional_matching':'not run; FairFace age 3-9 bin vs conditional GMDB years, compatible sex labels do not remove source confounding',
              'placebo':'method context only; no candidate-specific inferential p',
              'output_sha256':{p.name:core.baseline.sha(p) for p in output.rglob('*') if p.is_file() and p.name!='run_manifest.json'}}
    core.write_json(output/'run_manifest.json',manifest)
    print(ranked[['candidate_rank','feature','score_earned','score_possible','case_fraction_expected_side_of_reference_median','loo_sign_agreement']].to_string(index=False))


def report(ranked,output,plan):
    lines=['# Internal robustness and follow-up priority of ten locked associations','',
           'All ten baselines and the full 9,600-test family were replayed on authorized real GMDB/FairFace tables. '
           'These are internally assessed, selected associations, not independent clinical validation. No quantitative clinical corroboration '
           'was performed; no observation is externally validated. Catalog absence is not novelty.','',
           'The versioned analysis_plan.json froze breakpoints before new candidate outcomes. Ranking uses earned/possible score, '
           'then evidence completeness, then locked candidate order. Different denominators require caution; scores are descriptive, not p-values.','',
           'Reference percentiles use midranks: less-than plus half of ties. U is case-first, so U/(n*m) estimates '
           'P(case>reference)+0.5 P(tie). This is not diagnostic sensitivity. LOO reports raw p only, with no new q values. '
           'Bootstrap intervals use 2,000 case-only resamples with the reference fixed; they do not account for discovery selection, '
           'reference uncertainty or source confounding. Original q values always use all 9,600 primary tests.','',
           'European reference images cannot vary. Their fixed effects earn no comparator-stability score; q may move through other cells. '
           'Permutation effects contrast other GMDB syndromes, not FairFace. Point-estimate permutation q values must not be called '
           'Monte Carlo-certain; all conservative simultaneous upper-bound decisions here remain nonsignificant.','',
           '## Ranked ten','']
    for rank,row in enumerate(ranked.itertuples(),1):
        distinct=json.loads(row.other_syndrome_same_feature_summary)
        lines += [f'### {rank}. Candidate {row.candidate_rank}: {row.disease.split("__")[0]} / {row.canonical_race} / {row.feature}','',
                  f'Priority score {row.score_earned:g}/{row.score_possible:g}; completeness {row.completeness:.0%}. '
                  f'n={row.case_n}, reference n={row.reference_n}; g={row.hedges_g:+.3f}, original q={row.q_primary_global_9600:.5g}. '
                  f'{row.cases_expected_side}/{row.case_n} ({row.case_fraction_expected_side_of_reference_median:.1%}) lie strictly on the '
                  f'expected side of the reference median. Median reference percentile {row.median_case_reference_percentile:.1f} '
                  f'(IQR {row.case_percentile_q25:.1f}–{row.case_percentile_q75:.1f}).','',
                  f'LOO g range {row.loo_g_min:+.3f} to {row.loo_g_max:+.3f}; sign retention {row.loo_sign_agreement:.0%}; '
                  f'maximum |Δg|={row.loo_max_abs_delta_g:.3f}. Case-only conditional bootstrap 95% percentile interval '
                  f'[{row.case_bootstrap_g_lower:+.3f}, {row.case_bootstrap_g_upper:+.3f}].','',
                  f'Reference discovery frequency {row.reference_q_lt_0_05_frequency:.1%}; sign agreement '
                  f'{row.reference_g_sign_agreement:.1%}. {row.reference_stability_assessment}. '
                  f'Within-GMDB point permutation q: unblocked {row.permutation_unblocked_q:.5g}, blocked {row.permutation_blocked_q:.5g}. '
                  f'{row.permutation_mc_conservative_status}.','',
                  f'Other syndromes: {distinct["significant_same_direction"]}/{distinct["eligible_other_syndromes"]} significant in the '
                  f'same direction, {distinct["significant_opposite_direction"]} opposite. Effect rank '
                  f'{distinct["candidate_directional_g_rank"]}/{distinct["rank_denominator"]} is exploratory, without a rank uncertainty interval '
                  'or diagnostic specificity claim. See complete cross-category CSV, including nonsignificant and missing categories.','',
                  f'Clinical evidence: {row.clinical_evidence_status}. Additional targeted source review required before using '
                  f'“underreported”: {row.additional_targeted_review_required}. {row.key_cautions}.','']
    lines += ['## Best next targeted reviews','']
    for row in ranked.head(3).itertuples():
        lines += [f'- Candidate {row.candidate_rank}: {row.disease.split("__")[0]} / {row.canonical_race} / {row.feature}; '
                  f'{row.score_earned:g}/{row.score_possible:g}. Main caveat: '+
                  ('no reference-selection variation; ' if not row.reference_varied else 'overlapping reference resamples are not replication; ')+
                  'raw numerator/denominator unavailable, so ratio direction cannot establish absolute anatomy.']
    lines += ['', '## Missing analyses and interpretation limits','',
              '- Numerator/denominator decomposition: **not_evaluable for all ten**. Stored extraction results contain final normalized metrics, '
              'not aligned landmark coordinates or raw scale components. No patient images were accessed. No pixel/millimeter or absolute-size claims are possible.',
              '- Independent landmark/clinical QC and external replication: unavailable. Source QC was accepted, but pose, expression, hairline and eyebrow geometry remain measurement limits.',
              '- Age/sex matching: optional check not run. FairFace supplies age bins; GMDB conditional years are not independently established. No precise FairFace age was invented.',
              '- Same-syndrome category comparisons and bilateral metrics share sources; they are context, not independent replication. Different q values do not establish interaction.',
              '- Other-syndrome effect ranks omit rank uncertainty and covariate adjustment; treat them as exploratory. Anatomically redundant features belong to the same trait family.',
              '- No person-level distribution plot is released. The private figure is an aggregate rubric heatmap; NA remains visibly missing. Original identifiers occur only in restricted source inputs.',
              '- FairFace placebo provides methodological context only, not a candidate-specific p-value. No expensive permutation rerun was performed.',
              '- Neither a high score nor catalog absence establishes a new phenotype, diagnostic usefulness or an effect caused by ancestry.', '',
              'See candidate_validation_summary.csv (locked order), candidate_evidence_long.csv, cross_category_evidence.csv, '
              'related_feature_context.csv, measurement_definitions_reviewed.csv, and run_manifest.json. All real outputs remain ignored and unpublished.']
    (output/'candidate_validation_report.md').write_text('\n'.join(lines)+'\n')


def figure(summary,output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fields=['patient_consistency_score','loo_robustness_score','normalization_interpretability_score',
            'comparator_stability_score','cross_category_score','clinical_interpretability_score']
    values=summary[fields].to_numpy(float)
    fig,ax=plt.subplots(figsize=(8,5))
    cmap=plt.get_cmap('Blues',3).copy()
    cmap.set_bad('#eeeeee')
    ax.imshow(values,vmin=0,vmax=2,cmap=cmap,aspect='auto')
    for i in range(10):
        for j in range(6):
            ax.text(j,i,'NA' if np.isnan(values[i,j]) else str(int(values[i,j])),ha='center',va='center',color='white' if values[i,j]==2 else 'black')
    ax.set_yticks(range(10),[f'Candidate {i}' for i in range(1,11)])
    ax.set_xticks(range(6),['Patient\nconsistency','LOO\nrobustness','Ratio\ncomponents','Reference\nvariation','Other\ncategory','Measurement\ninterpretability'],fontsize=9)
    ax.set_title('Internal follow-up rubric: 0–2 per dimension; NA = unavailable',fontsize=11,pad=12)
    fig.text(.5,.02,'Descriptive priority only; selected findings; no independent clinical validation',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.05,1,1])
    target=output/'figures'
    target.mkdir(exist_ok=True)
    fig.savefig(target/'internal_robustness_rubric.png',dpi=300)
    fig.savefig(target/'internal_robustness_rubric.pdf')
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 5))
    y = np.arange(len(summary))
    ax.hlines(y, summary.case_bootstrap_g_lower, summary.case_bootstrap_g_upper,
              color='#315C7A', linewidth=1.5, label='95% case-only bootstrap percentile interval')
    ax.hlines(y, summary.loo_g_min, summary.loo_g_max,
              color='#BD782C', linewidth=5, label='Leave-one-patient-out g range')
    ax.scatter(summary.hedges_g, y, color='#18354A', s=24, zorder=3, label='Original primary g')
    ax.axvline(0, color='#888888', linewidth=.8, linestyle='--')
    ax.set_yticks(y, [f'{r.candidate_rank}: {r.feature}' for r in summary.itertuples()], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel('Hedges’ g (GMDB − same-category FairFace)')
    ax.set_title('Conditional uncertainty and single-patient influence', fontsize=12)
    fig.legend(*ax.get_legend_handles_labels(), fontsize=8, loc='lower center', bbox_to_anchor=(.65, .06))
    ax.spines[['top', 'right']].set_visible(False)
    fig.text(.5, .015, 'Reference fixed; 2,000 case-only resamples; intervals do not adjust for selection or multiple testing',
             ha='center', fontsize=8)
    fig.tight_layout(rect=[0, .22, 1, 1])
    fig.savefig(target / 'conditional_effect_uncertainty.png', dpi=300)
    fig.savefig(target / 'conditional_effect_uncertainty.pdf')
    plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--candidate-document',type=Path,required=True)
    parser.add_argument('--audit-directory',type=Path,default=ROOT/'local_review/race_stratified_literature_20261115')
    parser.add_argument('--output-dir',type=Path,default=ROOT/'local_review/undercharacterized_feature_validation')
    args=parser.parse_args()
    run(args.candidate_document,args.audit_directory,args.output_dir)

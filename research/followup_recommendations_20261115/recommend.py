import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from research.race_stratified_sensitivity import core
from research.race_stratified_literature_20261115.directions import markdown_table

HERE = Path(__file__).resolve().parent
BASE = core.baseline.PRIVATE
AUDIT = ROOT/'local_review/race_stratified_literature_20261115'
VALIDATION = ROOT/'local_review/undercharacterized_feature_validation'
KEYS = core.KEYS


def contrast(means, sds, ns):
    means,sds,ns = np.asarray(means,float),np.asarray(sds,float),np.asarray(ns,float)
    if means.shape!=(4,) or not np.isfinite(np.r_[means,sds,ns]).all() or (ns<2).any() or (sds<0).any():
        raise ValueError('Four independent supported groups required')
    v = sds**2/ns
    value = np.dot([1,-1,-1,1],means)
    if v.sum()==0:
        return value,0.,np.nan,np.nan,1.,'undefined_zero_variance'
    se = np.sqrt(v.sum())
    df = v.sum()**2/np.sum(v**2/(ns-1))
    stat = value/se
    return value,se,df,stat,2*t.sf(abs(stat),df),'tested'


def shortlist(frame, only_missing=False):
    keep = frame.loc[frame.q_heterogeneity_BH.lt(.05)&(frame.primary_q_1.lt(.05)|frame.primary_q_2.lt(.05))].copy()
    if only_missing:
        keep = keep.loc[keep.missing_catalog_annotation_either_side]
    keep = keep.sort_values(['BY_significant','q_heterogeneity_BH','abs_t','disease','feature','race_1','race_2'],
                            ascending=[False,True,False,True,True,True,True])
    return keep.drop_duplicates(['disease','trait_family']).head(10).copy()


def build(output):
    output = core.baseline.private_output(output)
    sources = [HERE/'analysis_plan.json',BASE/'primary_all_tests.csv',BASE/'statistics_manifest.json',
               AUDIT/'all_982_association_literature_audit.csv',AUDIT/'review_manifest.json',
               AUDIT/'feature_to_clinical_trait_mapping.csv',VALIDATION/'ranked_followup_priorities.csv',
               VALIDATION/'run_manifest.json']
    source_hashes = {str(p):core.baseline.sha(p) for p in sources}
    plan = json.loads(sources[0].read_text())
    core.write_json(output/'analysis_plan_frozen.json',plan)
    manifest = json.loads((BASE/'statistics_manifest.json').read_text())
    if manifest['study_id']!=plan['study'] or source_hashes[str(BASE/'primary_all_tests.csv')]!=manifest['output_sha256']['primary_all_tests.csv']:
        raise ValueError('Locked primary source mismatch')
    audit_manifest = json.loads((AUDIT/'review_manifest.json').read_text())
    if core.baseline.sha(AUDIT/'all_982_association_literature_audit.csv')!=audit_manifest['result_sha256']['all_982_association_literature_audit.csv']:
        raise ValueError('Literature audit source mismatch')
    validation_manifest = json.loads((VALIDATION/'run_manifest.json').read_text())
    if core.baseline.sha(VALIDATION/'ranked_followup_priorities.csv')!=validation_manifest['output_sha256']['ranked_followup_priorities.csv']:
        raise ValueError('Internal validation source mismatch')
    primary = pd.read_csv(BASE/'primary_all_tests.csv',float_precision='round_trip')
    audit = pd.read_csv(AUDIT/'all_982_association_literature_audit.csv',float_precision='round_trip')
    prior = pd.read_csv(VALIDATION/'ranked_followup_priorities.csv',float_precision='round_trip')
    mapping = pd.read_csv(AUDIT/'feature_to_clinical_trait_mapping.csv').set_index('feature')
    if len(primary)!=9600 or primary[KEYS].duplicated().any() or primary.q.lt(.05).sum()!=982:
        raise ValueError('Unexpected primary family')
    ten = prior.merge(audit[KEYS+['literature_status','targeted_narrative_review_performed','comparability_note']],on=KEYS,validate='one_to_one')
    if len(ten)!=10 or not ten.literature_status.eq('not_identified_in_reviewed_annotation_catalogs').all():
        raise ValueError('Ten prior shortlist audit labels do not match')
    joined = ten.merge(primary[KEYS+['hedges_g','q','case_n']],on=KEYS,validate='one_to_one',suffixes=('','_primary'))
    np.testing.assert_array_equal(joined.hedges_g,joined.hedges_g_primary)
    np.testing.assert_array_equal(joined.q_primary_global_9600,joined.q)
    np.testing.assert_array_equal(joined.case_n,joined.case_n_primary)
    if not joined.q.lt(.05).all():
        raise ValueError('Shortlist must retain original primary significance')
    ten.to_csv(output/'ten_undercharacterized_review_candidates.csv',index=False)
    audit_lookup = audit.set_index(KEYS).literature_status.to_dict()
    rows = []
    for (disease,feature), group in primary.groupby(['disease','feature'],sort=True):
        group = group.sort_values('canonical_race')
        for a,b in itertools.combinations(group.itertuples(),2):
            supported = a.supported and b.supported
            means=[a.case_mean,a.reference_mean,b.case_mean,b.reference_mean]
            sds=[a.case_sd,a.reference_sd,b.case_sd,b.reference_sd]
            ns=[a.case_n,a.reference_n,b.case_n,b.reference_n]
            value,se,df,stat,p,status = contrast(means,sds,ns) if supported else (np.nan,np.nan,np.nan,np.nan,1.,'unsupported')
            statuses=[audit_lookup.get((disease,race,feature),'not_in_significant_audit') for race in [a.canonical_race,b.canonical_race]]
            rows.append({'disease':disease,'feature':feature,'trait_family':mapping.loc[feature,'family'],
                         'race_1':a.canonical_race,'race_2':b.canonical_race,'case_n_1':a.case_n,'case_n_2':b.case_n,
                         'reference_n_1':a.reference_n,'reference_n_2':b.reference_n,'g_1':a.hedges_g,'g_2':b.hedges_g,
                         'primary_q_1':a.q,'primary_q_2':b.q,'mean_gap_1':a.case_mean-a.reference_mean,
                         'mean_gap_2':b.case_mean-b.reference_mean,'difference_of_mean_gaps':value,'contrast_se':se,
                         'df_Welch_Satterthwaite':df,'t_contrast':stat,'abs_t':abs(stat),'p_heterogeneity':p,'status':status,
                         'opposite_primary_g_signs':bool(a.hedges_g*b.hedges_g<0),'literature_status_1':statuses[0],
                         'literature_status_2':statuses[1],'missing_catalog_annotation_either_side':'not_identified_in_reviewed_annotation_catalogs' in statuses,
                         'formula':a.formula,'units':a.units})
    result = pd.DataFrame(rows)
    if len(result)!=1920:
        raise ValueError('Full expected 1,920 cross-category contrasts required')
    result['q_heterogeneity_BH']=core.bh(result.p_heterogeneity)
    result['q_heterogeneity_BY']=np.minimum(1,result.q_heterogeneity_BH*np.sum(1/np.arange(1,len(result)+1)))
    result['BY_significant']=result.q_heterogeneity_BY.lt(.05)
    result.to_csv(output/'all_1920_cross_category_mean_contrasts.csv',index=False)
    result.loc[result.q_heterogeneity_BH.lt(.05)].to_csv(output/'all_BH_significant_contrasts.csv',index=False)
    result.loc[result.BY_significant].to_csv(output/'all_BY_significant_contrasts.csv',index=False)
    selected = shortlist(result)
    missing = shortlist(result,True)
    selected.to_csv(output/'ten_cross_category_followup_candidates.csv',index=False)
    missing.to_csv(output/'cross_category_candidates_missing_catalog_annotations.csv',index=False)
    report(ten,selected,missing,result,output)
    for name,h in source_hashes.items():
        if core.baseline.sha(Path(name))!=h:
            raise ValueError('Source changed during analysis')
    core.write_json(output/'run_manifest.json',{'plan':plan,'source_hashes':source_hashes,
                    'code_sha256':core.baseline.sha(Path(__file__)),'n_heterogeneity_tests':len(result),
                    'BH_significant':int(result.q_heterogeneity_BH.lt(.05).sum()),'BY_significant':int(result.BY_significant.sum()),
                    'shortlisted_pairs':len(selected),'missing_annotation_shortlisted_pairs':len(missing),
                    'new_tests_are_secondary_approximate_mean_contrasts':True,'original_q_values_modified':False,
                    'output_sha256':{p.name:core.baseline.sha(p) for p in output.iterdir() if p.is_file() and p.name!='run_manifest.json'}})
    print(result[['q_heterogeneity_BH','q_heterogeneity_BY']].lt(.05).sum().to_string())
    print(selected[['disease','feature','race_1','race_2','g_1','g_2','q_heterogeneity_BH','q_heterogeneity_BY','missing_catalog_annotation_either_side']].to_string(index=False))


def report(ten,selected,missing,all_results,output):
    columns=['disease','canonical_race','feature','case_n','hedges_g','q_primary_global_9600','score_earned','score_possible']
    pair_columns=['disease','feature','race_1','race_2','case_n_1','case_n_2','g_1','g_2','primary_q_1','primary_q_2','q_heterogeneity_BH','q_heterogeneity_BY']
    by_rows = all_results.loc[all_results.BY_significant]
    lines=['# Significant undercharacterized candidates and exploratory category heterogeneity','',
           '**No finding is certified absent from all clinical literature.** The first list consists of ten significant '
           'associations missing from the reviewed HPO/Orphadata trait mappings, with targeted narrative checks where recorded. '
           'This is a bounded literature gap, not proof of novelty. Ratios and 2D proxies can operationalize related known traits. '
           'In particular, nasal tip/columella appearance, eyebrow configuration, chin/jaw gestalt and forehead appearance '
           'may already be clinically described without the exact FaceKit measurement.','',
           '## Ten internally assessed candidates for targeted clinical review','',markdown_table(ten[columns]),'',
           'These retain the prior locked ten, ordered by existing descriptive priority. All retained effect signs in every LOO run; '
           'European reference images cannot vary and were not scored as successful reference-selection replication. '
           'Raw numerator/denominator components were unavailable for all ten. KBG forehead width and Coffin–Siris chin height '
           'had catalog-only coverage and specifically require further targeted narrative review. Existing targeted syndrome reviews '
           'do not establish absence of a particular metric from all publications.','',
           '## New exploratory tests of differences between category-specific case–reference mean gaps','',
           f'All {len(all_results):,} eligible category-pair/feature contrasts were tested, before selecting rows for display. '
           f'{all_results.q_heterogeneity_BH.lt(.05).sum()} met BH q<.05 and {all_results.BY_significant.sum()} met the supplementary '
           'Benjamini–Yekutieli threshold for arbitrary dependence. These are a separate secondary family; original '
           'primary q values continue to refer to all 9,600 Mann–Whitney tests.','',
           'The tested contrast is (case mean A − reference mean A) − (case mean B − reference mean B), in feature units. '
           'For four independent groups, v_i=s_i²/n_i, SE²=Σv_i, df=(Σv_i)²/Σ[v_i²/(n_i−1)], '
           'and a two-sided t tail is used. This generalized Welch–Satterthwaite approximation accounts for four sample variances. '
           'It tests relative-to-reference mean gaps, not equality of standardized Hedges g or the original rank-based estimand. '
           'For the unequal-variance principle, see [NIST](https://itl.nist.gov/div898/software/dataplot/refman1/auxillar/diffmean.htm); '
           'the four-group extension and all formulas are explicit here and in code. '
           '[SciPy documents the BH and BY dependence assumptions](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.false_discovery_control.html).','',
           'Small cells, skew/outliers, correlated traits, shared references, source acquisition, age/sex/pose and broad recorded '
           'categories limit inference. This is an approximate, exploratory interaction contrast; no biological-race, causal, '
           'independent replication or clinical diagnostic conclusion follows. Different primary q values or opposite signs alone '
           'were never used as evidence of heterogeneity.','',
           '## All contrasts surviving the conservative dependence correction','',markdown_table(by_rows[pair_columns]),'',
           'Some interaction contrasts are significant even when neither individual primary test survives its separate 9,600-test correction. '
           'Those contrasts remain in the complete tables but are not included in the primary-significant follow-up shortlist.','',
           '## Ten heterogeneous comparisons for follow-up','',markdown_table(selected[pair_columns]),'',
           'Ranked by passing the conservative BY check, then BH q, absolute t and exact keys. At least one primary association '
           'must be significant. One entry per syndrome/trait-family prevents filling the list with left/right versions of one measurement. '
           'These are two-category comparisons; a single syndrome/feature/race row cannot itself demonstrate category variation.','',
           '## Heterogeneous comparisons also missing a reviewed catalog annotation','',markdown_table(missing[pair_columns]),'',
           'The same full-family q values apply to this filtered list; no favorable re-correction was performed after filtering.','',
           '## Individual measurement cautions','']
    for row in ten.itertuples():
        lines += [f'- {row.disease.split("__")[0]} / {row.canonical_race} / {row.feature}: {row.comparability_note}']
    lines += ['', '## Clinical source scope','',
              'The supporting audit/source registry remain in local_review/race_stratified_literature_20261115. '
              '[Kabuki GeneReviews](https://www.ncbi.nlm.nih.gov/books/NBK62111/) describes a recognizable facial gestalt, '
              'not validation of these particular normalized chin/hairline measures. '
              '[Nicolaides–Baraitser original phenotype paper](https://onlinelibrary.wiley.com/doi/full/10.1002/ajmg.a.32956) '
              'describes facial appearance; a missing exact transverse-forehead annotation is not proof of an unreported clinical trait. '
              'No source’s omission is treated as proof of clinical absence.','',
              'All aggregate exports are ignored/private. No baseline statistics, sensitivity artifacts, patient data or manuscript files were modified.']
    (output/'RECOMMENDATIONS_REPORT.md').write_text('\n'.join(lines)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output-dir',type=Path,default=ROOT/'local_review/followup_recommendations_20261115')
    build(parser.parse_args().output_dir)

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import rankdata

from .common import (ROOT, DEFAULT_OUT, KEY, SEED, canonicalize, features, load_private, offline,
                     private_output, protect_baseline, require, sha, write_aggregate, write_json)


def correlation(values, method='spearman'):
    values = np.asarray(values, float)
    if method == 'spearman':
        values = rankdata(values, axis=0)
    matrix = np.corrcoef(values, rowvar=False)
    require(np.isfinite(matrix).all(), 'Constant/nonfinite features prevent declared correlation clustering')
    return matrix


def clusters(matrix, threshold=.3):
    distance = np.clip(1 - np.abs(matrix), 0, 1)
    np.fill_diagonal(distance, 0)
    return fcluster(linkage(squareform(distance, checks=True), method='average'), threshold, criterion='distance')


def adjusted_rand(a, b):
    table = pd.crosstab(pd.Series(a), pd.Series(b)).to_numpy()
    choose = lambda v: np.sum(v * (v - 1) / 2)
    observed = choose(table)
    row, column, total = choose(table.sum(axis=1)), choose(table.sum(axis=0)), len(a) * (len(a) - 1) / 2
    expected = row * column / total
    denominator = .5 * (row + column) - expected
    return float((observed - expected) / denominator) if denominator else 1.


def dependency_edges(definitions):
    known = set(definitions.feature)
    edges = []
    for row in definitions.itertuples():
        parents = sorted(set(re.findall(r"f\[['\"]([^'\"]+)['\"]\]", row.formula)))
        require(set(parents).issubset(known), 'S6 refers to an unknown derived feature')
        for parent in parents:
            edges.append(dict(parent=parent, child=row.feature, relation='explicit_algebraic_reference', formula=row.formula))
        for denominator in ['bizygomatic_width', 'midface_height', 'face_height', 'mouth_width']:
            if re.search(r'\b' + re.escape(denominator) + r'\b', row.formula):
                edges.append(dict(parent='normalizer:' + denominator, child=row.feature,
                                  relation='shared_normalization_not_an_identity', formula=row.formula))
    return pd.DataFrame(edges)


def representatives(definitions, labels):
    frame = definitions[['feature', 'formula', 'anatomical_region', 'feature_family']].copy()
    frame['cluster_id'] = labels
    frame['derived_references'] = frame.formula.map(lambda s: len(re.findall(r"f\[['\"]", s)))
    frame['asymmetry'] = frame.feature.str.contains('asym|gaze|offset', regex=True)
    frame['formula_length'] = frame.formula.str.len()
    frame['representative'] = False
    indices = frame.sort_values(['cluster_id', 'asymmetry', 'derived_references', 'formula_length', 'feature']).groupby('cluster_id', sort=True).head(1).index
    frame.loc[indices, 'representative'] = True
    frame['selection_rule'] = 'source-reference clustering; prefer non-asymmetry, fewer derived refs, simpler formula, lexical ID; no outcome selection'
    return frame


def bootstrap_clusters(values, source, names, labels, out, repetitions):
    file = out / f'checkpoints/clusters_{source}_v2.npz'
    metadata = file.with_suffix('.json')
    signature = {'code_sha256': sha(Path(__file__)), 'source': source, 'shape': list(values.shape),
                 'seed': SEED, 'repetitions': repetitions, 'threshold': .3}
    if file.exists():
        require(json.loads(metadata.read_text()) == signature, 'Changed cluster-bootstrap signature; choose a new run')
        stored = np.load(file)
        coassign, ari = stored['coassign'], stored['ari']
    else:
        coassign, ari = np.zeros((120, 120)), []
        for rep in range(repetitions):
            rng = np.random.default_rng([SEED, 20, 0 if source == 'FairFace' else 1, rep])
            sample = values[rng.integers(len(values), size=len(values))]
            current = clusters(correlation(sample))
            coassign += current[:, None] == current[None, :]
            ari.append(adjusted_rand(labels, current))
            if (rep + 1) % 25 == 0:
                print(f'Cluster stability {source}: {rep + 1}/{repetitions}', flush=True)
        coassign /= repetitions
        file.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
        np.savez_compressed(file, coassign=coassign, ari=ari)
        write_json(signature, metadata)
    rows = []
    for label in np.unique(labels):
        indices = np.flatnonzero(labels == label)
        pairs = coassign[np.ix_(indices, indices)][np.triu_indices(len(indices), 1)]
        rows.append(dict(source=source, cluster_id=int(label), features=len(indices),
                         mean_within_cluster_coassignment=float(np.mean(pairs)) if len(pairs) else np.nan,
                         repetitions=repetitions, median_partition_adjusted_rand=float(np.median(ari))))
    return pd.DataFrame(rows), coassign


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUT)
    parser.add_argument('--bootstrap', type=int, default=100)
    args = parser.parse_args()
    offline()
    out = private_output(args.output_dir)
    protect_baseline(out)
    real, reference, raw, alternate, paths = load_private()
    definitions = features()
    names = definitions.feature.tolist()
    cohorts = [('FairFace', reference), ('GMDB_selected', real)] + [
        (label, part) for label, part in real.groupby('disease', sort=True) if len(part) >= 20]
    assignments, matrices, stability, source_labels = [], [], [], {}
    for source, frame in cohorts:
        values = frame[names].to_numpy(float)
        spearman, pearson = correlation(values), correlation(values, 'pearson')
        labels = clusters(spearman)
        source_labels[source] = labels
        for i, first in enumerate(names):
            for j in range(i, 120):
                matrices.append(dict(source=source, n=len(frame), feature_a=first, feature_b=names[j],
                                     spearman=spearman[i, j], pearson=pearson[i, j]))
        assignment = representatives(definitions, labels)
        assignment['source'], assignment['n'] = source, len(frame)
        assignments.append(assignment)
        if source in ['FairFace', 'GMDB_selected']:
            table, coassign = bootstrap_clusters(values, source, names, labels, out, args.bootstrap)
            stability.append(table)
            write_aggregate(pd.DataFrame(coassign, columns=names).assign(feature=names),
                            out / f'aggregate_candidates/P1_coclustering_{source}.csv')
    write_aggregate(pd.DataFrame(matrices), out / 'aggregate_candidates/P1_source_correlations.csv')
    assignment = pd.concat(assignments, ignore_index=True)
    write_aggregate(assignment, out / 'aggregate_candidates/P1_cluster_assignments.csv')
    write_aggregate(pd.concat(stability, ignore_index=True), out / 'aggregate_candidates/P1_cluster_stability.csv')
    write_aggregate(dependency_edges(definitions), out / 'aggregate_candidates/P1_feature_dependency_edges.csv')
    reference_groups = assignment.loc[assignment.source.eq('FairFace')]
    matrix = pd.read_csv(out / 'evidence/P0_exhaustive_evidence_matrix.csv', dtype={'canonical_syndrome_id': 'string'})
    matrix = matrix.merge(reference_groups[['feature', 'cluster_id', 'representative']], on='feature', validate='many_to_one')
    summaries = []
    for syndrome, part in matrix.groupby('syndrome', sort=True):
        summaries.append(dict(syndrome=syndrome, raw_features=120, reference_clusters=part.cluster_id.nunique(),
                              raw_primary_discoveries=int(part.ff_mannwhitney_q_fdr_global.lt(.05).sum()),
                              strict_representative_discoveries=int((part.representative & part.ff_mannwhitney_q_fdr_global.lt(.05)).sum()),
                              clusters_with_any_discovery=int(part.loc[part.ff_mannwhitney_q_fdr_global.lt(.05), 'cluster_id'].nunique()),
                              interpretation='descriptive correlated groups; original 5760-test BH unchanged; any-feature summary outcome-selected'))
    write_aggregate(canonicalize(pd.DataFrame(summaries)), out / 'aggregate_candidates/P1_cluster_discovery_summary.csv')
    subset = ['Cornelia de Lange syndrome', 'Noonan syndrome', 'WILLIAMS-BEUREN SYNDROME; WBS', 'Rubinstein-Taybi syndrome']
    signatures = matrix.loc[matrix.representative & matrix.syndrome.str.split('__').str[0].isin(subset)].copy()
    write_aggregate(signatures, out / 'aggregate_candidates/P1_signature_source.csv')
    comparison = []
    for source, labels in source_labels.items():
        comparison.append(dict(source=source, features=120, clusters=len(np.unique(labels)),
                               adjusted_rand_to_FairFace=adjusted_rand(source_labels['FairFace'], labels),
                               role='descriptive source/cohort partition; not independent biological dimensions'))
    write_aggregate(pd.DataFrame(comparison), out / 'aggregate_candidates/P1_cluster_source_comparison.csv')
    reference_n = len(np.unique(source_labels['FairFace']))
    gmdb_n = len(np.unique(source_labels['GMDB_selected']))
    report = ('# P1 feature dependence\n\n'
        f'All 120 features were correlated separately in FairFace, 1,845 selected GMDB patients and {len(cohorts) - 2} syndrome cohorts with n>=20. '
        f'Average linkage on 1−|Spearman rho|, cut .3, produces {reference_n} FairFace and {gmdb_n} GMDB groups. '
        f'Global source partition adjusted Rand agreement is {adjusted_rand(source_labels["FairFace"], source_labels["GMDB_selected"]):.6g}. '
        f'{args.bootstrap} independent row/patient bootstrap draws per global source assess co-clustering; FairFace path identity uniqueness remains unverified.\n\n'
        'Pearson matrices are supplied alongside Spearman. Source/cohort differences reflect both population and sampling differences. '
        'Clusters are descriptive and do not define independent biological dimensions or an effective hypothesis count. The original 5,760 BH family is unchanged.\n\n'
        'Dependency edges distinguish explicit f[feature] algebraic references from shared normalization: a common denominator is not itself an identity between outputs. '
        'Representatives are chosen without outcomes by non-asymmetry, fewer derived references, formula length and lexical ID. '
        'Signatures retain all such representatives for four predeclared adequately supported syndromes, including nonsignificant results; '
        'their source marginal confidence intervals are not multiplicity-adjusted. Cluster-any-discovery summaries are explicitly outcome-selected. '
        'Candidate rank sensitivity to representatives/families is evaluated separately and remains exploratory.\n')
    (out / 'P1_feature_dependence.md').write_text(report)
    protect_baseline(out)
    print(f'Dependence analyses complete: {reference_n} reference groups; {gmdb_n} GMDB groups')


if __name__ == '__main__':
    main()

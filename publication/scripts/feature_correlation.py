#!/usr/bin/env python3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, leaves_list, linkage
from scipy.spatial.distance import squareform

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "revision_results/final_20260914"
EVIDENCE = ROOT / "local_review/bibm2026_evidence"
OUTPUT = ROOT / "publication/feature_correlation"


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    features = pd.read_csv(RESULTS / "exact_feature_names.csv").feature.tolist()
    mapping = pd.read_csv(EVIDENCE / "feature_region_family_mapping.csv")
    reference = pd.read_csv(ROOT / "data/master_dataset_healthy.csv", low_memory=False)
    if len(features) != 125 or set(features) != set(mapping.feature): raise ValueError("Historical feature schema is not exactly 125 mapped features")
    if any(feature not in reference for feature in features): raise ValueError("Reference matrix does not contain the historical 125-feature schema")
    reference = reference[["file", *features]].copy()
    reference = reference.assign(reference_key=reference.file.astype("string").str.strip().str.replace("\\", "/", regex=False))
    reference = reference.drop_duplicates("reference_key")
    if len(reference) != 942: raise ValueError(f"Expected 942 distinct historical reference paths, found {len(reference)}")
    matrix = reference[features].apply(pd.to_numeric, errors="coerce")
    corr = matrix.corr(method="spearman").reindex(index=features, columns=features)
    if corr.isna().all(axis=0).any(): raise ValueError("At least one historical feature has no estimable reference correlation")
    corr.index.name = "feature"
    corr.reset_index().to_csv(OUTPUT / "feature_correlation_matrix.csv", index=False)

    values = corr.fillna(0).to_numpy()
    values = (values + values.T) / 2
    np.fill_diagonal(values, 1)
    distance = np.clip(1 - np.abs(values), 0, 1)
    np.fill_diagonal(distance, 0)
    tree = linkage(squareform(distance, checks=True), method="average")
    cluster_ids = fcluster(tree, t=.3, criterion="distance")
    order = leaves_list(tree)
    rows = pd.DataFrame({"feature": features, "cluster_id": cluster_ids}).merge(mapping[["feature", "anatomical_region", "feature_family"]], on="feature", validate="one_to_one")
    sizes = rows.groupby("cluster_id").size().rename("cluster_size")
    rows = rows.join(sizes, on="cluster_id")
    medians = {}
    for cluster_id, group in rows.groupby("cluster_id"):
        indices = [features.index(feature) for feature in group.feature]
        upper = np.abs(values[np.ix_(indices, indices)])[np.triu_indices(len(indices), 1)]
        medians[cluster_id] = float(np.median(upper)) if len(upper) else np.nan
    rows["median_absolute_within_cluster_rho"] = rows.cluster_id.map(medians)
    rows.sort_values(["cluster_size", "cluster_id", "feature"], ascending=[False, True, True]).to_csv(OUTPUT / "feature_clusters.csv", index=False)

    ordered = values[np.ix_(order, order)]
    labels = [features[i] for i in order]
    figure, ax = plt.subplots(figsize=(8.2, 7.4))
    image = ax.imshow(ordered, cmap="RdBu_r", vmin=-1, vmax=1, interpolation="nearest", rasterized=True)
    ticks = np.arange(0, 125, 10)
    ax.set_xticks(ticks, [labels[i].replace("_", " ") for i in ticks], rotation=55, ha="right", fontsize=6)
    ax.set_yticks(ticks, [labels[i].replace("_", " ") for i in ticks], fontsize=6)
    ax.set_title("Spearman correlation among 125 historical reference features", loc="left", fontsize=11)
    colorbar = figure.colorbar(image, ax=ax, pad=.02, shrink=.84)
    colorbar.set_label("Spearman rho")
    figure.tight_layout()
    figure.savefig(OUTPUT / "feature_correlation_heatmap.pdf", bbox_inches="tight")
    plt.close(figure)

    upper = np.abs(values[np.triu_indices(125, 1)])
    counts = {threshold: int((upper >= threshold).sum()) for threshold in [.7, .8, .9, .95]}
    mixed = rows.groupby("cluster_id").anatomical_region.nunique().gt(1).sum()
    report = f"""# Feature redundancy report

The reference-cohort correlation analysis uses 942 distinct FairFace image paths and the exact ordered 125-feature schema in `revision_results/final_20260914/exact_feature_names.csv`. No 120-feature output or regenerated FaceKit measurement is used.

## Results

- {counts[.7]:,} of 7,750 feature pairs have |Spearman rho| >= 0.70.
- {counts[.8]:,} pairs have |rho| >= 0.80, {counts[.9]:,} have |rho| >= 0.90, and {counts[.95]:,} have |rho| >= 0.95.
- Average-linkage clustering on `1 - |rho|`, cut at distance 0.30, yields {rows.cluster_id.nunique()} clusters. {int((sizes > 1).sum())} clusters contain more than one feature.
- {mixed} clusters span more than one anatomical region, showing that shared scale factors and geometry can link measurements that have different anatomical labels.

## Interpretation

The 125 tests within each syndrome are not biologically independent. Bilateral left/right measurements, mean and asymmetry derivatives, related distances, and features that share normalization denominators can move together. A cluster of significant results can therefore describe one underlying facial pattern rather than several independent discoveries. The matrix and cluster table characterize that dependence; they do not redefine the prespecified global Benjamini-Hochberg family and should not be used to retroactively reduce the multiplicity correction.

Spearman correlation captures monotonic dependence in the reference cohort only. It does not prove that two features use identical landmarks, and it can differ in a disease cohort. Exact landmark overlap would require a separate source-code-level dependency map of every feature formula.
"""
    (OUTPUT / "FEATURE_REDUNDANCY_REPORT.md").write_text(report)
    print(f"Generated reference correlation matrix and {rows.cluster_id.nunique()} descriptive clusters")


if __name__ == "__main__": main()

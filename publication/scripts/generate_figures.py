#!/usr/bin/env python3
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "revision_results/final_20260914"
EVIDENCE = ROOT / "local_review/bibm2026_evidence"
OUTPUT = ROOT / "publication/figures"
MANUSCRIPT = ROOT / "manuscript/figures"


def save(fig, stem):
    fig.savefig(OUTPUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUTPUT / f"{stem}.png", dpi=350, bbox_inches="tight")
    shutil.copy2(OUTPUT / f"{stem}.pdf", MANUSCRIPT / f"{stem}.pdf")
    plt.close(fig)


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    MANUSCRIPT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 10, "axes.labelsize": 9,
                         "xtick.labelsize": 7, "ytick.labelsize": 7, "pdf.fonttype": 42, "ps.fonttype": 42})

    fig1 = pd.read_csv(EVIDENCE / "figure_cross_syndrome_heatmap_source.csv")
    fig1.to_csv(OUTPUT / "fig1_source.csv", index=False)
    syndromes, features = list(dict.fromkeys(fig1.syndrome)), list(dict.fromkeys(fig1.feature))
    matrix = fig1.pivot(index="syndrome", columns="feature", values="hedges_g").reindex(index=syndromes, columns=features)
    limit = max(1, np.nanmax(np.abs(matrix.to_numpy())))
    figure, ax = plt.subplots(figsize=(7.2, 5.2))
    image = ax.imshow(matrix, cmap="RdBu_r", vmin=-limit, vmax=limit, aspect="auto", interpolation="nearest")
    ax.set_xticks(range(len(features)), [x.replace("_", " ") for x in features], rotation=45, ha="right")
    ax.set_yticks(range(len(syndromes)), [x.replace(";", ";\n", 1) for x in syndromes])
    ax.set_title("Cross-syndrome effect sizes for selected associations", loc="left")
    ax.set_xlabel("Feature")
    ax.set_ylabel("Syndrome")
    colorbar = figure.colorbar(image, ax=ax, pad=.02, shrink=.86)
    colorbar.set_label("Hedges' g")
    figure.tight_layout()
    save(figure, "fig1_cross_syndrome_heatmap")

    broad = pd.read_csv(RESULTS / "age_3_9/ethnicity_analysis/layer1_broad_ethnicity/analysis_A_healthy_kruskal_all_groups.csv.gz")
    detailed = pd.read_csv(RESULTS / "age_3_9/ethnicity_analysis/layer2_subcategory/analysis_A_healthy_kruskal_all_groups.csv.gz")
    fig2 = broad[["feature", "kruskal_H", "p_value", "fdr_q_value", "significant_fdr"]].merge(
        detailed[["feature", "kruskal_H", "p_value", "fdr_q_value", "significant_fdr"]], on="feature", suffixes=("_broad", "_detailed"), validate="one_to_one")
    fig2.to_csv(OUTPUT / "fig2_source.csv", index=False)
    x, y = -np.log10(fig2.fdr_q_value_broad), -np.log10(fig2.fdr_q_value_detailed)
    both = fig2.significant_fdr_broad.astype(bool) & fig2.significant_fdr_detailed.astype(bool)
    figure, ax = plt.subplots(figsize=(6.6, 4.3))
    ax.scatter(x[~both], y[~both], s=20, color="#8da0ae", alpha=.75, label="Not significant in both layers")
    ax.scatter(x[both], y[both], s=24, color="#2166ac", alpha=.85, label=f"Significant in both layers (n={both.sum()})")
    threshold = -np.log10(.05)
    ax.axvline(threshold, color="#555555", lw=.8, ls="--")
    ax.axhline(threshold, color="#555555", lw=.8, ls="--")
    ax.set(xlabel="Broad-group omnibus -log10(q)", ylabel="Detailed-group omnibus -log10(q)")
    ax.set_title("Reference-group variation across demographic layers", loc="left")
    ax.text(.98, .96, f"Broad: {fig2.significant_fdr_broad.astype(bool).sum()}/125\nDetailed: {fig2.significant_fdr_detailed.astype(bool).sum()}/125",
            transform=ax.transAxes, ha="right", va="top", fontsize=8)
    ax.legend(frameon=False, loc="upper left", fontsize=7)
    ax.spines[["top", "right"]].set_visible(False)
    figure.tight_layout()
    save(figure, "fig2_reference_variation")

    primary = pd.read_csv(RESULTS / "age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz")
    alternative = pd.read_csv(RESULTS / "age_3_9_id_sensitivity/significance_testing/stats_each_disease_vs_healthy.csv.gz")
    p = set(map(tuple, primary.loc[primary.mannwhitney_q_fdr_global.lt(.05), ["disease", "feature"]].to_numpy()))
    a = set(map(tuple, alternative.loc[alternative.mannwhitney_q_fdr_global.lt(.05), ["disease", "feature"]].to_numpy()))
    fig3 = pd.DataFrame({"category": ["Primary only", "Shared", "Alternative only"], "count": [len(p - a), len(p & a), len(a - p)]})
    fig3["primary_total"] = len(p)
    fig3["alternative_total"] = len(a)
    fig3["union"] = len(p | a)
    fig3["jaccard"] = len(p & a) / len(p | a)
    fig3.to_csv(OUTPUT / "fig3_source.csv", index=False)
    figure, ax = plt.subplots(figsize=(5.8, 3.5))
    bars = ax.bar(fig3.category, fig3["count"], color=["#b2182b", "#2166ac", "#ef8a62"], width=.62)
    ax.bar_label(bars, padding=3, fontsize=9)
    ax.set(ylabel="Globally significant comparisons", ylim=(0, max(fig3["count"]) * 1.15))
    ax.set_title("Photograph-selection sensitivity", loc="left")
    ax.text(.98, .92, f"Jaccard overlap: {len(p & a) / len(p | a):.1%}\nPrimary retained: {len(p & a) / len(p):.1%}", transform=ax.transAxes, ha="right", va="top")
    ax.spines[["top", "right"]].set_visible(False)
    figure.tight_layout()
    save(figure, "fig3_selection_overlap")
    print("Generated three source-data CSVs and PDF/PNG figure pairs")


if __name__ == "__main__": main()

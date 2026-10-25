#!/usr/bin/env python3
import argparse
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, FancyBboxPatch
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "revision_results/final_20260914"
EVIDENCE = ROOT / "local_review/bibm2026_evidence"
OUTPUT = ROOT / "publication/figures"
MANUSCRIPT = ROOT / "manuscript/figures"
LAYERS = {"broad": "layer1_broad_ethnicity", "detailed": "layer2_subcategory"}
BLUE, ORANGE, GRAY = "#2166ac", "#b65a18", "#777777"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def save(fig, stem, numbered=True):
    fig.savefig(OUTPUT / f"{stem}.pdf", bbox_inches="tight", metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(OUTPUT / f"{stem}.png", dpi=350, bbox_inches="tight")
    if numbered:
        shutil.copy2(OUTPUT / f"{stem}.pdf", MANUSCRIPT / f"{stem}.pdf")
    plt.close(fig)


def load_primary():
    primary = pd.read_csv(RESULTS / "age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz")
    schema = pd.read_csv(RESULTS / "exact_feature_names.csv")["feature"]
    require(len(schema) == 125 and set(primary.feature) == set(schema), "Historical feature schema changed")
    require(len(primary) == 6125 and primary.disease.nunique() == 49, "Primary test family changed")
    require(not primary.duplicated(["disease", "feature"]).any(), "Duplicate primary comparison")
    return primary


def figure1(primary):
    source = pd.read_csv(EVIDENCE / "figure_cross_syndrome_heatmap_source.csv")
    checked = source.merge(primary, left_on=["syndrome", "feature"], right_on=["disease", "feature"],
                           suffixes=("", "_frozen"), validate="one_to_one", indicator=True)
    require(checked._merge.eq("both").all(), "Unmatched heatmap cell")
    for column in ["hedges_g", "mannwhitney_q_fdr_global", "hedges_g_ci_95_low", "hedges_g_ci_95_high", "n_disease"]:
        require(np.allclose(checked[column], checked[f"{column}_frozen"], rtol=1e-12, atol=1e-14),
                f"Heatmap {column} differs from frozen results")
    sig = source.mannwhitney_q_fdr_global.lt(.05)
    require(source.cell_status.eq(np.where(sig, "significant", "tested_nonsignificant")).all(),
            "Incorrect heatmap cell status")
    require(len(source) == 120 and sig.sum() == 79, "Heatmap selection changed")
    source.to_csv(OUTPUT / "fig1_source.csv", index=False)
    syndromes, features = list(dict.fromkeys(source.syndrome)), list(dict.fromkeys(source.feature))
    matrix = source.pivot(index="syndrome", columns="feature", values="hedges_g").reindex(
        index=syndromes, columns=features)
    require(matrix.notna().all().all(), "Unexpected unavailable heatmap cell")
    limit = max(1, np.nanmax(np.abs(matrix.to_numpy())))
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    im = ax.pcolormesh(np.arange(len(features) + 1) - .5, np.arange(len(syndromes) + 1) - .5,
                       matrix, cmap="RdBu_r", vmin=-limit, vmax=limit, shading="flat", edgecolors="face", linewidth=.1, zorder=0)
    ax.set_ylim(len(syndromes) - .5, -.5)
    ax.set_xticks(range(len(features)), [x.replace("_", " ") for x in features], rotation=45, ha="right")
    ax.set_yticks(range(len(syndromes)), [x.replace(";", ";\n", 1) for x in syndromes])
    ax.set_title("Cross-syndrome effect sizes for selected associations", loc="left")
    ax.set(xlabel="Feature", ylabel="Syndrome")
    colorbar = fig.colorbar(im, ax=ax, pad=.02, shrink=.86)
    colorbar.set_label("Hedges' g (disease minus reference)")
    colorbar.solids.set_rasterized(False)
    colorbar.solids.set_edgecolor("face")
    fig.tight_layout(rect=(0, .085, 1, 1))
    fig.text(.5, .035, "Color shows Hedges' g for all 120 tested comparisons; significance is not encoded.\n"
             "No unavailable cells in this selection.", ha="center", fontsize=7)
    save(fig, "fig1_cross_syndrome_heatmap")
    print("Figure 1: 120 exact frozen matches; 79 significant, 41 tested nonsignificant, 0 unavailable")


def demographic_sources():
    omnibus, summaries = {}, []
    expected = {("A", "broad"): (72, 125), ("A", "detailed"): (89, 125),
                ("B", "broad"): (375, 36000), ("B", "detailed"): (452, 33875),
                ("C", "broad"): (479, 11500), ("C", "detailed"): (633, 15000)}
    files = {"A": "analysis_A_healthy_kruskal_all_groups.csv.gz",
             "B": "analysis_B_matched_disease_vs_healthy.csv.gz", "C": "analysis_C_all_interactions.csv.gz"}
    supplements = {"A": ("S2_reference_variation", "corrected_q"),
                   "B": ("S3_matched_comparisons", "global_q"), "C": ("S4_interactions", "global_q")}
    for analysis, filename in files.items():
        supplement, supplement_q = supplements[analysis]
        table = pd.read_csv(ROOT / f"publication/supplementary/Supplementary_Table_{supplement}.csv")
        for layer, directory in LAYERS.items():
            path = RESULTS / f"age_3_9/ethnicity_analysis/{directory}/{filename}"
            data = pd.read_csv(path)
            q_col = "fdr_q_value" if analysis == "A" else "global_fdr_q_value"
            flag = "significant_fdr" if analysis == "A" else "global_significant_fdr"
            significant = data[q_col].lt(.05)
            require(data[q_col].notna().all() and data[q_col].between(0, 1).all(), "Invalid demographic q")
            require(data[flag].equals(significant), "Significance flag differs from q < 0.05")
            count, total = int(significant.sum()), len(data)
            require((count, total) == expected[analysis, layer], f"Analysis {analysis}/{layer} changed")
            subset = table.loc[table.demographic_layer.eq(layer)]
            require(len(subset) == total and subset[supplement_q].lt(.05).sum() == count,
                    f"Supplement disagrees with Analysis {analysis}/{layer}")
            if analysis == "A":
                omnibus[layer] = data[["feature", "kruskal_H", "p_value", q_col, flag]]
            summaries.append({"analysis": analysis, "layer": layer, "significant": count, "tested": total,
                              "percent_significant": 100 * count / total, "q_column": q_col,
                              "bh_family": f"Analysis {analysis}, {layer} layer", "source": str(path.relative_to(ROOT))})
    source = omnibus["broad"].merge(omnibus["detailed"], on="feature", suffixes=("_broad", "_detailed"),
                                     validate="one_to_one")
    source.to_csv(OUTPUT / "fig2_source.csv", index=False)
    summary = pd.DataFrame(summaries)
    summary.to_csv(OUTPUT / "fig2_demographic_summary_source.csv", index=False)
    return source, summary


def figure2():
    source, summary = demographic_sources()
    fig = plt.figure(figsize=(7.2, 5.4))
    grid = fig.add_gridspec(2, 2, width_ratios=[1.25, 1], hspace=.86, wspace=.65)
    ax = fig.add_subplot(grid[:, 0])
    x, y = -np.log10(source.fdr_q_value_broad), -np.log10(source.fdr_q_value_detailed)
    broad, detailed = source.fdr_q_value_broad.lt(.05), source.fdr_q_value_detailed.lt(.05)
    categories = [(~broad & ~detailed, GRAY, "o", "Neither"), (broad & ~detailed, ORANGE, "^", "Broad only"),
                  (~broad & detailed, "#6c4991", "s", "Detailed only"), (broad & detailed, BLUE, "o", "Both")]
    handles = []
    for mask, color, marker, label in categories:
        ax.scatter(x[mask], y[mask], s=18, color=color, marker=marker, alpha=.8, linewidths=.3, edgecolors="white")
        handles.append(Line2D([], [], linestyle="", marker=marker, color=color, markersize=5,
                              label=f"{label} ({mask.sum()})"))
    threshold = -np.log10(.05)
    ax.axvline(threshold, color=GRAY, lw=.8, ls="--")
    ax.axhline(threshold, color=GRAY, lw=.8, ls="--")
    ax.set(xlabel="Broad omnibus −log₁₀(q)", ylabel="Detailed omnibus −log₁₀(q)", xlim=(-1, x.max() * 1.07),
           ylim=(-1, y.max() * 1.07))
    ax.set_title("A  Reference-only differences", loc="left", pad=23, fontsize=9, fontweight="bold")
    ax.text(0, 1.02, "Broad: 72/125; detailed: 89/125 features", transform=ax.transAxes, fontsize=7.2)
    fig.legend(handles=handles, title="Panel A: q < 0.05 in", frameon=False, loc="upper left",
               bbox_to_anchor=(.08, .925), ncol=4, fontsize=7, title_fontsize=7)
    ax.text(.98, .065, "Dashed lines: q = 0.05", transform=ax.transAxes, ha="right", fontsize=7)
    ax.spines[["top", "right"]].set_visible(False)
    for index, analysis, title in [(0, "B", "Matched syndrome–reference tests"),
                                    (1, "C", "Disease-by-group interactions")]:
        panel = fig.add_subplot(grid[index, 1])
        data = summary.loc[summary.analysis.eq(analysis)]
        bars = panel.barh([1, 0], data.percent_significant, color=[BLUE, ORANGE], height=.44)
        panel.set_yticks([1, 0], ["Broad", "Detailed"])
        panel.set(xlim=(0, 5), ylim=(-.55, 1.8), xlabel="Significant tests (%)")
        panel.set_xticks([0, 1, 2, 3, 4, 5])
        panel.set_title(f"{analysis}  {title}", loc="left", pad=12, fontsize=8.2, fontweight="bold")
        for bar, row in zip(bars, data.itertuples()):
            panel.text(.03, bar.get_y() + bar.get_height() + .09,
                       f"{row.significant:,}/{row.tested:,} ({row.percent_significant:.2f}%)", fontsize=7)
        panel.spines[["top", "right", "left"]].set_visible(False)
        panel.tick_params(axis="y", length=0)
    fig.subplots_adjust(left=.09, right=.98, bottom=.17, top=.75)
    fig.text(.09, .955, "Demographic analyses", fontsize=12, fontweight="bold")
    fig.text(.09, .04, "BH q < 0.05: A spans 125 features per layer; B and C each use a separate global family per layer.", fontsize=7)
    save(fig, "fig2_demographic_analyses")
    print(summary[["analysis", "layer", "significant", "tested", "percent_significant"]].to_string(index=False))


def figure3(primary):
    alternative = pd.read_csv(RESULTS / "age_3_9_id_sensitivity/significance_testing/stats_each_disease_vs_healthy.csv.gz")
    require(len(alternative) == 6125 and not alternative.duplicated(["disease", "feature"]).any(),
            "Alternative test family changed")
    require(set(map(tuple, primary[["disease", "feature"]].to_numpy())) ==
            set(map(tuple, alternative[["disease", "feature"]].to_numpy())), "Selections test different comparisons")
    p = set(map(tuple, primary.loc[primary.mannwhitney_q_fdr_global.lt(.05), ["disease", "feature"]].to_numpy()))
    a = set(map(tuple, alternative.loc[alternative.mannwhitney_q_fdr_global.lt(.05), ["disease", "feature"]].to_numpy()))
    counts = [len(p - a), len(p & a), len(a - p)]
    require(counts == [41, 770, 64] and (len(p), len(a), len(p | a)) == (811, 834, 875), "Overlap changed")
    source = pd.DataFrame({"category": ["Primary only", "Shared", "Alternative only"], "count": counts})
    for column, value in {"primary_total": len(p), "alternative_total": len(a), "union": len(p | a),
                          "jaccard": len(p & a) / len(p | a)}.items():
        source[column] = value
    source.to_csv(OUTPUT / "fig3_source.csv", index=False)
    photo = pd.read_csv(RESULTS / "selection_sensitivity_summary.csv").iloc[0]
    require((photo.patients, photo.selected_image_changed) == (1850, 94), "Photograph sensitivity changed")
    metrics = pd.DataFrame({"metric": ["primary_retained", "jaccard_overlap", "photographs_changed"],
                            "numerator": [len(p & a), len(p & a), photo.selected_image_changed],
                            "denominator": [len(p), len(p | a), photo.patients]})
    metrics["percent"] = 100 * metrics.numerator / metrics.denominator
    metrics.to_csv(OUTPUT / "fig3_metrics_source.csv", index=False)
    supplement = pd.read_csv(ROOT / "publication/supplementary/Supplementary_Table_S5_photo_sensitivity.csv")
    require(len(supplement) == len(primary) and supplement.primary_q.lt(.05).sum() == len(p) and
            supplement.alternative_q.lt(.05).sum() == len(a) and
            (supplement.primary_q.lt(.05) & supplement.alternative_q.lt(.05)).sum() == len(p & a),
            "S5 differs from overlap")
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    left = 0
    for count, color in zip(counts, [ORANGE, BLUE, "#6c4991"]):
        ax.barh(0, count, left=left, height=.42, color=color, edgecolor="white", linewidth=1.3)
        ax.text(left + count / 2, 0, str(count), ha="center", va="center", color="white", fontweight="bold")
        left += count
    for center, label in [(20.5, "Primary only"), (426, "Significant under both"), (843, "Alternative only")]:
        ax.annotate(label, xy=(center, .23), xytext=(center, .56), ha="center", fontsize=8,
                    arrowprops={"arrowstyle": "-", "color": GRAY, "lw": .7})
    ax.set(xlim=(-45, 920), ylim=(-.42, .87), yticks=[], xlabel="Significant syndrome–feature comparisons in union (875)")
    ax.set_xticks([0, 200, 400, 600, 800, 875])
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.set_title("Photograph-selection sensitivity", loc="left", pad=15, fontsize=11, fontweight="bold")
    fig.subplots_adjust(left=.065, right=.96, top=.83, bottom=.55)
    fig.text(.065, .34, f"Primary: {len(p)}     Alternative: {len(a)}     Shared: {len(p & a)}", fontsize=9)
    fig.text(.065, .23, f"Primary retained: {len(p & a) / len(p):.1%} (770/811)     "
             f"Jaccard overlap: {len(p & a) / len(p | a):.1%} (770/875)", fontsize=9)
    fig.text(.065, .12, f"Selected photographs changed: {int(photo.selected_image_changed)}/{int(photo.patients):,} "
             f"({photo.selected_image_changed / photo.patients:.1%}).", fontsize=8)
    save(fig, "fig3_selection_overlap")
    print("Figure 3: 41 / 770 / 64; 94.9% retained; 88.0% Jaccard; 94/1850 photographs changed")


def face(ax):
    ax.add_patch(Ellipse((0, 0), 1.65, 2.25, fill=False, edgecolor="#687581", lw=1.2))
    for x in [-.32, .32]:
        ax.add_patch(Ellipse((x, .28), .38, .15, fill=False, edgecolor="#87929a", lw=.8))
        ax.add_patch(Ellipse((x, .28), .075, .075, facecolor="#87929a", edgecolor="none"))
        ax.plot([x - .2, x, x + .18], [.52, .57, .52], color="#87929a", lw=2)
    ax.plot([-.08, -.13, .13, .08], [.25, -.15, -.15, .25], color="#87929a", lw=.8)
    ax.plot([-.27, -.08, 0, .08, .27], [-.45, -.42, -.45, -.42, -.45], color="#87929a", lw=.8)
    ax.plot([-.27, 0, .27], [-.45, -.55, -.45], color="#87929a", lw=.8)
    ax.set(xlim=(-1.05, 1.05), ylim=(-1.25, 1.25), aspect="equal")
    ax.axis("off")


def schematics(primary):
    formulas = pd.read_csv(ROOT / "publication/feature_definitions/historical_feature_formulas.csv").set_index("feature")
    features = ["inter_pupillary_distance", "philtrum_length", "philtrum_width"]
    require(formulas.loc[features, "extraction_status"].eq("EXACT_FROM_SOURCE").all(), "Missing exact formulas")
    formulas.loc[features].to_csv(OUTPUT / "schematic_workflow_formula_source.csv")
    fig = plt.figure(figsize=(7.2, 5.5))
    fig.text(.04, .95, "FaceKit extraction workflow | schematic", fontsize=12, fontweight="bold")
    flow = fig.add_axes([.02, .62, .96, .25])
    flow.axis("off")
    boxes = [("Facial image", "Input"), ("MediaPipe Face\nLandmarker", "Pretrained model"),
             ("478\nlandmarks", "Coordinates"), ("Historical FaceKit\nalignment", "Geometric step"),
             ("125 geometric\nmeasurements", "Historical schema")]
    for i, (title, subtitle) in enumerate(boxes):
        x = i * .202
        flow.add_patch(FancyBboxPatch((x + .005, .25), .175, .6, boxstyle="round,pad=0.005,rounding_size=0.025",
                                     facecolor="#eaf0f6" if i == 1 else "#f3f4f5", edgecolor="#98a4ae"))
        flow.text(x + .092, .60, title, ha="center", va="center", fontsize=7.5)
        flow.text(x + .092, .08, subtitle, ha="center", fontsize=6.7, color=GRAY)
        if i < 4:
            flow.annotate("", xy=(x + .205, .55), xytext=(x + .182, .55), arrowprops={"arrowstyle": "->", "lw": .9})
    flow.set_xlim(0, 1.01)
    points = [((-.32, .28), (.32, .28)), ((-.78, -.02), (.78, -.02)), ((0, -.15), (0, -.42))]
    titles = ["Interpupillary distance", "Bizygomatic width", "Philtrum length"]
    notes = ["468 ↔ 473; divided by\nbizygomatic width", "234 ↔ 454; a denominator\nfor many historical features",
             "2 ↔ 0; divided by\nmidface height = |y₁₆₈ − y₂|"]
    for i, (pair, title, note) in enumerate(zip(points, titles, notes)):
        ax = fig.add_axes([.045 + i * .325, .20, .26, .37])
        face(ax)
        ax.annotate("", xy=pair[1], xytext=pair[0], arrowprops={"arrowstyle": "<->", "color": BLUE, "lw": 1.4})
        ax.scatter(*np.array(pair).T, color=BLUE, s=8)
        ax.set_title(title, fontsize=8.5, pad=7)
        ax.text(.5, -.07, note, transform=ax.transAxes, ha="center", va="top", fontsize=7)
    fig.text(.04, .06, "Alignment centers on inner canthi 133/362 and makes their connecting line horizontal (in-plane roll).\n"
             "All drawn coordinates are schematic, not measured faces.", fontsize=7)
    save(fig, "schematic_A_facekit_workflow", numbered=False)
    mapping = pd.read_csv(EVIDENCE / "feature_region_family_mapping.csv")
    require(len(mapping) == 125 and mapping.feature.nunique() == 125, "Incomplete region mapping")
    regional = primary.merge(mapping[["feature", "anatomical_region"]], on="feature", validate="many_to_one")
    regional["significant"] = regional.mannwhitney_q_fdr_global.lt(.05)
    regional = regional.groupby("anatomical_region").agg(tested=("feature", "size"), significant=("significant", "sum"))
    regional["percent_significant"] = 100 * regional.significant / regional.tested
    require(regional.tested.sum() == 6125 and regional.significant.sum() == 811, "Region totals differ")
    for name, counts in {"forehead": (74, 245), "chin / jaw": (66, 245),
                         "global face shape / proportions": (50, 294)}.items():
        require(tuple(regional.loc[name, ["significant", "tested"]]) == counts, f"Region {name} changed")
    prior = pd.read_csv(EVIDENCE / "regional_trends.csv").set_index("anatomical_region").reindex(regional.index)
    require(np.array_equal(regional[["tested", "significant"]].to_numpy(),
                           prior[["total_tested_comparisons", "number_significant"]].to_numpy()),
            "Frozen region trends differ from recomputed region counts")
    regional = regional.sort_values("percent_significant", ascending=False)
    regional.to_csv(OUTPUT / "schematic_region_results_source.csv")
    fig = plt.figure(figsize=(7.2, 4.7))
    fig.text(.04, .94, "Historical feature regions", fontsize=12, fontweight="bold")
    ax = fig.add_axes([.02, .20, .39, .61])
    face(ax)
    labels = [("Forehead", (0, .81), (-.98, 1.16)), ("Brows", (-.32, .55), (-1, .70)),
              ("Eyes / iris", (.32, .28), (.50, .92)), ("Nose", (.02, .02), (.60, .43)),
              ("Philtrum /\nmouth", (0, -.40), (.58, -.40)),
              ("Cheeks / malar", (-.50, -.20), (-1.1, -.80)), ("Jaw / chin", (0, -.98), (.15, -1.23))]
    for label, point, target in labels:
        ax.annotate(label, xy=point, xytext=target, fontsize=7, ha="center",
                    arrowprops={"arrowstyle": "-", "lw": .7, "color": BLUE}, color=BLUE)
    ax.set_title("A  Region map | schematic", fontsize=8.5, loc="left")
    ax.set_xlim(-1.5, 1.45)
    ax = fig.add_axes([.66, .24, .30, .56])
    short = {"global face shape / proportions": "Global proportions", "philtrum / mouth / lips": "Philtrum / mouth",
             "ocular position / iris / gaze": "Iris / gaze", "eyes / orbital dimensions": "Eyes / orbital",
             "cheeks / malar region": "Cheeks / malar"}
    positions = np.arange(len(regional))
    ax.barh(positions, regional.percent_significant, color=BLUE, height=.63)
    ax.set_yticks(positions, [short.get(x, x.capitalize()) for x in regional.index])
    ax.invert_yaxis()
    ax.set(xlim=(0, 39), xlabel="Significant comparisons (%)")
    ax.set_xticks([0, 10, 20, 30])
    for i, row in enumerate(regional.itertuples()):
        ax.text(row.percent_significant + .6, i, f"{row.significant}/{row.tested}", va="center", fontsize=6.5)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    fig.text(.49, .835, "B  Region summary | measured results", fontsize=8.5)
    fig.text(.04, .06, "Regions describe the historical feature mapping, not inferred biological localization.\n"
             "Result panel: significant/tested comparisons, pooled global BH q < 0.05 across 6,125 tests.", fontsize=7)
    save(fig, "schematic_B_anatomical_regions", numbered=False)
    print("Schematics: exact historical formula references; all 125 features mapped; region totals 811/6125")


def main():
    parser = argparse.ArgumentParser(description="Generate publication figures from frozen historical aggregate outputs")
    parser.add_argument("--figures", nargs="+", choices=["1", "2", "3", "schematics"], default=["1", "2", "3", "schematics"])
    args = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    MANUSCRIPT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 10, "axes.labelsize": 8,
                         "xtick.labelsize": 7, "ytick.labelsize": 7, "pdf.fonttype": 42, "ps.fonttype": 42})
    primary = load_primary()
    for name in args.figures:
        if name == "1":
            figure1(primary)
        elif name == "2":
            figure2()
        elif name == "3":
            figure3(primary)
        else:
            schematics(primary)


if __name__ == "__main__":
    main()

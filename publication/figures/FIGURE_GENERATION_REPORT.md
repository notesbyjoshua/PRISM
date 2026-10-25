# Publication figure completion report

Completed October 18, 2026. All graphics use deterministic Matplotlib drawing and the frozen historical 125-feature results. No image-generation models, patient images, landmark extraction, new inferential analysis, manuscript edits, Google Docs edits, or table generation were used.

## Delivered files and reproduction

From the repository root:

```bash
.venv/bin/python publication/scripts/generate_figures.py
```

Individual outputs can be regenerated with `--figures 1`, `--figures 2`, `--figures 3`, or `--figures schematics`; multiple selections are accepted. Required dependencies remain NumPy, Pandas, and Matplotlib. Each graphic has a vector PDF with embedded TrueType text and a 350-dpi PNG in this directory. Numbered PDFs are copied into `manuscript/figures/`.

| Graphic | PDF/PNG stem | Source extract |
|---|---|---|
| Figure 1 | `fig1_cross_syndrome_heatmap` | `fig1_source.csv` (unchanged) |
| Figure 2 | `fig2_demographic_analyses` | `fig2_source.csv` (unchanged, panel A); `fig2_demographic_summary_source.csv` (all families) |
| Figure 3 | `fig3_selection_overlap` | `fig3_source.csv` (unchanged); `fig3_metrics_source.csv` |
| Schematic A | `schematic_A_facekit_workflow` | `schematic_workflow_formula_source.csv` |
| Schematic B | `schematic_B_anatomical_regions` | `schematic_region_results_source.csv` |

Changed existing files: `publication/scripts/generate_figures.py`; Figure 1 and Figure 3 PDF/PNG pairs; their two manuscript PDF copies; and `publication/README.md`. New files: this report; Figure 2 A/B/C PDF/PNG pair and manuscript PDF copy; two schematic PDF/PNG pairs; and four additional source CSVs listed above. No frozen input or supplementary table was modified.

Figure 1 retains the feature order, syndrome order, diverging effect-size scale, and all 120 effect sizes. At the author's request, the black significance dots were removed; the footer now states that significance is not encoded. No cell is unavailable in this selected grid. The Google Doc still contains its previously embedded version.

Figure 2 uses a new name because `manuscript/main.tex` has a reference-only caption for `fig2_reference_variation.pdf`. The old reference-only PDF/PNG and manuscript copy remain unchanged as legacy artifacts; the updated generator produces the new A/B/C graphic. The new numbered PDF is synchronized to `manuscript/figures/fig2_demographic_analyses.pdf`. Switching the LaTeX includegraphics path and caption together remains an author integration step; no manuscript reference or caption was silently changed.

Schematic C was omitted because the optional region map and workflow already illustrate measurement locations. A further syndrome-labeled schematic would substantially duplicate them, and the mean malar-bulge feature is a contour-derived quantity rather than a single cheek-width arrow. No disease-specific face shapes were drawn.

## Upstream sources

Paths below are relative to the repository root.

- All pooled effects: `revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz`; historical schema: `revision_results/final_20260914/exact_feature_names.csv`.
- Figure 1 selection: `local_review/bibm2026_evidence/figure_cross_syndrome_heatmap_source.csv`. Every selected g, global q, confidence bound, and disease sample size is joined back to the frozen pooled table and checked. The source CSV is retained without changing its contents.
- Figure 2: `revision_results/final_20260914/age_3_9/ethnicity_analysis/{layer1_broad_ethnicity,layer2_subcategory}/analysis_A_healthy_kruskal_all_groups.csv.gz`, `analysis_B_matched_disease_vs_healthy.csv.gz`, and `analysis_C_all_interactions.csv.gz`. The summary extract records each exact path, q-value column, and BH family. A uses `fdr_q_value`; B and C use `global_fdr_q_value`. Saved significance flags must agree with q < 0.05. All six totals are cross-checked against Supplementary Tables S2-S4.
- Figure 3: primary pooled table above and `revision_results/final_20260914/age_3_9_id_sensitivity/significance_testing/stats_each_disease_vs_healthy.csv.gz`; photograph-change counts: `revision_results/final_20260914/selection_sensitivity_summary.csv`. Shared/exclusive sets use `(disease, feature)` keys and global q < 0.05. The primary and alternative test universes must match; aggregate overlap is cross-checked against Supplementary Table S5. No patient rows are exported.
- Workflow: `publication/feature_definitions/historical_feature_formulas.csv`, its dependency map, and the historical source at FaceKit revision `e474628f3cb17a83edad59ac6e1d9e5438753032`, `src/facekit/core/geometric/extractor.py`. The 478-landmark Tasks image entry point and inner-canthal alignment are confirmed in that historical source. Source compatibility does not establish the exact model asset or dependency versions used for every historical record. This provenance limit remains documented here; the author requested removal of that sentence from the diagram. Shown coordinates are hand-specified schematic positions, never extracted image coordinates. Formula source extract also includes philtrum width for normalization context, although that measurement is not drawn.
- Region result panel: pooled table joined to all 125 rows of `local_review/bibm2026_evidence/feature_region_family_mapping.csv`. All nine region counts are independently recomputed and checked against `local_review/bibm2026_evidence/regional_trends.csv`. The generic outline is schematic; only the separately labeled bar panel displays results.

## Numerical checks passed

- Historical schema: 125 exact features; 49 syndromes; 6,125 unique pooled tests.
- Heatmap: 120 exact source matches; 79 significant and 41 tested nonsignificant cells; zero unavailable cells. Hedges' g uses the disease-minus-reference sign convention and symmetric limits around zero.
- Analysis A: broad 72/125 (57.6%); detailed 89/125 (71.2%). Across layers: 35 neither, 1 broad only, 18 detailed only, 71 both.
- Analysis B: broad 375/36,000 = 1.041667%; detailed 452/33,875 = 1.334317%.
- Analysis C: broad 479/11,500 = 4.165217%; detailed 633/15,000 = 4.220000%. These are formal interaction tests; B/C bars share a 0-5% axis, with no inferential comparison between families.
- Selection: 41 primary-only + 770 shared + 64 alternative-only = 875 union; primary 811; alternative 834. Retained 770/811 = 94.944513%; Jaccard 770/875 = 88.0%; changed photographs 94/1,850 = 5.081081%.
- Regions: all nine sum to 811/6,125; forehead 74/245; chin/jaw 66/245; global proportions 50/294. Mapping is complete and unique at the feature level.
- Measurements illustrated in the workflow: interpupillary distance 468-473 divided by 234-454 distance; bizygomatic width 234-454; philtrum length 2-0 divided by midface height, the absolute aligned y difference at 168 and 2. Source formulas retain their original denominator rules; arrows depict locations, not empirical magnitudes.

## Suggested captions

**Figure 1.** Hedges' g for 120 selected syndrome-feature comparisons from the historical 125-feature analysis. Positive values indicate higher disease-cohort measurements than reference measurements. Color encodes effect size rather than statistical significance; 79 displayed comparisons meet global BH q < 0.05 across 6,125 pooled tests, and 41 do not. No displayed cell is unavailable. Feature values are normalized according to their historical definitions; the selected grid does not establish diagnostic specificity.

**Figure 2.** Demographic analyses in broad and detailed annotation layers. (A) Reference-only omnibus q-values for 125 features per layer; dashed lines indicate q = 0.05 and symbols distinguish significance in neither, one, or both layers. (B) Significant matched syndrome-reference tests: 375/36,000 broad and 452/33,875 detailed. (C) Significant formal disease-by-group interaction tests: 479/11,500 broad and 633/15,000 detailed. Bars show percentages with numerators and denominators. Each analysis and layer retains its separate BH family. Layers reuse observations, and counts are not additive or independent replications.

**Figure 3.** Overlap of pooled syndrome-feature comparisons meeting global BH q < 0.05 under two within-patient photograph-selection rules. The proportional segments show 41 primary-only, 770 shared, and 64 alternative-only comparisons. Of 811 primary associations, 94.9% were retained; Jaccard overlap was 88.0%. The selected photograph changed for 94/1,850 patients (5.1%). This is a within-patient photograph-selection sensitivity check, not independent validation or replication.

**Schematic A.** Code-drawn illustration of the historical FaceKit workflow, separating MediaPipe's pretrained landmark detector from subsequent alignment and geometric calculations. Examples indicate landmark-defined measurement locations and denominators. Inner-canthal alignment corrects in-plane roll; it does not perform full 3D frontalization. Coordinates and face outlines are schematic. The historical image entry point expects 478 landmarks, but the exact detector asset used for each stored record is unverified.

**Schematic B.** (A) Generic schematic showing regions used to organize historical measurement families; global proportions span the outline. (B) Measured proportions of globally significant pooled comparisons within all nine mapped regions, with significant/tested counts. Region labels describe the feature mapping rather than biological localization, and correlated measurements are not independent discoveries.

## Execution and visual review

Executed the figure generator, targeted Figure 2 regeneration after layout adjustments, Python syntax compilation, and a second complete generation for determinism. No full reproduction run was needed: it would regenerate tables, audits, covariate models, provenance, and other outputs outside this task.

Rendered every delivered PDF through Poppler `pdftoppm -scale-to 1600 -png -singlefile` and visually inspected the resulting pages. Also inspected all five delivered PNGs. Corrected a scatter legend intersecting data and tight Figure 3 spacing, then re-rendered and inspected the affected outputs. Final checks cover labels, arrows, bar counts, zero-centered color scale, missing/significance states, legends, margins, mathematical glyphs, and source-versus-schematic labeling. PDF inspection confirms one page each, embedded TrueType fonts, extractable text, and no raster image objects. PNG inspection confirms approximately 350 dpi. Repeated generation produces identical PDF, PNG, and source-CSV SHA-256 hashes in this environment; synchronized manuscript PDFs match their publication originals.

## Author review still needed

- Replace the embedded Google Doc Figure 1 if adopting this revised export; no Google Doc was accessed or edited.
- Adopt the new Figure 2 file and expanded caption together in the manuscript. The legacy reference-only file remains available.
- Review journal sizing and the supplied captions at the final typeset width. Keep the explicit sensitivity and schematic qualifications.
- Public redistribution remains subject to the existing data-provider and institutional approvals. Only aggregate outputs were used, and no public-release decision was made.

## Requested figure edits

Removed the black dots and their legend from Figure 1, retaining all effect sizes and source checks. Removed the sentences beginning “No complete” and “478-point” from Schematic A. The schematic-coordinate label and alignment explanation remain. Regenerated both PDF/PNG pairs and synchronized the Figure 1 manuscript PDF.

Completed the remaining requested annotation removals: Schematic B no longer includes the sentences beginning “Global proportions” and “Features are correlated”; Figure 2 no longer includes “Layers reuse observations” or “Counts across”; Figure 3 no longer includes the sentence beginning “Within-patient”. Regenerated all five PDF/PNG pairs, visually inspected the PNGs and Poppler PDF renders, checked that the requested text was absent, and verified the three numbered manuscript PDF copies match. Numerical source checks still pass. These removals apply to figure artwork; the scientific context remains documented in the suggested captions and source notes above.

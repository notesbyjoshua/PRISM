# Conditional PDIDB–GMDB consistency results

This secondary analysis uses frozen current measurements without reprocessing images or refitting the primary study. The mentor reports that StyleGAN was trained on GMDB and PDIDB sampled that learned distribution. This statement is recorded separately from provider metadata: the identities, training membership and relationship of v1/v2 checkpoints to the published generator remain unverified. Synthetic observations do not increase the real-patient n or independently confirm the primary 955 associations.

## Coverage and methods

The actual files contain 1,727 synthetic records and 1,563 QC-passing rows. Nine exact-name matches support comparisons of 314 selected GMDB patients and 1424 synthetic images on 120 features each (1,080 pooled syndrome-feature cells). The 139 QC-passing synthetic 22q11.2 rows remain excluded because the generic label cannot establish equivalence to DiGeorge or distal deletion labels. All accepted comparisons have at least ten observations per source; available selected real n ranges from 16 to 67. Age_range, gender and race are missing for all 1,727 synthetic rows. No age-matched or demographic-adjusted synthetic sensitivity is possible, while real GMDB ages are 3 <= age < 10. The fixed reference is the same 3,409 FairFace paths in every unmodified profile.

The analysis plan was written after metadata inspection and before fidelity computation; it was not registered externally. It defines sample SD (ddof=1), synthetic-minus-real shifts, real-SD scaling, empirical Wasserstein distance, sample variance ratios and Hedges’ g against the fixed reference. The second profile scaling divides both mean contrasts by the identical reference SD. Winsorized sensitivity clips each source separately at its 1st/99th percentiles, including a once-fixed clipped reference. No p-values, inferential intervals, equivalence margins or generator-level uncertainty are computed.

## Main descriptive findings

Across the nine equally weighted syndrome summaries, median absolute standardized mean shift is 0.166 real SD, median variance ratio is 0.986, median effect-direction agreement is 85.0%, and median Pearson/Spearman profile correlation is 0.896/0.890. Across all 1,080 cells (equal 120-feature contribution per syndrome), 615 have lower synthetic sample variance. These medians summarize correlated measurements, not independent biological replications.

| syndrome | real_n | synthetic_n | median_absolute_standardized_shift | median_variance_ratio | fraction_synthetic_variance_lower | direction_agreement | pearson_g | spearman_g |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ANGELMAN SYNDROME; AS | 18 | 144 | 0.278958 | 1.06315 | 0.408333 | 0.783333 | 0.752653 | 0.645857 |
| Cornelia de Lange syndrome | 64 | 140 | 0.161994 | 0.838823 | 0.783333 | 0.758333 | 0.902967 | 0.89034 |
| Kabuki syndrome | 35 | 142 | 0.158046 | 0.896934 | 0.591667 | 0.858333 | 0.896026 | 0.872748 |
| KBG SYNDROME; KBGS | 32 | 173 | 0.21092 | 0.784903 | 0.75 | 0.716667 | 0.824132 | 0.786444 |
| NICOLAIDES-BARAITSER SYNDROME; NCBRS | 24 | 141 | 0.164809 | 1.0366 | 0.466667 | 0.891667 | 0.9171 | 0.90691 |
| Noonan syndrome | 34 | 178 | 0.12695 | 0.986306 | 0.533333 | 0.916667 | 0.961114 | 0.961754 |
| Rubinstein-Taybi syndrome | 24 | 167 | 0.166425 | 1.02383 | 0.458333 | 0.833333 | 0.893278 | 0.891895 |
| SMITH-MAGENIS SYNDROME; SMS | 16 | 156 | 0.22674 | 0.855689 | 0.625 | 0.85 | 0.859333 | 0.847288 |
| WILLIAMS-BEUREN SYNDROME; WBS | 67 | 183 | 0.201522 | 0.992331 | 0.508333 | 0.941667 | 0.96489 | 0.95999 |

Cornelia de Lange has a median variance ratio of 0.839 and lower synthetic variance in 78.3% of features; KBG has a median ratio of 0.785. Angelman, Nicolaides–Baraitser and Rubinstein–Taybi have syndrome-level median ratios above one. Thus compression is heterogeneous, and a pooled near-one dispersion summary does not establish fidelity for each syndrome or measurement.

## Five anchors, retained regardless of agreement

| syndrome | feature | real_g | synthetic_g | standardized_mean_shift | variance_ratio |
| --- | --- | --- | --- | --- | --- |
| Cornelia de Lange syndrome | philtrum_length | 1.16283 | 1.20899 | 0.0462958 | 1.07872 |
| Noonan syndrome | inter_pupillary_distance | 1.24369 | 1.15199 | -0.0523482 | 0.925288 |
| Rubinstein-Taybi syndrome | malar_bulge_mean | -1.00903 | -0.405529 | 0.47104 | 1.09548 |
| WILLIAMS-BEUREN SYNDROME; WBS | eb_thickness_mean | 1.26614 | 0.696375 | -0.448839 | 0.995194 |
| WILLIAMS-BEUREN SYNDROME; WBS | outer_canthal_distance | 1.03943 | 0.547738 | -0.434031 | 1.02135 |

Philtrum and pupil-spacing effects are close in sign and magnitude. The Williams–Beuren eyebrow-area and outer-canthal effects and the Rubinstein–Taybi malar effect are attenuated in these synthetic samples, although all five retain the primary effect sign. The malar index is a projected cheek/outline construct and the eyebrow index is polygon area; neither is a direct diagnosis of a clinical trait. No significance criterion was used to choose these anchors.

## Largest discrepancies (exploratory ranking)

| syndrome | feature | standardized_mean_shift | variance_ratio | real_n | synthetic_n |
| --- | --- | --- | --- | --- | --- |
| KBG SYNDROME; KBGS | chin_height | 1.24797 | 2.1095 | 32 | 173 |
| ANGELMAN SYNDROME; AS | inter_canthal_distance | -1.03328 | 1.1929 | 18 | 144 |
| ANGELMAN SYNDROME; AS | malar_bulge_mean | 1.02381 | 1.29266 | 18 | 144 |
| ANGELMAN SYNDROME; AS | face_triangularity | -0.977767 | 1.71407 | 18 | 144 |
| ANGELMAN SYNDROME; AS | chin_width | 0.952271 | 1.64427 | 18 | 144 |

## Sensitivities

| sensitivity | selection_status | syndromes | median_absolute_standardized_shift | median_variance_ratio | median_direction_agreement | median_pearson_g |
| --- | --- | --- | --- | --- | --- | --- |
| all_supported | predeclared support | 9 | 0.166425 | 0.986306 | 0.85 | 0.896026 |
| real_n_at_least_20 | exploratory cutoff after initial pooled results | 7 | 0.164809 | 0.986306 | 0.858333 | 0.902967 |
| real_n_at_least_30 | exploratory cutoff after initial pooled results | 5 | 0.161994 | 0.896934 | 0.858333 | 0.902967 |

| model_tag | sensitivity | syndromes | median_shift | median_variance_ratio | median_direction_agreement | median_pearson_g |
| --- | --- | --- | --- | --- | --- | --- |
| pooled | unmodified | 9 | 0.166425 | 0.986306 | 0.85 | 0.896026 |
| pooled | winsorized_01_99 | 9 | 0.175273 | 0.983217 | 0.841667 | 0.892617 |
| v1 | unmodified | 9 | 0.192788 | 0.915234 | 0.833333 | 0.8613 |
| v2 | unmodified | 9 | 0.202423 | 1.00138 | 0.858333 | 0.893007 |

Model-tag strata describe version-labeled release subsets, not independently trained generators. Every supported tag result is supplied in the feature and syndrome tables. Fixed-reference-SD profile correlations and contrasts are supplied alongside pooled-SD Hedges’ g. The frozen feature-correlation clusters have separate feature-within-cluster and equal-cluster-weighted summaries; they do not modify inferential denominators. The n >= 20 and n >= 30 real-cohort cutoffs were added after initial pooled results and are explicitly exploratory. Low-n exclusions and winsorization change descriptive summaries but do not resolve age/source mismatch.

| size_proxy | metric | syndromes | descriptive_spearman |
| --- | --- | --- | --- |
| age-restricted available selected GMDB n; NOT training n | median_absolute_standardized_shift | 9 | -0.577411 |
| age-restricted available selected GMDB n; NOT training n | median_variance_ratio | 9 | -0.309626 |
| age-restricted available selected GMDB n; NOT training n | direction_agreement | 9 | 0.276153 |
| age-restricted available selected GMDB n; NOT training n | pearson_g | 9 | 0.67783 |

The size variable here is the available age-restricted GMDB selection, not the actual number of StyleGAN training patients/images. These nine descriptive correlations cannot establish a training-size effect. No exact duplicate encoded-byte groups or exact synthetic/selected-real byte overlaps were found; that check cannot establish patient independence, absence of memorization or privacy.

## Relationship to FaceKit and limitations

[FaceKit §§2.6/3.6 and Fig. 9](https://arxiv.org/html/2610.09288v1#S3.SS6) already assess synthetic phenotype preservation. Its 4,559 real/synthetic images and HPO-mapped contrasts differ from this age-restricted selection and fixed FairFace baseline. Its reported median Pearson correlation (0.866) and variance ratio (0.94) are not directly comparable benchmarks. Exact generator-release equivalence is unverified. This audit extends a question already studied in FaceKit; it does not establish a new extraction method, independent replication or general generator validity.

All images, per-image measurements, hashes and sensitive linkage records remain local. The output directory contains only syndrome/feature statistics, source-level category counts, code/documentation and original aggregate plots. Existing PI publication review still applies; this pass does not perform privacy/IRB/legal release review.

## Figures and reproducibility

Figures 8–10 in manuscript/figures show only predeclared-anchor quartiles, aggregate effect profiles and syndrome-level fidelity. Captions and alt text are in publication/final_research_120/FIGURE_CAPTIONS.md. CSVs reproduce the plotted values exactly. Run stats_testing/synthetic_consistency.py with --active-config ACTIVE_ANALYSIS.json and a private --output-dir, followed by publication/scripts/final_research_figures.py and publication/scripts/final_research_text.py. Source hashes and the mentor/provider provenance distinction are in synthetic_consistency_manifest.json.

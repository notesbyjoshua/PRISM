# Final research figure captions and alt text

## fig5_facekit_prism_workflow

FaceKit measurement platform and PRISM study design. FaceKit provides landmark detection, pose correction and 120 measurements; PRISM supplies cohort selection and the statistical branches. Synthetic fidelity is separate and conditional on a GMDB-trained generator (mentor statement; checkpoint-specific provenance unverified). FaceKit: arXiv:2610.09288v1, §§2.1–2.2. This original schematic contains no patient images.

Alt text: Flow diagram separates the FaceKit measurement steps from PRISM real-patient testing and the synthetic-only fidelity branch.

## fig6_anatomical_definitions

Original anatomical schematic for five FaceKit indices. A: subnasale-to-upper-lip distance divided by vertical nasion-to-subnasale midface height. B: iris-center distance divided by bizygomatic width. C: outer-canthal distance divided by bizygomatic width. D: bilateral mean star-polygon eyebrow area divided by squared bizygomatic width; this measures geometry, not hair density. E: bilateral signed horizontal cheek-to-interpolated-face-oval distance divided by bizygomatic width, not three-dimensional malar projection. Blue indicates numerator geometry, orange the length scale. Points are illustrative artwork; exact landmark sets, signs and branches are defined in S6 and the pinned FaceKit source.

Alt text: Five original face-outline panels illustrate philtrum, iris and canthal distances, eyebrow polygons and cheek-to-oval geometry with their normalization lengths.

## fig7_regional_discovery_rates

Significant/tested syndrome-feature counts and rates in nine anatomical regions. All tests belong to the single 5,760-hypothesis primary global-BH family (q < .05). Region denominators equal 48 syndromes times the number of features in the region and therefore differ; correlated features are not independent anatomical discoveries. Source: regional_discovery_source.csv, independently cross-checked against S1 and regional_summary_120.csv.

Alt text: Horizontal bars show discovery rates for all nine regions, with counts and percentages. Chin/jaw and forehead have the highest rates; denominators differ.

## fig8_synthetic_anchor_distributions

Real and synthetic interquartile distributions for all five predeclared syndrome-feature anchors. Boxes span Q1–Q3 and horizontal lines show medians in actual normalized feature units; no whiskers, extrema, individual observations or photographs are displayed. GMDB contributes one selected age-eligible image per recorded patient. PDIDB model tags are pooled as an observed release mixture; all synthetic ages are missing. These are within-source distribution summaries, not clinical replication.

Alt text: Five paired interquartile boxes compare GMDB and PDIDB for the predeclared philtrum, pupil-spacing, malar, eyebrow-area and canthal-width anchors.

## fig9_synthetic_effect_profiles

Real-versus-synthetic effect profiles against the identical frozen 3,409-path FairFace reference. Each point is a syndrome-feature aggregate (120 per panel), not an individual; dashed lines indicate equality. Orange marks strictly opposite effect directions. Titles report within-syndrome Pearson correlation and direction agreement. Hedges’ g uses syndrome-minus-reference mean difference and pooled SD with the stated small-sample correction. Correlated features, a shared reference and GMDB-trained generation induce dependence; no independent confirmatory p-values or confidence intervals are inferred.

Alt text: Nine effect-profile scatter panels show real and synthetic syndrome-minus-reference Hedges g, an identity line and per-syndrome correlation and sign agreement.

## fig10_synthetic_fidelity

Per-syndrome descriptive fidelity across 120 correlated features. A: median absolute synthetic-minus-real mean shift divided by the real sample SD. B: median sample variance ratio; dashed line is equal dispersion. C: agreement of the signs of real/synthetic Hedges’ g against the same fixed FairFace reference. Model tags are pooled as the empirical release mixture; stratified and 1%/99% winsorized results are supplied separately. No bootstrap intervals or biological independence across features are assumed.

Alt text: Three aligned bar panels show standardized mean shift, variance preservation and effect-direction agreement for nine matched syndromes.

# Final 120-feature results insertions

Prepared November 9, 2026 for the Formal Paper tab of FaceKit Expansion Results. The Google Doc was inspected read-only; no tab or comment was edited. These replacements use the frozen current aggregate package. Exact machine-readable values and formulas are in publication/final_research_120; displayed numbers below are rounded to six significant digits. No primary models were rerun or replaced.

## Primary recurrence: replace the top-five-feature placeholder

Across the 48 eligible syndromes, the most recurrent globally significant measurements were bilateral mean lateral eyebrow thickness (32 syndromes), forehead height (29), chin height (27), left lateral eyebrow thickness (27) and right lateral eyebrow thickness (26). Eyebrow and forehead effects were positive; all 27 significant chin-height effects were negative. Left/right eyebrow indices share anatomy, so their recurrence is not evidence for independent phenotypes. Significance means q < 0.05 in the single 5,760-test global BH family, and all five measurements were tested in all 48 cohorts.

| feature | feature_name | anatomical_region | significant_syndromes | tested_syndromes | positive | negative |
| --- | --- | --- | --- | --- | --- | --- |
| eb_lateral_thickness_mean | Bilateral lateral eyebrow thickness index | eyebrows | 32 | 48 | 32 | 0 |
| forehead_height | Forehead height / bizygomatic width | forehead | 29 | 48 | 29 | 0 |
| chin_height | Chin height / bizygomatic width | chin / jaw | 27 | 48 | 0 | 27 |
| eb_lateral_thickness_l | Left lateral eyebrow thickness index | eyebrows | 27 | 48 | 27 | 0 |
| eb_lateral_thickness_r | Right lateral eyebrow thickness index | eyebrows | 26 | 48 | 26 | 0 |

Table caption: Five measurements ranked by number of primary significant syndrome cohorts, with lexical feature-ID tie breaking. Counts use global BH q < .05, not within-syndrome q. Forehead/chin height and lateral eyebrow thickness use the precise S6 normalized geometry; neither pixel distances nor clinical hair density are implied. Exact formulas are retained in the source CSV.

## Regional results: replace region-table/figure placeholders

Discovery rates varied across anatomical regions. Chin/jaw contributed 87/240 significant hypotheses (36.25%) and forehead 60/192 (31.25%), compared with 197/960 for eyebrows (20.52%) and 35/384 for cheeks/malar geometry (9.11%). Eyebrows generated more absolute discoveries while chin/jaw had a higher rate. Different regional feature counts produce unequal denominators, and shared landmarks make these correlated summaries rather than independent anatomical enrichment tests.

| anatomical_region | tested | significant | features | fraction_significant |
| --- | --- | --- | --- | --- |
| cheeks / malar region | 384 | 35 | 8 | 0.0911458 |
| chin / jaw | 240 | 87 | 5 | 0.3625 |
| eyebrows | 960 | 197 | 20 | 0.205208 |
| eyes / orbital dimensions | 1392 | 158 | 29 | 0.113506 |
| forehead | 192 | 60 | 4 | 0.3125 |
| global face shape / proportions | 288 | 59 | 6 | 0.204861 |
| nose | 816 | 133 | 17 | 0.16299 |
| ocular position / iris / gaze | 768 | 96 | 16 | 0.125 |
| philtrum / mouth / lips | 720 | 130 | 15 | 0.180556 |

Table caption: Significant and tested primary syndrome-feature hypotheses in all nine regions. Regional denominators are 48 times the region’s feature count; all cells remain in the same global correction family. Figure 7 uses this exact source table.

## Recorded-demographic examples: replace Analysis A/B/C placeholders

Selection was declared in SELECTION_RULES.md before extracting examples: A uses three fixed readable ratios; B/C prioritize supported cell counts and use the same feature sequence, without filtering by q or effect direction. Nonsignificant examples are retained. The selected effects describe source-labeled groups and cannot be interpreted as genetic ancestry, intrinsic group characteristics or causal demographic effects.

The examples use these exact S6 definitions (all unitless):

| feature | formula | units |
| --- | --- | --- |
| philtrum_length | _dist(lm[SUBNASALE], lm[UPPER_LIP_TOP]) / midface_height | unitless |
| inter_pupillary_distance | _dist(lm[R_IRIS_CENTER], lm[L_IRIS_CENTER]) / bizygomatic_width | unitless |
| forehead_height | abs(lm[TRICHION, 1] - lm[GLABELLA, 1]) / bizygomatic_width | unitless |

### Broad recorded groups

**A: reference-only differences.** philtrum_length: H=74.3519, q=1.49669e-15; inter_pupillary_distance: H=376.875, q=9.03182e-80; forehead_height: H=7.8944, q=0.0546174. Each omnibus q uses the layer’s 120-feature family. An omnibus association does not identify a particular group pair or a disease effect; no pairwise interval is inferred from this output. The following means/medians retain the six-significant-digit precision encoded in the source A table.

| feature | group | n | mean | median |
| --- | --- | --- | --- | --- |
| philtrum_length | African | 500 | 0.284012 | 0.281947 |
| philtrum_length | Asian | 1909 | 0.284933 | 0.283111 |
| philtrum_length | European | 500 | 0.306656 | 0.304248 |
| philtrum_length | Others | 500 | 0.290099 | 0.291628 |
| inter_pupillary_distance | African | 500 | 0.477039 | 0.477635 |
| inter_pupillary_distance | Asian | 1909 | 0.452612 | 0.452368 |
| inter_pupillary_distance | European | 500 | 0.454004 | 0.453248 |
| inter_pupillary_distance | Others | 500 | 0.455542 | 0.45483 |
| forehead_height | African | 500 | 0.199324 | 0.199978 |
| forehead_height | Asian | 1909 | 0.199869 | 0.201234 |
| forehead_height | European | 500 | 0.199962 | 0.20082 |
| forehead_height | Others | 500 | 0.202266 | 0.20342 |

**B: syndrome versus matched reference.**

Within the recorded European group, Coffin-Siris syndrome (33 disease, 500 reference) has philtrum_length Hedges’ g=+0.718637, marginal 95% bootstrap CI [0.265289, 1.20428], global layer-B q=0.126031 across 34,560 tests. The normalized means are 0.345906 versus 0.306656, a difference of +0.0392502 feature units.

Within the recorded European group, Cornelia de Lange syndrome (29 disease, 500 reference) has inter_pupillary_distance Hedges’ g=-0.145394, marginal 95% bootstrap CI [-0.575921, 0.307251], global layer-B q=0.743263 across 34,560 tests. The normalized means are 0.450638 versus 0.454004, a difference of -0.00336634 feature units.

Within the recorded Asian group, WILLIAMS-BEUREN SYNDROME; WBS (25 disease, 1909 reference) has forehead_height Hedges’ g=+1.45296, marginal 95% bootstrap CI [1.07485, 1.83247], global layer-B q=4.99523e-06 across 34,560 tests. The normalized means are 0.228429 versus 0.199869, a difference of +0.0285595 feature units.

**C: between-group differences in disease contrasts.**

For Cornelia de Lange syndrome, philtrum_length, group 1=Others (disease/reference 21/500) and group 2=European (29/500), the group-2-minus-group-1 disease/reference interaction is -0.00714586 normalized feature units (95% HC3 CI [-0.040461, 0.0261693]; global layer-C q=0.977567, 11,040 tests). This is a difference of contrasts, not a Hedges’ g or reference-only group difference.

For WILLIAMS-BEUREN SYNDROME; WBS, inter_pupillary_distance, group 1=Asian (disease/reference 25/1909) and group 2=Others (21/500), the group-2-minus-group-1 disease/reference interaction is -0.0114787 normalized feature units (95% HC3 CI [-0.026832, 0.00387456]; global layer-C q=0.633276, 11,040 tests). This is a difference of contrasts, not a Hedges’ g or reference-only group difference.

For Kabuki syndrome, forehead_height, group 1=Asian (disease/reference 19/1909) and group 2=European (13/500), the group-2-minus-group-1 disease/reference interaction is -0.0200829 normalized feature units (95% HC3 CI [-0.0361435, -0.00402229]; global layer-C q=0.160727, 11,040 tests). This is a difference of contrasts, not a Hedges’ g or reference-only group difference.

### Detailed recorded groups

**A: reference-only differences.** philtrum_length: H=127.622, q=9.04303e-25; inter_pupillary_distance: H=692.415, q=1.06221e-144; forehead_height: H=66.8028, q=2.91613e-12. Each omnibus q uses the layer’s 120-feature family. An omnibus association does not identify a particular group pair or a disease effect; no pairwise interval is inferred from this output. The following means/medians retain the six-significant-digit precision encoded in the source A table.

| feature | group | n | mean | median |
| --- | --- | --- | --- | --- |
| philtrum_length | African | 500 | 0.284012 | 0.281947 |
| philtrum_length | East Asian | 500 | 0.290228 | 0.286277 |
| philtrum_length | European | 500 | 0.306656 | 0.304248 |
| philtrum_length | Latino/Hispanic | 500 | 0.290099 | 0.291628 |
| philtrum_length | Middle Eastern | 409 | 0.296792 | 0.294443 |
| philtrum_length | SE Asian | 500 | 0.275243 | 0.273986 |
| philtrum_length | South Asian | 500 | 0.279629 | 0.278817 |
| inter_pupillary_distance | African | 500 | 0.477039 | 0.477635 |
| inter_pupillary_distance | East Asian | 500 | 0.438228 | 0.437306 |
| inter_pupillary_distance | European | 500 | 0.454004 | 0.453248 |
| inter_pupillary_distance | Latino/Hispanic | 500 | 0.455542 | 0.45483 |
| inter_pupillary_distance | Middle Eastern | 409 | 0.455743 | 0.454575 |
| inter_pupillary_distance | SE Asian | 500 | 0.451908 | 0.45045 |
| inter_pupillary_distance | South Asian | 500 | 0.465138 | 0.465253 |
| forehead_height | African | 500 | 0.199324 | 0.199978 |
| forehead_height | East Asian | 500 | 0.194103 | 0.195261 |
| forehead_height | European | 500 | 0.199962 | 0.20082 |
| forehead_height | Latino/Hispanic | 500 | 0.202266 | 0.20342 |
| forehead_height | Middle Eastern | 409 | 0.199598 | 0.20152 |
| forehead_height | SE Asian | 500 | 0.201978 | 0.202797 |
| forehead_height | South Asian | 500 | 0.203749 | 0.204119 |

**B: syndrome versus matched reference.**

Within the recorded European group, Coffin-Siris syndrome (33 disease, 500 reference) has philtrum_length Hedges’ g=+0.718637, marginal 95% bootstrap CI [0.298185, 1.19245], global layer-B q=0.130001 across 32,520 tests. The normalized means are 0.345906 versus 0.306656, a difference of +0.0392502 feature units.

Within the recorded European group, Cornelia de Lange syndrome (29 disease, 500 reference) has inter_pupillary_distance Hedges’ g=-0.145394, marginal 95% bootstrap CI [-0.569771, 0.281113], global layer-B q=0.740749 across 32,520 tests. The normalized means are 0.450638 versus 0.454004, a difference of -0.00336634 feature units.

Within the recorded European group, KBG SYNDROME; KBGS (24 disease, 500 reference) has forehead_height Hedges’ g=+0.654195, marginal 95% bootstrap CI [0.211008, 1.13715], global layer-B q=0.144568 across 32,520 tests. The normalized means are 0.21327 versus 0.199962, a difference of +0.0133082 feature units.

**C: between-group differences in disease contrasts.**

For Cornelia de Lange syndrome, philtrum_length, group 1=Latino/Hispanic (disease/reference 21/500) and group 2=European (29/500), the group-2-minus-group-1 disease/reference interaction is -0.00714586 normalized feature units (95% HC3 CI [-0.040461, 0.0261693]; global layer-C q=0.98628, 14,400 tests). This is a difference of contrasts, not a Hedges’ g or reference-only group difference.

For WILLIAMS-BEUREN SYNDROME; WBS, inter_pupillary_distance, group 1=Latino/Hispanic (disease/reference 15/500) and group 2=South Asian (15/500), the group-2-minus-group-1 disease/reference interaction is +0.00945443 normalized feature units (95% HC3 CI [-0.00963051, 0.0285394]; global layer-C q=0.867815, 14,400 tests). This is a difference of contrasts, not a Hedges’ g or reference-only group difference.

For Kabuki syndrome, forehead_height, group 1=East Asian (disease/reference 16/500) and group 2=European (13/500), the group-2-minus-group-1 disease/reference interaction is -0.0293148 normalized feature units (95% HC3 CI [-0.0461129, -0.0125166]; global layer-C q=0.0175854, 14,400 tests). This is a difference of contrasts, not a Hedges’ g or reference-only group difference.

These support-selected examples show what global counts conceal. Broad-group forehead height has no omnibus BH association at q=0.0546, although the broad Asian Williams–Beuren syndrome-reference contrast is large (g=1.453). The detailed Kabuki forehead interaction is significant (q=0.0176), while the selected detailed Cornelia de Lange philtrum and Williams–Beuren pupil-spacing interactions are not. A shared reference difference therefore cannot be substituted for evidence that a disease contrast differs across groups. All C models remain exploratory and source-confounded.

Table caption: Preselected A/B/C examples with source-label names, complete cell counts, exact feature formulas/units and the appropriate correction family. A reports omnibus H; B reports syndrome-minus-reference Hedges’ g and marginal bootstrap intervals; C reports a group-2-minus-group-1 interaction in normalized units and HC3 intervals. Global q is separate for every layer/analysis. Full estimates and formulas are in demographic_examples.csv; detailed A summaries are in demographic_A_group_summaries.csv.

## Adjusted results: replace Table 5 placeholder

| syndrome | feature | comparator | beta | ci_lower | ci_upper | n_complete | n_disease | n_comparator | q | family_hypotheses |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Cornelia de Lange syndrome | philtrum_length | all other age-eligible GMDB patients | 0.0366636 | 0.0225898 | 0.0507374 | 1810 | 63 | 1747 | 1.31693e-06 | 4 |
| Noonan syndrome | inter_pupillary_distance | all other age-eligible GMDB patients | 0.0210082 | 0.00927624 | 0.0327401 | 1810 | 34 | 1776 | 0.00059821 | 4 |
| WILLIAMS-BEUREN SYNDROME; WBS | eb_thickness_mean | all other age-eligible GMDB patients | 0.00119796 | 0.000719539 | 0.00167639 | 1810 | 67 | 1743 | 1.84287e-06 | 4 |
| Rubinstein-Taybi syndrome | malar_bulge_mean | all other age-eligible GMDB patients | -0.00220797 | -0.00563492 | 0.00121897 | 1810 | 24 | 1786 | 0.20666 | 4 |
| Cornelia de Lange syndrome | philtrum_length | frozen FairFace reference | 0.0451972 | 0.032053 | 0.0583415 | 3471 | 62 | 3409 | 4.76988e-11 | 3 |
| Noonan syndrome | inter_pupillary_distance | frozen FairFace reference | 0.0355073 | 0.0236378 | 0.0473767 | 3435 | 26 | 3409 | 6.81017e-09 | 3 |
| Rubinstein-Taybi syndrome | malar_bulge_mean | frozen FairFace reference | -0.00606285 | -0.0092521 | -0.0028736 | 3433 | 24 | 3409 | 0.000194586 | 3 |
| Cornelia de Lange syndrome | philtrum_length | all other age-eligible GMDB patients | 0.0366636 | 0.0225898 | 0.0507374 | 1810 | 63 | 1747 | 4.39973e-05 | 5880 |
| Noonan syndrome | inter_pupillary_distance | all other age-eligible GMDB patients | 0.0210082 | 0.00927624 | 0.0327401 | 1810 | 34 | 1776 | 0.00945558 | 5880 |
| WILLIAMS-BEUREN SYNDROME; WBS | eb_thickness_mean | all other age-eligible GMDB patients | 0.00119796 | 0.000719539 | 0.00167639 | 1810 | 67 | 1743 | 9.91512e-05 | 5880 |
| WILLIAMS-BEUREN SYNDROME; WBS | outer_canthal_distance | all other age-eligible GMDB patients | 0.0177318 | 0.00945907 | 0.0260044 | 1810 | 67 | 1743 | 0.00126008 | 5880 |
| Rubinstein-Taybi syndrome | malar_bulge_mean | all other age-eligible GMDB patients | -0.00220797 | -0.00563492 | 0.00121897 | 1810 | 24 | 1786 | 0.537409 | 5880 |

Table caption: Seven prespecified adjusted models and five anchor cells from the full adjusted screen. Coefficients and HC3 95% intervals are in the feature’s normalized units, not pooled-SD units. The first four rows use age/sex adjustment and a four-model BH family; the next three use recorded sex, broad group and yaw/pitch/roll and a three-model BH family. The remaining five retain the full 5,880-model family rather than being recorrected as a selected subset. Complete-case n and target/comparator counts are shown; all unsupported screen cells retain correction p=1 in the planned denominator. All exact family labels, formulas and source filenames are in adjusted_examples.csv.

Against FairFace, the covariate-adjusted philtrum_length coefficient for Cornelia de Lange syndrome is +0.045197208 (95% HC3 CI [0.032052959, 0.058341456]; three-model q=4.7698794e-11; complete-case n=3471, disease/reference=62/3409).

Against FairFace, the covariate-adjusted inter_pupillary_distance coefficient for Noonan syndrome is +0.035507268 (95% HC3 CI [0.023637832, 0.047376705]; three-model q=6.8101686e-09; complete-case n=3435, disease/reference=26/3409).

Against FairFace, the covariate-adjusted malar_bulge_mean coefficient for Rubinstein-Taybi syndrome is -0.0060628527 (95% HC3 CI [-0.0092521045, -0.0028736008]; three-model q=0.00019458582; complete-case n=3433, disease/reference=24/3409).

The age/sex-adjusted within-GMDB Rubinstein–Taybi malar coefficient is -0.0022079736 (95% HC3 CI [-0.0056349159, 0.0012189688]; four-model q=0.20666045). It is not significant at q < .05. The FairFace association must not be described as persisting under every comparator or adjustment.

### Why 48 primary cohorts but 49 adjusted targets

Primary eligibility is computed from the current selected-patient counts with n >= 10. The separately prespecified adjusted target list retains 49 labels, including FBXW7__cbba0000, which has nine current selected patients. All other 48 targets match the primary set. The preparation source fixed the main list from current counts and retained the earlier prespecified target list for extension bookkeeping; no additional current measurement or obsolete estimate is used to enlarge the primary family. extensions.py passes plan["targets"] to 48-cohort comparisons and plan["adjustment_targets"] to the full adjusted screen. Its support rules reject sparse target cells; the FBXW7 rows and their correction denominator remain rather than forcing the lists to coincide. The active ANALYSIS_PLAN.json, syndrome_eligibility.csv, extensions.py and full_adjusted_screen.csv establish this distinction. Archived preparation code was consulted only to audit how the target list was declared.

### Comparator nonagreement

The primary FairFace family has 955 associations and the all-other-GMDB family 222; 173 are shared, 782 FairFace-only and 49 GMDB-only. These comparisons estimate different contrasts. They do not imply that 782 primary results are false positives or estimate a clinical false-positive rate. Source effects, broad disease-shared morphology and comparator variance can all change which associations pass their separate global BH thresholds.

## Clinical literature reassessment

All 17 appendix questions were linked to a unique active syndrome-feature row and its exact S6 definition; their old estimates were not reused. All 17 active associations pass the current global threshold, but that does not validate their clinical meaning. Four have directly related trait support: Cornelia de Lange/philtrum length, Noonan/pupil spacing, Wiedemann–Steiner/inner-canthal spacing and Sotos/mid-forehead width. Even these do not validate the exact normalized ratios as clinical measurements. The other 13 remain related geometric observations, scope mismatches or unverified interpretations, and should not be retained as established clinical confirmations. The five original anchors may remain as geometric examples with their existing qualifications.

**A1. WILLIAMS-BEUREN SYNDROME; WBS / eb_thickness_mean.** n=67; g=+1.26614, 95% CI [0.98665, 1.58683], global q=6.24183e-13. Periorbital fullness does not establish eyebrow polygon area or hair density. Status: related construct only. [Clinical source: Clinical Description / facial features](https://www.ncbi.nlm.nih.gov/books/NBK1249/).

**A2. Cornelia de Lange syndrome / philtrum_length.** n=64; g=+1.16283, 95% CI [0.885793, 1.41151], global q=1.08731e-12. Long/smooth philtrum supports length-related interpretation; midface denominator not clinically validated. Status: direct related trait. [Clinical source: Suggestive Findings / craniofacial appearance](https://www.ncbi.nlm.nih.gov/books/NBK1104/).

**A3. WILLIAMS-BEUREN SYNDROME; WBS / outer_canthal_distance.** n=67; g=+1.03943, 95% CI [0.779987, 1.30702], global q=1.26103e-10. Periorbital fullness does not establish increased outer-canthal/bizygomatic ratio. Status: related construct only. [Clinical source: Clinical Description / facial features](https://www.ncbi.nlm.nih.gov/books/NBK1249/).

**A4. Cornelia de Lange syndrome / philtrum_width.** n=64; g=-0.399052, 95% CI [-0.593889, -0.217961], global q=0.0213725. Long/smooth philtrum does not establish philtrum width. Status: unsupported exact trait. [Clinical source: Suggestive Findings / craniofacial appearance](https://www.ncbi.nlm.nih.gov/books/NBK1104/).

**A5. Noonan syndrome / inter_pupillary_distance.** n=34; g=+1.24369, 95% CI [0.757464, 1.68374], global q=5.4355e-06. Widely spaced eyes supports spacing-related interpretation; normalized pupil ratio is a distinct measurement. Status: direct related trait. [Clinical source: Suggestive Findings / characteristic facies](https://www.ncbi.nlm.nih.gov/books/NBK1124/).

**A6. Mucopolysaccharidoses / chin_pointedness_angle.** n=37; g=+0.645967, 95% CI [0.405549, 0.877937], global q=0.00192438. An MPS-II description cannot validate a pooled Mucopolysaccharidoses cohort or this jaw angle. Status: scope and construct mismatch. [Clinical source: Clinical Description / facial features](https://www.ncbi.nlm.nih.gov/books/NBK1274/).

**A7. Noonan syndrome / nose_bridge_width.** n=34; g=+1.1297, 95% CI [0.698003, 1.56309], global q=1.73697e-05. Depressed bridge and wide nasal base are not lateral bridge width; age-dependent description. Status: related construct only. [Clinical source: Clinical Description / age-dependent facies](https://www.ncbi.nlm.nih.gov/books/NBK1124/).

**A8. KBG SYNDROME; KBGS / jaw_bigonial_width.** n=32; g=-1.24157, 95% CI [-1.52135, -0.936602], global q=2.7789e-08. Triangular face does not establish a particular bigonial/bizygomatic ratio. Status: related construct only. [Clinical source: Summary / characteristic facial features](https://www.ncbi.nlm.nih.gov/books/NBK487886/).

**A9. WIEDEMANN-STEINER SYNDROME; WDSTS / inter_canthal_distance.** n=14; g=+1.79782, 95% CI [1.25702, 2.31914], global q=1.79482e-05. Widely spaced eyes supports spacing-related interpretation; inner-canthal normalization remains distinct. Status: direct related trait. [Clinical source: Clinical Description / characteristic dysmorphic features](https://www.ncbi.nlm.nih.gov/books/NBK580718/).

**A10. Kabuki syndrome / outer_canthal_distance.** n=35; g=+0.978914, 95% CI [0.662411, 1.30012], global q=7.20609e-06. Long palpebral fissures do not directly establish total outer-canthal separation. Status: related construct only. [Clinical source: Summary / clinical characteristics](https://www.ncbi.nlm.nih.gov/books/NBK62111/).

**A11. Holoprosencephaly / mouth_opening.** n=11; g=+1.90458, 95% CI [1.04734, 2.7332], global q=0.00210709. Structural midline anomalies do not validate expression-dependent mouth opening in a photograph. Status: unsupported exact trait. [Clinical source: Clinical Characteristics / craniofacial spectrum](https://www.ncbi.nlm.nih.gov/books/NBK1530/).

**A12. Coffin-Siris syndrome / philtrum_length.** n=43; g=+0.834805, 95% CI [0.439058, 1.28472], global q=0.00422404. Retrieved clinical source does not verify the exact long-philtrum claim; retain as unresolved clinical interpretation. Status: not verified. [Clinical source: Clinical Description / facial features](https://www.ncbi.nlm.nih.gov/books/NBK131811/).

**A13. KBG SYNDROME; KBGS / jaw_to_forehead_ratio.** n=32; g=-1.29814, 95% CI [-1.61619, -0.954496], global q=3.10055e-08. Triangular face is related to outline, not proof of the exact jaw/forehead ratio. Status: related construct only. [Clinical source: Summary / characteristic facial features](https://www.ncbi.nlm.nih.gov/books/NBK487886/).

**A14. SOTOS SYNDROME; SOTOS / forehead_width_mid.** n=25; g=+1.31644, 95% CI [0.897713, 1.7262], global q=5.46577e-06. Broad prominent forehead supports a width-related interpretation; facial-width denominator remains distinct. Status: direct related trait. [Clinical source: Suggestive Findings / characteristic facial appearance](https://www.ncbi.nlm.nih.gov/books/NBK1479/).

**A15. Kabuki syndrome / inter_canthal_distance.** n=35; g=+0.902678, 95% CI [0.525977, 1.29085], global q=0.000241116. Long palpebral fissures and eyelid eversion do not establish inner-canthal distance. Status: related construct only. [Clinical source: Summary / clinical characteristics](https://www.ncbi.nlm.nih.gov/books/NBK62111/).

**A16. MANNOSIDOSIS, ALPHA B, LYSOSOMAL; MANSA / inter_pupillary_distance.** n=14; g=+1.45126, 95% CI [0.945106, 1.97278], global q=0.000365101. Coarse facies and widely spaced teeth do not establish pupil spacing. Status: unsupported exact trait. [Clinical source: Clinical Description / facial features](https://www.ncbi.nlm.nih.gov/books/NBK1396/).

**A17. Rubinstein-Taybi syndrome / malar_bulge_mean.** n=24; g=-1.00903, 95% CI [-1.53038, -0.513744], global q=0.00251863. Facial gestalt and grimacing smile do not clinically validate the projected cheek/oval index. Status: related construct only. [Clinical source: Clinical Description / craniofacial features](https://www.ncbi.nlm.nih.gov/books/NBK1526/).

## Separate synthetic consistency: ready-to-paste framing

The primary PRISM study analyzes real GMDB patient photographs against selected FairFace reference photographs. The mentor reports that PDIDB contains syndrome-conditioned images sampled from a StyleGAN model trained on GMDB. Because this learned-source relationship does not supply independent patients, PDIDB measurements were evaluated separately for conditional distributional fidelity. Checkpoint-specific generation lineage remains unverified. Nine matched syndromes cover 314 selected GMDB patients and 1,424 QC-passing synthetic images. Across syndrome summaries, median absolute mean shift was 0.166 real SD, median variance ratio 0.986, median effect-direction agreement 85.0% and median Pearson effect-profile correlation 0.896 against the shared frozen FairFace baseline. Age/sex/race are missing for all synthetic records, preventing equivalent age filtering and limiting interpretation of these similarities and discrepancies. No synthetic-derived hypothesis tests contribute support to the 955 primary associations. Full findings and sensitivities are in synthetic_consistency_results.md.

## Figure captions and placement


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

Place Figures 5/6 in Methods to replace workflow/anatomy placeholders, Figure 7 beside regional Results, and Figures 8–10 in the separate synthetic subsection or supplement. These repository numbers continue after existing Figures 1–4; renumber in the Google Doc only when inserted.

## Results-source map and unresolved items

| Formal Paper placeholder / claim | Exact current source and insertion |
| --- | --- |
| Top-five recurrent features | S1 → top_five_recurrent_features.csv; primary recurrence paragraph |
| Nine regions and regional figure | regional_summary_120.csv + independent S1 regrouping → regional_discovery_source.csv; Figure 7 |
| Demographic A/B/C examples | Complete layer1/layer2 A/B/C CSVs → demographic_examples.csv and demographic_A_group_summaries.csv; six separate families |
| Table 5 adjusted models | within_gmdb_age_sex_adjusted.csv; three_predefined_reference_adjusted_models.csv; full_adjusted_screen.csv → adjusted_examples.csv |
| Comparator overlap | S1 + within_gmdb_results.csv → comparator_overlap.csv |
| Appendix A1–A17 | S1 + S6 + cited clinical sections → clinical_literature_audit_17.csv |
| Extraction/workflow and geometry art | ACTIVE_ANALYSIS.json + S6 + exact pinned source → Figures 5/6 |
| New PDIDB fidelity | synthetic_consistency_manifest.json and its aggregate CSVs → Figures 8–10 and separate results report |

Pending scientific judgments: confirm the nine broad syndrome mappings clinically; establish checkpoint/generator release identity and exact training membership; obtain supported age metadata before claiming age-matched synthetic fidelity; validate the exact eyebrow/malar and other related geometric constructs clinically. The Coffin–Siris long-philtrum interpretation remains unverified in the retrieved source. No new clinical novelty is asserted from absence in a general description. The separate FaceKit overlap audit defines the limited statistical contribution. No external clinical validation, diagnostic accuracy, independent synthetic replication or new extractor claim is supported. Google Docs placeholders/comments still require insertion/resolution; BIBM and PRISM tabs were not edited. The handoff excludes privacy/IRB/legal release review, history remediation and full frontend testing; existing PI items remain pending.

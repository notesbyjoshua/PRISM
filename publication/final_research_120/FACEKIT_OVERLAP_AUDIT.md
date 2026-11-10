# Scientific overlap and differences

Source inspected: [FaceKit, arXiv:2610.09288v1](https://arxiv.org/html/2610.09288v1), November 9, 2026. This audit concerns scientific questions and figure content, not artwork reuse rights. All new schematics and plots are original.

| FaceKit overlap (precise location) | PRISM difference and exact output |
| --- | --- |
| Landmark/normalization platform: §§2.1–2.3, Fig. 1A | Reused platform, explicitly attributed; current source/model pins in ACTIVE_ANALYSIS.json and S6. Figure 5 adds patient selection and analytical branches. |
| HPO correspondence: §3.1, Fig. 3 | S1 tests all 120 features in each of 48 current age-restricted syndrome cohorts against 3,409 references, with one global 5,760-test family. Clinical descriptors are not novel. |
| Philtrum, iris/canthal, eyebrow and cheek geometry: §2.2; qualitative term mapping in §3.1 | Five clinical_highlights_120.csv anchors quantify selected contrasts. No claim that these indices or their general clinical interpretations are new. Exact pair-specific estimates require PRISM source rows. Figure 6 explains their definitions without reusing artwork. |
| Pose/repeat photographs: §3.2, Table 1, Figs. 4–5 | PRISM uses the pinned correction; alternate and random-photo outputs assess selection sensitivity rather than validating a new correction. |
| Age/reference sensitivity: §3.4, Fig. 7B–C | PRISM fixes childhood GMDB eligibility and compares expanded reference compositions; no age-conditional-reference superiority claim. |
| Recorded-group/reference effects and interaction: §3.5, Figs. 7A/8 | General demographic sensitivity and disease-by-group questions overlap. PRISM's two layers and separately corrected A/B/C families use different cohorts/contrasts, including sparse cells. demographic_examples.csv reports its actual estimates. |
| Synthetic fidelity, variance and effect profiles: §§2.6/3.6, Fig. 9 | Same broad question. synthetic_consistency/ adds nine explicitly mapped PDIDB labels, age-restricted one-patient selection, all 120 features, fixed-FairFace profiles and model-tag/winsorization checks. It does not replicate an independent clinical sample. |
| Cohort-versus-remaining contrasts: §3.6, Fig. 9A | PRISM's within_gmdb_results.csv applies a related comparator to selected real patients. Its 5,880-model age/sex screen and targeted families have explicitly different target/support rules; no generic contrast or regression-method novelty is claimed. |
| External dataset/privacy studies: §§3.3/3.7, Figs. 6/10, Table 2 | Not performed here. Encoded-byte checks and geometric fidelity do not establish privacy, identity independence or external clinical validity. |

The contribution that can be stated is a reproducible statistical investigation of specified current cohorts, contrast definitions, multiplicity families and sensitivity analyses using FaceKit. Novelty of any individual biological association requires a dedicated systematic literature review. The 17-entry clinical audit supports only trait-related interpretations, never the exact ratio's clinical validation. Generator-release equivalence to FaceKit's paper remains unverified.

## Five highlighted pairs

This maps each highlighted pair to the reused platform and its current result. The cited general experiment does not establish that an exact pair was reported previously, or that it is new.

| FaceKit overlap / cited experiment | Current PRISM pair and interpretation |
| --- | --- |
| Philtrum catalogue; §2.2 and §3.1/Fig. 3 | Cornelia de Lange / philtrum_length: clinical_highlights_120.csv, g=1.16283; related clinical length trait, exact ratio unvalidated. |
| Ocular spacing; §2.2 and §3.1/Fig. 3 | Noonan / inter_pupillary_distance: clinical_highlights_120.csv, g=1.24369; related spacing trait, exact ratio unvalidated. |
| Eyebrow geometry; §2.2 and §3.1/Fig. 3 | Williams–Beuren / eb_thickness_mean: clinical_highlights_120.csv, g=1.26614; polygon area, not confirmed hair density. |
| Canthal geometry; §2.2 and §3.1/Fig. 3 | Williams–Beuren / outer_canthal_distance: clinical_highlights_120.csv, g=1.03943; periorbital appearance is only a related construct. |
| Cheek geometry; §2.2 and §3.1/Fig. 3 | Rubinstein–Taybi / malar_bulge_mean: clinical_highlights_120.csv, g=−1.00903; projected geometry, not validated clinical malar prominence. |

Exact current rows and clinical citations for all 17 earlier questions are in clinical_literature_audit_17.csv. They reuse the §2.2 measurement catalogue and are interpreted within the §3.1 HPO-correspondence question; their PRISM effect estimates and cohort contrasts are newly computed current outputs. None is labeled a new biological finding on that basis.

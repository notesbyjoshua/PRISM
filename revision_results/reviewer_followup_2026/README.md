# Formal Paper reviewer follow-up code

Local analyses of the frozen historical 125-feature PRISM cohort. The code reconstructs the original cohorts, compares the 49 eligible syndrome groups against other GMDB syndromes, evaluates photograph-selection and subgroup sensitivity, inventories metadata, and produces aggregate reports and figures. It does not edit frozen results, existing figures, or manuscripts, and does not process patient photographs.

## Run locally

From the repository root, with a Python environment and the required local inputs (this workflow does not establish data-use authorization):

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python revision_results/reviewer_followup_2026/followup.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python revision_results/reviewer_followup_2026/audit_extensions.py
.venv/bin/python revision_results/reviewer_followup_2026/verify_hpo.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m pytest revision_results/reviewer_followup_2026/test_followup.py revision_results/reviewer_followup_2026/test_review_contracts.py -q --tb=short
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python revision_results/reviewer_followup_2026/build_reports.py
```

Dependencies are NumPy, Pandas, SciPy, statsmodels, Matplotlib, and pytest. `verify_hpo.py` downloads the public HPO ontology; no study metadata is transmitted. Statistical commands use local inputs only. The primary analysis defaults to 1,000 patient-row bootstrap samples and 100 random one-photograph-per-patient selections, with recorded seeds and separate analysis-specific BH families.

Generated results, reports, logs, figures, and the local metadata/permissions inventory are ignored by Git. They require data-owner redistribution review before public release. Repository visibility does not grant permission to use or redistribute patient data. Do not force-add local derivatives or patient records. The code does not establish clinical validity, genetic confirmation, or approved geometric-to-HPO equivalence.

## Independent statistical review

For an isolated dependency environment, create a new virtual environment and install `requirements.txt` from this directory. Run the five commands above from a fresh clone without initializing the FaceKit submodule, using that environment's Python. Preserve command logs, exit statuses, `pip freeze --all`, and the recorded provenance locally. The application checkout is informational; an uninitialized submodule must be recorded as unavailable, not as the parent repository's commit.

After the pipeline, run `.venv/bin/python revision_results/reviewer_followup_2026/independent_review.py`. This separately reconstructs patient/photo selection from the local masters, verifies all 12,250 point-estimate hypotheses across both comparator designs with rank-based calculations and a manual BH calculation, replays representative bootstrap intervals in every target syndrome using explicitly materialized resamples, checks four HC3 models by matrix calculations, and verifies the 16 restricted threshold families. Aggregate review files are written to ignored `local_only/independent_review/`; patient/image keys are kept in memory and never exported. Independent calculations verify software behavior and recorded-ID separation, not biological independence or clinical validity.

The eight HPO correspondences remain review candidates. `verify_hpo.py` retains the public ontology at ignored `local_only/hp.obo` with its hash/version. Pass `--ontology-file PATH` to replay that exact snapshot offline; the acquisition method and verification timestamp distinguish a download from local reuse. Report labels read the actual bootstrap and photograph-replicate counts from result files.

The original unadjusted main and eligible-49 sensitivity analyses each use a separate 6,125-hypothesis family. The original workflow adjusts four highlighted measurements in a separate four-model BH family. Subgroup thresholds restrict saved hypotheses and recompute BH; they do not refit models. Local planning documents are not an external preregistration. Generated review reports and results remain excluded from Git; only code, tests, requirements and documentation are published.

## New full-screen and pairwise extension

Read [EXTENSION_METHODS.md](EXTENSION_METHODS.md) for the analysis contract, support rules, fixed pair list construction, FDR denominators, bootstrap design and interpretation limits. The extension reuses the existing unadjusted tables and reference feature correlations rather than replacing or duplicating them. It adds a full age/sex HC3 screen, targeted pairwise contrasts, concordance uncertainty/attenuation by anatomy and feature family, cluster discovery concentration, and exploratory age strata.

After the original pipeline has supplied its local aggregate inputs, use a fresh ignored output directory:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python revision_results/reviewer_followup_2026/extended_analysis.py \
  --output-dir revision_results/reviewer_followup_2026/local_only/extension_new_run
PRISM_EXTENSION_RESULTS=revision_results/reviewer_followup_2026/local_only/extension_new_run \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m pytest \
  revision_results/reviewer_followup_2026/test_extended_analysis.py \
  revision_results/reviewer_followup_2026/test_followup.py \
  revision_results/reviewer_followup_2026/test_review_contracts.py -q --tb=short
```

Defaults are 1,000 patient bootstrap draws per pair/stratum and 1,000 joint concordance draws. `--bootstrap` and `--concordance-bootstrap` accept at least 100; actual counts are exported. Existing directories are refused, as are destinations outside this workflow's `local_only/`. Missing local inputs are fatal and recorded; the extension does not silently regenerate older results. `RESULTS_FOR_AUTHOR_REVIEW.md`, full and highlighted CSVs, sparse-model flags, five SVG/350-dpi PNG plots, execution status, numerical warnings and provenance are local-only. No patient/image keys are exported.

Tests without `PRISM_EXTENSION_RESULTS` run synthetic support, rank, HC3, resampling, independence and failure-recording checks; four completed-run checks explicitly skip. With that variable, the suite also checks every supported adjusted coefficient/HC3 SE/p by independent matrix calculation, all pair/age denominators and g/rank p-values, all three global BH families, representative explicit bootstrap intervals, concordance and cluster accounting. Bootstrap concordance intervals also export bias/SE and interval-exclusion flags; their nominal population coverage is unverified, especially when the distribution is displaced from the observed statistic. Do not use significance as clinical specificity or treat age-bin differences as tested developmental interactions.

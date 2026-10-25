# Formal Paper reviewer follow-up code

Local analyses of the frozen historical 125-feature PRISM cohort. The code reconstructs the original cohorts, compares the 49 eligible syndrome groups against other GMDB syndromes, evaluates photograph-selection and subgroup sensitivity, inventories metadata, and produces aggregate reports and figures. It does not edit frozen results, existing figures, or manuscripts, and does not process patient photographs.

## Run locally

From the repository root, with the existing Python environment and authorized local inputs:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python revision_results/reviewer_followup_2026/followup.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python revision_results/reviewer_followup_2026/audit_extensions.py
.venv/bin/python revision_results/reviewer_followup_2026/verify_hpo.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m pytest revision_results/reviewer_followup_2026/test_followup.py -q --tb=short
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python revision_results/reviewer_followup_2026/build_reports.py
```

Dependencies are NumPy, Pandas, SciPy, statsmodels, Matplotlib, and pytest. `verify_hpo.py` downloads the public HPO ontology; no study metadata is transmitted. Statistical commands use local inputs only. The primary analysis defaults to 1,000 patient-row bootstrap samples and 100 random one-photograph-per-patient selections, with recorded seeds and separate analysis-specific BH families.

Generated results, reports, logs, figures, and the local metadata/permissions inventory are ignored by Git. They require data-owner redistribution review before public release. Repository visibility does not grant permission to use or redistribute patient data. Do not force-add local derivatives or patient records. The code does not establish clinical validity, genetic confirmation, or approved geometric-to-HPO equivalence.

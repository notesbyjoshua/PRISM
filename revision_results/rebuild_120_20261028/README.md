# Sole active PRISM study

PRISM is a statistical investigation using 120 measurements supplied by FaceKit. The active manuscript is `manuscript/main.tex`; `manuscript_120.tex` here is its generation template. `ACTIVE_ANALYSIS.json` selects validated private inputs and the supported specifications. Main study: 1,845 recorded GMDB patients, 3,409 expanded FairFace reference paths, 48 primary cohorts, 5,760 tests and 955 global-BH associations. Synthetic PDIDB remains separate.

Active results include recorded-demographic A/B/C analyses, within-GMDB comparisons, prespecified adjusted/age/pairwise families and reference/photo robustness analyses. Planned correction denominators and unsupported cells are retained. Concordance is descriptive; its biased bootstrap uncertainty is excluded from active exports.

## Verification and reproduction

From the repository root:

```sh
.venv/bin/python stats_testing/statistical_analysis.py --output-dir local_review/NEW_PRIMARY_REPLAY
.venv/bin/python revision_results/rebuild_120_20261028/validate.py --scope full
.venv/bin/python revision_results/rebuild_120_20261028/build_materials.py
python3 publication/scripts/generate_figures.py
python3 publication/scripts/generate_supplementary.py
python3 publication/scripts/verify_reproducible_outputs.py
python3 -m pytest tests -q
python3 -m pytest revision_results/reviewer_followup_2026/test_extended_analysis.py -q
```

Install `publication/requirements.txt` in an isolated environment and initialize the pinned FaceKit submodule for package tests. Restricted inputs are not required for CI aggregate verification. Full numerical validation and extraction require the authorized local source tables and model. Image extraction is documented in `data/dataset_creation/CURRENT_120_EXTRACTION.md` and must remain private.

`prepare.py --dry-run`, `run_sensitivities.py --dry-run` and `extensions.py --dry-run` inspect the current-only plan. Extraction resumes with `--resume` under exact input/model/code signatures. Statistical jobs verify completed checkpoint hashes; interrupted jobs restart their specification. Changed numerical code or inputs require a new run. Source-location and export-policy changes retain the original execution signatures with explicit private postprocessing records; numerical utility function ASTs were checked unchanged.

The private source run is `local_only/run_20261028`. Current reporting, complete primary replay comparisons, build/test logs and `FINAL_ACTIVE_STUDY_AUDIT.json` are in `local_review/prism_120_primary_20261028`. Earlier empirical results, comparisons and superseded reports are preserved in private, checksummed, read-only archives. They are not manuscript-facing inputs. Git history was not rewritten.

`publication/PI_REVIEW.md` records the required IRB, image-use, synthetic-lineage and existing repository-history remediation review. No real GMDB photographs, per-person rows or landmarks are included in new manuscript exports. The dashboard returns aggregate histograms. The live Google Doc was not edited.

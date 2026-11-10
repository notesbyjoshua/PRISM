# Conditional synthetic fidelity supplement

These are aggregate publication candidates from the current 120-feature measurements. They describe the observed PDIDB release relative to selected age-eligible GMDB patients; they do not add clinical subjects or independent confirmation to the primary study. The mentor's GMDB-trained StyleGAN statement is distinguished from unresolved checkpoint-specific provider provenance in the manifest.

Read `predeclared_analysis_plan.md`, `synthetic_consistency_results.md` and `synthetic_consistency_manifest.json`. The eligibility audit includes all 48 primary cohorts and the excluded generic synthetic 22q11.2 label; a blank real n for that label means the mapping is unresolved, not that GMDB has no related records. Thirty-nine primary cohorts have no synthetic label. The explicit mapping is syndrome-level only and is not source-patient linkage. Clinical review of cohort equivalence remains pending.

`feature_fidelity_summary.csv` contains n >= 10 per-source aggregate statistics for pooled, model-tag and pooled winsorized strata. It contains no individual measurements or image/patient identifiers. `metadata_audit.csv` contains category counts, including missing values; ages/sex/race are absent for every synthetic record. Figures 8–10 are in `manuscript/figures/`; their source CSVs and captions are supplied. No photographs or generated image bytes are included. Existing export and PI publication review rules apply.

From the repository root, using an environment with `publication/requirements.txt` installed:

```sh
python stats_testing/synthetic_consistency.py --active-config ACTIVE_ANALYSIS.json --output-dir local_review/NEW_SYNTHETIC_FIDELITY
python publication/scripts/final_research_tables.py
python publication/scripts/final_research_figures.py
python publication/scripts/final_research_text.py
python stats_testing/validate_synthetic_consistency.py --report local_review/NEW_SYNTHETIC_VALIDATION.json
python publication/scripts/verify_final_research.py
python publication/scripts/verify_reproducible_outputs.py --write-manifest
python publication/scripts/verify_reproducible_outputs.py
```

The first and fifth commands require authorized private measurements. The first command prohibits network connections and never reprocesses images. Primary results remain frozen. Review intentional aggregate regeneration before updating the package manifest. The private full-validation command is documented in the completion report; public CI checks aggregate integrity only. All calculations are deterministic without resampling; the reserved seed is recorded for traceability. The n >= 20/30 real-cohort sensitivity cutoffs are explicitly exploratory additions after initial pooled results.

# Reference-selection stability diagnostic

This isolated secondary workstream preserves the November 15 primary study, selectors, case cohorts, crosswalk, original sensitivity checkpoints, manuscript and publication files. Real outputs and all detailed/small-cell tables stay in verified ignored `local_review/reference_stability_diagnostic/`. No images are processed, no clinical annotation study or alternative normalization is run, and nothing is committed or uploaded.

## Reproduce

From the repository root, using the existing environment:

```sh
.venv/bin/python -m pytest research/reference_stability_diagnostic/tests -q
.venv/bin/python research/reference_stability_diagnostic/run.py --dry-run
.venv/bin/python research/reference_stability_diagnostic/run.py --stage audit
.venv/bin/python research/reference_stability_diagnostic/run.py --stage B
.venv/bin/python research/reference_stability_diagnostic/run.py --stage C --repetitions 1000
.venv/bin/python research/reference_stability_diagnostic/run.py --stage all --resume
```

The completed output directory already has checkpoints; use `--resume` for existing B/C computations, or choose a fresh ignored directory such as `--output-dir local_review/reference_stability_diagnostic_new_run`. The default seed is 20261116 and repetitions are 1000. A run with fewer than 1000 repetitions is labeled PILOT and cannot receive a completed final-experiment status. The frozen plan records repetition counts, so changing them requires a fresh run rather than extending after observing favorable results.

`--dry-run` validates input/checkpoint provenance and feasibility, prints the plan and creates or overwrites no result files. `--stage A` reuses/audits the available original 1,000 draws; it does not regenerate them. Missing or corrupt original arrays block this run explicitly rather than substituting an unlabeled replica. In the completed study those arrays were available and correctly checksummed.

`engine.py` implements weighted shared-reference rank contributions, exact-U generating polynomials, tie corrections, centered numerical moment sums and effect calculations. This numerical centering is an implementation of ordinary sample variance, not a change to facial-feature normalization. Every draw constructs all 9,600 p values and one full-family BH correction. Exact/asymptotic MW auto rules and exact gamma Hedges correction are retained. Direct SciPy computations on a deterministic real subset and independent statsmodels BH checks validate the optimized calculations, including with-replacement duplicates/ties.

The exact-U lookup is the coefficient distribution of the Gaussian-binomial polynomial `product((1-z^(m+i))/(1-z^i), i=1..n) / choose(n+m,n)` for the smaller n. This avoids constructing huge artificial rank matrices at full reference sizes. Public SciPy tests check small and full-pool-size configurations.

Checkpoints reuse the existing atomic checksummed chunk implementation. Resume requires identical validated source, execution code (`run.py`, `engine.py`) and analysis-plan signatures. Original sensitivity inputs are read-only. Plotting/report postprocessing is isolated in `reporting.py`; its current hash is separately recorded. To regenerate only reports after a layout edit, without recomputing numerical arrays:

```sh
.venv/bin/python -c "import sys; from pathlib import Path; sys.path.insert(0,'research/reference_stability_diagnostic'); from reporting import report; report(Path('local_review/reference_stability_diagnostic'))"
```

## Designs and interpretation

- A: the exact original primary sample and original 1,000 reference-subset checkpoint draws. The selected sample, 982 primary discoveries and 76/193/580 original stable flags are verified before new work. Europe used its full pool; fixed European g is not selection stability.
- B: one deterministic analysis using every eligible deduplicated category reference; a separate 1,000-repetition, full-size with-replacement conditional reference bootstrap for each category, including Europe. Cases stay fixed. Duplicates/ties and truly distinct multiset draws are recorded.
- C: paired doubling stages, capped by the actual pools, with Europe fixed. The original quota stage reuses A; each of four new non-full stages receives 1,000 shared-reference subset draws. Full-pool selection has one unique subset and no recurrence/dispersion estimate. Its point values are separately labeled.

The original >=80% q<.05 recurrence and >=95% original-g-sign rule is fixed. Finite nonzero g is required for directional success. Null/zero/undefined values do not silently become stable. Reference quantiles are selection distributions, not clinical confidence intervals. A/C subsets and the B bootstrap quantify different conditional uncertainties; the finite-pool subset variance must vanish as a pool is exhausted. Neither scheme provides independent replication or removes the perfect GMDB-versus-FairFace source association.

## Output map

- `analysis_plan.json`, `source_and_execution_manifest.json`, `verification_report.json`: frozen design, source/commit/software/code/seed provenance, signatures, hashes and executed checks.
- `original_118_failure_reasons.csv`, `original_387_African_Asian_discovery_diagnosis.csv`, `failure_mechanism_counts.csv`: exact diagnosis keys and mutually exclusive operational failure labels, including zero categories.
- `original_reference_per_feature.csv` / `original_reference_summary.csv`: original 1,000-draw summaries and per-feature reference-moment variability.
- `A_vs_B_primary_effects.csv`, `A_vs_B_category_summary.csv`, `A_vs_B_region_summary.csv`, `A_vs_B_case_size_summary.csv`: deterministic full-reference changes, retention, gains/losses and descriptive strata.
- `reference_size_ladder_summary.csv`, `reference_size_ladder_per_feature.csv`, `incremental_ladder_changes.csv`: full size-dependent results. Deterministic full-size uncertainty fields are NA, not artificial 100% stability or zero-width confidence intervals.
- `B_reference_bootstrap_*`, `C_stage_*`, `full_reference_*`, `checkpoints/`: separate schemes, full-family numerical results and resumable checksummed aggregate arrays. No patient rows are exported.
- `followup_evidence_all_387_original_discoveries.csv`, `internally_consistent_followup_candidates.csv`: all original African/Asian findings with pass/fail evidence per scheme; the intersection is a descriptive follow-up filter, not a new p-value or clinical robustness criterion.
- `FAILURE_DIAGNOSIS.md`, `REFERENCE_DESIGN_COMPARISON.md`, `RESULTS_AND_LIMITATIONS.md`, `manuscript_ready_summary.md`: observed findings, stage status, cautious decision guidance and manuscript-ready draft text; manuscript files are unchanged.
- `COMPUTATIONAL_FORMULAS.md`, `FIGURE_CAPTIONS.txt`, `*_plot_values.csv`, `subset_stability_curve_plot_values.csv`: exact formulas, captions and plotting inputs. Other plotted columns are in the ladder summary.
- `figures/`: failure mechanisms, original q versus recurrence, stability versus reference size, g-selection dispersion versus size, and original-versus-full-reference g, each PNG 300 DPI plus vector SVG/PDF.

All counts remain correlated measurements rather than independent phenotypic discoveries. Larger references can reduce selection dispersion without preserving every original discovery. Any future primary replacement requires explicit scientific approval and a transparent new analysis plan.

Per-feature column convention: `p`, `q`, `hedges_g`, `case_n`, and `reference_n` preserve the original A values. `evaluated_reference_n`, `median_g`, `g_q025`, `g_q975`, `median_q`, and recurrence/sign fields describe the labeled resampling design. Full deterministic rows instead expose `deterministic_g`/`deterministic_q` and leave distribution/recurrence estimates NA; `A_vs_B_primary_effects.csv` contains the corresponding `B_p`, `B_q_global_9600`, `B_g` and `B_reference_n`. Pass flags are numeric 0/1 with explicit NA when not evaluable.

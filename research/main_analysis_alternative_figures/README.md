# Alternative main-analysis figures

From the repository root:

```sh
.venv/bin/python research/main_analysis_alternative_figures/generate_facekit_alternative_figures.py --figure forest
.venv/bin/python research/main_analysis_alternative_figures/generate_facekit_alternative_figures.py --figure regions
.venv/bin/python research/main_analysis_alternative_figures/generate_facekit_alternative_figures.py --figure both
.venv/bin/python -m pytest research/main_analysis_alternative_figures/tests -q
```

Defaults: 2,000 conditional case-only bootstrap repetitions, seed 20261115, private output `local_review/main_analysis_alternative_figures/`. The region-only branch never calls the patient loader. Each branch writes its own plot/aggregate CSV/caption plus shared selection and manifest; it does not regenerate the other figure.

Selection reuses the exact count/category-only selector from Panels A/B. The original 9,600-test primary output must match its manifest hash, study dimensions, original 982 discoveries, feature definitions and exact per-cell reference allocation. No original U/p/q values are recalculated. Forest point g and group means are checked against the authorized source observations; the reference sample is replayed from the manifest's quotas/seed and compared to every stored selected path, image ID, content hash and category. Patient-level arrays remain in memory.

The companion script was unavailable. This implementation freezes the four requested features and an explicit nine-to-seven region map in code: merge orbital dimensions with ocular/iris/gaze; merge global face shape with cheeks; retain forehead, eyebrows, nose, mouth/lips and chin/jaw separately. The resulting feature denominators are 4, 20, 45, 17, 15, 5 and 14. These are planned feature counts; every exported regional row also carries its actual supported denominator. Missing/ineligible rows have NA values, not zero discoveries.

Forest uncertainty uses joint case-row resampling across four features and NumPy SeedSequence([seed, zero-based selected-syndrome index, zero-based category index]). References never resample. These conditional percentile intervals omit reference uncertainty, selection, multiplicity, acquisition effects and dependencies; marker fill alone represents the original global BH decision. The region chart applies the original q flags without any regional correction. Neither figure demonstrates clinical validation, causal ancestry effects or independent phenotypic discoveries.

Outputs are 300-DPI PNG and vector PDF/SVG with aggregate tables, figure-specific captions, selection and manifest. All real outputs remain private and ignored. No photos, patient identifiers, per-person measurements, alternative normalization, model fitting or manuscript changes are produced. Code and synthetic tests are not committed automatically.

# Updated sensitivity figure

The composite compares the original reference allocation with larger resampleable subsets in Panel C. Panels A and B retain their original results. All figures are available as 300 DPI PNG and vector PDF/SVG.

- [Composite PNG](sensitivity_ABC_composite.png), [PDF](sensitivity_ABC_composite.pdf), [SVG](sensitivity_ABC_composite.svg)
- [Panel A PNG](sensitivity_panel_A.png), [PDF](sensitivity_panel_A.pdf), [SVG](sensitivity_panel_A.svg)
- [Panel B PNG](sensitivity_panel_B.png), [PDF](sensitivity_panel_B.pdf), [SVG](sensitivity_panel_B.svg)
- [Panel C PNG](sensitivity_panel_C.png), [PDF](sensitivity_panel_C.pdf), [SVG](sensitivity_panel_C.svg)
- [Caption](FIGURE_CAPTION.txt), [exact Panel C values](panel_C_reference_comparison_values.csv), [plotting inputs](reference_comparison_inputs.json), [release provenance](release_manifest.json)

| Recorded category | Original reference n | Original stable | Larger reference n | Larger-subset stable |
|---|---:|---:|---:|---:|
| African | 52 | 76/135 (56.3%) | 832 | 112/135 (83.0%) |
| Asian | 190 | 193/252 (76.6%) | 3,040 | 233/252 (92.5%) |

Combined stability changes from 269/387 (69.5%) to 345/387 (89.1%). Both designs use 1,000 without-replacement reference draws with fixed cases, retaining the same original discovery denominators. Stable means global BH q < 0.05 in at least 80% of draws and original g-sign agreement in at least 95%. Each draw uses the full 9,600-test family. The larger allocation is selected by reference size, not by stability outcome.

European reference photographs remain fixed at n = 647 and are omitted from the stability bars. The deterministic full-pool diagnostic and separate with-replacement bootstrap are not plotted. Near-full subsets reduce finite-pool selection variation; greater recurrence does not establish independent replication, diagnostic accuracy or elimination of image-source confounding. Recorded demographic categories are dataset annotations, not verified genetic ancestry.

The released inputs contain only category-level counts and the allocation plan. They contain no patient, image, syndrome/feature or small-cell records. Private diagnostic checkpoints, detailed tables and original figure archives remain under ignored `local_review/`.

Reproduce from the repository root using the existing project environment:

```sh
.venv/bin/python research/sensitivity_figure/generate_prism_sensitivity_figure.py --panel all --output-dir local_review/sensitivity_summary_figure_reproduction
.venv/bin/python -m pytest research/sensitivity_figure/tests -q
```

Reproduction uses the pinned released inputs by default and writes only to the ignored directory. It does not rerun inference. `release_manifest.json` records exact output hashes, input provenance and the reproduction checks. PNG comparisons use the recorded rendering environment; PDF timestamps and SVG identifiers can differ between runs, and fonts or library versions can change rendering.

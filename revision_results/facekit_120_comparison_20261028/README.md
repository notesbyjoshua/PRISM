# FaceKit 120 versus historical 125: local comparison audit

**Historical audit stage, superseded.** The later authorized extraction and validated
120-feature rebuild are complete and active. See
[`../rebuild_120_20261028/README.md`](../rebuild_120_20261028/README.md).
The blocked-status text below records the evidence available at the earlier stage;
it does not describe the current analysis or permissions.

**Decision: INDETERMINATE.** The latest FaceKit submodule is pinned to
`5764bdd4af9360b55a046a16579fb5a890a8182f` (upstream main checked October 28, 2026).
Original selected photographs and compatible saved landmarks were not found in
the audited repository locations, and image-reprocessing permission is unconfirmed.
The 120-feature patient reanalysis is blocked. Historical 125-feature results remain primary.

**Later offline cache discovery:** the author identified two GMDB Arrow shards outside
the repository. The follow-up input audit in
`local_only/hf_cache_audit_20261028/audit_results/README.md` confirms 10,980 embedded,
decodable restored images, including all 1,850 selected historical childhood images,
and matching historical metadata. A verified owner-only backup is in the adjacent
`private_backup/` directory. This supersedes the earlier report specifically on GMDB
image availability. Historical image-byte provenance, FairFace images and separate
reprocessing approval remain unresolved. No extraction pilot or new statistical
analysis was run; the author requested an input audit only at this stage.

The completed package is in `local_only/run_20261028_final/`. Start with
`PRISM_FaceKit_120_vs_125_Comparison_Report.pdf` or its Markdown equivalent.
It contains the full feature crosswalk/formula audit, historical numerical checks,
explicit unavailable 120-feature metrics, clinical highlights, source-linked mapping
graphic, privacy inventory, decision and blockers, claim-to-source table and manifests.
All generated results are Git-ignored and require review before redistribution.

## What ran

- Reconstructed the frozen selection: 1,850 patients, 942 reference images, 49 eligible syndromes.
- Recomputed all 6,125 historical U/g/p/q point estimates and 811 global-BH discoveries.
- Independently verified rank sums, tie-corrected p-values and BH; bootstrap intervals were not rerun.
- Confirmed 120 finite base outputs on synthetic landmarks, without patient-image extraction.
- Audited every historical feature and preserved the hairline-to-forehead many-to-one relation.
- Strict full-pipeline equivalences certified: 0/120 current and 0/125 historical; changed alignment prevents unconditional certification.
- Three 90% agreement fractions: unavailable, with null numerators and denominators.
- Updated FaceKit, refreshed its editable installation without dependency upgrades, and checked dependency consistency and relevant tests.

The application submodule was updated from `63a35875bd7e633d67a495af5f7248cc74b82667`.
Its old commit is retained in the submodule branch `codex/backup-before-update-20261028`.
The upstream update changes packaged feature-group defaults and adds related tests;
the feature extractor and canonical-depth numeric array remain unchanged relative to that prior application checkout.
The audited source is not claimed to be the exact code revision used for the paper.

## Reproduce

Run from the PRISM root. The existing `.venv` supplies NumPy, pandas, SciPy,
statsmodels, Matplotlib and pytest; `environment_freeze.txt` records installed versions.
The PDF builder uses the bundled document Python with ReportLab. No new analysis
dependencies were installed. Commands only read local study inputs; the public-source
download sends no study data.

The source archive is pinned by `local_only/source_acquisition.json`. To recover it
when unavailable, download the GitHub source archive at the exact commit above
and extract it to `local_only/sources/FaceKit-<commit>/`. Verify its recorded SHA-256.
The checked-out FaceKit submodule alone is sufficient to inspect the same source,
but `audit.py` deliberately requires that isolated source snapshot.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python \
  revision_results/facekit_120_comparison_20261028/audit.py \
  --output revision_results/facekit_120_comparison_20261028/local_only/fresh_run

.venv/bin/python revision_results/facekit_120_comparison_20261028/freeze_selection.py \
  revision_results/facekit_120_comparison_20261028/local_only/fresh_run

PRISM_FACEKIT_AUDIT=revision_results/facekit_120_comparison_20261028/local_only/fresh_run \
  .venv/bin/python -m pytest revision_results/facekit_120_comparison_20261028/test_audit.py -q

.venv/bin/python revision_results/facekit_120_comparison_20261028/plot_mapping.py \
  revision_results/facekit_120_comparison_20261028/local_only/fresh_run

/Users/joshua/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 \
  revision_results/facekit_120_comparison_20261028/build_report.py \
  revision_results/facekit_120_comparison_20261028/local_only/fresh_run
```

`audit.py` refuses an existing output directory and any destination outside this
workflow's `local_only/`. Never force-add private outputs. PDF rendering and inspection
are required before delivering a rebuilt report. Current results, input hashes,
scripts, environment and frozen-output integrity are covered by the audit tests.

Relevant update checks:

```sh
.venv/bin/python -m pip check
.venv/bin/facekit --help
.venv/bin/python -m pytest FaceKit/tests/test_canonicalize_pose.py \
  FaceKit/tests/test_api_plugin_registry.py FaceKit/tests/test_extract_features_custom.py \
  FaceKit/tests/test_score.py FaceKit/tests/test_feature_groups.py tests -q
```

No original patient images, 120-feature study results, live manuscript patches,
external document edits, commits, pushes or publication are part of this completed
feasibility package. Recovering authorized inputs and confirming a paper-specific
source/model pin are the next required steps.

## Subsequent extraction-only run

After the initial feasibility and cache audits, the author authorized local extraction
from GMDB and downloading the official FairFace margin025 inputs. Both full extraction
scans completed. See `local_review/feature_extraction_120/20261028_current120/README.md`
and `data/dataset_creation/CURRENT_120_EXTRACTION.md` for the new measurements, QC,
sampling manifests and restart commands. The older feasibility package above remains
an unchanged record of that earlier stage. No statistical comparison, switch decision,
live manuscript edit or public release was performed in the extraction stage.

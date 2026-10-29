# PRISM statistical study using FaceKit

Statistical investigation using FaceKit measurements, with the validated 48-syndrome, 5,760-test current analysis. Primary reference: 3,409 FairFace image paths. No synthetic observations in primary inference. Feature definitions and figures are current; older source material is archived separately. No patient photographs or per-patient tables are included.

GitHub Actions verifies the committed aggregate package against the current pinned
FaceKit source, exact 120-feature schema, primary BH/effect arithmetic, matching
manuscript exports and `ARTIFACT_MANIFEST.json`. It runs without private inputs,
landmark models or dataset credentials. It does not re-extract patient images or
claim full patient-level replication on the hosted runner.

Run `python publication/scripts/verify_reproducible_outputs.py` after installing
`publication/requirements.txt` and initializing the pinned FaceKit submodule.
The workflow's temporary Markdown report is added to the Actions summary; missing
reports cannot obscure an earlier failure. Only use `--write-manifest` after
intentionally regenerating and reviewing the aggregate exports locally. Full
extraction and analysis require the authorized private run described in
`revision_results/rebuild_120_20261028/README.md`.

The active manuscript, tables and figures contain only the current study. Concordance points are descriptive and its assessed bootstrap uncertainty is excluded. See `PI_REVIEW.md` before submission or any image reuse.

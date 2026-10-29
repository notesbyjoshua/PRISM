# Private current FaceKit extraction

`extract_current_120.py` produces local 120-feature measurements from official
FairFace margin025 images or the verified private GMDB Arrow backup. It performs
extraction and QC bookkeeping only. It never runs PRISM hypothesis tests, trains
models, uploads data or changes source measurement tables.

The completed October 28 run is documented in
`local_review/feature_extraction_120/20261028_current120/README.md`.
All inputs, person-level measurements, checkpoints, logs and manifests in that
directory are private and Git-ignored. Public redistribution requires separate review.

Use the existing system Python, which has PyArrow, MediaPipe and Pillow. The script
loads the local FaceKit source directly. It requires commit
`5764bdd4af9360b55a046a16579fb5a890a8182f`, the recorded cached landmark asset checksum,
exactly 120 base features, and no feature plugins. Extraction blocks network access
and cannot silently download new weights. FairFace downloads are a separate,
explicitly authorized preparation step.

## New run

Choose unused directories under `local_review/feature_extraction_120/`. For FairFace:

```sh
/Library/Frameworks/Python.framework/Versions/3.11/bin/python3 \
  data/dataset_creation/extract_current_120.py fairface \
  --output local_review/feature_extraction_120/NEW_RUN/fairface \
  --archive PATH/TO/fairface-img-margin025-trainval.zip \
  --train-labels PATH/TO/fairface_label_train.csv \
  --val-labels PATH/TO/fairface_label_val.csv
```

For GMDB, use a backup directory with the existing verified `backup_manifest.json`
in its parent directory:

```sh
/Library/Frameworks/Python.framework/Versions/3.11/bin/python3 \
  data/dataset_creation/extract_current_120.py gmdb \
  --output local_review/feature_extraction_120/NEW_RUN/gmdb \
  --gmdb-backup PATH/TO/private_backup
```

`--pilot-only` runs demographic/selection representatives and repeat checks first.
Then repeat the command with `--resume` and without `--pilot-only` to finish the
entire candidate pool. Without `--pilot-only`, the pilot and full scan run in sequence.
Inputs, model, code, environment and feature schema must match before any restart.
The completed run exercised the actual pilot-to-full restart for both datasets.

## Output and selection policy

FairFace keeps age as source category `3-9`. Every unique candidate is screened
before sampling. A selects up to 500 per source race category without replacement;
B is the full accepted pool. C selects the same count in every race x Female/Male
cell, using `min(250, minimum cell capacity)`. Unknown sex is preserved in A/B and
excluded from C. These samples overlap. In the shortfall CSV, C's `target` is the
effective uniform target and `shortfall` refers to the nominal 250-per-cell goal.

The shared FairFace master contains accepted measurements and A/B/C flags. The
three sampling manifests reference those same records. A separate historical
942-key cross-reference is provenance bookkeeping and does not force A/B/C membership.

GMDB exports every cache record, including explicit failures, and a separate
quality-passed subset. Finite measurements from pose-failed images remain in the
all-records table. Failed detection/ambiguity rows have missing vectors, never values
inferred from historical columns. Historical primary membership flags are preserved.

The detector screens up to two faces for ambiguity, then uses the standard one-face
detector's 478 floating-point pixel landmarks and matrix with current FaceKit
geometry. Coordinates, original photographs and embeddings are never exported.
The frontal gate is 15/15/10 degrees; accepted vectors must contain all 120 finite
measurements and use full matrix-based pose correction. Missing occlusion,
expression, identity and clinical-quality assessments are not invented.

## Validation

```sh
/Library/Frameworks/Python.framework/Versions/3.11/bin/python3 -m pytest \
  tests/test_current_120_extraction.py -q

/Library/Frameworks/Python.framework/Versions/3.11/bin/python3 \
  data/dataset_creation/validate_current_120.py \
  local_review/feature_extraction_120/20261028_current120
```

The validator checks every serialized measurement against the checkpoint, accepted
feature finiteness, source labels, all candidate keys, QC counts, manifests, pinned
assets, original cache integrity and historical CSV hashes. It does not compute
disease effects or significance. Original image-byte equivalence to the historical
extraction and exact paper execution provenance remain unverified.

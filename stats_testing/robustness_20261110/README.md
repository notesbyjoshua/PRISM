# Private robustness extension, November 10, 2026

Implements the high-value research handout without modifying/promoting PRISM_120_20261028. Code contains no patient observations. All empirical outputs stay in Git-ignored local_review/prism_robustness_20261110/. Candidate and low-support aggregates require PI/disclosure review.

Read ANALYSIS_PLAN.md and the private 00_preflight.md. common.py verifies schema, unique recorded patient selection, canonical taxonomy, measured-source/model pins, and frozen files/private input hashes. Missing inputs fail closed. Protected computations block network connections. No credentials, pixels, extractor or synthetic records are needed.

## Commands used

From the repository root, use the existing authorized local environment:

    export OPENBLAS_NUM_THREADS=1
    export OMP_NUM_THREADS=1
    export PYTHONPATH=.:FaceKit/src
    export MPLCONFIGDIR=local_review/prism_robustness_20261110/matplotlib

Run stages sequentially with milestone review:

    local_review/ci_120_fix/venv/bin/python -m stats_testing.robustness_20261110.evidence

Review Milestone A before continuing:

    local_review/ci_120_fix/venv/bin/python -m stats_testing.robustness_20261110.source_bias --repetitions 100
    local_review/ci_120_fix/venv/bin/python -m stats_testing.robustness_20261110.dependence --bootstrap 100
    local_review/ci_120_fix/venv/bin/python -m stats_testing.robustness_20261110.stability
    local_review/ci_120_fix/venv/bin/python -m stats_testing.robustness_20261110.geometry
    local_review/ci_120_fix/venv/bin/python -m stats_testing.robustness_20261110.review
    local_review/ci_120_fix/venv/bin/python -m stats_testing.robustness_20261110.sensitivity_summary
    local_review/ci_120_fix/venv/bin/python -m stats_testing.robustness_20261110.literature
    local_review/ci_120_fix/venv/bin/python -m stats_testing.robustness_20261110.plots
    local_review/ci_120_fix/venv/bin/python -m stats_testing.robustness_20261110.reporting
    local_review/ci_120_fix/venv/bin/python -m stats_testing.robustness_20261110.validate --full-private-baseline

Literature retrieval was a separate public-source review. literature.py renders its reviewed registry offline; it does not claim systematic coverage of all 955 associations. Full-text, indexed-passage, indirect and abstract-only access states remain explicit.

The read-only full-validation wrapper executes the unchanged frozen validator, routing only its two expected final report writes into the new private directory. It checks the earlier reports remain byte-identical and does not modify numerical functions. New validation independently recalculates BH, point estimates and two control repetitions, and checks scopes, schemas, figure hashes and claim sources.

Public-fixture tests and aggregate verification:

    local_review/ci_120_fix/venv/bin/python -m pytest tests -q --tb=short
    local_review/ci_120_fix/venv/bin/python -m pytest revision_results/reviewer_followup_2026/test_extended_analysis.py -q --tb=short
    local_review/ci_120_fix/venv/bin/python publication/scripts/verify_reproducible_outputs.py --report local_review/prism_robustness_20261110/validation/public_aggregate_verification.md

No manifest-writing flag or new remote CI run is used. After review, with the required completion report present:

    local_review/ci_120_fix/venv/bin/python -m stats_testing.robustness_20261110.provenance
    local_review/ci_120_fix/venv/bin/python -m stats_testing.robustness_20261110.provenance --check

The private manifest hashes final sources/drafts, records environment/seeds/input digests, and states blockers. Rerendering requires revalidation and regeneration of that private manifest.

## Audit and interpretation

The first negative-control allocation left pitch-sorted case pools ordered. That v1 checkpoint is invalid and retained solely for audit; corrected v2 randomizes allocation within pools. Final control tables derive exclusively from v2. A dependency-parser substring bug was corrected with exact normalizer tokens; final cluster caches use v2 signatures. Both defects have regression fixtures. Superseded checkpoints are not evidence.

Case-only bootstrap ranges are conditional on fixed empirical comparators and remain descriptive after undercoverage checks. No concordance-bootstrap interval is computed. New geometry/quality families never replace primary q. Rankings are post hoc, without a validated score or held-out advantage.

Raw decomposition is BLOCKED_PRIVATE_INPUT; real-image transformations are BLOCKED_IMAGE_REPROCESSING_APPROVAL. New method development is NOT_JUSTIFIED, a conditional scientific no-go. Public release and manuscript integration require the reviews specified in the handout and publication/PI_REVIEW.md.

# Covariate feasibility

This audit evaluates whether fields recorded in both source datasets can support a narrow sensitivity analysis.
It does not treat harmonized dataset labels as genetic ancestry and does not treat the FairFace `3-9` category as
an exact age.

## Analytic cohorts

- Disease: 1,850 selected patients aged 3 to under 10 under the frozen selection rule.
- Reference: 942 deduplicated FairFace images with a passing stored frontal-quality flag.

## Decision

A limited sensitivity analysis is feasible for the three predefined highlighted candidates. It uses yaw, pitch, roll, sex / gender label, broad demographic group. Stored frontal quality is an inclusion rule rather than a covariate because all retained
records pass it. Detailed demographic labels are not used because syndrome-specific cells are sparse.

Exact-age adjustment is not feasible. Disease ages can be reconstructed conditionally from year and month fields,
whereas the reference dataset supplies the categorical label `3-9`. The sensitivity models therefore retain the
frozen childhood cohort but do not include age as a numeric or matched covariate.

## Field decisions

| Covariate | Comparable semantics | Use in model | Reason |
|---|---|---|---|
| stored frontal-quality information | yes as a common eligibility rule | no as a covariate | Both analytic cohorts require a passing stored flag, leaving no modeled variation. |
| yaw | yes | yes | Numeric degrees from the same FaceKit pose fields; modeled continuously. |
| pitch | yes | yes | Numeric degrees from the same FaceKit pose fields; modeled continuously. |
| roll | yes | yes | Numeric degrees from the same FaceKit pose fields; modeled continuously. |
| sex / gender label | limited; harmonized dataset labels | yes | Only common Male/Female labels are used; labels are dataset annotations. |
| broad demographic group | limited; established publication mapping | yes | Uses the existing broad harmonization and excludes unmapped labels. |
| detailed demographic group | limited; established publication mapping | no | Comparable mapped labels exist, but sparse syndrome-specific cells make the broad layer more stable for three models. |
| age information / age category | no for exact age | no | Disease age can be conditionally reconstructed; FairFace supplies the 3-9 category, not exact age. |

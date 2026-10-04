# Cohort balance report

This report describes recorded covariates in the frozen childhood disease cohort and the deduplicated FairFace
reference cohort. It does not report a single balance score because disease status and source dataset are perfectly
confounded, and exact age is not comparable across sources.

## Pose

| Variable | Cohort | N available | Mean | SD | Median | Q1 | Q3 |
|---|---|---:|---:|---:|---:|---:|---:|
| pose_yaw | disease | 1850 | -0.037 | 5.107 | -0.104 | -3.019 | 2.879 |
| pose_yaw | reference | 942 | -0.634 | 7.032 | -0.436 | -5.428 | 4.040 |
| pose_pitch | disease | 1850 | 1.014 | 6.876 | 1.307 | -3.682 | 6.134 |
| pose_pitch | reference | 942 | 3.697 | 7.149 | 4.583 | -0.944 | 9.330 |
| pose_roll | disease | 1850 | -0.180 | 1.437 | -0.142 | -1.029 | 0.689 |
| pose_roll | reference | 942 | 0.174 | 2.058 | 0.142 | -0.947 | 1.492 |

## Categorical distributions

Sex-label, broad demographic-label, and available age-category counts are recorded in `cohort_balance.csv`.
These are dataset annotations. The disease exact-age field and reference age category are shown separately
and are not presented as directly balanced measurements.

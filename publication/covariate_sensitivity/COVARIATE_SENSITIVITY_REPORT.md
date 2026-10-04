# Covariate sensitivity report

## Predefined analysis

Before fitting any adjusted model, the analysis set was fixed to three manuscript-highlighted examples:

- Cornelia de Lange syndrome / `philtrum_length`
- Noonan syndrome / `inter_pupillary_distance`
- Rubinstein-Taybi syndrome / `malar_bulge_mean`

Each exploratory model is ordinary least squares with HC3 robust standard errors. The estimand is the adjusted
disease-versus-reference mean difference in the historical normalized feature. Models include one selected image
per age-eligible disease patient and deduplicated reference images with complete sex-label, broad demographic-label,
yaw, pitch, and roll fields. Missing modeled covariates use complete-case exclusion. The three disease coefficients
form one predefined exploratory family and receive Benjamini–Hochberg correction across three tests.

The original Mann–Whitney analyses remain primary. These models neither replace them nor change their correction family.

## Results

| Syndrome | Feature | Disease N | Reference N | Adjusted coefficient | HC3 95% CI | Raw p | Exploratory q |
|---|---|---:|---:|---:|---|---:|---:|
| Cornelia de Lange syndrome | `philtrum_length` | 62 | 942 | 0.052347 | [0.0382491, 0.0664448] | 3.399e-13 | 1.02e-12 |
| Noonan syndrome | `inter_pupillary_distance` | 26 | 942 | 0.0371447 | [0.0248052, 0.0494842] | 3.637e-09 | 5.455e-09 |
| Rubinstein-Taybi syndrome | `malar_bulge_mean` | 24 | 942 | -0.00641233 | [-0.0096, -0.00322467] | 8.058e-05 | 8.058e-05 |

These analyses assess sensitivity to recorded covariates but cannot separate disease effects from unmeasured
acquisition, curation, or source-dataset effects because disease status and dataset source remain perfectly confounded.

The adjusted coefficients should be read only as a recorded-covariate sensitivity check. They do not establish
replication, eliminate dataset effects, or support new clinical phenotype claims.

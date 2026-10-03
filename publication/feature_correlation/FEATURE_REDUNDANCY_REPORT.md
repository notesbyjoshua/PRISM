# Feature redundancy report

The reference-cohort correlation analysis uses 942 distinct FairFace image paths and the exact ordered 125-feature schema in `revision_results/final_20260914/exact_feature_names.csv`. No 120-feature output or regenerated FaceKit measurement is used.

## Results

- 268 of 7,750 feature pairs have |Spearman rho| >= 0.70.
- 212 pairs have |rho| >= 0.80, 131 have |rho| >= 0.90, and 68 have |rho| >= 0.95.
- Average-linkage clustering on `1 - |rho|`, cut at distance 0.30, yields 59 clusters. 22 clusters contain more than one feature.
- 7 clusters span more than one anatomical region, showing that shared scale factors and geometry can link measurements that have different anatomical labels.

## Interpretation

The 125 tests within each syndrome are not biologically independent. Bilateral left/right measurements, mean and asymmetry derivatives, related distances, and features that share normalization denominators can move together. A cluster of significant results can therefore describe one underlying facial pattern rather than several independent discoveries. The matrix and cluster table characterize that dependence; they do not redefine the prespecified global Benjamini-Hochberg family and should not be used to retroactively reduce the multiplicity correction.

Spearman correlation captures monotonic dependence in the reference cohort only. It does not prove that two features use identical landmarks, and it can differ in a disease cohort. Exact landmark overlap would require a separate source-code-level dependency map of every feature formula.

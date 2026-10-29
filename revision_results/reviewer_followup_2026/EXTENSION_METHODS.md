# Computational extension contract

Reusable numerical functions support the current FaceKit-based PRISM statistical study through `revision_results/rebuild_120_20261028/extensions.py`. Direct report-building entry points in this utility module are disabled. The active correction-family definitions are in `stats_results/extensions_120/FAMILY_DEFINITIONS.json`.

Within-GMDB specifications retain separate 5,760-test families. The adjusted, pairwise and age families retain 5,880, 22,320 and 11,760 planned hypotheses, including unsupported cells with correction p=1 and missing reported estimates/q-values. Adjustment uses HC3 with the declared complete-case covariates; effect and bootstrap definitions remain fixed. Synthetic tests check these numerical contracts without patient inputs.

Concordance is descriptive. Its assessed resampling method is biased and supplies no active confidence intervals or uncertainty summaries. Reference-feature correlation clustering is descriptive and does not reduce the correction denominator. Synthetic PDIDB observations are excluded from real-patient inference.

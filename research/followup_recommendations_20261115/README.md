# Significant follow-up candidates and exploratory category heterogeneity

```sh
.venv/bin/python research/followup_recommendations_20261115/recommend.py
.venv/bin/python -m pytest research/followup_recommendations_20261115/tests -q
```

`analysis_plan.json` fixed the new mean-contrast procedures before examining their results. The first shortlist reuses the ten previously locked/validated candidates; it does not certify literature-wide absence. It requires exact primary and audit keys and unchanged values. Literature/source boundaries remain explicit in the report.

The new secondary family contains all 1,920 eligible same-syndrome, same-feature category pairs. The estimand is the difference between two case-minus-reference mean gaps in original feature units. A four-independent-group Welch–Satterthwaite contrast uses the stored primary means/SDs/counts. Original Mann–Whitney p/q/g values are neither substituted nor modified. BH and supplementary BY apply to the complete new family; filters do not produce a new correction.

The follow-up display requires at least one significant original primary association and keeps one comparison per syndrome/measurement family. A separate list additionally requires a missing reviewed catalog annotation. All 40 BH and all 11 BY contrasts remain in separate complete tables, including rows where both individual primary tests were nonsignificant. No claim of independent traits follows from correlated bilateral/derived metrics. Mean-gap contrasts are not tests of Hedges g equality or of the original rank-based estimand.

All data-derived outputs are private, Git-ignored under `local_review/followup_recommendations_20261115/`, including the report, ten literature-review candidates, cross-category lists, every planned contrast, frozen plan and provenance. The implementation reads aggregate tables only; source/run hashes are checked and inputs remain unchanged.

These are approximate exploratory tests with sparse groups, distribution/independence assumptions, source acquisition and unadjusted age/sex/pose limitations. BY addresses dependence among valid p-values; it does not fix model assumptions or confounding. Recorded categories are not genetically verified ancestry. Missing annotations do not prove novelty, absence of a related clinical descriptor, or diagnostic usefulness. No manuscript or public release is made.

# Three sensitivities for the race-stratified PRISM redesign

This replaces the old pooled-reference sensitivity entrypoint with three experiments using the
completed November 15 primary baseline. `SENSITIVITY_ANALYSIS.json` selects this workflow;
`ACTIVE_ANALYSIS.json`, frozen data, historical numerical outputs and manuscript artifacts stay unchanged.
The original untracked November 10 extension is unrelated and preserved.

## Designs fixed before outcomes

The loader verifies the new primary manifest, source/code hashes, original reference sample,
eligible syndrome/category cells and full testing family. It reuses the exact primary source
loader and recorded-category crosswalk, including Middle Eastern/West Asian within broad Asian.
Each eligible cell has >=5 distinct GMDB patients with one accepted photograph, conditional age
3 to <10 and exactly 120 finite pinned FaceKit measurements. Unsupported source categories are excluded.
No old pooled, synthetic PDIDB or newly extracted measurements can substitute for this baseline.

**A — same-source placebo:** retain the original case counts and reference allocations. Within each
race, first attempt one shared FairFace reference disjoint from all pseudo-cases; when impossible,
try per-pair disjoint draws with cross-pair reuse explicitly allowed. If a contrast remains impossible,
retain all its hypotheses as p=1 in the full family. Never replace observations or shrink groups.
The fixed European allocation exhausts its FairFace pool, making its placebo contrasts unsupported.
This calibrates a same-source workflow and does not supply a p-value for GMDB/FairFace effects.

**B — GMDB label permutations:** compare each eligible syndrome with every other available GMDB
syndrome patient in that race, including patients from statistically ineligible cells. Shuffle one
multiclass label assignment across whole patients per race/repetition, shared across features and
syndromes; preserve every label count. Use |U−n_case*n_other/2| and empirical (1+exceedances)/(B+1),
including ties. An independent blocked family runs when any cell has age |SMD|>=.5 or a recorded-sex
fraction gap>=.20. Blocks are sex (including Unknown) × conditional age [3,5), [5,7), [7,10).
Only blocks with multiple labels shuffle. The blocked statistic is centered at its exact conditional
expected U, calculated from each label's block counts and each feature's block mean ranks.

The permutation schedule is fixed at 1,000 exploration, 10,000 and then 100,000 final draws per enabled
null, with no outcome-dependent stopping. Default final B is 100,000 (cap 100,000). The p floor is
1/(B+1). Exact binomial Monte Carlo tail bounds with simultaneous Bonferroni coverage across the
planned family are propagated through BH to flag uncertain Monte Carlo FDR decisions; these are not
clinical effect intervals. Different nulls have separate correction families.

**C — reference stability:** hold GMDB fixed, and draw one new shared reference per race without
replacement in each repetition, retaining the primary X and quotas. Repeat the complete primary
test and BH family. Report recurrence, sign consistency and g selection quantiles. The descriptive
stability rule is >=80% q<.05 recurrence and >=95% sign agreement. A cross-design shortlist additionally
requires main q<.05, empirical q<.05 in every enabled GMDB null and concordant g direction. Sort by
worst empirical q, then main q and exact identifiers. Neither shortlist nor recurrence is replication.
The European sample is its entire pool, so its g cannot vary; global q can vary through other groups.

All families retain the same fixed planned denominator. Two-sided MW follows primary per-feature
auto semantics (exact for min n<=8 without ties; asymptotic with ties/continuity correction otherwise).
Hedges' g is syndrome/pseudo-case minus comparator, pooled sample SD, exact gamma J correction;
zero pooled variance gives flagged undefined g. Complete finite features are mandatory. Mean and
rank calculations are optimized without changing tests. Exact-U lookup tables are constructed with
public SciPy calls on artificial rank configurations, avoiding private SciPy APIs. Fixed-case/reference
rank contributions and full-population GMDB ranks are precomputed, then validated against direct tests.

## Run from repository root

Use the existing `.venv`; no dependency changes or installs are needed. Versions match
`revision_results/race_stratified_20261115/requirements-ml.txt`; this work uses NumPy, pandas,
SciPy, matplotlib and PyArrow, without fitting ML. Synthetic tests additionally use pytest/statsmodels.
Runtime is measured in execution_commands.json; reference checkpoints use about 160 MB for 1,000
draws, with additional modest rank caches and permutation-count checkpoints. No images are processed.

```sh
.venv/bin/python stats_testing/reference_sensitivity.py --dry-run
.venv/bin/python -m pytest research/race_stratified_sensitivity/tests -q
.venv/bin/python stats_testing/reference_sensitivity.py --analysis all --repetitions 20 \
  --permutation-repetitions 20 --output-dir local_review/race_stratified_sensitivity_smoke
.venv/bin/python stats_testing/reference_sensitivity.py --analysis placebo --repetitions 1000
.venv/bin/python stats_testing/reference_sensitivity.py --analysis permutation --permutation-repetitions 1000
.venv/bin/python stats_testing/reference_sensitivity.py --analysis permutation --permutation-repetitions 10000 --resume
.venv/bin/python stats_testing/reference_sensitivity.py --analysis permutation --permutation-repetitions 100000 --resume
.venv/bin/python stats_testing/reference_sensitivity.py --analysis reference-resampling --repetitions 1000
.venv/bin/python research/race_stratified_sensitivity/validate.py
```

The implementation can also be invoked directly with
`.venv/bin/python research/race_stratified_sensitivity/run.py`. CLI options are `--analysis
placebo|permutation|reference-resampling|all`, `--repetitions`, `--permutation-repetitions`, `--seed`,
`--output-dir`, `--dry-run`, and `--resume`. Output must be inside verified Git-ignored local_review.
Use a fresh directory after numerical code, source or design changes. Same-design extension preserves
earlier iterations. Each repetition uses `SeedSequence([seed*100+analysis_index, variant, repetition,
race_index])`, independently of scheduling, with default analysis seeds 2026111501/02/03. Atomic
checksummed chunks record contiguous completed ranges. Existing checkpoints require explicit resume.

Plans, source/exclusion/capacity audits and exact eligible cells are written before random experiments.
The run is BLOCKED if the new primary or private sources are missing or hashes/schema fail; it never
falls back to old tables. Short runs are explicitly exploratory, not claimed as full execution.
The final private report is `local_review/race_stratified_sensitivity/RESULTS_AND_LIMITATIONS.md`.
The requested aggregate CSVs, figures, data dictionary and manifests are beside it; all small-cell
candidates remain private pending disclosure review. No public-export or cloud-upload action is included.

Recorded categories are imperfect annotations, not verified ancestry. Source acquisition remains
confounded in the main comparison; placebo support is incomplete. GMDB labels may not be exchangeable
even within the recorded age/sex blocks, and FairFace paths/content do not identify unique people.
Overlapping reference draws do not establish independent validation. Clinical/HPO, ML and manuscript
work are outside this replacement task.

Method references: [SciPy Mann–Whitney](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html)
and [SciPy permutation conventions](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.permutation_test.html).

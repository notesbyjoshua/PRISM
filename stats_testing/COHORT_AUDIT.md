# Current 120-feature cohort audit

The active configuration is `ACTIVE_ANALYSIS.json`. Run
`.venv/bin/python stats_testing/cohort_audit.py --dry-run` to verify the exact
120-column schema, finite vectors, current quality flags, unique recorded GMDB
patient IDs, unique FairFace paths and actual eligible-syndrome counts.

Current: 1,845 selected GMDB patients, 3,409 reference paths, 48 main syndromes.
Reference identity uniqueness and clinical health are unverified. PDIDB is excluded.

The historical audit and source were frozen privately with the original analysis.
Historical reproduction uses the frozen source and matching archived inputs and
results. The active CLI rejects `--historical` to prevent mixing old data with
current result tables. Historical helper functions remain for frozen workflows.
See `revision_results/rebuild_120_20261028/README.md` for the full validation package.

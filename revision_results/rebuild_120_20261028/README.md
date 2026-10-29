# Current PRISM 120-feature rebuild

The validated new analysis is active. It contains 120 features, 1,845 selected
GMDB patients, 3,409 FairFace reference paths, 48 eligible syndromes, 5,760 main
hypotheses and 955 global-BH-significant associations. Synthetic PDIDB is kept
separate. See repository `ACTIVE_ANALYSIS.json` for the active configuration.

All protected inputs, row-level results, read-only historical archives and the
full PDF/report package are under ignored `local_only/run_20261028/`. Its README
documents provenance, commands, dry runs, checkpoint resumption and validation.
The historical 125-feature archive is preserved; active exports have been replaced.
The validation-stage report records the state before the user requested a Git
commit and push. Protected inputs, row-level measurements and private reports
remain local; the repository changes contain code and aggregate working exports.
No live Google Doc edit or deployment was performed.

Scripts here prepare authorized data, extract synthetic PDIDB, build the separate
source catalog, rerun statistics/extensions, independently validate results,
generate manuscript materials, promote once and produce the private report.
Use the private run README before executing them. The exact paper execution pin
and image-publication permissions remain distinct from technical validation.

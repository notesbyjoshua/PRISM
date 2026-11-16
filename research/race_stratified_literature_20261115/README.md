# Main-analysis clinical-source audit and cross-category direction report

This code reads the verified `PRISM_RACE_STRATIFIED_20261115` main-analysis output. It does not rerun statistics or alter the manuscript. All detailed outputs, including small aggregate cells, are restricted to git-ignored `local_review/race_stratified_literature_20261115/` through the existing private-output guard.

The literature inventory covers every original significant main-analysis result, not every nonsignificant planned test. It is exhaustive over result rows and the downloaded annotation corpus, with targeted narrative checks. It cannot establish absence from all published literature or a novel phenotype. Clinical mappings are explicit, conservative proxies; no exact FaceKit measurement has been validated against a clinical source here. HPO and Orphadata can share evidence.

The separate direction report compares every eligible same-syndrome, same-feature category pair, including nonsignificant rows. Its sign reversals are descriptive and are not interaction tests.

Run from the repository root:

```sh
.venv/bin/python research/race_stratified_literature_20261115/directions.py
.venv/bin/python research/race_stratified_literature_20261115/audit.py
.venv/bin/python research/race_stratified_literature_20261115/reports.py
.venv/bin/python -m pytest tests/test_race_stratified_literature_20261115.py
```

The audit requires the frozen source files under the output directory's `public_sources/`. `download_manifest.json` records source URLs, hashes and sizes. `audit.py` refuses changed source hashes. Sources are HPO ontology v2026-09-01, HPO annotations 2026-09-02, Orphadata nomenclature/phenotypes 2026-06-23. Current Orphadata URLs are mutable, so retain the frozen copies; a later release is a new review, not an exact rerun. The nomenclature XML is parsed directly; fuzzy candidate-name output is not used.

`curation.py` records reviewed diagnosis identities and 120-feature mappings. `narrative.py` records targeted source passages, access limits, paraphrased notes and the positive traits actually added. Some NCBI direct requests returned reCAPTCHA; indexed clinical passages were used where available. A retained PMID from HPO does not assert full-text review. Never treat failed retrieval as absent clinical evidence.

`reports.py` writes the clinical comparison report, all 982 rows, per-classification lists, targeted source registry and checksummed manifest. `directions.py` writes all cross-category comparisons, every reversal, and significant subsets. No raw patient records or photographs are exported. Detailed reports remain local; this task does not authorize publishing them.

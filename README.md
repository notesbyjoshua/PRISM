# PRISM

**Phenotype Representation for Interpretable Syndrome Mapping**

PRISM is a local research application for exploring facial phenotype measurements
across rare-disease cohorts. It combines FaceKit feature extraction, a
pandas/SciPy statistical pipeline, a FastAPI and SQLite backend, and a
Next.js/Recharts dashboard.

> PRISM is an exploratory research tool. Its results are not medical advice,
> a diagnosis, or a substitute for clinical evaluation.

## Contents

- [Architecture](#architecture)
- [Requirements](#requirements)
- [First-time setup](#first-time-setup)
- [Preparing input data](#preparing-input-data)
- [Running the statistical analysis](#running-the-statistical-analysis)
- [Statistical methods](#statistical-methods)
- [Analysis outputs](#analysis-outputs)
- [Running the website](#running-the-website)
- [Real data versus demo data](#real-data-versus-demo-data)
- [Running the frontend and backend separately](#running-the-frontend-and-backend-separately)
- [API](#api)
- [Configuration](#configuration)
- [Project structure](#project-structure)
- [Troubleshooting](#troubleshooting)
- [Database roadmap](#database-roadmap)
- [Acknowledgements](#acknowledgements)

## Architecture

```text
Face images / landmark data
          │
          ▼
      FaceKit
          │ phenotype CSVs
          ▼
  pandas + SciPy analysis
          │ stats_results/*.csv
          ▼
  FastAPI ─── SQLite
          │ JSON API
          ▼
 Next.js + React + Recharts
```

The application currently runs entirely on the local computer. No patient
images or analysis data are uploaded by PRISM.

## Requirements

- macOS, Linux, or another Unix-like development environment
- Python 3.10 or newer
- Node.js 20 or newer
- npm
- Git

The current project has been tested with Python 3.11 and Node.js 24.

## First-time setup

Run these commands from the PRISM directory:

```bash
cd /Users/joshua/Documents/PRISM

git submodule update --init --recursive

python3 -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip setuptools wheel
python -m pip install -r backend/requirements.txt
python -m pip install -e './FaceKit[all]'

cd frontend
npm install
cd ..
```

The Python virtual environment lives in `.venv/`, and frontend packages live in
`frontend/node_modules/`. Both directories are ignored by Git.

### Using the virtual environment

Activate it in each new terminal before running Python commands:

```bash
source .venv/bin/activate
```

Exit it with:

```bash
deactivate
```

Activation is optional if the environment's Python is called directly:

```bash
.venv/bin/python stats_testing/statistical_analysis.py
```

For an editor's Run button, configure its Python interpreter as:

```text
/Users/joshua/Documents/PRISM/.venv/bin/python
```

## Preparing input data

The analysis supports either two separate CSVs or one combined CSV.

### Separate healthy and disease CSVs

Both files should have:

- A `frontal_ok` column. Only rows equal to `True` are analyzed.
- The same numeric FaceKit phenotype columns.
- At least the configured minimum number of valid observations per group.
- An optional `image_id` column.

The disease CSV may have a `disease` column containing individual cohort names.
If it does not, every row is assigned to one generic `disease` group. Any
existing disease label in the healthy CSV is replaced with `healthy`.

### Combined CSV

A combined file must contain a `disease` column, and healthy rows must use the
exact label `healthy`.

### Producing phenotype CSVs with FaceKit

FaceKit is installed as the local `FaceKit/` submodule. See its available
commands with:

```bash
source .venv/bin/activate
facekit --help
facekit extract-features --help
```

See [FaceKit/README.md](FaceKit/README.md) for its full extraction workflow.

## Running the statistical analysis

The simplest workflow is to edit the two constants near the top of
`stats_testing/statistical_analysis.py`:

```python
DISEASE_FILE = "data/phenotypes_all_disease.csv"
HEALTHY_FILE = "data/phenotypes_all_healthy.csv"
```

Paths may be absolute or relative to the PRISM directory. Then run:

```bash
cd /Users/joshua/Documents/PRISM
source .venv/bin/activate
python stats_testing/statistical_analysis.py
```

This also works with an editor's Run button when the editor uses the `.venv`
interpreter.

### Command-line input overrides

Command-line arguments override the constants in the file:

```bash
python stats_testing/statistical_analysis.py \
  --healthy-csv data/phenotypes_all_healthy.csv \
  --disease-csv data/phenotypes_all_disease.csv
```

For a combined file:

```bash
python stats_testing/statistical_analysis.py \
  --input-csv data/phenotypes_all.csv
```

### Adjusting analysis settings

```bash
python stats_testing/statistical_analysis.py \
  --healthy-csv data/phenotypes_all_healthy.csv \
  --disease-csv data/phenotypes_all_disease.csv \
  --bootstrap-iterations 1000 \
  --rank-bootstrap-iterations 500 \
  --max-missing-fraction 0.20 \
  --min-group-n 15 \
  --effect-threshold 0.50 \
  --kw-effect-threshold 0.06 \
  --mad-threshold 5 \
  --top-k 10
```

Defaults:

| Setting | Default | Meaning |
| --- | ---: | --- |
| `bootstrap-iterations` | 500 | Resamples used for effect-size confidence intervals |
| `rank-bootstrap-iterations` | 200 | Resamples used to estimate ranking stability |
| `bootstrap-seed` | 2026 | Random seed for reproducible bootstraps |
| `max-missing-fraction` | 0.30 | Maximum allowed missing fraction in either group |
| `min-group-n` | 10 | Minimum valid observations per group |
| `effect-threshold` | 0.50 | Minimum absolute Hedges' g for high priority |
| `kw-effect-threshold` | 0.06 | Minimum Kruskal-Wallis epsilon-squared |
| `mad-threshold` | 5 | Scaled MAD cutoff for the robust analysis |
| `top-k` | 10 | Number of positions used for rank stability |

Larger bootstrap settings improve precision but increase runtime substantially,
especially when the disease CSV contains many cohorts.

## Statistical methods

The script performs three main analyses:

1. All disease rows combined versus healthy rows.
2. Every eligible disease cohort versus healthy rows.
3. Kruskal-Wallis comparison across all sufficiently large groups.

The result tables include:

- Group counts, means, medians, and missingness.
- Mann-Whitney U tests and Benjamini-Hochberg FDR correction.
- Welch's t-test for the combined disease comparison.
- Cohen's d and small-sample-corrected Hedges' g.
- Cliff's delta and rank-biserial correlation.
- Bootstrap 95% confidence intervals for Hedges' g and Cliff's delta.
- Directional and direction-independent univariate ROC AUC.
- A robust reanalysis after removing values beyond the configured scaled-MAD
  threshold independently within each group.
- Bootstrap ranking stability: the fraction of bootstrap runs in which a
  feature appears in the top K by absolute Hedges' g.

### Phenotype Difference Score (PDF)

Each two-group result also receives a Phenotype Difference Score from 0 to 100:

```text
PDF = 100 × (0.45E + 0.20S + 0.25B + 0.10N)
```

Its components are:

```text
E = min(|Hedges' g| / 1.5, 1)
S = min(-log10(FDR q) / 5, 1)
B = fraction of effect-size bootstrap samples with the observed direction
N = min(min(n_disease, n_healthy) / 50, 1)
```

Effect magnitude receives the greatest weight. The output includes the four
individual `pdf_*_component` columns, the bootstrap direction-stability value,
and the final `phenotype_difference_score_pdf`, making the calculation fully
auditable. PDF is not calculated when bootstrapping is disabled because its
direction-stability component is unavailable.

A two-group result is marked `high_priority` when:

```text
global FDR q < 0.05 and |Hedges' g| >= 0.5
```

The Kruskal-Wallis high-priority rule is:

```text
FDR q < 0.05 and epsilon-squared >= 0.06
```

Feature specificity is calculated as:

```text
1 - (diseases with a high-priority effect / diseases tested)
```

The specificity-weighted effect score multiplies that specificity by the median
absolute Hedges' g among qualifying diseases. This prevents a feature with no
supported effects from appearing valuable merely because its raw specificity is
1.

## Analysis outputs

The analysis creates `stats_results/` automatically and writes:

| File | Contents |
| --- | --- |
| `stats_healthy_vs_disease.csv` | Complete combined disease-versus-healthy results |
| `stats_each_disease_vs_healthy.csv` | Complete disease-specific results |
| `stats_kruskal_all_groups.csv` | Complete multi-group results |
| `significant_healthy_vs_disease.csv` | Combined results with FDR q below 0.05 |
| `significant_each_disease_vs_healthy.csv` | Disease-specific results with global FDR q below 0.05 |
| `significant_kruskal_all_groups.csv` | Multi-group results with FDR q below 0.05 |
| `high_priority_healthy_vs_disease.csv` | Significant combined results meeting the effect threshold |
| `high_priority_each_disease_vs_healthy.csv` | Significant disease-specific results meeting the effect threshold |
| `high_priority_kruskal_all_groups.csv` | Significant multi-group results meeting the effect threshold |
| `feature_reliability.csv` | Valid counts, missing fractions, and eligibility by feature and group |
| `feature_rank_stability.csv` | Top-K bootstrap frequencies |
| `feature_specificity.csv` | Disease counts, specificity, and specificity-weighted effects |
| `analysis_configuration.csv` | Exact settings used for the run |

## Running the website

After generating statistics, start both services with:

```bash
cd /Users/joshua/Documents/PRISM
./scripts/dev.sh
```

Open:

- Dashboard: <http://localhost:3000>
- FastAPI documentation: <http://localhost:8000/docs>
- API health check: <http://localhost:8000/api/health>
- Data-source summary: <http://localhost:8000/api/summary>

Stop both services with `Control-C` in the terminal running `dev.sh`.

At startup, the backend imports available `stats_results` CSVs into the ignored
SQLite file `backend/data/prism.db`.

## Real data versus demo data

If no analysis CSVs are available, the backend inserts a small illustrative
dataset so the interface can be previewed. **Those records are fake and must not
be interpreted as scientific results.** The dashboard displays a yellow demo
banner in this state.

Check the current source at <http://localhost:8000/api/summary>:

```json
{
  "data_source": "demo",
  "is_demo": true
}
```

After real results are imported, those fields become:

```json
{
  "data_source": "analysis",
  "is_demo": false
}
```

To replace demo data:

1. Configure the healthy and disease files.
2. Run `stats_testing/statistical_analysis.py`.
3. Confirm CSVs exist in `stats_results/`.
4. Restart `./scripts/dev.sh`.

The backend clears demo records and imports the analysis automatically. While
the backend is already running, results can instead be reloaded with:

```bash
curl -X POST http://localhost:8000/api/admin/reload
```

## Running the frontend and backend separately

Use two terminals when separate logs or debugging are helpful.

### Terminal 1: backend

```bash
cd /Users/joshua/Documents/PRISM
source .venv/bin/activate
cd backend
uvicorn app.main:app --reload --port 8000
```

### Terminal 2: frontend

```bash
cd /Users/joshua/Documents/PRISM/frontend
npm run dev
```

The frontend expects the backend at `http://localhost:8000` by default.

## API

FastAPI provides interactive OpenAPI documentation at
<http://localhost:8000/docs>.

Main endpoints:

| Method | Endpoint | Purpose |
| --- | --- | --- |
| `GET` | `/api/health` | Server health check |
| `GET` | `/api/summary` | Counts, average AUROC, and data-source status |
| `GET` | `/api/diseases` | Available disease labels |
| `GET` | `/api/results` | Filtered and sorted statistical results |
| `GET` | `/api/specificity` | Ranked feature-specificity results |
| `GET` | `/api/features/{feature}` | All results for one feature |
| `POST` | `/api/admin/reload` | Re-import CSVs from `stats_results/` |

Example queries:

```bash
curl 'http://localhost:8000/api/results?high_priority=true&limit=20'

curl --get \
  --data-urlencode 'disease=Noonan syndrome' \
  --data-urlencode 'sort_by=roc_auc' \
  --data-urlencode 'descending=true' \
  'http://localhost:8000/api/results'
```

## Configuration

### Frontend API URL

For a backend running somewhere other than port 8000:

```bash
cd frontend
cp .env.example .env.local
```

Then edit `.env.local`:

```text
NEXT_PUBLIC_API_URL=http://localhost:8000
```

Restart the frontend after changing environment variables.

### CORS

The local backend currently permits frontend requests from:

- `http://localhost:3000`
- `http://127.0.0.1:3000`

Update `allow_origins` in `backend/app/main.py` before hosting the frontend at a
different origin.

## Project structure

```text
PRISM/
├── FaceKit/                         # Git submodule and feature extractor
├── backend/
│   ├── app/
│   │   ├── database.py              # SQLite connection and schema
│   │   ├── import_results.py        # CSV normalization/import and demo seed
│   │   └── main.py                  # FastAPI routes
│   ├── data/prism.db                # Generated local DB; ignored by Git
│   └── requirements.txt
├── data/                             # Input phenotype/metadata CSVs
├── frontend/
│   ├── app/                          # Next.js app and global styles
│   ├── components/                   # Dashboard, cards, and charts
│   ├── lib/api.ts                    # Typed API client
│   └── package.json
├── scripts/dev.sh                    # Starts frontend and backend together
├── stats_results/                    # Generated analysis CSVs
└── stats_testing/
    └── statistical_analysis.py       # Statistical pipeline
```

## Troubleshooting

### `No input files configured`

Set both `DISEASE_FILE` and `HEALTHY_FILE` near the top of
`stats_testing/statistical_analysis.py`, or provide both command-line flags.

### `The configured ... CSV does not exist`

Check the spelling and location. Relative paths are resolved from the PRISM
project directory, not from `stats_testing/`.

### `ModuleNotFoundError`

Activate the environment and reinstall Python requirements:

```bash
source .venv/bin/activate
python -m pip install -r backend/requirements.txt
python -m pip install -e './FaceKit[all]'
```

### `npm` or Next.js package errors

```bash
cd frontend
npm install
```

If the shared npm cache has a permissions problem, use a temporary cache:

```bash
npm install --cache /tmp/prism-npm-cache
```

### Ports 3000 or 8000 are already in use

Stop the existing development process with `Control-C`, or start services on
different ports. If the backend port changes, update `frontend/.env.local`.

### The dashboard still shows demo data

Confirm that `stats_results/stats_each_disease_vs_healthy.csv` exists, restart
the backend, and inspect <http://localhost:8000/api/summary>. You can also call
the reload endpoint described above.

### The analysis is slow

Bootstrap work scales with the number of cohorts, features, and iterations. For
a quick validation run, lower the iteration counts:

```bash
python stats_testing/statistical_analysis.py \
  --bootstrap-iterations 50 \
  --rank-bootstrap-iterations 25
```

Use the larger defaults or higher values for final research results.

## Database roadmap

SQLite is appropriate for the current local, single-user application. The API
and CSV import layer keep storage concerns isolated so the project can later move
to PostgreSQL or Supabase. That migration should add schema migrations,
environment-based credentials, authentication, row-level security, and a clear
policy for sensitive facial or medical data.

## Future direction

- Add feature and disease detail pages.
- Add downloadable filtered result tables.
- Add analysis-run versioning and provenance.
- Add authenticated multi-user support with PostgreSQL/Supabase.
- Consider an opt-in FaceKit image workflow only after privacy, consent,
  retention, and clinical-safety requirements are defined.
- Train a neural network using PDS score/implement on new input images? (Each disease might have a PDS 
  vector signature (signed))
- add an example photo fro each disease as well as where the top features are on that photo

## Acknowledgements

PRISM builds upon [FaceKit](https://github.com/HongzhuoChen/FaceKit), an
open-source toolkit for rare-disease facial phenotype analysis developed by
Hongzhuo Chen at the University of Pennsylvania.

I gratefully acknowledge Hongzhuo Chen and Dr. Kai Wang of the
Children's Hospital of Philadelphia for their guidance and support throughout
the development of this project.

I acknowledge the GestaltMatcher Database (GMDB) and FairFace as the image data
resources used in this work.

I extend my sincere appreciation to the individuals, patients, and families
whose contributions make this research possible. Use of these resources remains
subject to their respective licenses, consent terms, data-use agreements, and
citation requirements.

I also thank the maintainers and contributors of the open-source software on
which PRISM depends, including FaceKit, Python, pandas, NumPy, SciPy,
statsmodels, FastAPI, SQLite, React, Next.js, and Recharts.

# PRISM (Phenotype Representation for Interpretable Syndrome Mapping)

PRISM combines FaceKit-derived facial measurements with interpretable statistical
analysis and a local research dashboard.

## Local web app

The application has three layers:

- `stats_testing/`: pandas/SciPy statistical analysis.
- `backend/`: FastAPI and SQLite API, with a CSV import layer designed to be
  replaceable by PostgreSQL or Supabase later.
- `frontend/`: Next.js/React dashboard with Recharts visualizations.

### First-time setup

```bash
source .venv/bin/activate
pip install -r backend/requirements.txt
cd frontend && npm install && cd ..
```

### Run locally

```bash
chmod +x scripts/dev.sh
./scripts/dev.sh
```

Then open <http://localhost:3000>. FastAPI documentation is available at
<http://localhost:8000/docs>.

When statistical result CSVs exist in `stats_results/`, the API imports them
into SQLite automatically. Otherwise, the dashboard displays an explicitly
labeled illustrative dataset so the interface can still be previewed.

## Future direction

- Allow users to process a face image through FaceKit and explore candidate
  phenotype matches. This will require careful privacy, consent, and clinical
  safety design before implementation.

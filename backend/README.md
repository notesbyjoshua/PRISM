# PRISM API

The API imports analysis CSVs from `../stats_results` into SQLite at startup.
When no result CSVs exist, it serves a small, clearly marked illustrative dataset.

```bash
source ../.venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --reload --port 8000
```

OpenAPI documentation is available at <http://localhost:8000/docs>.


#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
export PYTHONDONTWRITEBYTECODE=1

cd "$PROJECT_DIR/backend"
"$PROJECT_DIR/.venv/bin/uvicorn" app.main:app --reload --port 8000 &
BACKEND_PID=$!

cd "$PROJECT_DIR/frontend"
npm run dev &
FRONTEND_PID=$!

cleanup() {
  kill "$BACKEND_PID" "$FRONTEND_PID" 2>/dev/null || true
}
trap cleanup EXIT INT TERM

wait

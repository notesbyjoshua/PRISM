#!/usr/bin/env python3
import io
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from pandas.api.types import is_float_dtype

ROOT = Path(__file__).resolve().parents[2]
FIGURE_ROOTS = ("publication/figures/", "publication/feature_correlation/", "manuscript/figures/")


def git(*args): return subprocess.check_output(["git", *args], cwd=ROOT)


def changed_files(): return git("diff", "--name-only", "--diff-filter=ACMRTUXB").decode().splitlines()


def committed_csv(path): return pd.read_csv(io.BytesIO(git("show", f"HEAD:{path}")), low_memory=False)


def sort_frame(frame):
    frame = frame.copy()
    keys = []
    for column in frame:
        series = frame[column]
        key = f"__key_{len(keys)}"
        if is_float_dtype(series): frame[key] = series.round(10).fillna(np.inf)
        else: frame[key] = series.astype("string").fillna("<NA>")
        keys.append(key)
    return frame.sort_values(keys, kind="mergesort").drop(columns=keys).reset_index(drop=True)


def compare_csv(path):
    expected, actual = committed_csv(path), pd.read_csv(ROOT / path, low_memory=False)
    if list(expected.columns) != list(actual.columns): raise AssertionError(f"{path}: columns changed")
    if len(expected) != len(actual): raise AssertionError(f"{path}: row count changed from {len(expected)} to {len(actual)}")
    expected, actual = sort_frame(expected), sort_frame(actual)
    for column in expected:
        left, right = expected[column], actual[column]
        if is_float_dtype(left) and is_float_dtype(right):
            if not np.allclose(left, right, rtol=1e-10, atol=1e-12, equal_nan=True):
                delta = np.nanmax(np.abs(left.to_numpy() - right.to_numpy()))
                raise AssertionError(f"{path}: numeric column {column!r} differs (maximum absolute delta {delta})")
        elif not left.astype("string").fillna("<NA>").equals(right.astype("string").fillna("<NA>")):
            raise AssertionError(f"{path}: nonnumeric column {column!r} differs")


def main():
    changed = changed_files()
    csvs, figures, unexpected = [], [], []
    for path in changed:
        if path.endswith(".csv"): csvs.append(path)
        elif path.endswith((".pdf", ".png")) and path.startswith(FIGURE_ROOTS): figures.append(path)
        else: unexpected.append(path)
    if unexpected: raise AssertionError("Unexpected regenerated changes:\n" + "\n".join(unexpected))
    for path in csvs: compare_csv(path)
    print(f"Reproducibility check passed: {len(csvs)} CSVs are structurally and numerically equivalent; "
          f"{len(figures)} platform-rendered figure files were regenerated successfully.")


if __name__ == "__main__": main()

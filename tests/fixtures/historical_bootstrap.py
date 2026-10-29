import numpy as np
import pandas as pd

SOURCE_COMMIT = "befee5c"
SOURCE_SHA256 = "3934d76a65d3ca8fa60ff83ad3775141bae0ae58f8bab49195e1ed96a74bc14c"


def cohens_d(x, y):
    """
    Cohen's d for two independent groups.
    Positive value means x has a larger mean than y.
    """
    x = np.asarray(x.dropna(), dtype=float)
    y = np.asarray(y.dropna(), dtype=float)
    nx = len(x)
    ny = len(y)
    if nx < 2 or ny < 2:
        return np.nan
    sx = np.var(x, ddof=1)
    sy = np.var(y, ddof=1)
    pooled_sd = np.sqrt(((nx - 1) * sx + (ny - 1) * sy) / (nx + ny - 2))
    if pooled_sd == 0:
        return np.nan
    return (np.mean(x) - np.mean(y)) / pooled_sd


def hedges_g(x, y):
    """
    Small-sample corrected Cohen's d.
    Usually better than plain Cohen's d for small disease groups.
    """
    d = cohens_d(x, y)
    x = np.asarray(x.dropna(), dtype=float)
    y = np.asarray(y.dropna(), dtype=float)
    nx = len(x)
    ny = len(y)
    if np.isnan(d) or nx + ny <= 3:
        return np.nan
    correction = 1 - (3 / (4 * (nx + ny) - 9))
    return d * correction


def cliffs_delta_fast(x, y):
    """Cliff's delta without constructing the full pairwise comparison matrix."""
    x = np.asarray(x, dtype=float)
    y = np.sort(np.asarray(y, dtype=float))
    if len(x) == 0 or len(y) == 0:
        return np.nan
    lesser = np.searchsorted(y, x, side="left").sum()
    greater = (len(y) - np.searchsorted(y, x, side="right")).sum()
    return (lesser - greater) / (len(x) * len(y))


def bootstrap_effect_ci(x, y, iterations, rng):
    """Effect-size CIs and the fraction of Hedges' g samples with the same direction."""
    x = np.asarray(pd.Series(x).dropna(), dtype=float)
    y = np.asarray(pd.Series(y).dropna(), dtype=float)
    if iterations == 0 or len(x) < 2 or len(y) < 2:
        return np.nan, np.nan, np.nan, np.nan, np.nan
    g_samples = np.empty(iterations)
    delta_samples = np.empty(iterations)
    for i in range(iterations):
        xb = x[rng.integers(0, len(x), len(x))]
        yb = y[rng.integers(0, len(y), len(y))]
        g_samples[i] = hedges_g(pd.Series(xb), pd.Series(yb))
        delta_samples[i] = cliffs_delta_fast(xb, yb)
    g_valid = g_samples[np.isfinite(g_samples)]
    delta_valid = delta_samples[np.isfinite(delta_samples)]
    g_ci = np.quantile(g_valid, [0.025, 0.975]) if len(g_valid) else [np.nan, np.nan]
    d_ci = np.quantile(delta_valid, [0.025, 0.975]) if len(delta_valid) else [np.nan, np.nan]
    observed_g = hedges_g(pd.Series(x), pd.Series(y))
    if not len(g_valid) or not np.isfinite(observed_g) or observed_g == 0:
        direction_stability = np.nan
    elif observed_g > 0:
        direction_stability = np.mean(g_valid > 0)
    else:
        direction_stability = np.mean(g_valid < 0)
    return g_ci[0], g_ci[1], d_ci[0], d_ci[1], direction_stability

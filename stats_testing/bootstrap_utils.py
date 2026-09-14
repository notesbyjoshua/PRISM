import numpy as np


def bootstrap_g(x, y, iterations, rng, batch_size=100):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    x, y = x[np.isfinite(x)], y[np.isfinite(y)]
    if len(x) < 2 or len(y) < 2 or iterations <= 0:
        return np.array([])
    df = len(x) + len(y) - 2
    correction = 1 - 3 / (4 * df - 1)
    samples = []
    for start in range(0, iterations, batch_size):
        size = min(batch_size, iterations - start)
        xb = x[rng.integers(len(x), size=(size, len(x)))]
        yb = y[rng.integers(len(y), size=(size, len(y)))]
        pooled = ((len(x) - 1) * xb.var(axis=1, ddof=1) + (len(y) - 1) * yb.var(axis=1, ddof=1)) / df
        with np.errstate(divide='ignore', invalid='ignore'):
            samples.append(correction * (xb.mean(axis=1) - yb.mean(axis=1)) / np.sqrt(pooled))
    values = np.concatenate(samples)
    return values[np.isfinite(values)]

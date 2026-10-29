import numpy as np


def bootstrap_pair_effects(x, y, iterations, rng, batch_size=100):
    x, y = np.asarray(x, float), np.asarray(y, float)
    nx, ny = len(x), len(y)
    correction = 1-3/(4*(nx+ny)-9)
    order = np.argsort(y)
    left, right = np.searchsorted(y[order], x, 'left'), np.searchsorted(y[order], x, 'right')
    effects, deltas = [], []
    for start in range(0, iterations, batch_size):
        size = min(batch_size, iterations-start)
        ix, iy = np.empty((size, nx), dtype=int), np.empty((size, ny), dtype=int)
        for index in range(size):
            ix[index] = rng.integers(0, nx, nx)
            iy[index] = rng.integers(0, ny, ny)
        xb, yb = x[ix], y[iy]
        pooled = ((nx-1)*xb.var(axis=1, ddof=1)+(ny-1)*yb.var(axis=1, ddof=1))/(nx+ny-2)
        with np.errstate(divide='ignore', invalid='ignore'):
            effects.append(correction*(xb.mean(axis=1)-yb.mean(axis=1))/np.sqrt(pooled))
        offsets = np.arange(size)[:, None]*ny
        weights = np.bincount((iy+offsets).ravel(), minlength=size*ny).reshape(size, ny)
        cumulative = np.pad(np.cumsum(weights[:, order], axis=1), ((0, 0), (1, 0)))
        differences = cumulative[:, left]-(ny-cumulative[:, right])
        deltas.append(np.take_along_axis(differences, ix, axis=1).sum(axis=1)/(nx*ny))
    return np.concatenate(effects), np.concatenate(deltas)


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

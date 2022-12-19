import numpy as np
import math

def LpNorm_1D(x, f, p = 2.0):
    # Check numpy arrays
    if isinstance(x, np.ndarray) and isinstance(f, np.ndarray):
        pass
    else:
        raise Exception("x and f must be ndarrays!")
    # Check for size consistency
    if x.size == f.size:
        numPoints = x.size
        pass
    else:
        raise Exception("x and f must have the same dimension!")
    # Check x to be sorted
    if np.all(np.diff(x) >= 0.0):
        pass
    else:
        raise Exception("x must be sorted!")
    # Copute integral: int_\Omega (f)^p dx; with \Omega = [min(x),
    # max(x)]
    norm = 0.0
    for k1 in range(numPoints - 1):
        x1 = x[k1]
        x2 = x[k1 + 1]
        f1 = f[k1]
        f2 = f[k1 + 1]
        norm += (((f2 + f1)/2.0)**p)*(x2 - x1)
    # Compute norm
    norm = norm**(1.0/p)
    return norm

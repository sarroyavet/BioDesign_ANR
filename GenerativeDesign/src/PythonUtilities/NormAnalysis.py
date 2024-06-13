# Libraries {{{
import numpy as np
from scipy import integrate
import unittest
# }}}

# Functions {{{
# LpNorm_1D {{{
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
# }}}

# LpNorm_ndim {{{
def LpNorm_ndim(f, x, p = 2):
    # Documentation {{{
    """
    Function to compute the Lp norm of a n--dimensional function.
    f must be a numpy array of size (m1 x m2 x ... x mn).
    x must be a list of numpy 1D arrays of size m1, m2, ... mn.
    """
    # }}}
    # Test f and x
    n = f.ndim
    if len(x) != n:
        raise("Error: see doc." + LpNorm_ndim.__doc__)
    ms = f.shape
    for k1 in range(n):
        if ms[k1] != x[k1].size:
            raise("Error: see doc." + LpNorm_ndim.__doc__)
    # Compute the integrals
    lpNorm = (RecursiveIntegral(f**p, x))**(1.0/p)
    return lpNorm
# }}}

# Recursive integral {{{
def RecursiveIntegral(f, x):
    integral = 0.0
    n = f.ndim
    ms = f.shape
    if n == 1:
        xi = x[0]
        yi = f
    else:
        xi = x[0]
        yi = []
        for k1 in range(ms[0]):
            ff = f[k1]
            xx = x[1:]
            yi.append(RecursiveIntegral(ff, xx))
        yi = np.array(yi)
    return integrate.simps(yi, xi)
# }}}
# }}}

# Test {{{
class Test_NormAnalysis(unittest.TestCase):
    # test_ values {{{
    def test_values(self):
        def hertz(x):
            return (1.0 - x*x)**0.5
        def func1(x, y, z):
            return x + y + z
        def func2(x, y, z):
            return np.exp(x + y + z)
        def func3(x, y, z):
            return x - y*y - z*z*z
        m1 = 81
        m2 = 81
        m3 = 81
        x = np.linspace(-1.0, 1.0, m1)
        y = np.linspace(-3.0, 2.0, m2)
        z = np.linspace(-2.0, 3.0, m3)
        f1 = np.zeros((m1, m2, m3))
        f2 = np.zeros((m1, m2, m3))
        f3 = np.zeros((m1, m2, m3))
        for k1 in range(m1):
            for k2 in range(m2):
                for k3 in range(m3):
                    f1[k1, k2, k3] =  func1(x[k1], y[k2], z[k3])
                    f2[k1, k2, k3] =  func2(x[k1], y[k2], z[k3])
                    f3[k1, k2, k3] =  func3(x[k1], y[k2], z[k3])
        h = hertz(x)
        return (x, y, z, f1, f2, f3, h)
    # }}}
    # Test RecursiveIntegral {{{
    def test_RecursiveIntegral(self):
        (x, y, z, f1, f2, f3, h) = self.test_values()
        I1 = RecursiveIntegral(f1, [x, y, z])
        I2 = RecursiveIntegral(f2, [x, y, z])
        I3 = RecursiveIntegral(f3, [x, y, z])
        I4 = RecursiveIntegral(h, [x])
        self.assertAlmostEqual(I1, 0.0    , places = 2)
        self.assertAlmostEqual(I2, 344.146, places = 2)
        self.assertAlmostEqual(I3, -279.17, places = 2)
        self.assertAlmostEqual(I4, np.pi/2, places = 2)
        return
    # }}}
    # Test LpNorm_ndim {{{
    def test_LpNorm_ndim(self):
        (x, y, z, f1, f2, f3, h) = self.test_values()
        I1 = LpNorm_ndim(f1, [x, y, z])
        I2 = LpNorm_ndim(f2, [x, y, z], p = 5)
        I3 = LpNorm_ndim(f3, [x, y, z], p = 10)
        I4 = LpNorm_ndim(h, [x])
        self.assertAlmostEqual(I1, 15.0   , places = 2)
        self.assertAlmostEqual(I2, 153.596, places = 2)
        self.assertAlmostEqual(I3, 30.8181, places = 2)
        self.assertAlmostEqual(I4, 1.15470, places = 2)
        return
    # }}}
# }}}

# Run tests {{{
if __name__ == '__main__':
    unittest.main()
# }}}

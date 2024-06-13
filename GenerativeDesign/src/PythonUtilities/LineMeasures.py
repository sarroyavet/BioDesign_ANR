# Libraries {{{
import os
import sys
import numpy as np
from scipy import interpolate
from scipy import integrate
import pandas as pd
from statsmodels.nonparametric.smoothers_lowess import lowess

# In-house modules
sys.path.append(os.path.dirname(os.path.abspath(__file__))[:-4] + 'PythonUtilities/')
# }}}

# Class discrete line {{{
class DiscreteLine(object):
    # Description {{{
    """
    This class represents a finite element line where a set of variables
    are defined in a dictionary.
    """
    # }}}
    # Properties {{{
    @property
    def x_array(self):
        return self._x_array
    @property
    def xdom(self):
        return self._xdom
    @property
    def numPoints(self):
        return self._numPoints
    @property
    def values(self):
        return self._values
    # }}}
    # __init__ {{{
    def __init__(self, x_array, **kwargs):
        # Initialisation {{{
        numPoints = x_array.shape[0]
        xmin = x_array.min()
        #if abs(xmin) > 1.0e-6:
        #    raise ValueError("'xmin' must be zero, tolerance = 1.0e-6.")
        xmax = x_array.max()
        xdom = np.array([xmin, xmax])
        # }}}
        # Assign values {{{
        # Obligatory
        self._x_array = x_array
        self._xdom = xdom
        self._numPoints = numPoints
        # Optional
        self._values = kwargs.get("values", {})
        # Test values {{{
        if not isinstance(self.values, dict):
            raise TypeError("'values' must be a dictionary")
        for key, arg in self.values.items():
            if arg.shape[0] != numPoints:
                message  = "The argument '" + key + "' must be a numpy "
                message += "array with the same number of points of "
                message += "x_array."
                raise ValueError(message)
        # }}}
        # }}}
        return
    # }}}
    # add_value {{{
    def add_value(self, key, array):
        self._values[key] = array
        return
    # }}}
    # Compute curvature {{{
    def Curvature(self, ykey = "y", **kwargs):
        # Make spline
        x = self.x_array
        y = self.values[ykey]
        cs = interpolate.CubicSpline(x, y)
        # Compute spline derivatives
        yprim     = cs.derivative(nu = 1)
        yprimprim = cs.derivative(nu = 2)
        # Compute kappa
        kappa = np.zeros(self.numPoints)
        for k1 in range(self.numPoints):
            xi = x[k1]
            kappa[k1] = yprimprim(xi)/((1.0 + yprim(xi)**2.0)**1.5)
        return kappa
    # }}}
    # Compute derivative {{{
    def Derivative(self, ykey, order, **kwargs):
        # Make spline
        x = self.x_array
        y = self.values[ykey]
        cs = interpolate.CubicSpline(x, y)
        # Make spline derivatives as function
        yprim_func = cs.derivative(nu = order)
        # Compute derivative
        array = np.zeros(self.numPoints)
        for k1 in range(self.numPoints):
            xi = x[k1]
            array[k1] = yprim_func(xi)
        return array
    # }}}
    # Compute lp-norm of a field {{{
    def LpNorm(self, key, p = 2.0):
        # Get data
        x = self.x_array
        arg = self.values[key]
        norm = integrate.simps(arg**p, x)**(1.0/p)
        return norm
    # }}}
    # Map from values from discrete line {{{
    def MapValues(self, line, keys_in, keys_out = None):
        # Set up keys
        if keys_out == None:
            keys_out = keys_in
        numKeys = len(keys_in)
        # Get variables
        numPoints = self.numPoints
        x = self.x_array
        xline = line.x_array
        # Interpolate each value
        for k1 in range(numKeys):
            key_in  = keys_in[k1]
            key_out = keys_out[k1]
            # Make interpolation function
            vline = line.values[key_in]
            int_func = interpolate.interp1d(xline, vline,
                    kind = "cubic", fill_value = "extrapolate")
            # Initialisation of array
            arr = np.zeros(numPoints)
            for k2 in range(numPoints):
                arr[k2] = int_func(x[k2])
            # Add to values
            self.add_value(key_out, arr)
        return
    # }}}
# }}}

# Evolving discrete line from data frame {{{
def EvolvingDiscreteLine(df, tkey = "Time", xkey = "x_array",
        smoothing = None, sample_at = 1):
    # Get dependent variables {{{
    y_names = df.columns.values.tolist()
    y_names.remove(tkey)
    y_names.remove(xkey)
    # }}}
    # Get times {{{
    times = df[tkey].unique()
    numTimes = len(times)
    # }}}
    # Make a dictionary for the evolving discrete line {{{
    # Initialisation
    evolvingLine = dict()
    for k1 in range(numTimes):
        # Get line data frame at time k1 {{{
        df_k1 = df[df[tkey] == times[k1]]
        df_k1 = df_k1.drop(tkey, axis = 1)
        # }}}
        # Make line at time k1 a new (expanded) data frame {{{
        x_array = df_k1[xkey][k1][::sample_at]
        values = dict()
        for colName in y_names:
            if smoothing == None:
                values[colName] = df_k1[colName][k1][::sample_at]
            else:
                values[colName] = lowess(df_k1[colName][k1][::sample_at],
                        x_array, frac = smoothing)[:, 1]
        dl_k1 = DiscreteLine(x_array, values = values)
        # }}}
        # Add to evolving line
        evolvingLine[times[k1]] = dl_k1
    # }}}
    return evolvingLine
# }}}

# Test {{{
if __name__ == '__main__':
    # Define mesh
    x1 = np.linspace(0.0, 1.0, 401)
    y1 = x1**2.0
    line1 = DiscreteLine(x1)
    line1.add_value("yesc", y1)

    x2 = np.linspace(0.0, 1.0, 201)
    y2 = 1.0 - np.sqrt(1.0 - x2**2.0)
    line2 = DiscreteLine(x2)
    line2.add_value("ymai", y2)
    # Map slave and master into reference coordinates
    xfine = np.linspace(0.0, 1.0, 1001)
    fineLine = DiscreteLine(xfine)
    fineLine.MapValues(line1, ["yesc"])
    fineLine.MapValues(line2, ["ymai"])
    # Compute curvatures
    kesc = fineLine.Curvature(ykey = "yesc")
    kmai = fineLine.Curvature(ykey = "ymai")
    fineLine.add_value("kesc", kesc)
    fineLine.add_value("kmai", kmai)
    # Compute curvature diference
    kdif = kesc - kmai
    fineLine.add_value("kdif", kdif)
    # Compute measures
    print(fineLine.LpNorm("kesc", p = 2.0))
    print(fineLine.LpNorm("kmai", p = 2.0))
    print(fineLine.LpNorm("kdif", p = 2.0))

# }}}

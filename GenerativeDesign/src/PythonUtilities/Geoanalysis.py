# Libraries {{{
import os
import sys
import numpy as np

# }}}

# Classes {{{
# Straight line {{{
class StraightLine(object):
    """This class defines a straight line object.
    Initialisation with two arguments in the following possible
    combinations:
        1. (m, b):  slope (scalar) and y--intercept (scalar).
        2. (m, p1): slope (scalar) and one point (vector).
        3. (p1, p2): two points (vector)."""
    # Properties {{{
    @property
    def m(self):
        return self._m
    @property
    def b(self):
        return self._b
    # }}}

    # Initialiser {{{
    def __init__(self, **kwargs):
        # Check items {{{
        items = kwargs.items()
        if len(items) != 2:
            print("Error: StraightLine initialisation requires 2 arguments.")
            raise
        items = dict(items)
        if 'm' in items:
            if 'b' in items:
                case = 1
            else:
                if 'p1' in items:
                    case = 2
                else:
                    print("Error: wrong combination of arguments.")
        else:
            if (('p1' in items) and ('p2' in items)):
                case = 3
            else:
                print("Error: wrong combination of arguments.")
        # }}}
        # Define slope and y--intercept {{{
        if case == 1:
            m = items["m"]
            b = items["b"]
        elif case == 2:
            m = items["m"]
            p = items["p1"]
            px = p[0]
            py = p[1]
            b = p[1] - m*p[0]
        elif case == 3:
            p1 = items["p1"]
            p2 = items["p2"]
            m = (p2[1] - p1[1])/(p2[0] - p1[0])
            b = p1[1] - m*p1[0]
        # }}}
        self._m = m
        self._b = b
        return
    # }}}

    # Call {{{
    def __call__(self, x):
        return self.m*x + self.b
    # }}}

    # String form {{{
    def __str__(self):
        if self.b > 0.0:
            return "{:6.3f}x + {:6.3f}".format(self.m, abs(self.b))
        else:
            return "{:6.3f}x - {:6.3f}".format(self.m, abs(self.b))
    # }}}

    # Closest point to a reference point r {{{
    def ClosestPointToPoint(self, r):
        rx = r[0]
        ry = r[1]
        m = self._m
        b = self._b
        px = (rx + m*(ry - b))/(1.0 + m**2.0)
        py = self(px)
        return np.array([px, py])
    # }}}

    # Distance to a point r {{{
    def DistanceToPoint(self, r):
        p = self.ClosestPointToPoint(r)
        dx = r[0] - p[0]
        dy = r[1] - p[1]
        return np.sqrt(dx**2.0 + dy**2.0)
    # }}}

    # Closest point from a set of points {{{
    def ClosestPointFromSetOfPoints(self, points):
        numPoints = len(points)
        dists = np.zeros(numPoints)
        for k1 in range(numPoints):
            dists[k1] = self.DistanceToPoint(points[k1])
        minDistId = np.argmin(dists)
        minDist = dists[minDistId]
        minPoint = points[minDistId]
        return minPoint, minDist, minDistId
    # }}}

    # Intersection between lines {{{
    def Intersection(self, sl2):
        if abs(self.m - sl2.m) < 0.000001:
            raise("Error: the lines are parallel between each other.")
        x = (sl2.b - self.b)/(self.m - sl2.m)
        if self.m < sl2.m:
            y = self(x)
        else:
            y = sl2(x)
        return np.array([x, y])
    # }}}
# }}}

# Vertical line {{{
class VerticalLine(object):
    # Properties {{{
    @property
    def xref(self):
        return self._xref
    # }}}

    # Initialiser {{{
    def __init__(self, xref):
        self._xref = xref
        return
    # }}}

    # String form {{{
    def __str__(self):
        return "x = {:6.3f}".format(self.xref)
    # }}}

    # Closest point to a reference point r {{{
    def ClosestPointToPoint(self, r):
        px = self.xref
        py = r[1]
        return np.array([px, py])
    # }}}

    # Distance to a point r {{{
    def DistanceToPoint(self, r):
        dx = r[0] - self.xref
        return abs(dx)
    # }}}

    # Closest point from a set of points (*** repeated ***) {{{
    def ClosestPointFromSetOfPoints(self, points):
        numPoints = len(points)
        dists = np.zeros(numPoints)
        for k1 in range(numPoints):
            dists[k1] = self.DistanceToPoint(points[k1])
        minDistId = np.argmin(dists)
        minDist = dists[minDistId]
        minPoint = points[minDistId]
        return minPoint, minDist, minDistId
    # }}}

    # Intersection between straight lines {{{
    def Intersection(self, sl2):
        try:
            m = sl2.m
            b = sl2.b
        except AttributeError:
            raise AttributeError("sl2 must be a straight line.")
        x = self.xref
        y = sl2(x)
        return np.array([x, y])
    # }}}
# }}}
# }}}

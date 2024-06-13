# Libraries {{{
import os
import sys
import gmsh
import getopt
import json
import math
import numpy as np
from scipy.spatial import KDTree
from multiprocessing import Pool, cpu_count
from shapely.geometry import LineString
from shapely.ops import polygonize, unary_union
from scipy import integrate
from scipy.interpolate import interp1d

# In-house modules
sys.path.append(os.path.dirname(os.path.abspath(__file__))[:-4] + 'PythonUtilities/')
from meshMaker import AsterMesh
from Geoanalysis import StraightLine
# }}}

# Geometry class {{{
class Geometry(object):
    # Properties {{{
    @property
    def points(self):
        return self._points
    @property
    def lines(self):
        return self._lines
    @property
    def ndim(self):
        return self._ndim
    @property
    def npoints(self):
        return self._npoints
    @property
    def nlines(self):
        return self._nlines
    @property
    def name(self):
        return self._name
    @property
    def pointIds(self):
        return self._pointIds
    # }}}
    # __init__ {{{
    def __init__(self, name, points, lines, **kwargs):
        if isinstance(points, list):
            points = np.array(points)
        shape = points.shape
        self._npoints = shape[0]
        self._ndim    = shape[1]
        self._points = points*kwargs.get("scale", 1.0)
        self._nlines = len(lines)
        self._lines = lines
        self._name  = name
        self._pointIds = []
        return
    # }}}
    # Make geometry {{{
    def MakeGeometry(self, makeFace = True, clearAll = True):
        # GMSH Initialisation {{{
        if clearAll:
            gmsh.initialize()
            gmsh.clear()
        model = gmsh.model
        geo = model.occ
        mesh = model.mesh
        # }}}
        # Make points {{{
        for k1 in range(self.npoints):
            coord_x = self.points[k1, 0]
            coord_y = self.points[k1, 1]
            if self.ndim == 2:
                coord_z = 0.0
            else:
                coord_z = self.points[k1, 2]
            self._pointIds.append(geo.addPoint(coord_x, coord_y, coord_z))
        # }}}
        # Make lines {{{
        pointIds = self.pointIds
        for line in self.lines:
            ltype = line["type"]
            ordered_points = line["ordered_points"]
            list_points = [pointIds[k1] for k1 in ordered_points]
            if ltype == "straight":
                #straight lines {{{
                line["id"] = geo.addLine(*list_points)
                # }}}
            elif ltype == "spline":
                # Splines {{{
                line["id"] = geo.addSpline(*list_points)
                # }}}
            elif ltype == "bspline":
                # B-splines {{{
                line["id"] = geo.addBSpline(list_points)
                # }}}
            elif ltype == "bezier":
                # Bézier {{{
                line["id"] = geo.addBezier(list_points)
                # }}}
            elif ltype == "circlearc":
                # Circle arc {{{
                line["id"] = geo.addCircleArc(*list_points)
                # }}}
            elif ltype == "ellipsearc":
                # Ellipse arc {{{
                line["id"] = geo.addEllipseArc(*list_points)
                # }}}
            else:
                raise ValueError("Error: unavailable type of line.")
        # }}}
        # Make face {{{
        if makeFace:
            loop = geo.addCurveLoop([k1 + 1 for k1 in range(self.nlines)])
            face = geo.addPlaneSurface([loop])
        # }}}
        # Synchronize
        geo.synchronize()
        # Add Physical groups {{{
        if makeFace:
            for line in self.lines:
                physical_name = line["physical_name"]
                line_id = line["id"]
                self.PhysicalGroup(1, [line_id], physical_name)
            self.PhysicalGroup(2, [face], self.name)
        # }}}
        return
    # }}}
    # Physical group {{{
    def PhysicalGroup(self, dim, lst_entities, gr_name):
        # self variables {{{
        model = gmsh.model
        # }}}
        # Get a list of the existing physical groups
        physical_groups = model.getPhysicalGroups(dim)
        physical_dictionary = {}
        if len(physical_groups) > 0:
            for gr_dim, gr_tag in physical_groups:
                name = model.getPhysicalName(gr_dim, gr_tag)
                physical_dictionary[name] = (gr_dim, gr_tag)
        lst_tags = lst_entities
        if gr_name in physical_dictionary:
            gr_dim, gr_tag = physical_dictionary[gr_name]
            for tag in model.getEntitiesForPhysicalGroup(gr_dim, gr_tag):
                lst_tags.append(tag)
            model.removePhysicalGroups([(gr_dim, gr_tag)])
        gr_id = model.addPhysicalGroup(dim, lst_tags)
        model.setPhysicalName(dim, gr_id, gr_name)
        return
    # }}}
    # Popup {{{
    def Popup(self):
        gmsh.fltk.run()
        return
    # }}}
    # Make mesh {{{
    def MakeMesh(self, meshAlgo = 8, smoothItes = 10, recombine = True,
            **kwargs):
        # GMSH environment parameters {{{
        model = gmsh.model
        geo = model.occ
        mesh = model.mesh
        # }}}
        gmsh.option.setNumber("Mesh.Algorithm", meshAlgo)
        gmsh.option.setNumber("Mesh.Smoothing", smoothItes)
        mesh.generate(1)
        mesh.generate(2)
        if recombine:
            mesh.recombine()
        if "setOrder" in kwargs:
            mesh.setOrder(kwargs["setOrder"])
        mesh.removeDuplicateNodes()
        return
    # }}}
    # Export mesh {{{
    def Export(self, fmt = 'unv', folder = ''):
        if folder == '':
            folder = os.getcwd()
        pathname = folder + '/' + self.name + '.' + fmt
        gmsh.write(pathname)
        return
    # }}}
    # Cut line {{{
    def CutLine(self, lineName, xmin = -10.0e+9, xmax = 10.0e+9,
            ymin = -10.0e+9, ymax = 10.0e+9):
        # Get nodes and line {{{
        nodes = self.points
        lines = self.lines
        for l_i in lines:
            if l_i["physical_name"] == lineName:
                line = l_i
                break
        # }}}
        # Add node to points if it is between the limits {{{
        # Initialisation
        points = []
        ordered_points = []
        node_Ids = line["ordered_points"]
        id_i = 0
        prev_in = False
        prev_out = False
        for nodeId in node_Ids:
            add = False
            node = nodes[nodeId]
            xi = node[0]
            yi = node[1]
            # Test node coordinates
            if (xmin <= xi) and (xi <= xmax):
                if (ymin <= yi) and (yi <= ymax):
                    add = True
            if add:
                if prev_out:
                    # Interpolate point {{{
                    node_prev = nodes[nodeId - 1]
                    xi_prev = node_prev[0]
                    yi_prev = node_prev[1]
                    # Cut at x?
                    if xi_prev < xmin:
                        cut_x = xmin
                    elif xi_prev > xmax:
                        cut_x = xmax
                    else:
                        cut_x = xi_prev
                    # Cut at y?
                    if yi_prev < ymin:
                        cut_y = ymin
                    elif yi_prev > ymax:
                        cut_y = ymax
                    else:
                        cut_y = yi_prev
                    # For vertical lines
                    if xi == xi_prev:
                        # For vertical lines
                        new_x = xi_prev
                        if ymax < yi_prev:
                            new_y = ymax
                        else:
                            new_y = ymin
                    elif yi == yi_prev:
                        # For horizontal lines
                        new_y = yi_prev
                        if xmax < xi_prev:
                            new_x = xmax
                        else:
                            new_x = xmin
                    else:
                        # For other lines
                        straight_interp_f_of_x = lambda x: StraightLine(p1 = [xi, yi],
                                p2 = [xi_prev, yi_prev])(x)
                        straight_interp_f_of_y = lambda y: StraightLine(p1 = [yi, xi],
                                p2 = [yi_prev, xi_prev])(y)
                        # Interpolate at cut_x and cut_y
                        y_of_cut_x = straight_interp_f_of_x(cut_x)
                        x_of_cut_y = straight_interp_f_of_y(cut_y)
                        if (xmin <= x_of_cut_y) and (x_of_cut_y <= xmax):
                            new_x = x_of_cut_y
                            new_y = cut_y
                        else:
                            new_y = y_of_cut_x
                            new_x = cut_x
                    # }}}
                else:
                    new_x = xi
                    new_y = yi
                points.append([new_x, new_y])
                ordered_points.append(id_i)
                id_i += 1
                prev_in = True
                prev_out = False
            else:
                prev_out = True
                if prev_in:
                    # Compute new cut point {{{
                    # Cut at x?
                    if xi < xmin:
                        cut_x = xmin
                    elif xi > xmax:
                        cut_x = xmax
                    else:
                        cut_x = xi
                    # Cut at y?
                    if yi < ymin:
                        cut_y = ymin
                    elif yi > ymax:
                        cut_y = ymax
                    else:
                        cut_y = yi
                    # Make interpolation line
                    node_prev = nodes[nodeId - 1]
                    xi_prev = node_prev[0]
                    yi_prev = node_prev[1]
                    if xi == xi_prev:
                        # For vertical lines
                        new_x = xi
                        if ymax < yi:
                            new_y = ymax
                        else:
                            new_y = ymin
                    elif yi == yi_prev:
                        # For horizontal lines
                        new_y = yi
                        if xmax < xi:
                            new_x = xmax
                        else:
                            new_x = xmin
                    else:
                        # For other lines
                        straight_interp_f_of_x = lambda x: StraightLine(p1 = [xi, yi],
                                p2 = [xi_prev, yi_prev])(x)
                        straight_interp_f_of_y = lambda y: StraightLine(p1 = [yi, xi],
                                p2 = [yi_prev, xi_prev])(y)
                        # Interpolate at cut_x and cut_y
                        y_of_cut_x = straight_interp_f_of_x(cut_x)
                        x_of_cut_y = straight_interp_f_of_y(cut_y)
                        if (xmin <= x_of_cut_y) and (x_of_cut_y <= xmax):
                            new_x = x_of_cut_y
                            new_y = cut_y
                        else:
                            new_y = y_of_cut_x
                            new_x = cut_x
                    # }}}
                    # Add new cut point
                    points.append([new_x, new_y])
                    ordered_points.append(id_i)
                    id_i += 1
                    prev_in = False
        # }}}
        points = np.array(points)
        lines = [{"physical_name" : lineName,
                  "ordered_points" : ordered_points,
                  "type" : line["type"]}]
        return points, lines
    # }}}
    # Cut geometry {{{
    def CutGeometry(self, xmin = -10.0e+9, xmax = 10.0e+9,
            ymin = -10.0e+9, ymax = 10.0e+9, flip_boundaries = []):
        # Get cut lines {{{
        cutLines = []
        for line_i in self.lines:
            lineName_i = line_i["physical_name"]
            cutLines.append(self.CutLine(lineName_i, xmin = xmin,
                xmax = xmax, ymin = ymin, ymax = ymax))
        # }}}
        # Set new points and lines {{{
        for k1 in range(self.nlines):
            cutLine_i = cutLines[k1]
            points_i = cutLine_i[0]
            line_i = cutLine_i[1][0]
            if line_i["physical_name"] in flip_boundaries:
                points_i = np.flip(points_i, axis = 0)
            # Add first line {{{
            if k1 == 0:
                points = points_i
                # Raise error if points is empty {{{
                if len(points) == 0:
                    message = "The code assumes that at least two nodes"
                    message += " of the first line remain. So that it "
                    message += "does not disappears. This is for the "
                    message += "initialisation; it can be modified. "
                    message += "You can also modify the order of the lines."
                    raise ValueError(message)
                # }}}
                lines = [line_i]
            # }}}
            # Add other lines {{{
            else:
                # Add non empty lines {{{
                if len(line_i["ordered_points"]) > 0:
                    # Last node id
                    last_id = len(points) - 1
                    # Add new points
                    points = np.concatenate([points, points_i[1:, :]])
                    # Add new lines
                    ordered_points_i = line_i["ordered_points"]
                    for k2 in range(len(ordered_points_i)):
                        ordered_points_i[k2] += last_id
                    # If this is the last line
                    if k1 == self.nlines - 1:
                        ordered_points_i[-1] = 0
                    line_i["ordered_points"] = ordered_points_i
                # }}}
                # Add empty lines {{{
                else:
                    # Last node id and line type as straight
                    last_id = len(points) - 1
                    line_i["type"] = "straight"
                    # If this is the last line {{{
                    if k1 == self.nlines - 1:
                        # Add line order
                        line_i["ordered_points"] = [last_id, 0]
                        # Add new line
                    # }}}
                    # For other lines {{{
                    else:
                        # Get next line points
                        points_iplus1 = cutLines[k1 + 1][0]
                        points = np.concatenate([points,
                            np.array([points_iplus1[0, :]])])
                        # Add line order
                        line_i["ordered_points"] = [last_id, last_id + 1]
                    # }}}
                # }}}
                # Add new line
                lines.append(line_i)
            # }}}
        # }}}
        return points, lines
    # }}}
    # Scale {{{
    def Scale(self, axis, func, line_names = []):
        """
        Scale the points at the indicated axis with func:
            self.points[line_Set, axis] = func(self.points[lineSet, axis])
            lineSet is defined by the lines in line_names
        """
        # Test line_names:
        if not isinstance(line_names, list):
            raise TypeError("line_names must be a list (of strings).")
        # Modify points
        if line_names == []:
            self._points[:, axis] = func(self.points[:, axis])
        else:
            # Create a set of nodes to be modified
            point_set = set()
            for line in self.lines:
                line_name = line["physical_name"]
                if line_name in line_names:
                    for node_id in line["ordered_points"]:
                        point_set.add(node_id)
            point_set = list(point_set)
            # Modify only this points
            for node_id in point_set:
                self._points[node_id, axis] = func(self.points[node_id,
                    axis])
        return self.points
    # }}}
    # __add__ {{{
    def __add__(self, geom):
        # New points
        newPoints = np.concatenate((self._points, geom._points),
                axis = 0)
        # New lines
        geomLines = geom.lines
        for k1 in range(len(geomLines)):
            line = geomLines[k1]
            for k2 in range(len(line["ordered_points"])):
                line["ordered_points"][k2] += self.npoints
        newLines = self.lines + geomLines
        sumGeom = Geometry(self.name, newPoints, newLines)
        return sumGeom
    # }}}
    # __str__ {{{
    def __str__(self):
        message = np.array_str(self.points) + "\n"
        message = str(self.lines) + "\n"
        return message
    # }}}
    # Change name {{{
    def change_name(self, name):
        self._name = name
        return
    # }}}
    # Get ordered node coordinates {{{
    def GetOrderedNodeCoords(self, line_index):
        # Get coordinates and line ordered points
        coords = self.points
        ordered_points = self.lines[line_index]["ordered_points"]
        # Make array of ordered coordinates
        ord_coords = []
        for k1 in ordered_points:
            ord_coords.append(coords[k1, :])
        ord_coords = np.array(ord_coords)
        return ord_coords
    # }}}
# }}}

# Reconstruct class {{{
class Reconstruct(Geometry):
    # __init__ {{{
    def __init__(self, name, oldpoints, oldlines, asterMeshFile,
            **kwargs):
        lineType = kwargs.get("lineType", "bspline")
        # Get new list of points {{{
        # Read old mesh
        oldmesh = AsterMesh(asterMeshFile)
        # Set up new points {{{
        # Get the list of nodes based on the physical names of the
        # original mesh
        nodeList = list()
        for line in oldlines:
            for node in oldmesh.nodeSets[line["physical_name"]]:
                if not node in nodeList:
                    nodeList.append(node)
        # Get the list of new nodes and a dictionary as a function
        # relating the old node name and the new node id
        nodeDict = {}
        newNodes = []
        nodeId = -1
        for node in nodeList:
            nodeId += 1
            nodeDict[node] = {"id" : nodeId}
            newNodes.append(oldmesh.GetNodeCoordinates(node))
        points = np.array(newNodes)
        # }}}
        # Update the lines with the new nodes {{{
        lines = []
        for line in oldlines:
            # Change line type
            newLine = {}
            newLine["type"] = lineType
            newLine["physical_name"] = line["physical_name"]
            # Set the new ordered list of points
            ordered_points = list()
            lineNodes = oldmesh.nodeSets[line["physical_name"]]
            for node in lineNodes:
                ordered_points.append(nodeDict[node]["id"])
            # Test/set non self intersection spline {{{
            new_ordered_points = SplineNoLoops(points, ordered_points)
            # }}}
            newLine["ordered_points"] = new_ordered_points
            lines.append(newLine)
        # }}}
        # }}}
        # Call geometry
        super().__init__(name, points, lines, **kwargs)
        return
    # }}}
# }}}

# Functions {{{
# Mesh field distance--threshold to curve {{{
def MeshFieldDistanceCurve(f0, curveList, sizeMin, distMin,
        sizeGrwRate, distGrwRate, sizeMax, numFields = 0,
        growthStyle = 'linear', sampling = 10000, at = "CurvesList"):
    # Initialisation {{{
    mesh = gmsh.model.mesh
    fieldList = []
    ff = f0
    # }}}
    # Make base distance field {{{
    di_fi = ff
    mesh.field.add("Distance", ff)
    mesh.field.setNumbers(ff, at, curveList)
    mesh.field.setNumber(ff, "Sampling", sampling)
    # }}}
    # Thresholds {{{
    size = sizeMin
    dist = distMin
    addField = True
    k1 = 0
    while addField:
        if (numFields > 0) and (numFields < k1) :
            addField = False
        else:
            k1 += 1
        # Set mesh size field as threshold
        ff += 1
        fieldList.append(ff)
        mesh.field.add("Threshold", ff)
        mesh.field.setNumber(ff, "InField", di_fi)
        mesh.field.setNumber(ff, "SizeMin", size)
        mesh.field.setNumber(ff, "SizeMax", sizeMax)
        mesh.field.setNumber(ff, "DistMin", dist)
        # Update sizes
        if growthStyle == 'linear':
            size += sizeMin*sizeGrwRate
            dist += distMin*distGrwRate
            mesh.field.setNumber(ff, "DistMax", dist)
        elif growthStyle == 'geometric':
            size *= sizeGrwRate
            dist *= distGrwRate
            mesh.field.setNumber(ff, "DistMax", dist)
        else:
            raise("Error: growth style not available")
        if size > sizeMax:
            addField = False
    # }}}
    # Background mesh {{{
    bacFieldNum = ff + 1
    mesh.field.add("Min", bacFieldNum)
    mesh.field.setNumbers(bacFieldNum, "FieldsList", fieldList)
    mesh.field.setAsBackgroundMesh(bacFieldNum)
    fieldList.append(bacFieldNum)
    # }}}
    return fieldList
# }}}

# Spline with no loops {{{
def SplineNoLoops(nodes, ordered_points):
    # List of nodes that compose the spline
    s_points = []
    for k1 in ordered_points:
        xx = nodes[k1, 0]
        yy = nodes[k1, 1]
        s_points.append((xx, yy))
    # Spline, test simplicity
    spline = LineString(s_points)
    if spline.is_simple:
        return ordered_points
    # Separate spline loops
    splines = unary_union(spline)
    # get splines that do not form a loop
    splines_no_loop = []
    for s_i in splines.geoms:
        if not s_i.is_ring:
            splines_no_loop.append(s_i)
    # Merge the new splines
    new_coords = []
    for line_i in splines_no_loop:
        new_coords += line_i.coords[:]
    # Set new ordered points
    old_coords = list(spline.coords)
    new_ordered_points = []
    for p_i in new_coords:
        for k1 in range(len(ordered_points)):
            if old_coords[k1] == p_i:
                new_ordered_points.append(ordered_points[k1])
    # Recursive function
    return_order = SplineNoLoops(nodes, new_ordered_points)
    return return_order
# }}}

# Spline with no sharp angles {{{
def SplineNoSharpAngles(nodes, ordered_points, min_angle_deg = 20.0):
    # Set coordinates {{{
    coords = []
    for k1 in ordered_points:
        xx = nodes[k1, 0]
        yy = nodes[k1, 1]
        coords.append([xx, yy])
    coords = np.array(coords)
    # }}}
    # Compute angles {{{
    angles = [np.pi]
    for k1 in range(1, len(coords) - 1):
        x1 = coords[k1 - 1, 0]
        x2 = coords[k1    , 0]
        x3 = coords[k1 + 1, 0]
        y1 = coords[k1 - 1, 1]
        y2 = coords[k1    , 1]
        y3 = coords[k1 + 1, 1]
        v21 = np.array([x1 - x2, y1 - y2])
        v23 = np.array([x3 - x2, y3 - y2])
        uni21 = v21/np.linalg.norm(v21)
        uni23 = v23/np.linalg.norm(v23)
        angles.append(np.arccos(np.dot(uni21, uni23)))
    angles.append(np.pi)
    angles = np.array(angles)*180.0/np.pi
    # }}}
    min_angle = min(angles)
    if min_angle > min_angle_deg:
        return ordered_points
    else:
        new_ordered_points = []
        for k1 in range(len(ordered_points)):
            if angles[k1] >= min_angle_deg:
                new_ordered_points.append(ordered_points[k1])
    return_order = SplineNoSharpAngles(nodes, new_ordered_points,
            min_angle_deg = min_angle_deg)
    return return_order
# }}}

# Circle bottom y {{{
def CircleBottomY(x, r, cx, cy):
    return cy - np.sqrt(r**2.0 - (x - cx)**2.0)
# }}}

# Make points from dimensionless geometry circle and rectangle {{{
def MakePoints(r, Ge, Be, h_fact, w_fact):
    # Compute length
    lfstar = 1.0/Ge
    l0star = Be/Ge
    # Get circle point based on its size and on the length
    if r >= lfstar:
        circ_x = lfstar
        circ_y = CircleBottomY(circ_x, r, 0.0, r)
    else:
        circ_x = r
        circ_y = r
    # Make points
    cir_points = [[0.0, 0.0],
                  [circ_x, circ_y],
                  [circ_x, h_fact*lfstar],
                  [0.0, h_fact*lfstar],
                  [0.0, r]]
    rec_points = [[0.0, 0.0],
                  [w_fact*lfstar, 0.0],
                  [w_fact*lfstar, -w_fact*lfstar],
                  [0.0, -w_fact*lfstar]]
    aux_points = [[0.0, 0.0],
                  [1.1*l0star, 0.0]]
    dict_points = {"esc" : cir_points,
                   "mai" : rec_points,
                   "aux" : aux_points}
    return dict_points
# }}}
# }}}

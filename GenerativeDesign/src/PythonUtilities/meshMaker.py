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

# In-house modules
sys.path.append(os.path.dirname(os.path.abspath(__file__))[:-4] + 'PythonUtilities/')
import ReadGmsh as rgmsh
# }}}

# Class mesh maker size gradient {{{
class MeshMakerSizeGradient(object):
    # Properties {{{
    @property
    def nonPhysicalNodeGroups(self):
        return self._nonPhysicalNodeGroups
    @property
    def name_ref_point(self):
        return self._name_ref_point
    @property
    def name_main_line(self):
        return self._name_main_line
    @property
    def list_name_points(self):
        return self._list_name_points
    @property
    def maxLength(self):
        return self._maxLength
    @property
    def cutAt_y(self):
        return self._cutAt_y
    @property
    def mirrorData(self):
        return self._mirrorData
    @property
    def name(self):
        return self._name
    @property
    def lcMin(self):
        return self._lcMin
    @property
    def lcMax(self):
        return self._lcMax
    @property
    def distMin(self):
        return self._distMin
    @property
    def sizeGrwRate(self):
        return self._sizeGrwRate
    @property
    def distGrwRate(self):
        return self._distGrwRate
    @property
    def numFields(self):
        return self._numFields
    @property
    def growthStyle(self):
        return self._growthStyle
    @property
    def sampling(self):
        return self._sampling
    @property
    def at(self):
        return self._at
    @property
    def fullGeo(self):
        return self._fullGeo
    @property
    def aux_lines(self):
        return self._aux_lines
    # }}}
    # __init__ {{{
    def __init__(self, lcMin, lcMax, **kwargs):
        # Parameters {{{
        self._lcMin = lcMin
        self._lcMax = lcMax
        self._name = kwargs.get("name", "geo")
        if len(self.name) > 4:
            raise("Error: the mesh name cannot be longer than four characters.")
        distMinFactor = kwargs.get("distMinFactor", 1.0)
        self._distMin = kwargs.get("distMin", lcMin*distMinFactor)
        self._sizeGrwRate = kwargs.get("sizeGrwRate", 1.0)
        self._distGrwRate = kwargs.get("distGrwRate", 1.0)
        self._numFields = kwargs.get("numFields",
                int(self.maxLength/self._distMin))
        self._growthStyle = kwargs.get("growthStyle", "linear")
        self._sampling = kwargs.get("sampling", 10000)
        self._at = kwargs.get("at", "CurvesList")
        self._fullGeo = kwargs.get("fullGeo", False)
        self._cutAt_y = kwargs.get("cutAt_y", None)
        self._nonPhysicalNodeGroups = {}
        # }}}
        return
    # }}}
    # set auxiliary line {{{
    def set_aux_lines(self, linelist):
        self._aux_lines = linelist
        return
    # }}}
    # Set mesh fields {{{
    def SetMeshFields(self, f0 = 1):
        # GMSH environment parameters {{{
        model = gmsh.model
        geo = model.occ
        mesh = model.mesh
        # }}}
        # Parameters {{{
        lcMin       = self.lcMin
        lcMax       = self.lcMax
        aux_lines   = self.aux_lines
        distMin     = self.distMin
        sizeGrwRate = self.sizeGrwRate
        distGrwRate = self.distGrwRate
        numFields   = self.numFields
        growthStyle = self.growthStyle
        sampling    = self.sampling
        at          = self.at
        maxLength   = self.maxLength
        # }}}
        numFields = int(maxLength/(lcMin*distMin))
        fieldList = MeshFieldDistanceCurve(mesh, f0, aux_lines, lcMin,
                distMin, sizeGrwRate, distGrwRate, lcMax, numFields,
                growthStyle = growthStyle, sampling = sampling, at = at)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
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
        # Correct groups of nodes if there are {{{
        if hasattr(self, 'list_name_points'):
            list_name_points = self.list_name_points
            name = self.name
            if len(list_name_points) > 0:
                physical_groups = gmsh.model.getPhysicalGroups(0)
                # Change entity type code: from 8 to 7
                with open(pathname, 'r') as fle:
                    lines = fle.readlines()
                    # Find each node group line
                    for group in list_name_points:
                        # Get the tag of the group
                        group_name = name + group
                        for _, gr_tag in physical_groups:
                            if group_name == gmsh.model.getPhysicalName(0, gr_tag):
                                tag = gr_tag
                                break
                        gr_nodes = gmsh.model.mesh.getNodesForPhysicalGroup(0, tag)
                        gr_nodes = gr_nodes[0]
                        index = lines.index(group_name + '\n')
                        infoline = lines[index - 1].split()
                        numNodes = eval(infoline[-1])
                        numLines = int(numNodes/2) + numNodes%2
                        for k1 in range(numLines):
                            nd1 = k1*2
                            nd2 = k1*2 + 1
                            line = lines[index + 1 + k1].split()
                            line = [eval(var) for var in line]
                            line[0] = 7
                            line[1] = gr_nodes[nd1]
                            if len(line) == 8:
                                line[4] = 7
                                line[5] = gr_nodes[nd2]
                            newline = ''
                            for var in line:
                                newline += "{:10d}".format(var)
                            lines[index + 1 + k1] = newline + '\n'
                with open(pathname, 'w') as fle:
                    for line in lines:
                        fle.write(line)
        # }}}
        # Add non-physical group of nodes {{{
        nonPhysicalNodeGroups = self.nonPhysicalNodeGroups
        if len(nonPhysicalNodeGroups) > 0:
            # Read file
            with open(pathname, 'r') as fle:
                lines = fle.readlines()
            # Find the location of the groups and set the number of the
            # non physical groups
            index = lines.index("  2477\n")
            index += 1
            nextline = lines[index]
            numGroups = []
            read = True
            while read:
                groupinfo = [eval(val) for val in nextline.split()]
                numGroups.append(groupinfo[0])
                numNodes = groupinfo[-1]
                numLines = int(numNodes/2) + numNodes%2
                index += 1 + 1 + numLines
                nextline = lines[index]
                if nextline == '    -1\n':
                    read = False
            numGroup = max(numGroups)
            # Set lines to be added
            newLines = ['    -1\n', '  2477\n']
            for name in nonPhysicalNodeGroups:
                numGroup += 1
                nodeIds = nonPhysicalNodeGroups[name]
                numNodes = len(nodeIds)
                lineToAdd = [numGroup, 0, 0, 0, 0, 0, 0, numNodes]
                newLine = ''
                for var in lineToAdd:
                    newLine += "{:10d}".format(var)
                newLines.append(newLine + '\n')
                newLines.append(name + '\n')
                numLines = int(numNodes/2) + numNodes%2
                for k1 in range(numLines):
                    nd1 = nodeIds[k1*2]
                    lineToAdd = [7, nd1, 0, 0]
                    try:
                        nd2 = nodeIds[k1*2 + 1]
                        lineToAdd += [7, nd2, 0, 0]
                    except IndexError:
                        pass
                    newLine = ''
                    for var in lineToAdd:
                        newLine += "{:10d}".format(var)
                    newLines.append(newLine + '\n')
            newLines.append('    -1')
            lines = lines + newLines
            with open(pathname, 'w') as fle:
                for line in lines:
                    fle.write(line)
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
    # Get list of main lines to remesh {{{
    def GetListMainLines(self, prevMeshFile):
        # Object parameters {{{
        name_main_line    = self.name_main_line
        name_ref_point    = self.name_ref_point
        name        = self.name
        # }}}
        # Get information from previous mesh {{{
        prevMesh = AsterMesh(prevMeshFile)
        # Get the node names of the curve to be remeshed (and split)
        nodeNames = prevMesh.nodeSets[name + name_main_line]
        node_1 = nodeNames[0]
        node_f = nodeNames[-1]
        # Get the reference points
        refPointNames = prevMesh.nodeSets[name + name_ref_point]
        # Remove the end points from refPointNames
        if node_1 in refPointNames:
            refPointNames.remove(node_1)
        if node_f in refPointNames:
            refPointNames.remove(node_f)
        # Make list of reference lines
        if len(refPointNames) == 0:
            refLines = [nodeNames]
        else:
            # Find the indices of each reference point
            indices = [0]
            for refPointName in refPointNames:
                indices.append(nodeNames.index(refPointName))
            indices.sort()
            refLines = []
            for k1 in range(len(indices) - 1):
                refLines.append(nodeNames[indices[k1]: indices[k1 + 1] + 1])
            refLines.append(nodeNames[indices[-1]:])
        # }}}
        return refLines
    # }}}
    # Add node group from list if coordinates {{{
    def Add_nonPhysical_node_group_from_coordinates(self, coordList,
            name):
        """
        Find the node id's of the closest node to each coordinate.
        """
        if len(name) > 4:
            raise("The length of name must be less or equal to 4.")
        # Get node coordinates
        nodeTags, nodeCoords, _ = gmsh.model.mesh.getNodes()
        numNodes = int(len(nodeCoords)/3)
        nodeCoords = np.array([nodeCoords[k1*3: k1*3 + 3] for k1 in range(numNodes)])
        # Get the node ids related to the closest nodes to the indicated
        # coordinates
        nodeIds = []
        for coord in coordList:
            cx = coord[0]
            cy = coord[1]
            distances = []
            for node in nodeCoords:
                nx = node[0]
                ny = node[1]
                distances.append(np.sqrt((cx - nx)**2.0 + (cy - ny)**2.0))
            distances = np.array(distances)
            id_ = np.argmin(distances)
            nodeIds.append(nodeTags[id_])
        # Save information related to the specified non physical group
        self._nonPhysicalNodeGroups[self.name + name] = nodeIds
        return
    # }}}
# }}}

# Sub class of MeshMakerSizeGradient: Circle {{{
class Circle(MeshMakerSizeGradient):
    # Properties {{{
    @property
    def R(self):
        return self._R
    @property
    def p_aux_1(self):
        return self._p_aux_1
    @property
    def p_aux_2(self):
        return self._p_aux_2
    @property
    def name_top(self):
        return self._name_top
    @property
    def name_middle(self):
        return self._name_middle
    # }}}
    # __init__ {{{
    def __init__(self, R, lcMin, lcMax, **kwargs):
        def ylow_circ(x, r, cx, cy):
            return -np.sqrt(r**2.0 - (x - cx)**2.0) + cy
        self._R = R
        self._maxLength = R
        super().__init__(lcMin, lcMax, **kwargs)
        if "fine_x" in kwargs:
            fine_x = kwargs["fine_x"]
            p_aux_1 = [ fine_x, ylow_circ(fine_x, R, 0.0, R)]
            p_aux_2 = [-fine_x, ylow_circ(fine_x, R, 0.0, R)]
        elif ("fine_x1" in kwargs) and ("fine_x2" in kwargs):
            fine_x1 = kwargs["fine_x1"]
            fine_x2 = kwargs["fine_x2"]
            p_aux_1 = [fine_x1, ylow_circ(fine_x1, R, 0.0, R)]
            p_aux_2 = [fine_x2, ylow_circ(fine_x2, R, 0.0, R)]
        elif ("p_aux_1" in kwargs) and ("p_aux_2" in kwargs):
            p_aux_1 = kwargs["p_aux_1"]
            p_aux_2 = kwargs["p_aux_2"]
        else:
            raise("Error: (fine_x) or (fine_x1, fine_x2) or (p_aux_1, p_aux_2) is not defined.")
        self._p_aux_1 = p_aux_1
        self._p_aux_2 = p_aux_2
        # Set group names from kwargs {{{
        name_ref_point = kwargs.get("name_ref_point", "_pcb")
        if len(name_ref_point) > 4:
            raise("The length of name_ref_point must be less or equal to 4.")
        name_main_line = kwargs.get("name_main_line", "_cir")
        if len(name_main_line) > 4:
            raise("The length of name_main_line must be less or equal to 4.")
        name_top = kwargs.get("name_top", "_top")
        if len(name_top) > 4:
            raise("The length of name_top must be less or equal to 4.")
        name_middle = kwargs.get("name_middle", "_mid")
        if len(name_middle) > 4:
            raise("The length of name_middle must be less or equal to 4.")
        self._name_ref_point = name_ref_point
        self._name_main_line = name_main_line
        self._name_middle = name_middle
        self._name_top    = name_top
        self._list_name_points = [name_ref_point]
        # }}}
        return
    # }}}
    # Make geometry {{{
    def MakeGeometry(self, **kwargs):
        # Set group names from kwargs {{{
        name_ref_point = self.name_ref_point
        name_main_line = self.name_main_line
        name_middle = self.name_middle
        name_top    = self.name_top
        name_ref_point = kwargs.get("name_ref_point", name_ref_point)
        if len(name_ref_point) > 4:
            raise("The length of name_ref_point must be less or equal to 4.")
        name_main_line = kwargs.get("name_main_line", name_main_line)
        if len(name_main_line) > 4:
            raise("The length of name_main_line must be less or equal to 4.")
        name_top = kwargs.get("name_top", name_top)
        if len(name_top) > 4:
            raise("The length of name_top must be less or equal to 4.")
        name_middle = kwargs.get("name_middle", name_middle)
        if len(name_middle) > 4:
            raise("The length of name_middle must be less or equal to 4.")
        self._name_ref_point = name_ref_point
        self._name_main_line = name_main_line
        self._name_middle = name_middle
        self._name_top    = name_top
        self._list_name_points = [name_ref_point]
        # }}}
        # Parameters {{{
        R       = self.R
        p_aux_1 = self.p_aux_1
        p_aux_2 = self.p_aux_2
        fullGeo = self.fullGeo
        name = self.name
        cutAt_y = self.cutAt_y
        # }}}
        # GMSH Initialisation {{{
        gmsh.initialize()
        gmsh.clear()
        gmsh.model.add(name)
        model = gmsh.model
        geo = model.occ
        mesh = model.mesh
        # }}}
        # Make base face {{{
        if cutAt_y:
            if cutAt_y < R:
                Tx = np.sqrt(R**2.0 - (cutAt_y - R)**2.0)
                Ty = cutAt_y
            else:
                Tx = R
                Ty = R
        else:
            Tx = R
            Ty = R
        # Points
        p_cc = geo.addPoint(0.0,   R, 0.0, R)
        p_cb = geo.addPoint(0.0, 0.0, 0.0, R)
        p_rt = geo.addPoint( Tx,  Ty, 0.0, R)
        p_ct = geo.addPoint(0.0,  Ty, 0.0, R)
        # Lines
        l_circle = geo.addCircleArc(p_cb, p_cc, p_rt)
        l_top    = geo.addLine(p_rt, p_ct)
        l_middle = geo.addLine(p_ct, p_cb)
        l_loop   = geo.addCurveLoop([l_circle, l_top, l_middle])
        # Face
        face = geo.addPlaneSurface([l_loop])
        # }}}
        # For full geometry {{{
        if fullGeo:
            # Points
            mp_lc = geo.addPoint(-Tx,  Ty, 0.0, R)
            # Lines
            ml_circle = geo.addCircleArc(mp_lc, p_cc, p_cb)
            ml_top    = geo.addLine(p_ct, mp_lc)
            ml_loop   = geo.addCurveLoop([ml_circle, -l_middle, ml_top])
            # Face
            mface = geo.addPlaneSurface([ml_loop])
        # }}}
        # Create auxiliary arc where the mesh will be finer {{{
        aux_p_1 = geo.addPoint(p_aux_1[0], p_aux_1[1], 0.0, R)
        aux_p_2 = geo.addPoint(p_aux_2[0], p_aux_2[1], 0.0, R)
        aux_c = geo.addCircleArc(aux_p_1, p_cc, aux_p_2)
        self.set_aux_lines([aux_c])
        # }}}
        # Synchronize
        geo.synchronize()
        # Make physical groups {{{
        # Points
        self.PhysicalGroup(0, [p_cb], name + name_ref_point)
        # Lines
        lst_middle = [l_middle]
        if fullGeo:
            lst_main = [l_circle, ml_circle]
            lst_top    = [l_top   , ml_top   ]
            lst_face   = [face    , mface    ]
        else:
            lst_main = [l_circle]
            lst_top    = [l_top   ]
            lst_face   = [face]
        self.PhysicalGroup(1, lst_main, name + name_main_line)
        self.PhysicalGroup(1, lst_top   , name + name_top )
        self.PhysicalGroup(1, lst_middle, name + name_middle )
        self._mainLine = lst_main
        # Faces
        self.PhysicalGroup(2, lst_face, name)
        # }}}
        return
    # }}}
    # Remesh {{{
    def Remesh(self, prevMeshFile, **kwargs):
        # Object parameters {{{
        name_ref_point = self.name_ref_point
        name_main_line = self.name_main_line
        name_middle = self.name_middle
        name_top    = self.name_top
        lcMin       = self.lcMin
        lcMax       = self.lcMax
        name        = self.name
        fullGeo     = self.fullGeo
        cutAt_y     = self.cutAt_y
        R           = self.R
        p_aux_1     = self.p_aux_1
        p_aux_2     = self.p_aux_2
        maxLength   = self.maxLength
        # }}}
        # Get kwargs {{{
        x_aux_1 = kwargs.get("x_aux_1", p_aux_1[0])
        x_aux_2 = kwargs.get("x_aux_2", p_aux_2[0])
        # }}}
        # GMSH Initialisation {{{
        gmsh.finalize()
        gmsh.initialize()
        gmsh.clear()
        gmsh.model.add(name)
        model = gmsh.model
        geo = model.occ
        mesh = model.mesh
        # }}}
        # Get information from previous mesh {{{
        prevMesh = AsterMesh(prevMeshFile)
        refPointNames = prevMesh.nodeSets[name + name_ref_point]
        refLines = self.GetListMainLines(prevMeshFile)
        refLineNodes = []
        for line in refLines:
            refLineNodes.append(prevMesh.GetNodeSetCoordinates(nodeSet = line))
        # }}}
        # Geometry reconstruction {{{
        if cutAt_y:
            # Recover of _cir points {{{
            # Select only the nodes below cutAt_y
            cutRefLines = []
            for line in refLines:
                cutLine = []
                for nodeName in line:
                    node = prevMesh.GetNodeCoordinates(nodeName)
                    yi = node[1]
                    # Create node if the y coordinate is below cutAt_y
                    if yi < cutAt_y:
                        cutLine.append(nodeName)
                cutRefLines.append(cutLine)
            cutRefLineNodes = []
            for line in cutRefLines:
                cutRefLineNodes.append(prevMesh.GetNodeSetCoordinates(nodeSet = line))
            # }}}
            # Add necessary points to complete the spline to make _cir {{{
            if fullGeo:
                # Add two points to ensure that the y0 = cutAt_y at each
                # side of the _cir
                # First point
                points = cutRefLineNodes[0]
                p0 = points[0]
                p1 = points[1]
                m = (p1[1] - p0[1])/(p1[0] - p0[0])
                b = p1[1] - m*p1[0]
                xi = (cutAt_y - b)/m
                cutRefLineNodes[0] = np.vstack((np.array([xi, cutAt_y]),
                    cutRefLineNodes[0]))
                # Last point
                points = cutRefLineNodes[-1]
                p0 = points[-2]
                p1 = points[-1]
                m = (p1[1] - p0[1])/(p1[0] - p0[0])
                b = p1[1] - m*p1[0]
                xi = (cutAt_y - b)/m
                cutRefLineNodes[1] = np.vstack((cutRefLineNodes[1],
                    np.array([xi, cutAt_y])))
            else:
                # Add one point to the finish _cir
                if len(cutRefLineNodes) > 1:
                    raise("Error: too many cutRefLineNodes")
                points = cutRefLineNodes[0]
                # Find the correct side
                p0 = points[0]
                pf = points[-1]
                if p0[1] < pf[1]:
                    pa = points[-2]
                    m = (pf[1] - pa[1])/(pf[0] - pa[0])
                    b = pf[1] - m*pf[0]
                    xi = (cutAt_y - b)/m
                    cutRefLineNodes[0] = np.vstack((cutRefLineNodes[0],
                        np.array([xi, cutAt_y])))
                    cutRefLineNodes[0][0, 0] = 0.0
                    p_rt = len(cutRefLineNodes[0])
                    p_cb = 1
                else:
                    pa = points[1]
                    m = (p0[1] - pa[1])/(p0[0] - pa[0])
                    b = p0[1] - m*p0[0]
                    xi = (cutAt_y - b)/m
                    cutRefLineNodes[0] = np.vstack((np.array([xi, cutAt_y],
                        cutRefLineNodes[0])))
                    cutRefLineNodes[0][-1, 0] = 0.0
                    p_rt = 1
                    p_cb = len(cutRefLineNodes[0])
            # }}}
            # Make reference point and main lines {{{
            # Make reference point
            xi, yi = prevMesh.GetNodeCoordinates(refPointNames[0])
            p_cb = geo.addPoint(xi, yi, 0.0, maxLength)
            # Make first line points
            firstLineIds = []
            for coords in cutRefLineNodes[0]:
                xi = coords[0]
                yi = coords[1]
                firstLineIds.append(geo.addPoint(xi, yi, 0.0, maxLength))
            l_main1 = geo.addBSpline(firstLineIds)
            if fullGeo:
                secondLineIds = [firstLineIds[-1]]
                for coords in cutRefLineNodes[1][1 :]:
                    xi = coords[0]
                    yi = coords[1]
                    secondLineIds.append(geo.addPoint(xi, yi, 0.0, maxLength))
                l_main2 = geo.addBSpline(secondLineIds)
                l_main_list = [l_main1, l_main2]
            else:
                l_circle = l_main1
                l_main_list = [l_circle]
            # }}}
            # Make the other lines and face {{{
            if fullGeo:
                l_top = geo.addLine(firstLineIds[0], secondLineIds[-1])
                l_loop = geo.addCurveLoop(l_main_list + [l_top])
            else:
                p_cb = firstLineIds[0]
                p_cc = geo.addPoint(0.0, cutAt_y, 0.0, R)
                l_top = geo.addLine(firstLineIds[-1], p_cc)
                l_middle = geo.addLine(p_cc, p_cb)
                l_loop = geo.addCurveLoop([l_circle, l_top, l_middle])
            face = geo.addPlaneSurface([l_loop])
            # }}}
        else:
            model = rgmsh.Remesh([gr_names], [file_names], [lst_lcMax],
                    [name], model, geoType = 'occ')
        # }}}
        # Make auxiliary spline {{{
        nds = prevMesh.GetNodeSetCoordinates(setName = name
                + name_main_line)
        auxPointIds = []
        for node in nds:
            xi = node[0]
            yi = node[1]
            # Create node if x_aux_1 < xi < x_aux_2
            if (x_aux_1 < xi) and (xi < x_aux_2):
                # Create node if x_aux_1 < yi < x_aux_2
                if (x_aux_1 < yi) and (yi < x_aux_2):
                    auxPointIds.append(geo.addPoint(xi, yi, 0.0, R))
        aux_c = geo.addBSpline(auxPointIds)
        self.set_aux_lines([aux_c])
        # }}}
        # Synchronize
        geo.synchronize()
        # Make physical groups {{{
        # Points
        self.PhysicalGroup(0, [p_cb], name + name_ref_point)
        # Lines
        #lst_main = [l_main1, l_main2]
        lst_main = l_main_list
        lst_top    = [l_top   ]
        lst_face   = [face]
        if not fullGeo:
            lst_middle = [l_middle]
            self.PhysicalGroup(1, lst_middle, name + name_middle )
        self.PhysicalGroup(1, lst_main, name + name_main_line)
        self.PhysicalGroup(1, lst_top   , name + name_top )
        # Faces
        self.PhysicalGroup(2, lst_face, name)
        # }}}
        return
    # }}}
# }}}

# Sub class of MeshMakerSizeGradient: Ellipse {{{
class Ellipse(MeshMakerSizeGradient):
    # Properties {{{
    @property
    def H(self):
        return self._H
    @property
    def x_func(self):
        return self._x_func
    @property
    def y_func(self):
        return self._y_func
    @property
    def sign(self):
        return self._sign
    @property
    def Cx(self):
        return self._Cx
    @property
    def Cy(self):
        return self._Cy
    @property
    def Rx(self):
        return self._Rx
    @property
    def Ry(self):
        return self._Ry
    @property
    def p_aux_1(self):
        return self._p_aux_1
    @property
    def p_aux_2(self):
        return self._p_aux_2
    @property
    def name_dy(self):
        return self._name_dy
    @property
    def name_dxM(self):
        return self._name_dxM
    @property
    def name_dxR(self):
        return self._name_dxR
    # }}}
    # __init__ {{{
    def __init__(self, Rx, Ry, Cx, Cy, H, lcMin, lcMax, sign = -1.0, **kwargs):
        def x_ellipse(y):
            return Cx + Rx*np.sqrt(1.0 - ((y - Cy)/Ry)**2.0)
        def y_ellipse(x):
            return Cy + sign*Ry*np.sqrt(1.0 - ((x - Cx)/Rx)**2.0)
        self._H = H
        self._Rx = Rx
        self._Ry = Ry
        self._Cx = Cx
        self._Cy = Cy
        self._sign = sign
        self._x_func = x_ellipse
        self._y_func = y_ellipse
        self._maxLength = max(Rx, Ry)
        super().__init__(lcMin, lcMax, **kwargs)
        self.Update_p_aux(**kwargs)
        # Set group names from kwargs {{{
        name_ref_point = kwargs.get("name_ref_point", "_pcb")
        if len(name_ref_point) > 4:
            raise("The length of name_ref_point must be less or equal to 4.")
        name_main_line = kwargs.get("name_main_line", "_cir")
        if len(name_main_line) > 4:
            raise("The length of name_main_line must be less or equal to 4.")
        name_dy = kwargs.get("name_dy", "_dy")
        if len(name_dy) > 4:
            raise("The length of name_dy must be less or equal to 4.")
        name_dxR = kwargs.get("name_dxR", "_dxR")
        if len(name_dxR) > 4:
            raise("The length of name_dxR must be less or equal to 4.")
        name_dxM = kwargs.get("name_dxM", "_dxM")
        if len(name_dxM) > 4:
            raise("The length of name_dxM must be less or equal to 4.")
        self._name_ref_point = name_ref_point
        self._name_main_line = name_main_line
        self._name_dxM = name_dxM
        self._name_dxR = name_dxR
        self._name_dy    = name_dy
        self._list_name_points = [name_ref_point]
        # }}}
        return
    # }}}
    # Make geometry {{{
    def MakeGeometry(self, **kwargs):
        # Set group names from kwargs {{{
        name_ref_point = self.name_ref_point
        name_main_line = self.name_main_line
        name_dxM = self.name_dxM
        name_dxR = self.name_dxR
        name_dy    = self.name_dy
        name_ref_point = kwargs.get("name_ref_point", name_ref_point)
        if len(name_ref_point) > 4:
            raise("The length of name_ref_point must be less or equal to 4.")
        name_main_line = kwargs.get("name_main_line", name_main_line)
        if len(name_main_line) > 4:
            raise("The length of name_main_line must be less or equal to 4.")
        name_dy = kwargs.get("name_dy", name_dy)
        if len(name_dy) > 4:
            raise("The length of name_dy must be less or equal to 4.")
        name_dxR = kwargs.get("name_dxR", name_dxR)
        if len(name_dxR) > 4:
            raise("The length of name_dxR must be less or equal to 4.")
        name_dxM = kwargs.get("name_dxM", name_dxM)
        if len(name_dxM) > 4:
            raise("The length of name_dxM must be less or equal to 4.")
        self._name_ref_point = name_ref_point
        self._name_main_line = name_main_line
        self._name_dxM = name_dxM
        self._name_dxR = name_dxR
        self._name_dy    = name_dy
        self._list_name_points = [name_ref_point]
        # }}}
        # Parameters {{{
        H       = self.H
        Rx      = self.Rx
        Ry      = self.Ry
        Cx      = self.Cx
        Cy      = self.Cy
        p_aux_1 = self.p_aux_1
        p_aux_2 = self.p_aux_2
        fullGeo = self.fullGeo
        name = self.name
        cutAt_y = self.cutAt_y
        x_func = self.x_func
        y_func = self.y_func
        maxLength = self.maxLength
        sign = self.sign
        lcMax = self.lcMax
        # }}}
        # GMSH Initialisation {{{
        gmsh.initialize()
        gmsh.clear()
        gmsh.model.add(name)
        model = gmsh.model
        geo = model.occ
        mesh = model.mesh
        # }}}
        # Make base face {{{
        # Major axis reference coordinates
        if Rx > Ry:
            mx = Cx + Rx/2.0
            my = Cy
        else:
            mx = Cy
            my = Cy + Ry/2.0
        # End points with and without cut at y
        Tx = Cx + Rx
        Ty = Cy
        if cutAt_y:
            ylims = [Cy, Cy + Ry*sign]
            ymin = min(ylims)
            ymax = max(ylims)
            if ymin < cutAt_y and cutAt_y < ymax:
                Tx = x_func(cutAt_y)
                Ty = cutAt_y
        # Points
        p_e_c = geo.addPoint(Cx, Cy, 0.0, lcMax)
        p_e_l = geo.addPoint(Cx, Cy + sign*Ry, 0.0, lcMax)
        p_e_r = geo.addPoint(Tx, Ty, 0.0, lcMax)
        p_v = geo.addPoint(Cx, Ty, 0.0, lcMax)

        p_e_major = geo.addPoint(mx, my, 0.0, lcMax)

        p_l = geo.addPoint(Cx, Cy - sign*(H - Ry), 0.0, lcMax)
        p_r = geo.addPoint(Tx, Cy - sign*(H - Ry), 0.0, lcMax)

        # Lines
        l_arc = geo.addEllipseArc(p_e_l, p_e_c, p_e_major, p_e_r)
        l_r    = geo.addLine(p_e_r, p_r)
        l_x    = geo.addLine(p_r, p_l)
        l_middle = geo.addLine(p_l, p_e_l)
        l_loop   = geo.addCurveLoop([l_arc, l_r, l_x, l_middle])
        # Face
        face = geo.addPlaneSurface([l_loop])
        # }}}
        # For full geometry {{{
        if fullGeo:
            raise ValueError("The code is not ready to use fullGeo.")
            # Points
            mp_lc = geo.addPoint(-Tx,  Ty, 0.0, R)
            # Lines
            ml_circle = geo.addCircleArc(mp_lc, p_cc, p_cb)
            ml_top    = geo.addLine(p_ct, mp_lc)
            ml_loop   = geo.addCurveLoop([ml_circle, -l_middle, ml_top])
            # Face
            mface = geo.addPlaneSurface([ml_loop])
        # }}}
        # Create auxiliary arc where the mesh will be finer {{{
        aux_p_1 = geo.addPoint(p_aux_1[0], p_aux_1[1], 0.0, maxLength)
        aux_p_2 = geo.addPoint(p_aux_2[0], p_aux_2[1], 0.0, maxLength)
        aux_c1 = geo.addEllipseArc(p_e_l, p_e_c, p_e_major, aux_p_1)
        aux_c2 = geo.addEllipseArc(p_e_l, p_e_c, p_e_major, aux_p_2)
        self.set_aux_lines([aux_c1, aux_c2])
        # }}}
        # Synchronize
        geo.synchronize()
        # Make physical groups {{{
        # Points
        self.PhysicalGroup(0, [p_e_l], name + name_ref_point)
        # Lines
        lst_dxM = [l_middle]
        if fullGeo:
            raise ValueError("The code is not ready to use fullGeo.")
        else:
            lst_main = [l_arc]
            lst_dy    = [l_x]
            lst_face   = [face]
            lst_dxR = [l_r]
        self.PhysicalGroup(1, lst_main, name + name_main_line)
        self.PhysicalGroup(1, lst_dy   , name + name_dy )
        self.PhysicalGroup(1, lst_dxM, name + name_dxM )
        self.PhysicalGroup(1, lst_dxR, name + name_dxR )
        self._mainLine = lst_main
        # Faces
        self.PhysicalGroup(2, lst_face, name)
        # }}}
        return
    # }}}
    # Remesh {{{
    def Remesh(self, prevMeshFile, **kwargs):
        # Object parameters {{{
        name_ref_point = self.name_ref_point
        name_main_line = self.name_main_line
        name_dxM = self.name_dxM
        name_dxR = self.name_dxR
        name_dy    = self.name_dy
        H       = self.H
        Rx      = self.Rx
        Ry      = self.Ry
        Cx      = self.Cx
        Cy      = self.Cy
        p_aux_1 = self.p_aux_1
        p_aux_2 = self.p_aux_2
        fullGeo = self.fullGeo
        name = self.name
        cutAt_y = self.cutAt_y
        x_func = self.x_func
        y_func = self.y_func
        maxLength = self.maxLength
        sign = self.sign
        # }}}
        # Get kwargs {{{
        x_aux_1 = kwargs.get("x_aux_1", p_aux_1[0])
        x_aux_2 = kwargs.get("x_aux_2", p_aux_2[0])
        # }}}
        # GMSH Initialisation {{{
        gmsh.finalize()
        gmsh.initialize()
        gmsh.clear()
        gmsh.model.add(name)
        model = gmsh.model
        geo = model.occ
        mesh = model.mesh
        # }}}
        # Get information from previous mesh {{{
        prevMesh = AsterMesh(prevMeshFile)
        refPointNames = prevMesh.nodeSets[name + name_ref_point]
        refLines = self.GetListMainLines(prevMeshFile)
        refLineNodes = []
        for line in refLines:
            linePoints = prevMesh.GetNodeSetCoordinates(nodeSet = line)
            refLineNodes.append(linePoints)
        # }}}
        # Geometry reconstruction {{{
        # Cut node list if necessary {{{
        if cutAt_y:
            # Cut list of points in the limits provided by Ty
            cutRefLines = []
            for line in refLines:
                cutLine = []
                for nodeName in line:
                    node = prevMesh.GetNodeCoordinates(nodeName)
                    yi = node[1]
                    # Create node if the y coordinate is below cutAt_y
                    if abs(yi) < abs(cutAt_y):
                        cutLine.append(nodeName)
                cutRefLines.append(cutLine)
            cutRefLineNodes = []
            for line in cutRefLines:
                cutRefLineNodes.append(prevMesh.GetNodeSetCoordinates(nodeSet = line))
            refLineNodes = cutRefLineNodes
            # }}}
        # Find the correct side and flip if necessary {{{
        for k1 in range(len(refLineNodes)):
            p0 = refLineNodes[k1][0]
            pf = refLineNodes[k1][-1]
            if abs(Cx - p0[0]) > abs(Cx - pf[0]):
                refLineNodes[k1] = np.flip(refLineNodes[k1], 0)
        # }}}
        # Make main line {{{
        if len(refLineNodes) > 1:
            message = "Too many main lines to reconstruct ({}).\n".format(len(refLineNodes))
            message += "Currently only one line can be used.\n"
            raise ValueError(message)
        refLineNodes = refLineNodes[0]
        firstLineIds = []
        for coords in refLineNodes:
            xi = coords[0]
            yi = coords[1]
            firstLineIds.append(geo.addPoint(xi, yi, 0.0, maxLength))
        l_main1 = geo.addBSpline(firstLineIds)
        # }}}
        # Make other necessary points and lines and construct the face {{{
        # Get coordinates and ids from spline
        p_e_l = firstLineIds[0]
        p_e_r = firstLineIds[-1]
        p_e_l_coords = refLineNodes[0]
        p_e_r_coords = refLineNodes[-1]
        # Make  points
        xi = p_e_l_coords[0]
        yi = Cy - sign*(H - Ry)
        p_l = geo.addPoint(xi, yi, 0.0, maxLength)
        xi = p_e_r_coords[0]
        p_r = geo.addPoint(xi, yi, 0.0, maxLength)
        # Make lines
        l_r    = geo.addLine(p_e_r, p_r)
        l_x    = geo.addLine(p_r, p_l)
        l_middle = geo.addLine(p_l, p_e_l)
        l_loop   = geo.addCurveLoop([l_main1, l_r, l_x, l_middle])
        # Face
        face = geo.addPlaneSurface([l_loop])
        # }}}
        # }}}
        # Make auxiliary spline {{{
        auxPointIds = []
        for node in refLineNodes:
            xi = node[0]
            yi = node[1]
            # Create node if x_aux_1 < xi < x_aux_2
            if (x_aux_1 < xi) and (xi < x_aux_2):
                # Create node if x_aux_1 < yi < x_aux_2
                if (x_aux_1 < yi) and (yi < x_aux_2):
                    auxPointIds.append(geo.addPoint(xi, yi, 0.0, maxLength))
        aux_c = geo.addBSpline(auxPointIds)
        self.set_aux_lines([aux_c])
        # }}}
        # Synchronize
        geo.synchronize()
        # Make physical groups {{{
        # Points
        self.PhysicalGroup(0, [p_e_l], name + name_ref_point)
        # Lines
        lst_dxM = [l_middle]
        if fullGeo:
            raise ValueError("The code is not ready to use fullGeo.")
        else:
            lst_main = [l_main1]
            lst_dy    = [l_x]
            lst_face   = [face]
            lst_dxR = [l_r]
        self.PhysicalGroup(1, lst_main, name + name_main_line)
        self.PhysicalGroup(1, lst_dy   , name + name_dy )
        self.PhysicalGroup(1, lst_dxM, name + name_dxM )
        self.PhysicalGroup(1, lst_dxR, name + name_dxR )
        self._mainLine = lst_main
        # Faces
        self.PhysicalGroup(2, lst_face, name)
        # }}}
        return
    # }}}
    # Update p_aux {{{
    def Update_p_aux(self, **kwargs):
        # Get parameters
        x_ellipse = self.x_func
        y_ellipse = self.y_func
        Rx = self.Rx
        Ry = self.Ry
        Cx = self.Cx
        Cy = self.Cy
        # Update p_aux
        if "fine_x" in kwargs:
            fine_x = kwargs["fine_x"]
            if abs(fine_x) < Rx:
                p_aux_1 = [-fine_x, y_ellipse(fine_x)]
                p_aux_2 = [ fine_x, y_ellipse(fine_x)]
            else:
                p_aux_1 = [Cx - Rx, Cy]
                p_aux_2 = [Cx + Rx, Cy]
        elif ("fine_x1" in kwargs) and ("fine_x2" in kwargs):
            fine_x1 = kwargs["fine_x1"]
            if abs(fine_x1) < Rx:
                p_aux_1 = [fine_x1, y_ellipse(fine_x1)]
            else:
                p_aux_1 = [Cx - Rx, Cy]
            fine_x2 = kwargs["fine_x2"]
            if abs(fine_x2) < Rx:
                p_aux_2 = [fine_x2, y_ellipse(fine_x2)]
            else:
                p_aux_2 = [Cx + Rx, Cy]
        elif ("p_aux_1" in kwargs) and ("p_aux_2" in kwargs):
            p_aux_1 = kwargs["p_aux_1"]
            p_aux_2 = kwargs["p_aux_2"]
        else:
            message = "(fine_x) or (fine_x1, fine_x2) or (p_aux_1, p_aux_2) is not defined."
            raise ValueError(message)
        self._p_aux_1 = p_aux_1
        self._p_aux_2 = p_aux_2
        return
    # }}}
# }}}

# Sub class of MeshMakerSizeGradient: Rectangle {{{
class Rectangle(MeshMakerSizeGradient):
    # Properties {{{
    @property
    def W(self):
        return self._W
    @property
    def H(self):
        return self._H
    @property
    def p_aux_1(self):
        return self._p_aux_1
    @property
    def p_aux_2(self):
        return self._p_aux_2
    @property
    def name_middle(self):
        return self._name_middle
    @property
    def name_bottom(self):
        return self._name_bottom
    @property
    def name_side(self):
        return self._name_side
    # }}}
    # __init__ {{{
    def __init__(self, W, H, lcMin, lcMax, **kwargs):
        self._W = W
        self._H = H
        self._maxLength = max(W, H)
        super().__init__(lcMin, lcMax, **kwargs)
        if "fine_x" in kwargs:
            fine_x = kwargs["fine_x"]
            p_aux_1 = [ fine_x, 0.0]
            p_aux_2 = [-fine_x, 0.0]
        elif ("fine_x1" in kwargs) and ("fine_x2" in kwargs):
            fine_x = kwargs["fine_x"]
            fine_x1 = kwargs["fine_x1"]
            fine_x2 = kwargs["fine_x2"]
            p_aux_1 = [fine_x1, 0.0]
            p_aux_2 = [fine_x2, 0.0]
        elif ("p_aux_1" in kwargs) and ("p_aux_2" in kwargs):
            p_aux_1 = kwargs["p_aux_1"]
            p_aux_2 = kwargs["p_aux_2"]
        else:
            raise("Error: (fine_x) or (fine_x1, fine_x2) or (p_aux_1, p_aux_2) is not defined.")
        self._p_aux_1 = p_aux_1
        self._p_aux_2 = p_aux_2
        # Set group names from kwargs {{{
        name_ref_point = kwargs.get("name_ref_point", "_pct")
        if len(name_ref_point) > 4:
            raise("The length of name_ref_point must be less or equal to 4.")
        name_bottom = kwargs.get("name_bottom", "_bot")
        if len(name_bottom) > 4:
            raise("The length of name_bottom must be less or equal to 4.")
        name_side = kwargs.get("name_side", "_sid")
        if len(name_side) > 4:
            raise("The length of name_side must be less or equal to 4.")
        name_main_line = kwargs.get("name_main_line", "_top")
        if len(name_main_line) > 4:
            raise("The length of name_main_line must be less or equal to 4.")
        name_middle = kwargs.get("name_middle", "_mid")
        if len(name_middle) > 4:
            raise("The length of name_middle must be less or equal to 4.")
        self._name_middle = name_middle
        self._name_main_line    = name_main_line
        self._name_side   = name_side
        self._name_bottom = name_bottom
        self._name_ref_point = name_ref_point
        self._list_name_points = [name_ref_point]
        # }}}
        return
    # }}}
    # Make geometry {{{
    def MakeGeometry(self, **kwargs):
        # Set group names from kwargs {{{
        name_middle = self.name_middle
        name_main_line    = self.name_main_line
        name_side   = self.name_side
        name_bottom = self.name_bottom
        name_ref_point = self.name_ref_point
        name_ref_point = kwargs.get("name_ref_point", name_ref_point)
        if len(name_ref_point) > 4:
            raise("The length of name_ref_point must be less or equal to 4.")
        name_bottom = kwargs.get("name_bottom", name_bottom)
        if len(name_bottom) > 4:
            raise("The length of name_bottom must be less or equal to 4.")
        name_side = kwargs.get("name_side", name_side)
        if len(name_side) > 4:
            raise("The length of name_side must be less or equal to 4.")
        name_main_line = kwargs.get("name_main_line", name_main_line)
        if len(name_main_line) > 4:
            raise("The length of name_main_line must be less or equal to 4.")
        name_middle = kwargs.get("name_middle", name_middle)
        if len(name_middle) > 4:
            raise("The length of name_middle must be less or equal to 4.")
        self._name_middle = name_middle
        self._name_main_line    = name_main_line
        self._name_side   = name_side
        self._name_bottom = name_bottom
        self._name_ref_point = name_ref_point
        # }}}
        # Parameters {{{
        W       = self.W
        H       = self.H
        p_aux_1 = self.p_aux_1
        p_aux_2 = self.p_aux_2
        fullGeo = self.fullGeo
        name = self.name
        cutAt_y = self.cutAt_y
        # }}}
        # GMSH Initialisation {{{
        gmsh.initialize()
        gmsh.clear()
        gmsh.model.add(name)
        model = gmsh.model
        geo = model.occ
        mesh = model.mesh
        # }}}
        # Make base face {{{
        R = max(H, W)
        if cutAt_y:
            H = min(cutAt_y, H)
        else:
            H = H
        # Points
        p_cb = geo.addPoint(0.0,  -H, 0.0, R)
        p_rb = geo.addPoint(  W,  -H, 0.0, R)
        p_rt = geo.addPoint(  W, 0.0, 0.0, R)
        p_ct = geo.addPoint(0.0, 0.0, 0.0, R)
        # Lines
        l_bottom = geo.addLine(p_cb, p_rb)
        l_right  = geo.addLine(p_rb, p_rt)
        l_top    = geo.addLine(p_rt, p_ct)
        l_middle = geo.addLine(p_ct, p_cb)
        l_loop   = geo.addCurveLoop([l_bottom, l_right, l_top, l_middle])
        # Face
        face = geo.addPlaneSurface([l_loop])
        # }}}
        # For full geometry {{{
        if fullGeo:
            # Points
            mp_lb = geo.addPoint(-W,  -H, 0.0, R)
            mp_lt = geo.addPoint(-W, 0.0, 0.0, R)
            # Lines
            ml_top    = geo.addLine( p_ct, mp_lt)
            ml_left   = geo.addLine(mp_lt, mp_lb)
            ml_bottom = geo.addLine(mp_lb,  p_cb)
            ml_loop   = geo.addCurveLoop([ml_top, ml_left, ml_bottom, - l_middle])
            # Face
            mface = geo.addPlaneSurface([ml_loop])
        # }}}
        # Create auxiliary arc where the mesh will be finer {{{
        aux_p_1 = geo.addPoint(p_aux_1[0], p_aux_1[1], 0.0, R)
        aux_p_2 = geo.addPoint(p_aux_2[0], p_aux_2[1], 0.0, R)
        aux_c = geo.addLine(aux_p_1, aux_p_2)
        self.set_aux_lines([aux_c])
        # }}}
        # Synchronize
        geo.synchronize()
        # Make physical groups {{{
        # Points
        self.PhysicalGroup(0, [p_ct], name + name_ref_point)
        # Lines
        lst_middle = [l_middle]
        if fullGeo:
            lst_bottom = [l_bottom, ml_bottom]
            lst_main    = [l_top   , ml_top   ]
            lst_side   = [l_right , ml_left  ]
            lst_face   = [face    , mface    ]
            self.PhysicalGroup(1, [l_right]  , name + name_side + "R" )
            self.PhysicalGroup(1, [ml_left]  , name + name_side + "L" )
        else:
            lst_bottom = [l_bottom]
            lst_main    = [l_top   ]
            lst_side   = [l_right ]
            lst_face   = [face    ]
            self.PhysicalGroup(1, lst_side  , name + name_side  )
        self.PhysicalGroup(1, lst_bottom, name + name_bottom)
        self.PhysicalGroup(1, lst_main   , name + name_main_line   )
        self.PhysicalGroup(1, lst_middle, name + name_middle)
        # Faces
        self.PhysicalGroup(2, lst_face  , name  )
        # }}}
        return
    # }}}
    # Remesh {{{
    def Remesh(self, prevMeshFile, **kwargs):
        # Object parameters {{{
        name_middle = self.name_middle
        name_main_line    = self.name_main_line
        name_side   = self.name_side
        name_bottom = self.name_bottom
        name_ref_point = self.name_ref_point
        lcMin       = self.lcMin
        lcMax       = self.lcMax
        name        = self.name
        fullGeo     = self.fullGeo
        cutAt_y     = self.cutAt_y
        W           = self.W
        H           = self.H
        p_aux_1     = self.p_aux_1
        p_aux_2     = self.p_aux_2
        maxLength   = self.maxLength
        # }}}
        # Get kwargs {{{
        x_aux_1 = kwargs.get("x_aux_1", p_aux_1[0])
        x_aux_2 = kwargs.get("x_aux_2", p_aux_2[0])
        # }}}
        # GMSH Initialisation {{{
        gmsh.finalize()
        gmsh.initialize()
        gmsh.clear()
        gmsh.model.add(name)
        model = gmsh.model
        geo = model.occ
        mesh = model.mesh
        # }}}
        # Get information from previous mesh {{{
        prevMesh = AsterMesh(prevMeshFile)
        refPointNames = prevMesh.nodeSets[name + name_ref_point]
        refLines = self.GetListMainLines(prevMeshFile)
        refLineNodes = []
        for line in refLines:
            refLineNodes.append(prevMesh.GetNodeSetCoordinates(nodeSet = line))
        # }}}
        # Geometry reconstruction {{{
        if cutAt_y:
            # Recover of _top points {{{
            # Make reference point
            xi, yi = prevMesh.GetNodeCoordinates(refPointNames[0])
            p_ct = geo.addPoint(xi, yi, 0.0, maxLength)
            # Make first line points
            firstLineIds = []
            for coords in refLineNodes[0]:
                xi = coords[0]
                yi = coords[1]
                if (xi < 0.99*W) and (xi > -0.99*W):
                    firstLineIds.append(geo.addPoint(xi, yi, 0.0, maxLength))
            # Complete line
            p0 = refLineNodes[0][0]
            pf = refLineNodes[0][-1]
            if pf[0] - p0[0] > 0:
                # The line is set from left to right
                # Add point to the right
                firstLineIds.append(geo.addPoint(W, pf[1], 0.0, maxLength))
            else:
                # The line is set from right to left
                # Add point to the right
                firstLineIds.insert(0, geo.addPoint(W, p0[1], 0.0, maxLength))
            l_main1 = geo.addBSpline(firstLineIds)
            if fullGeo:
                secondLineIds = [firstLineIds[-1]]
                for coords in refLineNodes[1][1 :]:
                    xi = coords[0]
                    yi = coords[1]
                    if (xi < 0.99*W) and (xi > -0.99*W):
                        secondLineIds.append(geo.addPoint(xi, yi, 0.0, maxLength))
                # Complete line
                p0 = refLineNodes[1][0]
                pf = refLineNodes[1][-1]
                if pf[0] - p0[0] > 0:
                    # The line is set from left to right
                    # Add point to the left
                    secondLineIds.insert(0, geo.addPoint(-W, p0[1], 0.0, maxLength))
                else:
                    # The line is set from right to left
                    # Add point to the left
                    secondLineIds.append(geo.addPoint(-W, pf[1], 0.0, maxLength))
                l_main2 = geo.addBSpline(secondLineIds)
                l_main_list = [l_main1, l_main2]
            else:
                l_top = l_main1
                l_main_list = [l_top]
            # }}}
            # Make the other lines and face {{{
            points = prevMesh.GetNodeSetCoordinates(setName = name
                    + name_main_line)
            if fullGeo:
                if points[0][0] < points[-1][0]:
                    p_lt = firstLineIds[ 0]
                    p_rt = secondLineIds[-1]
                    sign = -1
                else:
                    p_lt = secondLineIds[-1]
                    p_rt = firstLineIds[ 0]
                    sign = 1
                p_rb = geo.addPoint(  W, -H, 0.0, maxLength)
                p_lb = geo.addPoint( -W, -H, 0.0, maxLength)
                l_left = geo.addLine(p_lt, p_lb)
                l_bottom = geo.addLine(p_lb, p_rb)
                l_right = geo.addLine(p_rb, p_rt)
                l_loop = geo.addCurveLoop([l_left, l_bottom, l_right, sign*l_main1, sign*l_main2])
            else:
                p_rb = geo.addPoint(0.0, -cutAt_y, 0.0, maxLength)
                p_lb = geo.addPoint(W, -cutAt_y, 0.0, maxLength)
                p_rt = firstLineIds[-1]
                p_lt = firstLineIds[0]
                l_middle    = geo.addLine(p_rt, p_rb)
                l_bottom = geo.addLine(p_rb, p_lb)
                l_left = geo.addLine(p_lb, p_lt)
                l_loop = geo.addCurveLoop([l_top, l_middle, l_bottom,
                    l_left])
            face = geo.addPlaneSurface([l_loop])
            # }}}
        else:
            model = rgmsh.Remesh([gr_names], [file_names], [lst_lcMax],
                    [name], model, geoType = 'occ')
        # }}}
        # Make auxiliary spline {{{
        auxPointIds = []
        for node in points:
            xi = node[0]
            yi = node[1]
            # Create node if the y coordinate is below cutAt_y
            if (x_aux_1 < xi) and (xi < x_aux_2):
                auxPointIds.append(geo.addPoint(xi, yi, 0.0, maxLength))
        aux_c = geo.addBSpline(auxPointIds)
        self.set_aux_lines([aux_c])
        # }}}
        # Synchronize
        geo.synchronize()
        # Make physical groups {{{
        # Points
        self.PhysicalGroup(0, [p_ct], name + name_ref_point)
        # Lines
        if fullGeo:
            lst_main    = [l_main1, l_main2]
            self.PhysicalGroup(1, [l_right]  , name + name_side + "R" )
            self.PhysicalGroup(1, [l_left]  , name + name_side + "L" )
        else:
            lst_side   = [l_left]
            lst_main    = [l_top]
            lst_middle = [l_middle]
            self.PhysicalGroup(1, lst_middle, name + name_middle)
            self.PhysicalGroup(1, lst_side  , name + name_side  )
        lst_bottom = [l_bottom]
        lst_face   = [face]
        self.PhysicalGroup(1, lst_bottom, name + name_bottom)
        self.PhysicalGroup(1, lst_main   , name + name_main_line   )
        # Faces
        self.PhysicalGroup(2, lst_face, name)
        # }}}
        return
    # }}}
# }}}

# Class to read ASTER mesh {{{
class AsterMesh(object):
    # Properties {{{
    @property
    def coorType(self):
        return self._coorType
    @property
    def filelines(self):
        return self._filelines
    @property
    def nodes(self):
        return self._nodes
    @property
    def elements(self):
        return self._elements
    @property
    def nodeSets(self):
        return self._nodeSets
    @property
    def elementSets(self):
        return self._elementSets
    @property
    def neighbourNodes(self):
        return self._neighbourNodes
    # }}}
    # __init__ {{{
    def __init__(self, filename):
        # Read file
        with open(filename, 'r') as fle:
            filelines = fle.readlines()
            filelines = [line[:-1] for line in filelines]
        self._filelines = filelines
        # Read coordinates
        if len(self.FindIndices('COOR_3D')) == 0:
            self._coorType = 'COOR_2D'
        else:
            self._coorType = 'COOR_3D'
        self.read_coor()
        # Read elements
        self._elements = {}
        # Read elements: SEG2
        self.read_elements('SEG2')
        # Read elements: TRIA3
        self.read_elements('TRIA3')
        # Read elements: QUAD4
        self.read_elements('QUAD4')
        # Read groups: node groups
        self._nodeSets = self.read_groups('GROUP_NO')
        # Read groups: element groups
        self._elementSets = self.read_groups('GROUP_MA')
        return
    # }}}
    # read_coor {{{
    def read_coor(self):
        # Find the location of the nodes
        filelines = self.filelines
        coorType = self.coorType
        indices = self.FindIndices(coorType)
        if len(indices) == 0:
            print("Error: there is no " + coorType + ".")
            raise
        elif len(indices) > 1:
            raise("Error: there are too many " + coorType + " lines.")
        elif len(indices) == 1:
            # Load node names and coordinates
            lineNum = indices[0]
            read = True
            nodes = {}
            neighbourNodes = {}
            while read:
                lineNum += 1
                line = filelines[lineNum]
                if 'FINSF' in line:
                    read = False
                else:
                    line = line.split()
                    nName = line[0]
                    neighbourNodes[nName] = set()
                    if self.coorType == 'COOR_2D':
                        nodes[nName] = np.array([eval(line[1]), eval(line[2])])
                    else:
                        nodes[nName] = np.array([eval(line[1]), eval(line[2]), eval(line[3])])
        self._nodes = nodes
        self._neighbourNodes = neighbourNodes
        return
    # }}}
    # read_elements {{{
    def read_elements(self, etype):
        # Find element type in file lines
        filelines = self.filelines
        indices = self.FindIndices(etype)
        eles = {}
        if len(indices) > 0:
            for index in indices:
                lineNum = index
                read = True
                while read:
                    lineNum += 1
                    line = filelines[lineNum]
                    if 'FINSF' in line:
                        read = False
                    else:
                        line = line.split()
                        #eles.append(line)
                        eles[line[0]] = line[1:]
                        eleNods = line[1:]
                        nnod = len(eleNods)
                        for k1 in range(nnod):
                            nodeName = eleNods[k1]
                            for k2 in range(nnod):
                                if k1 != k2:
                                    neighName = eleNods[k2]
                                    self._neighbourNodes[nodeName].add(neighName)
        self._elements[etype] = eles
        return
    # }}}
    # read_groups {{{
    def read_groups(self, gtype):
        # Find group type in file lines
        filelines = self.filelines
        indices = self.FindIndices(gtype)
        dic = {}
        if len(indices) > 0:
            for index in indices:
                lineNum = index
                read = True
                ginfo = []
                while read:
                    lineNum += 1
                    line = filelines[lineNum]
                    if 'FINSF' in line:
                        read = False
                    else:
                        ginfo.append(line)
                gname = ginfo[0].replace(" ", "")
                gobjects = []
                for k1 in range(1, len(ginfo)):
                    line = ginfo[k1].split()
                    for obj in line:
                        gobjects.append(obj)
                dic[gname] = gobjects
        return dic
    # }}}
    # Find string indices {{{
    def FindIndices(self, string):
        filelines = self.filelines
        indices = []
        for k1 in range(len(filelines)):
            line = filelines[k1].replace(' ', '')
            if line == string:
                indices.append(k1)
        return indices
    # }}}
    # Get node coordinates {{{
    def GetNodeCoordinates(self, nodeName):
        nodes = self.nodes
        coords = nodes[nodeName]
        return coords
    # }}}
    # Get node set coordinates {{{
    def GetNodeSetCoordinates(self, setName = None, nodeSet = None):
        if setName and nodeSet:
            raise("Error: you should indicate either setName or nodeSet, not both.")
        if not(setName or nodeSet):
            raise("Error: you should indicate either setName or nodeSet.")
        if setName:
            group = self.nodeSets[setName]
        else:
            group = nodeSet
        coords = []
        for nodeName in group:
            coords.append(self.GetNodeCoordinates(nodeName))
        return np.array(coords)
    # }}}
    # Write mesh {{{
    def Write(self, filename):
        # Update filelines
        self._Update_filelines()
        # Write filelines
        with open(filename, 'w') as fle:
            for line in self.filelines:
                fle.write(line)
        return
    # }}}
    # Add element based on connectivities {{{
    def add_element(self, etype, connectivities):
        # Check input {{{
        if etype == 'SEG2':
            if not len(connectivities) == 2:
                raise("Error: the length of connectivities should be 2 for 'SEG2' elements.")
            nnod = 2
        elif etype == 'TRIA3':
            if not len(connectivities) == 3:
                raise("Error: the length of connectivities should be 3 for 'TRIA3' elements.")
            nnod = 3
        elif etype == 'QUAD4':
            if not len(connectivities) == 4:
                raise("Error: the length of connectivities should be 4 for 'QUAD4' elements.")
            nnod = 4
        else:
            raise("Error: etype must belong to ('SEG2', 'TRIA3', 'QUAD4'). Insted '" + etype + "' was provided.")
        # }}}
        # Check that the indicated nodes exist {{{
        for nodeName in connectivities:
            if not isinstance(nodeName, str):
                raise("Error: the elements of 'connectivities' must be strings, the name of the nodes.")
            if not nodeName in self.nodes:
                print("**************************************************")
                print("**************************************************")
                print("Error: the indicated node(" + nodeName + ") has not been created yet.")
                print("**************************************************")
                print("**************************************************")
                raise("Error: the indicated node(" + nodeName + ") has not been created yet.")
        # }}}
        # Find the current number of elements {{{
        numEles = 0
        for key in self.elements:
            numEles += len(self.elements[key])
        # }}}
        # Add new element
        eleId = numEles + 1
        eleId = 'M' + str(eleId)
        ele = [eleId] + connectivities
        #self.elements[etype].append(ele)
        self.elements[etype][eleId] = connectivities
        # Update neighbour nodes
        eleNods = connectivities
        nnod = len(eleNods)
        for k1 in range(nnod):
            nodeName = eleNods[k1]
            for k2 in range(nnod):
                if k1 != k2:
                    neighName = eleNods[k2]
                    self._neighbourNodes[nodeName].add(neighName)
        return eleId
    # }}}
    # Update filelines {{{
    def _Update_filelines(self):
        # Add title {{{
        newLines = [' TITRE\n']
        newLines.append('ASTER 14.06.00 CONCEPT mesh\n')
        newLines.append('MAILLAGE_SDASTER\n')
        newLines.append('FINSF\n')
        newLines.append('%\n')
        # }}}
        # Add coordinates {{{
        newLines.append(' ' + self.coorType + '\n')
        numNods = len(self.nodes)
        for nodeName in self.nodes:
            coord = self.nodes[nodeName]
            node_x = '{:1.14E}'.format(coord[0])
            if coord[0] >= 0.0:
                node_x = ' ' + node_x
            node_y = '{:1.14E}'.format(coord[1])
            if coord[1] >= 0.0:
                node_y = ' ' + node_y
            newLines.append(' {:8s} {:21s} {:21s}'.format(
                nodeName, node_x, node_y))
            if self.coorType == 'COOR_3D':
                node_z = '{:1.14E}'.format(coord[2])
                if coord[2] >= 0.0:
                    node_z = ' ' + node_z
                newLines.append(' {:21s}'.format(node_z))
            newLines.append('\n')
        newLines.append('FINSF\n')
        newLines.append('%\n')
        # }}}
        # Add elements {{{
        for etype in self.elements:
            connectivities = self.elements[etype]
            newLines.append(' ' + etype + '\n')
            for eleId in connectivities:
                line = ' {:8s}'.format(eleId)
                for obj in connectivities[eleId]:
                    line += ' {:8s}'.format(obj)
                line += '\n'
                newLines.append(line)
            newLines.append('FINSF\n')
            newLines.append('%\n')
        # }}}
        # Add node sets (group of nodes) {{{
        dictionary = self.nodeSets
        heading = 'GROUP_NO'
        for key in dictionary:
            array = dictionary[key]
            newLines.append(' ' + heading + '\n')
            newLines.append(key + '\n')
            numObjs = len(array)
            numLines = int((numObjs - 1)/7) + 1
            for k1 in range(numLines):
                numObjInLine = min(numObjs - k1*7, 7)
                line = ''
                for k2 in range(numObjInLine):
                    k3 = (k1*7) + k2
                    obj = array[k3]
                    line += ' {:8s}'.format(obj)
                line += '\n'
                newLines.append(line)
            newLines.append('FINSF\n')
            newLines.append('%\n')
        # }}}
        # Add element sets (group of elements) {{{
        dictionary = self.elementSets
        heading = 'GROUP_MA'
        for key in dictionary:
            array = dictionary[key]
            newLines.append(' ' + heading + '\n')
            newLines.append(key + '\n')
            numObjs = len(array)
            numLines = int((numObjs - 1)/7) + 1
            for k1 in range(numLines):
                numObjInLine = min(numObjs - k1*7, 7)
                line = ''
                for k2 in range(numObjInLine):
                    k3 = (k1*7) + k2
                    obj = array[k3]
                    line += ' {:8s}'.format(obj)
                line += '\n'
                newLines.append(line)
            newLines.append('FINSF\n')
            newLines.append('%\n')
        # }}}
        newLines.append('FIN\n')
        self._filelines = newLines
        return
    # }}}
    # Add elements to connect surfaces {{{
    def Add_elements_to_connect_surfaces(self, masterSetName,
            slaveSetName, groupName, oneElementPerSlaveNode = 1,
            xlim = [], rewrite = False):
        # Check name {{{
        if not rewrite:
            if groupName in self.elementSets:
                raise("Error: the indicated group name already exist.")
        # }}}
        # Get the sets {{{
        nodes = self.nodes
        masterSet = self.nodeSets[masterSetName]
        slaveSet  = self.nodeSets[slaveSetName]
        newSet = []
        # Make master and slave dictionaries
        master = {}
        for nodeName in masterSet:
            master[nodeName] = nodes[nodeName]
        slave = {}
        for nodeName in slaveSet:
            slave[nodeName] = nodes[nodeName]
        # }}}
        # Make elements {{{
        if oneElementPerSlaveNode == 1:
            # Create one element for each slave node {{{
            # Find the closest node on the master surface for each node on
            # the slave surface.
            # Using the KD-Tree is a data structure
            tree = KDTree(np.array(list(master.values())))
            pool = Pool(processes=cpu_count())
            arguments = []
            for slaveName, slaveCoord in slave.items():
                arguments.append((slaveName, slaveCoord, KDTree(np.array(list(master.values())))))
            chunksize = int(len(slave)/cpu_count())
            connectivities = pool.starmap(find_closest_node, arguments, chunksize = chunksize)
            if len(xlim) == 0:
                for connect in connectivities:
                    ele = [connect[0], list(master.keys())[connect[1]]]
                    eleId = self.add_element('SEG2', ele)
                    newSet.append(eleId)
            elif len(xlim) == 2:
                xmin = xlim[0]
                xmax = xlim[1]
                for connect in connectivities:
                    # Add element [slaveName, masterName]
                    node = nodes[connect[0]]
                    if xmin <= node[0] and node[0] <= xmax:
                        ele = [connect[0], list(master.keys())[connect[1]]]
                        eleId = self.add_element('SEG2', ele)
                        newSet.append(eleId)
            else:
                raise("Error: len(xlim) can be either 0 or 2")
            # }}}
        elif oneElementPerSlaveNode == 2:
            # Connect all the slave nodes to all the master nodes {{{
            if len(xlim) == 0:
                for slaveName in slave:
                    for masterName in master:
                        ele = [slaveName, masterName]
                        eleId = self.add_element('SEG2', ele)
                        newSet.append(eleId)
            elif len(xlim) == 2:
                xmin = xlim[0]
                xmax = xlim[1]
                for slaveName in slave:
                    node = nodes[slaveName]
                    if xmin <= node[0] and node[0] <= xmax:
                        for masterName in master:
                            ele = [slaveName, masterName]
                            eleId = self.add_element('SEG2', ele)
                            newSet.append(eleId)
            else:
                raise("Error: len(xlim) can be either 0 or 2")
            # }}}
        # }}}
        # Add element set {{{
        self._elementSets[groupName] = newSet
        # }}}
        return
    # }}}
    # Delete node set {{{
    def DeleteNodeSet(self, setName):
        del self.nodeSets[setName]
        return
    # }}}
    # Add neighbour group of nodes {{{
    def Add_neighbour_node_set(self, refSetName, setName):
        # Information {{{
        """
        Create a node set with the union of the sets of neighbour nodes
        of the nodes of refSetName.
        In addition, this neighbour set has no intersection with
        refSetName.
        """
        # }}}
        # Check that setName does not exist
        if setName in self.nodeSets:
            raise("Error: the given setName already exists.")
        # Get the union of neighbouring nodes and the reference set
        unionSet = set()
        refSet = set(self.nodeSets[refSetName])
        for nodeName in refSet:
            neighSet = self.neighbourNodes[nodeName]
            unionSet = unionSet | neighSet
        # Compute the new set as the difference between unionSet and
        # refSet
        self.nodeSets[setName] = list(unionSet.difference(refSet))
        return
    # }}}
# }}}

# Class mesh maker uniform size {{{
class UnionOfLines(object):
    # Documentation {{{
    """
    - points: a numpy array of points with coordinates x, y (and z).
    - lines: a list of dictionaries; each dictionary refers to a line
        - {"physical_name" : name,
           "type" : straight/spline/...,
           "ordered_points" : list of point Ids.}
    """
    # }}}
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
    def lcMin(self):
        return self._lcMin
    @property
    def order(self):
        return self._order
    @property
    def reference_curve(self):
        return self._reference_curve
    # }}}
    # __init__ {{{
    def __init__(self, points, lines, **kwargs):
        if isinstance(points, list):
            points = np.array(points)
        shape = points.shape
        self._npoints = shape[0]
        self._ndim    = shape[1]
        self._points = points*kwargs.get("scale", 1.0)
        self._nlines = len(lines)
        self._lines = lines
        self._name  = kwargs.get("name", "face")
        self._lcMin = kwargs.get("lcMin", None)
        self._order = kwargs.get("order", 1)
        self._reference_curve = kwargs.get("reference_curve", dict())
        return
    # }}}
    # Make geometry {{{
    def MakeGeometry(self):
        # GMSH Initialisation {{{
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
            if self.lcMin is None:
                geo.addPoint(coord_x, coord_y, coord_z)
            else:
                geo.addPoint(coord_x, coord_y, coord_z, self.lcMin)
        # }}}
        # Make lines {{{
        for line in self.lines:
            ltype = line["type"]
            ordered_points = line["ordered_points"]
            if ltype == "straight":
                #straight lines {{{
                p_1 = ordered_points[0] + 1
                p_2 = ordered_points[1] + 1
                line["id"] = geo.addLine(p_1, p_2)
                # }}}
            elif ltype == "spline":
                # Splines {{{
                list_points = [pi + 1 for pi in ordered_points]
                line["id"] = geo.addSpline(list_points)
                # }}}
            elif ltype == "bspline":
                # B-splines {{{
                list_points = [pi + 1 for pi in ordered_points]
                line["id"] = geo.addBSpline(list_points)
                # }}}
            elif ltype == "bezier":
                # Bézier {{{
                list_points = [pi + 1 for pi in ordered_points]
                line["id"] = geo.addBezier(list_points)
                # }}}
            elif ltype == "circlearc":
                # Circle arc {{{
                list_points = [pi + 1 for pi in ordered_points]
                line["id"] = geo.addCircleArc(*list_points)
                # }}}
            else:
                raise("Error: unavailable type of line.")
        # }}}
        # Make face {{{
        loop = geo.addCurveLoop([k1 + 1 for k1 in range(self.nlines)])
        face = geo.addPlaneSurface([loop])
        # }}}
        # Synchronize
        geo.synchronize()
        # Add Physical groups {{{
        for line in self.lines:
            physical_name = line["physical_name"]
            line_id = line["id"]
            self.PhysicalGroup(1, [line_id], physical_name)
        self.PhysicalGroup(2, [face], self.name)
        # }}}
        return
    # }}}
    # Set mesh size {{{
    def SetMeshSize(self, lcMin):
        mesh = gmsh.model.mesh
        mesh.setSize(gmsh.model.getEntities(0), lcMin)
        bg_field = mesh.field.add("MathEval")
        mesh.field.setString(bg_field, "F", "{}".format(lcMin))
        mesh.field.setAsBackgroundMesh(bg_field)
        return
    # }}}
    # Set reference curve mesh size {{{
    def SetReferenceCurveMeshSize(self, **kwargs):
        # Get curve ids {{{
        # Set line id to line name
        line_name_id = {}
        for line in self.lines:
            line_name_id[line["physical_name"]] = line["id"]
        # List of curve ids
        reference_curve = self.reference_curve
        curveNames = reference_curve["curveList"]
        curveList = []
        for curveName in curveNames:
            curveList.append(line_name_id[curveName])
        # }}}
        field_parameters = self.reference_curve.get("field_parameters", {})
        self.SetMeshField_distance_to_curves(0, curveList, **field_parameters)
        return
    # }}}
    # Make mesh / reapeted !!!{{{
    def MakeMesh(self, meshAlgo = 8, smoothItes = 10, recombine = True,
            **kwargs):
        # GMSH environment parameters {{{
        model = gmsh.model
        geo = model.occ
        mesh = model.mesh
        # }}}
        gmsh.option.setNumber("Mesh.Algorithm", meshAlgo)
        gmsh.option.setNumber("Mesh.Smoothing", smoothItes)
        gmsh.option.setNumber('Mesh.SecondOrderIncomplete', 1)
        mesh.generate(2)
        if recombine:
            mesh.recombine()
        self._order = kwargs.get("order", self.order)
        mesh.setOrder(self.order)
        mesh.removeDuplicateNodes()
        return
    # }}}
    # Popup / repeated !!! {{{
    def Popup(self):
        gmsh.fltk.run()
        return
    # }}}
    # Physical group / repeated !!! {{{
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
    # Export mesh {{{
    def Export(self, fmt = 'unv', folder = ''):
        if folder == '':
            folder = os.getcwd()
        pathname = folder + '/' + self.name + '.' + fmt
        gmsh.write(pathname)
        return
    # }}}
    # Remesh {{{
    def Remesh(self, filename, lineType = "bspline",
            min_angle_deg = 20.0):
        # Read the old mesh in ASTER format {{{
        oldmesh = AsterMesh(filename)
        # }}}
        # Set up new points {{{
        # Get the list of nodes based on the physical names of the
        # original mesh
        nodeList = list()
        for line in self.lines:
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
        self._points = np.array(newNodes)
        shape = self.points.shape
        self._npoints = shape[0]
        # }}}
        # Update the lines with the new nodes {{{
        for line in self._lines:
            # Change straight lines to b-splines
            if line["type"] == "straight":
                line["type"] = lineType
            # Set the new ordered list of points
            ordered_points = list()
            lineNodes = oldmesh.nodeSets[line["physical_name"]]
            for node in lineNodes:
                ordered_points.append(nodeDict[node]["id"])
            # Test/set non self intersection spline {{{
            points = self.points
            new_ordered_points = SplineNoLoops(points, ordered_points)
            new_ordered_points = SplineNoSharpAngles(points,
                    new_ordered_points, min_angle_deg = min_angle_deg)
            # }}}
            line["ordered_points"] = new_ordered_points
        # }}}
        # Make new geometry and mesh {{{
        self.MakeGeometry()
        # }}}
        return newNodes
    # }}}
    # Mesh volume {{{
    def MeshVolume(self, dim):
        gmsh.plugin.setNumber("MeshVolume", "Dimension", dim)
        gmsh.plugin.run("MeshVolume")
        _, _, data = gmsh.view.getListData(0)
        volume = data[0][-1]
        return volume
    # }}}
    # Set mesh size by curvature {{{
    def SetCurvatureMeshSize(self,
            rho = lambda kappa: abs(kappa)**0.5 + 1.0, **kwargs):
        # Compute curvature, monitor function and size function {{{
        s_points = self.GetFullContour()
        kappa = Curvature(s_points[:, 0], s_points[:, 1], closed = True)
        rho_x = rho(kappa)
        size = self.lcMin/rho_x
        # }}}
        # Make curvature dependent mesh sizes {{{
        model = gmsh.model
        mesh = model.mesh
        max_lcMin = kwargs.get("lcMax", self.lcMin)
        fields = []
        for k1 in range(len(size)):
            size_k1 = size[k1]
            x = s_points[k1, 0]
            y = s_points[k1, 1]
            Fi = mesh.field.add("Ball")
            mesh.field.setNumber(Fi, "XCenter", x)
            mesh.field.setNumber(Fi, "YCenter", y)
            mesh.field.setNumber(Fi, "Radius", size_k1)
            mesh.field.setNumber(Fi, "VIn", size_k1)
            mesh.field.setNumber(Fi, "Thickness", size_k1)
            mesh.field.setNumber(Fi, "VOut", max_lcMin)
            fields.append(Fi)
        bg = mesh.field.add("Min")
        mesh.field.setNumbers(bg, "FieldsList", fields)
        mesh.field.setAsBackgroundMesh(bg)
        # }}}
        return
    # }}}
    # Get full contour {{{
    def GetFullContour(self):
        points = self.points
        lines = self.lines
        s_points = []
        for line in lines:
            ordered_points = line["ordered_points"]
            for k1 in ordered_points[:-1]:
                xx = points[k1, 0]
                yy = points[k1, 1]
                s_points.append([xx, yy])
        return np.array(s_points)
    # }}}
    # Set mesh fields : distance to a set of curves {{{
    def SetMeshField_distance_to_curves(self, f0, curveList, **kwargs):
        # Set parameters {{{
        mesh = gmsh.model.mesh
        sizeMin = kwargs.get("sizeMin", self.lcMin)
        sizeMax = kwargs.get("sizeMax", 10.0*sizeMin)
        distMinFactor = kwargs.get("distMinFactor", 2.0)
        distMin = kwargs.get("distMin", sizeMin*distMinFactor)
        sizeGrwRate = kwargs.get("sizeGrwRate", 1.1)
        distGrwRate = kwargs.get("distGrwRate", 1.1)
        numFields = kwargs.get("numFields", 1000)
        growthStyle = kwargs.get("growthStyle", "linear")
        # }}}
        fieldList = MeshFieldDistanceCurve(mesh, f0, curveList, sizeMin,
                distMin, sizeGrwRate, distGrwRate, sizeMax, numFields,
                growthStyle = growthStyle)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        return
    # }}}
# }}}

# Functions {{{
# Mesh field distance--threshold to curve {{{
def MeshFieldDistanceCurve(mesh, f0, curveList, sizeMin, distMin,
        sizeGrwRate, distGrwRate, sizeMax, numFields,
        growthStyle = 'linear', sampling = 10000, at = "CurvesList"):
    # Initialisation {{{
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
    for k1 in range(numFields):
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
            break
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

# Mesh from parameters {{{
def MeshFromParams(meshParams):
    geoType = meshParams["geoType"].lower()
    if geoType in ['circ', 'circle']:
        newMesh = Circle(**meshParams)
    elif geoType in ['rect', 'rectangle']:
        newMesh = Rectangle(**meshParams)
    elif geoType in ['elli', 'ellipse']:
        newMesh = Ellipse(**meshParams)
    else:
        raise("Error: unknown geometry type.")
    try:
        boundary_names = meshParams["boundary_names"]
    except KeyError:
        boundary_names = {}
    newMesh.MakeGeometry(**boundary_names)
    return newMesh
# }}}

# Find closest node {{{
def find_closest_node(name, coord, tree):
    dist, index = tree.query(coord, k = 1)
    return [name, index]
# }}}

# Compute curvature {{{
def Curvature(x, y, **kwargs):
    kappa = SignedCurvature(x, y, **kwargs)
    return np.abs(kappa)
# }}}

# Compute signed curvature {{{
def SignedCurvature(x, y, **kwargs):
    # Make equivalent closed array if necessary {{{
    closed = kwargs.get("closed", False)
    if closed:
        xx = np.zeros(x.size + 4)
        yy = np.zeros(y.size + 4)
        xx[0]  = x[-2]
        yy[0]  = y[-2]
        xx[1]  = x[-1]
        yy[1]  = y[-1]
        xx[-2] = x[0]
        yy[-2] = y[0]
        xx[-1] = x[1]
        yy[-1] = y[1]
        for k1 in range(x.size):
            xx[k1 + 2] = x[k1]
            yy[k1 + 2] = y[k1]
    else:
        xx = x
        yy = y
    # }}}
    dxdt   = np.gradient(xx, edge_order = 2)
    dydt   = np.gradient(yy, edge_order = 2)
    d2xdt2 = np.gradient(dxdt, edge_order = 2)
    d2ydt2 = np.gradient(dydt, edge_order = 2)
    kappa_num = np.abs(dxdt*d2ydt2 - dydt*d2xdt2)
    kappa_sign = np.sign(dxdt*d2ydt2 - dydt*d2xdt2)
    kappa_det = (dxdt**2.0 + dydt**2.0)**(3.0/2.0)
    kappa_det[kappa_det < 1.0e-6] = 1.0e-6
    kappa = kappa_num/kappa_det
    kappa[kappa < 1.0e-6] = 0.0
    kappa = kappa*kappa_sign
    if closed:
        kappa = kappa[2:-2]
    return kappa
# }}}

# Compute normal {{{
def Normal(x, y, **kwargs):
    # Make equivalent closed array if necessary {{{
    closed = kwargs.get("closed", False)
    if closed:
        xx = np.zeros(x.size + 4)
        yy = np.zeros(y.size + 4)
        xx[0]  = x[-2]
        yy[0]  = y[-2]
        xx[1]  = x[-1]
        yy[1]  = y[-1]
        xx[-2] = x[0]
        yy[-2] = y[0]
        xx[-1] = x[1]
        yy[-1] = y[1]
        for k1 in range(x.size):
            xx[k1 + 2] = x[k1]
            yy[k1 + 2] = y[k1]
    else:
        xx = x
        yy = y
    # }}}
    dxdt   = np.gradient(xx, edge_order = 2)
    dydt   = np.gradient(yy, edge_order = 2)
    tan_mag = np.sqrt(dxdt**2.0 + dydt**2.0)
    nx =  dydt/tan_mag
    ny = -dxdt/tan_mag
    normal = np.column_stack((nx, ny))
    return normal
# }}}

# Compute size by Mackenzie monitor function {{{
def Size_Mackenzie(points, **kwargs):
    numPoints, dim = points.shape
    # Density by curvature function
    beta = kwargs.get("beta", 0.5)
    f_kappa = lambda kappa: abs(kappa)**beta
    # Get arc length from 0 to k1 point
    acc_length = [0.0]
    for k1 in range(1, numPoints):
        arc = LineString(points[:k1 + 1])
        acc_length.append(arc.length)
    acc_length = np.array(acc_length)
    totLength = acc_length[-1]
    # Curvature
    kappa = Curvature(points[:, 0], points[:, 1])
    # rho_floor
    integral = integrate.simps(f_kappa(kappa), acc_length)
    rho_floor = integral/totLength
    # Mackenzie monitor function
    rho = []
    for k1 in range(numPoints):
        rho.append(0.5*(rho_floor + f_kappa(kappa[k1])))
    rho = np.array(rho)
    # Accumulated monitor function
    acc_rho = [0.0]
    for k1 in range(1, numPoints):
        l_k1   = acc_length[:k1 + 1]
        rho_k1 = rho[:k1 + 1]
        acc_rho.append(integrate.simps(rho_k1, l_k1))
    acc_rho = np.array(acc_rho)
    # Mesh size
    if "numEles" in kwargs:
        numEles = kwargs["numEles"]
    elif "lcMin" in kwargs:
        lcMin = kwargs["lcMin"]
        numEles = int(totLength/lcMin)
    else:
        raise("Error: either numEles or lcMin must be provided.")
    optLength = acc_rho[-1]/numEles
    nodeSize = []
    for rho_i in rho:
        nodeSize.append(optLength/rho_i)
    nodeSize = np.array(nodeSize)
    return nodeSize
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
# }}}

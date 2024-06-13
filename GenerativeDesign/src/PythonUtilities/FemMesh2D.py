import numpy as np
import math

class FemMesh2D(object):
    """"Class to store and manipulate a finite-element mesh."""

    # Properties {{{
    @property
    def nodes(self):
        return self._nodes
    @property
    def elementSets(self):
        return self._elementSets
    def add_elementSets(self, *eleSets):
        if (len(eleSets) == 0):
            raise("Error: no element data has been provided.")
        for eleSet, setName in eleSets:
            # Check the type of data
            try:
                if not(np.issubdtype(eleSet.dtype, np.integer)):
                    raise("Error: the elements must be numpy" +
                            " arrays of integers.")
            except AttributeError:
                raise("Error: the elements must be numpy arrays.")
            # Get the type of element set and the number of elements
            self._elementSets[setName] = FemMesh2D.Elements(eleSet,
                setName)
        return
        # }}}

    # Initialiser and printer {{{
    def __init__(self, nodes, *eleSets):
        # Load the nodes
        try:
            if not nodes.ndim == 2:
                raise("Error: nodes must be a 2D-numpy array.")
        except AttributeError:
            raise("Error: nodes must be a 2D-numpy array.")
        numNods, dim = nodes.shape
        if not dim == 2:
            raise("Error: the input 'nodes' have less or more than 2 "
                    + "columns. They must be exactly 2 for a 2D mesh.")
        self._nodes = FemMesh2D.Nodes(nodes)
        # Load element sets
        if len(eleSets) == 0:
            raise("Error: no element data has been provided.")
        self._elementSets = {}
        self.add_elementSets(*eleSets)
        return
    def __str__(self):
        output = 'Mesh information\n'
        output += self.nodes.__str__()
        listEleSetNames = list(self.elementSets.keys())
        for eleSetName in listEleSetNames:
            output += self.elementSets[eleSetName].__str__()
        return output
    #}}}

    # Mesh from MAIL_PY {{{
    @staticmethod
    def MeshFromMail_Py(mail_py, *eleSetNames):
        # If no element set is specified then all are taken
        if len(eleSetNames) == 0:
            eleSetNames = mail_py.gma.keys()
        # Get nodes
        nodes = mail_py.cn[:,:2]
        # Get elements
        elementList = mail_py.co
        inv_dic = {value : key for key, value in mail_py.dic.items()}
        # Get element sets
        eleSets = []
        for name in eleSetNames:
            points = []
            lines = []
            triangles = []
            quadrangles = []
            lines_3 = []
            triangles_6 = []
            quadrangles_8 = []
            eleIds = mail_py.gma[name]
            for eleId in eleIds:
                ele = elementList[eleId]
                eleType = inv_dic[mail_py.tm[eleId]]
                if eleType == 'POI1':
                    points.append(ele)
                elif eleType == 'SEG2':
                    lines.append(ele)
                elif eleType == 'TRIA3':
                    triangles.append(ele)
                elif eleType == 'QUAD4':
                    quadrangles.append(ele)
                elif eleType == 'SEG3':
                    lines_3.append(ele)
                elif eleType == 'TRIA6':
                    triangles_6.append(ele)
                elif eleType == 'QUAD8':
                    quadrangles_8.append(ele)
                else:
                    print("Element set name: " + name)
                    print("Element size: " + str(ele.size))
                    print("Error: the type of element is not recognisable. " +
                        "Currently we have: linear edges, linear " +
                        "triangles and linear quadrangles. You can " +
                        " program other sort of elements though.")
                    raise
            if len(points) > 0:
                eleSets.append((np.array(points), name + '_points'))
            if len(lines) > 0:
                eleSets.append((np.array(lines), name + '_lines'))
            if len(triangles) > 0:
                eleSets.append((np.array(triangles), name + '_triangles'))
            if len(quadrangles) > 0:
                eleSets.append((np.array(quadrangles), name + '_quadrangles'))
            if len(lines_3) > 0:
                eleSets.append((np.array(lines_3), name + '_lines_3'))
            if len(triangles_6) > 0:
                eleSets.append((np.array(triangles_6), name + '_triangles_6'))
            if len(quadrangles_8) > 0:
                eleSets.append((np.array(quadrangles_8), name + '_quadrangles_8'))
        mesh = FemMesh2D(nodes, *eleSets)
        return mesh
        # }}}

    # Element set nodes {{{
    def ElementSetNodes(self, key):
        # Test input {{{
        if not isinstance(key, str):
            raise("Error: the key must be a string.")
        if not key in list(self.elementSets.keys()):
            raise("Error: key is not in elementSets.")
        # }}}
        # Get the list of nodes {{{
        elementSet = self.elementSets[key].conncts
        nodIds = []
        for ele in elementSet:
            nodIds += list(ele)
        nodIds = set(nodIds)
        nodIds = list(nodIds)
        # }}}
        # Make array with node ids and coordinates {{{
        numNods = len(nodIds)
        nodes = self.nodes.coords
        nods = [None]*numNods
        for k1 in range(numNods):
            nod = nodIds[k1]
            nods[k1] = [nodes[nod, 0], nodes[nod, 1]]
        nods = np.array(nods)
        # }}}
        return nodIds, nods
    # }}}

    # Distance between two element sets {{{
    def DistanceTwoElementSets(self, maitSetName, esclSetName,
            direction = None):
        # Description {{{
        """
        Compute the distance between to element sets as the minimum
        distance between the nodes.
        If direction is None, the distance is compute by the minimum
        euclidean distance between the two sets.
        If direction is a vector: (i) a straight line from each point if
        esclSetName is defined with the direction of the vector, (ii)
        the closest point to from maitSetName is found for each line,
        thus, a set of point pairs is stablished, (iii) the distance
        of each set of point pairs is computed; (iv) the minimum
        distance is selected.
        When direction is a vector, the return distance includes a sign.
        The sign indicates whether the direction vector and the distance
        vector (a vector from the slave point to the master point) go in
        the same direction.
        """
        # }}}
        # Get the two set of points {{{
        maitEleSet = self.elementSets[maitSetName].conncts
        esclEleSet = self.elementSets[esclSetName].conncts
        maitNodSet = set([])
        esclNodSet = set([])
        for ele in maitEleSet:
            for nod in ele:
                maitNodSet.add(nod)
        for ele in esclEleSet:
            for nod in ele:
                esclNodSet.add(nod)
        maitNodSet = list(maitNodSet)
        esclNodSet = list(esclNodSet)
        # }}}
        # Get node coordinates {{{
        coords = self.nodes.coords
        maitCoords = []
        esclCoords = []
        for nod in maitNodSet:
            maitCoords.append([coords[nod, 0], coords[nod, 1]])
        for nod in esclNodSet:
            esclCoords.append([coords[nod, 0], coords[nod, 1]])
        maitCoords = np.array(maitCoords)
        esclCoords = np.array(esclCoords)
        # }}}
        # Compute distance from each slave point to each master point {{{
        if direction is None:
            from scipy.spatial.distance import cdist as cdist
            distanceMatrix = cdist(esclCoords, maitCoords)
            distance = distanceMatrix.min()
        else:
            from Geoanalysis import StraightLine as sl
            from Geoanalysis import VerticalLine as vl
            # Test vector {{{
            if not isinstance(direction, np.ndarray):
                raise("Error: direction must be a numpy array.")
            if not (direction.shape == (2,)):
                raise("Error: direction must have a (2,) shape.")
            uniDir = direction/np.linalg.norm(direction)
            # }}}
            # Compute slope
            if uniDir[0] == 0.0:
                m = None
                m1 = 0.0
            else:
                m = uniDir[1]/uniDir[0]
                if abs(m) > 0.0:
                    m1 = -1.0/m
                else:
                    raise ValueError("You need to program a vertical line for intersection.")
            # Compute the distance with uniDir for each slave node
            distanceToMait = []
            maitIds = []
            for esclCoord in esclCoords:
                if m == None:
                    line = vl(esclCoord[0])
                else:
                    line = sl(p1 = esclCoord, m =m)
                # Find the closest two points
                p1, _, maitId1 = \
                    line.ClosestPointFromSetOfPoints(maitCoords)
                # Make a line with this points
                auxLine = sl(p1 = p1, m = m1)
                intersec = line.Intersection(auxLine)
                distanceVector = intersec - esclCoord
                signDistance = np.dot(distanceVector, uniDir)
                distanceToMait.append(signDistance)
            # Get the minimum distance and the id of the nodes
            distanceToMait = np.array(distanceToMait)
            distance = distanceToMait.min()
        # }}}
        return distance
    # }}}

    # Compute mesh area {{{
    def ComputeMeshArea(self):
        # Get a list of the two--dimensional elements {{{
        twoDimSets = []
        for key in self.elementSets:
            eleset = self.elementSets[key]
            if eleset.dim == 2:
                twoDimSets.append(key)
        # }}}
        # Compute area as the sum of the area of each element set {{{
        area = 0.0
        coords = self.nodes.coords
        for key in twoDimSets:
            area += self.elementSets[key].compute_area(coords)
        # }}}
        return area
    # }}}

    # Compute mesh quality {{{
    def Quality(self, noSet = []):
        # Initialisation of the quality array
        quality = {}
        minQ = 1.0
        # Get node coordinates
        coords = self.nodes.coords
        # Get each element set
        eleSets = self.elementSets
        meanQ = 0.0
        totNumEles = 0
        for eleSetName in eleSets:
            if eleSetName in noSet:
                continue
            eleSet = eleSets[eleSetName]
            # Compute the quality of the elements of dimension 2 or 3
            if eleSet.dim != 1:
                numEles = eleSet.numEles
                totNumEles += numEles
                # Add mesh-quality report of the element set
                meanQua, minQua, maxQua = eleSet.compute_quality(coords)
                meanQ += meanQua*numEles
                quality[eleSetName] = {"meanQua" : meanQua,
                                       "maxQua" : maxQua,
                                       "minQua" : minQua}
                # Update the minimum mesh quality
                if minQua < minQ:
                    minQ = minQua
        # Set minimum and mean quality
        quality["minQua"] = minQ
        quality["meanQua"] = meanQ/totNumEles
        return quality
    # }}}

    # Class nodes{{{
    class Nodes(object):
        """Class for the nodes of a finite element mesh"""
        # Initialiser and printer {{{
        def __init__(self, coords):
            numNods, dim = coords.shape
            self._coords = coords
            self._numNods = numNods
            self._dim = dim
        def __str__(self):
            output = "Nodes:\n"
            output += f"    Number of nodes: {self.numNods}\n"
            output += f"    Dimensions: {self.dim}\n"
            output += "    Coordinates: \n"
            output += self.coords.__str__()
            output += "\n"
            return output
        # }}}

        # Properties {{{
        @property
        def coords(self):
            return self._coords
        @property
        def numNods(self):
            return self._numNods
        @property
        def dim(self):
            return self._dim
            # }}}
    # }}}

    # Class Elements {{{
    class Elements(object):
        """Class for the elements of finite-element mesh"""
    # Properties {{{
        @property
        def numEles(self):
            return self._numEles
        @property
        def nnod(self):
            return self._nnod
        @property
        def name(self):
            return self._name
        @name.setter
        def name(self, nameStr):
            self._name = nameStr
            return
        @property
        def conncts(self):
            return self._conncts
        @property
        def type_string(self):
            return self._type_string
        @property
        def type_id(self):
            return self._type_id
        @property
        def area(self):
            return self._area
        @property
        def dim(self):
            return self._dim
        # }}}

        # Initialiser and printer {{{
        def __init__(self, conncts, name):
            numEles, nnod = conncts.shape
            self._numEles = numEles
            self._nnod = nnod
            self._name = name
            self._conncts = conncts
            self._area = None
            # Set the type of element
            if nnod == 1:
                self._dim = 0
                self._type_string = "VTK_VERTEX"
                self._type_id = 1
            elif nnod == 2:
                self._dim = 1
                self._type_string = "VTK_LINE"
                self._type_id = 3
            elif nnod == 3:
                if '_lines_3' in name:
                    self._dim = 1
                    self._type_string = "VTK_QUADRATIC_EDGE"
                    self._type_id = 21
                else:
                    self._dim = 2
                    self._type_string = "VTK_TRIANGLE"
                    self._type_id = 5
            elif nnod == 4:
                self._dim = 2
                self._type_string = "VTK_QUAD"
                self._type_id = 9
            elif nnod == 6:
                self._dim = 2
                self._type_string = "VTK_QUADRATIC_TRIANGLE"
                self._type_id = 22
            elif nnod == 8:
                self._dim = 2
                self._type_string = "VTK_QUADRATIC_QUAD"
                self._type_id = 23
            else:
                raise("Error: the type of element is not recognisable. " +
                        "Currently we have: linear edges, linear " +
                        "triangles and linear quadrangles. You can " +
                        " program other sort of elements though.")
            return
        def __str__(self):
            output = "Set of elements:\n"
            output += f"    Name: {self.name}\n"
            output += f"    Type: {self.type_string}\n"
            output += f"          {self.type_id}\n"
            output += f"    Number of elements: {self.numEles}\n"
            output += f"    Number of nodes per element: {self.nnod}\n"
            output += f"    Elements:\n"
            output += self.conncts.__str__()
            output += "\n"
            return output
        # }}}

        # Compute area {{{
        def compute_area(self, coords):
            # Test coords {{{
            if not isinstance(coords, np.ndarray):
                raise("Error: coords must be a numpy array.")
            if not coords.ndim == 2:
                raise("Error: coords.ndim must be 2.")
            if not coords.shape[0] >= self.conncts.max():
                raise("Error: coords.shape[1] must be greater or equal to conncts.max().")
            # }}}
            # Compute area based on dimension {{{
            area = 0.0
            if self.dim == 1:
                func = FemMesh2D.Elements.EdgeArea
            else:
                func = FemMesh2D.Elements.FaceArea
            for ele in self.conncts:
                x = []
                y = []
                for nod in ele:
                    x.append(coords[nod, 0])
                    y.append(coords[nod, 1])
                area += func(np.array(x), np.array(y))
            self._area = area
            # }}}
            return area
        # Edge area {{{
        @staticmethod
        def EdgeArea(x, y):
            # Tests arrays
            if not isinstance(x, np.ndarray):
                raise("Error: x must be a numpy array.")
            if not isinstance(y, np.ndarray):
                raise("Error: y must be a numpy array.")
            if not x.ndim == 1:
                raise("Error: x.ndim must be 1.")
            if not y.ndim == 1:
                raise("Error: y.ndim must be 1.")
            if not x.size == 2:
                raise("Error: x.size must be 2.")
            if not y.size == 2:
                raise("Error: y.size must be 2.")
            return math.sqrt((x[1] - x[0])**2.0 +(y[1] - y[0])**2.0)
        # }}}

        # Face area {{{
        @staticmethod
        def FaceArea(x, y):
            # Tests arrays
            if not isinstance(x, np.ndarray):
                raise("Error: x must be a numpy array.")
            if not isinstance(y, np.ndarray):
                raise("Error: y must be a numpy array.")
            if not x.ndim == 1:
                raise("Error: x.ndim must be 1.")
            if not y.ndim == 1:
                raise("Error: y.ndim must be 1.")
            if not x.size == y.size:
                raise("Error: x.size must be y.size.")
            # Compute area as half the determinant
            area = 0.0
            xx = np.array(list(x) + [x[0]])
            yy = np.array(list(y) + [y[0]])
            for k1 in range(x.size):
                area += xx[k1]*yy[k1 + 1] - yy[k1]*xx[k1 + 1]
            area *= 0.5
            return area
        # }}}
        # }}}

        # Compute quality {{{
        def compute_quality(self, coords, beta = 2.0):
            # Test coords {{{
            if not isinstance(coords, np.ndarray):
                raise("Error: coords must be a numpy array.")
            if not coords.ndim == 2:
                raise("Error: coords.ndim must be 2.")
            if not coords.shape[0] >= self.conncts.max():
                raise("Error: coords.shape[1] must be greater or equal to conncts.max().")
            # }}}
            # Quality based on dimension {{{
            meanQua = 0.0
            maxQua = 0.0
            minQua = 2.0
            for ele in self.conncts:
                x = []
                y = []
                if self.type_id in [1, 3]:
                    continue
                elif self.type_id in [5, 9]: # Linear elements
                    for nod in ele:
                        x.append(coords[nod, 0])
                        y.append(coords[nod, 1])
                elif self.type_id == 22:
                    nodOrder = [0, 3, 1, 4, 2, 5]
                    orederedEle = [ele[ii] for ii in nodOrder]
                    for nod in orederedEle:
                        x.append(coords[nod, 0])
                        y.append(coords[nod, 1])
                elif self.type_id == 23:
                    nodOrder = [0, 4, 1, 5, 2, 6, 3, 7]
                    orederedEle = [ele[ii] for ii in nodOrder]
                    for nod in orederedEle:
                        x.append(coords[nod, 0])
                        y.append(coords[nod, 1])
                else:
                    message = 'Invalid element type_id ' + str(self.type_id)
                    raise ValueError(message)
                eleQuality = FemMesh2D.Elements.FaceQuality(np.array(x),
                        np.array(y), beta)
                meanQua += eleQuality
                if eleQuality > maxQua:
                    maxQua = eleQuality
                if eleQuality < minQua:
                    minQua = eleQuality
            meanQua = meanQua/float(self.numEles)
            # }}}
            return meanQua, minQua, maxQua
        # Quality for linear faces {{{
        @staticmethod
        def FaceQuality(x, y, beta = 2):
            # Tests arrays {{{
            if not isinstance(x, np.ndarray):
                raise("Error: x must be a numpy array.")
            if not isinstance(y, np.ndarray):
                raise("Error: y must be a numpy array.")
            if not x.ndim == 1:
                raise("Error: x.ndim must be 1.")
            if not y.ndim == 1:
                raise("Error: y.ndim must be 1.")
            if not x.size == y.size:
                raise("Error: x.size must be y.size.")
            # }}}
            # Area and perimeter {{{
            area = FemMesh2D.Elements.FaceArea(x, y)
            if area <= 0.0:
                return -1.0
            perimeter = 0.0
            for k1 in range(x.size - 1):
                x1 = x[k1]
                x2 = x[k1 + 1]
                y1 = y[k1]
                y2 = y[k1 + 1]
                perimeter += np.sqrt((x2 - x1)**2.0 + (y2 - y1)**2.0)
            x1 = x[-1]
            x2 = x[0]
            y1 = y[-1]
            y2 = y[0]
            perimeter += np.sqrt((x2 - x1)**2.0 + (y2 - y1)**2.0)
            # }}}
            # Compute perimeter of reference element {{{
            if x.size == 3 or x.size == 6:
                l_r = 2.0*np.sqrt(area/np.sqrt(3.0))
                num_sides = 3
            elif x.size == 4 or x.size == 8:
                l_r = np.sqrt(area)
                num_sides = 4
            else:
                raise("Error: the quality for faces is only available"\
                        + " quadrangles and triangles.")
            p_ref = num_sides*l_r
            quality = 1.0 - abs(1.0 - (p_ref/perimeter)**beta)
            # }}}
            return quality
        # }}}
        # }}}
    # }}}

# Test {{{
if __name__ == '__main__':
    # Define mesh
    nodes = np.array([[0.0, 0.0],
                      [0.5, 0.0],
                      [1.0, 0.0],
                      [0.0, 0.5],
                      [1.0, 0.5],
                      [0.0, 1.0],
                      [0.5, 1.0],
                      [1.0, 1.0]])
    edges = np.array([[0, 1],
                      [1, 2],
                      [2, 3]])
    triangles = np.array([[0, 1, 3],
                          [1, 2, 4],
                          [3, 6, 5],
                          [4, 7, 6]])
    quadrangles = np.array([[1, 4, 6, 3]])
    mesh = FemMesh2D(nodes, (edges, 'edges'), (triangles, 'triangles'),
            (quadrangles, 'quadrangles'))
    print(mesh)
# }}}

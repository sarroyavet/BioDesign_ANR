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
        # Get element sets
        eleSets = []
        for name in eleSetNames:
            points = []
            lines = []
            triangles = []
            quadrangles = []
            eleIds = mail_py.gma[name]
            for eleId in eleIds:
                ele = elementList[eleId]
                if ele.size == 1:
                    points.append(ele)
                elif ele.size == 2:
                    lines.append(ele)
                elif ele.size == 3:
                    triangles.append(ele)
                elif ele.size == 4:
                    quadrangles.append(ele)
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
                self._dim = 2
                self._type_string = "VTK_TRIANGLE"
                self._type_id = 5
            elif nnod == 4:
                self._dim = 2
                self._type_string = "VTK_QUAD"
                self._type_id = 9
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
                for nod in ele:
                    x.append(coords[nod, 0])
                    y.append(coords[nod, 1])
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
            if x.size == 3:
                l_r = 2.0*np.sqrt(area/np.sqrt(3.0))
            elif x.size == 4:
                l_r = np.sqrt(area)
            else:
                raise("Error: the quality for faces is only available"\
                        + " quadrangles and triangles.")
            p_ref = float(x.size)*l_r
            quality = (p_ref/perimeter)**beta
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

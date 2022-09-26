import numpy as np
from FemMesh2D import FemMesh2D as fm2

class FemVtk(object):
    """Class to organise and export finite-element data."""
    # Initialiser{{{
    def __init__(self, mesh):
        if not(isinstance(mesh, fm2)):
            raise("Error: the mesh argument must be an instance of " +
                    "fm2.")
        self._mesh = mesh
        self._data = []
        return
    # }}}

    # Properties {{{
    @property
    def mesh(self):
        return self._mesh
    @property
    def data(self):
        return self._data
        # }}}

    # Add data {{{
    def add_data(self, dataSet, dataName, dataType = 'float'):
        # Check that the number of nodes for the data is the same as
        # for the mesh
        mesh = self.mesh
        if not(self.mesh.nodes.numNods == dataSet.shape[0]):
            raise("Error: the number of nodes of dataSet must be equal " +
                    "to the number of nodes of the mesh.")
        self.data.append(FemVtk.VtkData(dataSet, dataName, dataType))
        return
    # }}}

    # From Voigt to tensor {{{
    @staticmethod
    def FromVoigtToTensor(voigt, model):
        """Convert a Voigt-like array into a symmetric matrix.
        Depending on the model the user should provide:
            - for plane-stress: array([SIXX, SIYY, SIXY], ...]),
            - for plane-strain: array([SIXX, SIYY, SIZZ, SIXY], ...]), and
            - for 3D: array([SIXX, SIYY, SIZZ, SIXY, SIYZ, SIZX], ...])."""
        # Voigt must be a 1D numpy array
        try:
            if not voigt.ndim == 1:
                raise("Error: voigt must be a 1D-numpy array")
        except:
            raise("Error: voigt must be a 1D-numpy array")
        # Create the 3-d array for each case
        numVals = voigt.size
        if model == 'plane-stress':
            numComps = 3
        elif model == 'plane-strain':
            numComps = 4
        elif model == '3D':
            numComps = 6
        else:
            raise("Error: the indicated model is not available. You can use: 'plane-stress', 'plane-strain' or '3D'")
        # Assign the values
        numNods = int(numVals/numComps)
        array = []
        for k1 in range(numNods):
            tensor = np.zeros((3, 3))
            row0 = k1*numComps
            tensor[0, 0] = voigt[row0    ]
            tensor[1, 1] = voigt[row0 + 1]
            if numComps == 3:
                tensor[0, 1] = voigt[row0 + 2]
                tensor[1, 0] = voigt[row0 + 2]
            else:
                tensor[2, 2] = voigt[row0 + 2]
                tensor[0, 1] = voigt[row0 + 3]
                tensor[1, 0] = voigt[row0 + 3]
                if numComps == 6:
                    tensor[0, 2] = voigt[row0 + 4]
                    tensor[2, 0] = voigt[row0 + 4]
                    tensor[1, 2] = voigt[row0 + 5]
                    tensor[2, 1] = voigt[row0 + 5]
            array.append(tensor)
        return np.array(array)
        # }}}

    # Write to vtk functions {{{
    # Write vtk file {{{
    def WriteVtk(self, fileName, title = 'Results', eleSetList = [],
            dataList = []):
        # Get outputs
        output = self.WriteHeader(title)
        output += self.WriteNodes()
        output += self.WriteElements(eleSetList)
        output += self.WriteData(dataList)
        try:
            with open(fileName, 'w') as File:
                File.write(output)
        except OSError:
            raise("Error: fileName must be a string.")
        except PermissionError:
            raise("Error: if fileName starts with / be sure you " +
                    "use the full path")
        return
        # }}}

    # Write header {{{
    def WriteHeader(self, title):
        output  = '# vtk DataFile Version 2.0\n'
        output += title + '\n'
        output += 'ASCII\n'
        output += 'DATASET UNSTRUCTURED_GRID\n'
        return output
        # }}}

    # Write nodes {{{
    def WriteNodes(self):
        coords = self.mesh.nodes.coords
        numNods = self.mesh.nodes.numNods
        output  = 'POINTS ' + str(numNods) + ' float\n'
        for k1 in range(numNods):
            point = (coords[k1, 0], coords[k1, 1], 0.0)
            output += ('%.8e  %.8e  %.8e\n' % point)
        return output
        # }}}

    # Write elements {{{
    def WriteElements(self, eleSetList):
        elementSets = self.mesh.elementSets
        if len(eleSetList) == 0:
            eleSetList = list(elementSets.keys())
        # Get the number cell list size
        size = 0
        numCells = 0
        for k1 in eleSetList:
            eleSet = elementSets[k1]
            eleSetSize = eleSet.numEles*(eleSet.nnod + 1)
            size += eleSetSize
            numCells += eleSet.numEles
        # Write the CELLS block
        output  = 'CELLS %d %d\n' % (numCells, size)
        for k1 in eleSetList:
            eleSet = elementSets[k1]
            nnod = eleSet.nnod
            numEles = eleSet.numEles
            for k2 in range(numEles):
                info = [nnod]
                info += list(eleSet.conncts[k2, :])
                fmt= '%d '*(nnod + 1) + '\n'
                output += fmt % tuple(info)
        # Write the CELL_TYPES block
        output += 'CELL_TYPES %g \n' % (numCells)
        for k1 in eleSetList:
            eleSet = elementSets[k1]
            type_id = eleSet.type_id
            output += ('%d\n' % (type_id))*eleSet.numEles
        return output
        # }}}

    # Write data {{{
    def WriteData(self, dataList):
        data = self.data
        output = ''
        if len(data) == 0:
            return output
        if len(dataList) == 0:
            dataList = range(len(data))
        # Write the data blocks
        output  = 'POINT_DATA %d\n' % (self.mesh.nodes.numNods)
        for k1 in dataList:
            numNods  = data[k1].numNods
            dataSet  = data[k1].dataSet
            dataName = data[k1].dataName
            dataType = data[k1].dataType
            nature   = data[k1].nature
            if nature == 'scalar':
                output += 'SCALARS %s %s 1\n' % (dataName, dataType)
                output += 'LOOKUP_TABLE default\n'
                for k2 in range(numNods):
                    output += '%.8e\n' % dataSet[k2]
            elif nature == 'vector':
                output += 'VECTORS %s %s \n' % (dataName, dataType)
                for k2 in range(numNods):
                    output += ('%.8e %.8e %.8e\n' %
                        tuple(dataSet[k2, :]))
            elif nature == 'tensor':
                output += 'TENSORS %s %s \n' % (dataName, dataType)
                for k2 in range(numNods):
                    for k3 in range(3):
                        output += ('%.8e %.8e %.8e\n' %
                            tuple(dataSet[k2, k3, :]))
        return output
        # }}}
    # }}}

    # Norm of a scalar field {{{
    def LpNormScalar(self, ID, p = 2.0):
        # Get nodal values
        data = self.data[ID]
        # Test nature
        if not data.nature == 'scalar':
            raise("Error: the indicated data (ID) does not have scalar nature.")
        # Initialisation
        norm = 0.0
        eleSetKeys = self.mesh.elementSets.keys()
        coords = self.mesh.nodes.coords
        dataSet = data.dataSet
        for key in eleSetKeys:
            eleSet = self.mesh.elementSets[key]
            numEles = eleSet.numEles
            for ele in eleSet.conncts:
                # Get area and mean value
                x = []
                y = []
                u = 0.0
                for nod in ele:
                    u += dataSet[nod]
                    x.append(coords[nod, 0])
                    y.append(coords[nod, 1])
                u /= ele.size
                area = fm2.Elements.FaceArea(np.array(x), np.array(y))
                # Add to norm
                norm += u**p*area
        norm = norm**(1.0/p)
        return norm
    # }}}

    # VtkData class {{{
    class VtkData(object):
        # Initialiser {{{
        def __init__(self, dataSet, dataName, dataType = 'float'):
            # Test the input
            if (not(isinstance(dataSet, np.ndarray))):
                raise("Error: the dataSet must be a numpy array.")
            if not(isinstance(dataName, str)):
                raise("Error: the dataName must be an string.")
            if not(isinstance(dataType, str)):
                raise("Error: the dataType must be an string.")
            # Infer dataSet nature
            if dataSet.ndim == 1:
                # Scalar
                self._nature = 'scalar'
            elif dataSet.ndim == 2:
                # Vector
                # Test number of columns
                if not dataSet.shape[1] == 3:
                    raise("Error: the vtk format only allows " +
                            "three-component vectors. So dataSet shape " +
                            "must be: (NumNods x 3).")
                self._nature = 'vector'
            elif dataSet.ndim == 3:
                # Tensor
                # Test the second and third dimensions
                if not((dataSet.shape[1] == 3) and
                        (dataSet.shape[2] == 3)):
                    raise("Error: the vtk format only allows " +
                            "three-component square tensors. So " +
                            "dataSet shape must be (NumNods x 3 x 3).")
                self._nature = 'tensor'
            else:
                raise("Error: the vtk format only allows scalar, " +
                        "3-vector or 3x3-tensor data. In addition, " +
                        "the first dimension of data must be the " +
                        "number of nodes.")
            self._dataSet = dataSet
            self._dataName = dataName
            self._dataType = dataType
            self._numNods = dataSet.shape[0]
            return
            # }}}

        # Properties {{{
        @property
        def dataSet(self):
            return self._dataSet
        @property
        def dataName(self):
            return self._dataName
        @property
        def dataType(self):
            return self._dataType
        @property
        def numNods(self):
            return self._numNods
        @property
        def nature(self):
            return self._nature
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
                      [1, 2]])
    triangles = np.array([[0, 1, 3],
                          [1, 2, 4],
                          [3, 6, 5],
                          [4, 7, 6]])
    quadrangles = np.array([[1, 4, 6, 3]])
    scalarData = np.random.rand(8)
    vectorData = np.random.rand(8, 3)
    tensorData = np.random.rand(8, 3, 3)
    mesh = fm2(nodes, (edges, 'edges'), (triangles, 'triangles'),
            (quadrangles, 'quadrangles'))
    vtk = FemVtk(mesh)
    vtk.add_data(scalarData, 's')
    vtk.add_data(vectorData, 'v')
    vtk.add_data(tensorData, 't')
    vtk.WriteVtk("testAll.vtk")
    vtk.WriteVtk("testEleSet.vtk", eleSetList = [1])
    vtk.WriteVtk("testDataSet.vtk", dataList = [1])
# }}}

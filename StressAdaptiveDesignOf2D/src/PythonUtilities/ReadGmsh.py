# Module to read some GMSH files
# Libraries {{{
import gmsh
from scipy.signal import savgol_filter
import numpy as np
import math
# }}}

# GmshOutput class  {{{
# Class to load $View information from a gmsh output file
# $PostFormat
# 1.2 0 8
class GmshOutput(object):
    # Constructor from an input file with a GMSH output {{{
    def __init__(self, fileName):
        #########################################################
        # Load POS ASCII file format (Legacy) (see documentation)
        #########################################################
        # Read file
        with open(fileName, 'r') as inFile:
            lines = inFile.read().splitlines()
        # Test existence and format
        try:
            postFormat = lines.index('$PostFormat')
        except ValueError:
            print('Error: the file: ' + fileName + ' has not the ' +
                    'heading $PostFormat')
            raise
        frmt = lines[postFormat + 1].split()[0]
        if not(eval(frmt) == 1.2):
            print('Error: the file: ' + fileName + ' has not the ' +
                    'format 1.2. You can rewrite this function to be '
                    + 'able to read the new format (' + frmt + ').')
            raise
        ###########
        # Get Views
        ###########
        numViews = lines.count("$View")
        if numViews == 0:
            print('Error: the file ' + fileName + ' has not the ' +
                    'has no "$View" headings.')
        views = [None]*numViews
        lineNum = lines.index('$View')
        for k1 in range(numViews):
            # Get a set of view lines
            viewLines = []
            readOpt = True
            while readOpt:
                lineNum += 1
                newLine = lines[lineNum]
                if newLine == '$EndView':
                    # Find next $View
                    readOpt = False
                    if k1<numViews - 1:
                        lines = lines[lineNum:]
                        lineNum = lines.index('$View')
                else:
                    viewLines.append(newLine)
            # Create View object
            views[k1] = self.View(viewLines)
        self.numViews = numViews
        self.views = views
        return
    # }}}

    # Class View {{{
    # No documentation for format 1.2 has been found, the new
    # documentation has been interpreted!!!!!
    class View(object):
        # ViewError {{{
        # To handle errors
        class ViewError(Exception):
            pass
            # Messages
            def notReady(option):
                print("Error: this option is not yet available. "
                            + "You can add it, though. ")
                print("Option: " + option)
            def notInView(option):
                print("Error: this View does not have any " + option +
                        ".")
            def notKeyWord(option):
                print("Error: this keyword is not available, check" +
                        " it. Keyword: " + option + ".")
        # }}}
        # This part is not complete yet! It is only available for
        # lines
        # Initialiser {{{
        def __init__(self, viewLines):
            # Get heading information
            line = viewLines[0].split()
            view_name = line[0]
            nb_time_steps = eval(line[1])
            self.view_name     = view_name
            self.nb_time_steps = nb_time_steps
            # Get points
            line = viewLines[1].split()
            line = [eval(obj) for obj in line]
            nb_scalar_points, nb_vector_points, nb_tensor_points = line
            self.nb_scalar_points = nb_scalar_points
            self.nb_vector_points = nb_vector_points
            self.nb_tensor_points = nb_tensor_points
            # Get lines
            line = viewLines[2].split()
            line = [eval(obj) for obj in line]
            nb_scalar_lines, nb_vector_lines, nb_tensor_lines = line
            self.nb_scalar_lines = nb_scalar_lines
            self.nb_vector_lines = nb_vector_lines
            self.nb_tensor_lines = nb_tensor_lines
            # Get triangles
            line = viewLines[3].split()
            line = [eval(obj) for obj in line]
            nb_scalar_triangles, nb_vector_triangles, nb_tensor_triangles = line
            self.nb_scalar_triangles = nb_scalar_triangles
            self.nb_vector_triangles = nb_vector_triangles
            self.nb_tensor_triangles = nb_tensor_triangles
            # Get quadrangles
            line = viewLines[4].split()
            line = [eval(obj) for obj in line]
            nb_scalar_quadrangles, nb_vector_quadrangles, nb_tensor_quadrangles = line
            self.nb_scalar_quadrangles = nb_scalar_quadrangles
            self.nb_vector_quadrangles = nb_vector_quadrangles
            self.nb_tensor_quadrangles = nb_tensor_quadrangles
            # Get tetrahedra
            line = viewLines[5].split()
            line = [eval(obj) for obj in line]
            nb_scalar_tetrahedra, nb_vector_tetrahedra, nb_tensor_tetrahedra = line
            self.nb_scalar_tetrahedra = nb_scalar_tetrahedra
            self.nb_vector_tetrahedra = nb_vector_tetrahedra
            self.nb_tensor_tetrahedra = nb_tensor_tetrahedra
            # Get hexahedra
            line = viewLines[6].split()
            line = [eval(obj) for obj in line]
            nb_scalar_hexahedra, nb_vector_hexahedra, nb_tensor_hexahedra = line
            self.nb_scalar_hexahedra = nb_scalar_hexahedra
            self.nb_vector_hexahedra = nb_vector_hexahedra
            self.nb_tensor_hexahedra = nb_tensor_hexahedra
            # Get prisms
            line = viewLines[7].split()
            line = [eval(obj) for obj in line]
            nb_scalar_prisms, nb_vector_prisms, nb_tensor_prisms = line
            self.nb_scalar_prisms = nb_scalar_prisms
            self.nb_vector_prisms = nb_vector_prisms
            self.nb_tensor_prisms = nb_tensor_prisms
            # Get pyramids
            line = viewLines[8].split()
            line = [eval(obj) for obj in line]
            nb_scalar_pyramids, nb_vector_pyramids, nb_tensor_pyramids = line
            self.nb_scalar_pyramids = nb_scalar_pyramids
            self.nb_vector_pyramids = nb_vector_pyramids
            self.nb_tensor_pyramids = nb_tensor_pyramids
            # Get texts
            line = viewLines[9].split()
            line = [eval(obj) for obj in line]
            nb_text2d, nb_text2d_chars, nb_text3d, nb_text3d_chars = line
            self.nb_text2d       = nb_text2d
            self.nb_text2d_chars = nb_text2d_chars
            self.nb_text3d       = nb_text3d
            self.nb_text3d_chars = nb_text3d_chars
            ############
            # Read field
            ############
            numLine = 9
            # Initialisation of variables
            time_steps = []
            scalar_lines = []
            for _ in range(nb_time_steps):
                numLine += 1
                time_step = eval(viewLines[numLine])
                time_steps.append(time_step)
                # Points
                #  Scalar
                if nb_scalar_points > 0:
                    self.ViewError.notReady("scalar_points")
                    raise self.ViewError
                #  Vector
                if nb_vector_points > 0:
                    self.ViewError.notReady("vector_points")
                    raise self.ViewError
                #  Tensor
                if nb_tensor_points > 0:
                    self.ViewError.notReady("tensor_points")
                    raise self.ViewError
                # Lines
                #  Scalar
                if nb_scalar_lines > 0:
                    # Create list
                    scalar_lines_t = []
                    numLine += 1
                    # Load each line: two points (p1, p2) each with
                    # the following information [xi, yi, zi, vi] where
                    # vi is the value of the field
                    for _ in range(nb_scalar_lines):
                        p1 = [None]*4
                        p2 = [None]*4
                        for k1 in range(3):
                            line = viewLines[numLine + k1].split()
                            line = [eval(obj) for obj in line]
                            p1[k1], p2[k1] = line
                        p1[3] = eval(viewLines[numLine + 3])
                        p2[3] = eval(viewLines[numLine + 4])
                        scalar_lines_t.append((p1, p2))
                        numLine += 5
                    # Add the set of lines to the whole set of scalar
                    # lines for each time step
                    scalar_lines.append(scalar_lines_t)
                #  Vector
                if nb_vector_lines > 0:
                    self.ViewError.notReady("vector_lines")
                    raise self.ViewError
                #  Tensor
                if nb_tensor_lines > 0:
                    self.ViewError.notReady("tensor_lines")
                    raise self.ViewError
                # Triangles
                #  Scalar
                if nb_scalar_triangles > 0:
                    self.ViewError.notReady("scalar_triangles")
                    raise self.ViewError
                #  Vector
                if nb_vector_triangles > 0:
                    self.ViewError.notReady("vector_triangles")
                    raise self.ViewError
                #  Tensor
                if nb_tensor_triangles > 0:
                    self.ViewError.notReady("tensor_triangles")
                    raise self.ViewError
                # Quadrangles
                #  Scalar
                if nb_scalar_quadrangles > 0:
                    self.ViewError.notReady("scalar_quadrangles")
                    raise self.ViewError
                #  Vector
                if nb_vector_quadrangles > 0:
                    self.ViewError.notReady("vector_quadrangles")
                    raise self.ViewError
                #  Tensor
                if nb_tensor_quadrangles > 0:
                    self.ViewError.notReady("tensor_quadrangles")
                    raise self.ViewError
                # Tetrahedra
                #  Scalar
                if nb_scalar_tetrahedra > 0:
                    self.ViewError.notReady("scalar_tetahedra")
                    raise self.ViewError
                #  Vector
                if nb_vector_tetrahedra > 0:
                    self.ViewError.notReady("vector_tetahedra")
                    raise self.ViewError
                #  Tensor
                if nb_tensor_tetrahedra > 0:
                    self.ViewError.notReady("tensor_tetahedra")
                    raise self.ViewError
                # Hexahedra
                #  Scalar
                if nb_scalar_hexahedra > 0:
                    print("Option: scalar_hexahedra")
                    raise self.ViewError
                #  Vector
                if nb_vector_hexahedra > 0:
                    print("Option: vector_hexahedra")
                    raise self.ViewError
                #  Tensor
                if nb_tensor_hexahedra > 0:
                    print("Option: tensor_hexahedra")
                    raise self.ViewError
                # Prisms
                #  Scalar
                if nb_scalar_prisms > 0:
                    self.ViewError.notReady("scalar_prisms")
                    raise self.ViewError
                #  Vector
                if nb_vector_prisms > 0:
                    self.ViewError.notReady("vector_prisms")
                    raise self.ViewError
                #  Tensor
                if nb_tensor_prisms > 0:
                    self.ViewError.notReady("tensor_prisms")
                    raise self.ViewError
                # Pyramids
                #  Scalar
                if nb_scalar_pyramids > 0:
                    self.ViewError.notReady("scalar_pyramids")
                    raise self.ViewError
                #  Vector
                if nb_vector_pyramids > 0:
                    self.ViewError.notReady("vector_pyramids")
                    raise self.ViewError
                #  Tensor
                if nb_tensor_pyramids > 0:
                    self.ViewError.notReady("tensor_pyramids")
                    raise self.ViewError
            # Add to self
            self.time_steps = time_steps
            if nb_scalar_lines > 0:
                self.scalar_lines = scalar_lines
            return
        # }}}

        # Get mesh from lines {{{
        def GetMeshFromLines(self, fieldType, step, deleteField = True):
            # Test and select field
            if fieldType.lower() == 'scalar':
                if self.nb_scalar_lines == 0:
                    self.ViewError.notInView("scalar lines")
                    raise self.ViewError
                else:
                    field = self.scalar_lines[step]
            elif fieldType.lower() == 'vector':
                if self.nb_vector_lines == 0:
                    self.ViewError.notInView("vector lines")
                    raise self.ViewError
                else:
                    field = self.vector_lines[step]
            elif fieldType.lower() == 'tensor':
                if self.nb_tensor_lines == 0:
                    self.ViewError.notInView("tensor lines")
                    raise self.ViewError
                else:
                    field = self.tensor_lines[step]
            else:
                slf.ViewError.notKeyWord(fieldType)
                raise self.ViewError
            # Initialisation
            nodes = []
            elements = []
            readPoints = []
            # Read first line and save it as first nodes and element
            li = field[0]
            numNod = 0
            nodes.append(li[numNod])
            numNod += 1
            nodes.append(li[numNod])
            elements.append((0, 1))
            for li in field[1:]:
                # Initialisation
                ele = []
                # Locate each point depending on whether they have
                # been already read in another element
                for pi in li:
                    # If read in another element get the index
                    if pi in nodes:
                        ele.append(nodes.index(pi))
                    else:
                        numNod += 1
                        ele.append(numNod)
                        nodes.append(pi)
                # Save element
                ele = tuple(ele)
                elements.append(ele)
            # Delete field if indicated
            if deleteField:
                nodes = [tuple(obj[:3]) for obj in nodes]
            else:
                nodes = [tuple(obj) for obj in nodes]
            return (nodes, elements)
        # }}}
        # }}}
# }}}

# Functions {{{
# Smooth nodes {{{
def SmoothNodes(nodes, factor = 1, order = 5, fix_contour = False):
    # Convert nodes to numpy array
    nodArray = np.array([list(node) for node in nodes])
    if not nodArray.ndim == 2:
        raise("Error: nodes must be a list of tuples or lists so that it can be converted into a 2D-numpy array.")
    rows, cols = nodArray.shape
    # Set filter window width
    length = int(rows/factor)
    length = length - 1 if length%2 == 0 else length
    # Smooth each column
    columns = [None]*cols
    for k1 in range(cols):
        columns[k1] = savgol_filter(nodArray[:, k1], length, order)
    smoothArray = np.column_stack(tuple(columns))
    # Convert to list of tuple format
    array = [tuple(smoothArray[k1, :]) for k1 in range(rows)]
    if fix_contour:
        nod0 = nodes[0]
        arr0 = array[0]
        diff = tuple([nod0[k1] - arr0[k1] for k1 in range(cols)])
        array = [(array[k1][0] + diff[0],
                  array[k1][1] + diff[1],
                  array[k1][2] + diff[2]) for k1 in range(rows)]
        array[0] = nodes[0]
        array.append(nodes[-1])
    return array
# }}}

# Organise nodes in 1D mesh {{{
# This function takes a 1D element list and returns the order in which
# the nodes appear. It is useful to reconstruct a meshed contour. If
# the mesh is closed, the organised mesh starts at the first node of
# elements and if it is open it starts at the first open node in
# elements
def OrganiseNodesIn1DMesh(nodes, elements):
    # Sort each element
    numEles = len(elements)
    sortEles = [list(tpl) for tpl in elements]
    for k1 in range(numEles):
        sortEles[k1].sort()
    # Closed or open contour
    allNodes = [nod for ele in elements for nod in ele]
    closed = numEles == len(set(allNodes))
    # Set first node
    if closed:
        nodei = firstNodes[0]
    else:
        firstNodes = [ele[0] for ele in sortEles]
        nodei = LeastFrequentElementsInList(firstNodes)
        if len(nodei) == 1:
            nodei = nodei[0]
        else:
            nodeif = nodei
            seconNodes = [ele[1] for ele in sortEles]
            nodei = LeastFrequentElementsInList(seconNodes)
            if (len(nodeif) == len(nodei)):
                nodei = [nod for nod in firstNodes if nod not in seconNodes]
            elif (len(nodeif) > len(nodei)):
                nodei = [nod for nod in firstNodes if nod not in seconNodes]
            elif (len(nodeif) < len(nodei)):
                nodei = [nod for nod in seconNodes if nod not in firstNodes]
            if len(nodei) == 1:
                nodei = nodei[0]
            elif len(nodei) == 2:
                nodei = nodei[0]
            else:
                raise('Error: There are probably more than one contour.')
    # Find the other nodes in order of appearance
    ordNodes = [nodei]
    cont = True
    while cont:
        found = False
        for ele in sortEles:
            if nodei in ele:
                if nodei == ele[0]:
                    nodei = ele[1]
                else:
                    nodei = ele[0]
                ordNodes.append(nodei)
                found = True
                break
        if found:
            sortEles.remove(ele)
            if len(sortEles) == 0:
                cont = False
        else:
            raise('Error: There are probably more than one contour.')
    # Construct organise nodes and elements
    newNodes = [nodes[nodId] for nodId in ordNodes]
    if closed:
        newNodes = newNodes[:-1]
    newElements = [None]*numEles
    for k1 in range(numEles):
        if k1 < numEles - 1:
            newElements[k1] = (k1 + 1, k1 + 2)
        else:
            if closed:
                newElements[k1] =  (k1 + 1, 1)
            else:
                newElements[k1] = (k1 + 1, k1 + 2)
    return (newNodes, newElements, closed)
# }}}

# Find the least frequent elements in a list {{{
def LeastFrequentElementsInList(inList):
    from collections import Counter
    count = Counter(inList)
    values = [key for key, value in count.items()
            if value == min(count.values())]
    return values
# }}}

# Remesh faces {{{
def Remesh(grpNamList, filNamList, lcList, faceNams, modelIn,
        recombine = False, smooth = False, factor = 1, order = 5):
    model = modelIn
    # Get the number of faces
    numFaces = len(faceNams)
    if numFaces > 1:
        numFaces -= 1
    if not isinstance(smooth, list):
        if smooth:
            smooth = [True]*numFaces
        else:
            smooth = [False]*numFaces
    if not isinstance(recombine, list):
        if recombine:
            recombine = [True]*numFaces
        else:
            recombine = [False]*numFaces
    # Test that the length of every list is numFaces
    if not((len(grpNamList) == numFaces) and (len(filNamList) ==
        numFaces) and (len(lcList) == numFaces) and
        (len(smooth) == numFaces) and (len(recombine) == numFaces)):
        raise("Error: all the list must have the same length. " +
        "In other words, they must have the same number of meshes. " +
        "In addition, take into account that the list faceNams has "
        + "length 1 for 1-face meshes and length n+1 for n-face " +
        "meshes (with n>1). This since faceNams = [namFac1, namFac2,"
        + " ..., namFacn, namWholeMesh]. Finally, if a given face " +
        "should not have a name, add None to the list.")
    # Setup of the mesh parameters
    numContours = [None]*numFaces
    nodes       = [None]*numFaces
    elements    = [None]*numFaces
    closed      = [None]*numFaces
    for k1 in range(numFaces):
        # Get the number of contours
        numConts = len(grpNamList[k1])
        if not isinstance(smooth[k1], list):
            if smooth[k1]:
                smooth[k1] = [True]*numConts
            else:
                smooth[k1] = [False]*numConts
        # Test that the length of every list k1 is numConts
        if not((len(filNamList[k1]) == numConts) and
                (len(lcList[k1]) == numConts) and
                (len(smooth[k1]) == numConts)):
            raise("Error: all the list must have the same length. " +
                    "In other words, they must have the same " +
                    "number pf contours since each index points to " +
                    "a specific contour. Thus, the number of indices "
                    + "must be equal.")
        numContours[k1]  = numConts
        # Get nodes, elements and closeness for each contour
        filNams = filNamList[k1]
        nds     = [None]*numConts
        elmnts  = [None]*numConts
        clsd    = [None]*numConts
        for k2 in range(numConts):
            resu = GmshOutput(filNams[k2])
            nds[k2], elmnts[k2] = resu.views[0].GetMeshFromLines(
                    'scalar', 0)
            nds[k2], elmnts[k2], clsd[k2] = OrganiseNodesIn1DMesh(
                    nds[k2], elmnts[k2])
            if smooth[k1][k2]:
                nds[k2] = SmoothNodes(nds[k2], fix_contour = True, order = 5)
        nodes[k1]       = nds
        elements[k1]    = elmnts
    # Generate unique ids for each point and merge the reapeted
    uniqueNodes = [None]*numFaces
    contours    = [None]*numFaces
    uniqueID = 0
    for k1 in range(numFaces):
        unqNds = []
        cntrs  = []
        nds = nodes[k1]
        lc = lcList[k1]
        for k2 in range(numContours[k1]):
            nodeSet = nds[k2]
            lci = lc[k2]
            contour = []
            for node in nodeSet:
                # test whether the node has already been defined
                undefined = True
                for uniqueNode in unqNds:
                    if node == uniqueNode[0]:
                        undefined = False
                        ID = uniqueNode[1]
                        break
                # If the node has not been defined yet, add it to the list
                # of nodes to be drawn and assign it a unique id
                if undefined:
                    uniqueID += 1
                    unqNds.append((node, uniqueID, lci))
                    ID = uniqueID
                # Indicate the ID of the node to the contour to define the
                # spline
                contour.append(ID)
            cntrs.append(contour)
        uniqueNodes[k1] = unqNds
        contours[k1] = cntrs
    # Create nodes
    for unqNds in uniqueNodes:
        for node, ID, lci in unqNds:
            xi = node[0]
            yi = node[1]
            zi = node[2]
            model.geo.addPoint(xi, yi, zi, lci, ID)
    # Create contours
    splines = [None]*numFaces
    for k1 in range(numFaces):
        splns = [None]*numContours[k1]
        cntrs = contours[k1]
        for k2 in range(numContours[k1]):
            splns[k2] = model.geo.addSpline(cntrs[k2])
        splines[k1] = splns
    # Closed curve loops
    closedContours = [None]*numFaces
    for k1 in range(numFaces):
        splns = splines[k1]
        closedContours[k1] = model.geo.addCurveLoop(splns)
    # Faces
    faces = [None]*numFaces
    for k1 in range(numFaces):
        face = faces[k1]
        faces[k1] = model.geo.addPlaneSurface([closedContours[k1]])
        if recombine[k1]:
            model.geo.mesh.setRecombine(2, faces[k1])
    # Synchronize
    model.geo.synchronize()
    # Create physical groups
    for k1 in range(numFaces):
        splGroups = [None]*numContours[k1]
        splns = splines[k1]
        grpNams = grpNamList[k1]
        face = faces[k1]
        # Add 1D groups
        for k2 in range(numContours[k1]):
            if not(grpNams[k2] == None):
                splGroups[k2] = model.addPhysicalGroup(1,
                        [splns[k2]])
                model.setPhysicalName(1, splGroups[k2],
                        grpNams[k2])
        # Add 2D groups
        if not(faceNams[k1] == None):
            faceGr = model.addPhysicalGroup(2, [face])
            model.setPhysicalName(2, faceGr, faceNams[k1])
    if numFaces > 1:
        if not(faceNams[numFaces] ==None):
            faceGr = model.addPhysicalGroup(2, faces)
            model.setPhysicalName(2, faceGr, faceNams[numFaces])
    return model
# }}}

# Remesh line {{{
def RemeshLine(model, recoFile, lc):
    # Read reconstruction file
    resu = GmshOutput(recoFile)
    nds, eles = resu.views[0].GetMeshFromLines('scalar', 0)
    nds, eles, _ = OrganiseNodesIn1DMesh(nds, eles)
    # Create points
    points = []
    for node in nds:
        xi = node[0]
        yi = node[1]
        zi = node[2]
        points.append(model.geo.addPoint(xi, yi, zi, lc))
    # Make line
    line = model.geo.addSpline(points)
    return points, line
# }}}

# Define box fields with increasing mesh size {{{
def BoxFieldsIncreasingSize(mesh, f0, x0, y0, minWidth, minHeight,
        minSize, maxSize, grwRatioX, grwRatioY, grwRatioS, numFields,
        fix):
    """Set of box fields with increasing element size. Returns a list
    with the Ids of the fields."""
    # Functions to define the limits based on fix {{{
    func_fix_low_min = lambda c0, l, r, n : c0
    func_fix_low_max = lambda c0, l, r, n : c0 + l*r**n
    func_fix_upp_min = lambda c0, l, r, n : c0 - l*r**n
    func_fix_upp_max = lambda c0, l, r, n : c0
    if fix == 'bl' or fix == 'lb':
        func_xmin = func_fix_low_min
        func_xmax = func_fix_low_max
        func_ymin = func_fix_low_min
        func_ymax = func_fix_low_max
    elif fix == 'tl' or fix == 'lt':
        func_xmin = func_fix_low_min
        func_xmax = func_fix_low_max
        func_ymin = func_fix_upp_min
        func_ymax = func_fix_upp_max
    elif fix == 'tr' or fix == 'rt':
        func_xmin = func_fix_upp_min
        func_xmax = func_fix_upp_max
        func_ymin = func_fix_upp_min
        func_ymax = func_fix_upp_max
    elif fix == 'br' or fix == 'bt':
        func_xmin = func_fix_upp_min
        func_xmax = func_fix_upp_max
        func_ymin = func_fix_low_min
        func_ymax = func_fix_low_max
    elif fix == 'cc' or fix == 'c':
        func_xmin = func_fix_upp_min
        func_xmax = func_fix_low_max
        func_ymin = func_fix_upp_min
        func_ymax = func_fix_low_max
    else:
        raise("Error: fix must be a string with a tuple from the sets: {l, r} (x direction) and {b, t} (y direction) or 'c' and 'cc' to fix the centre.")
    # }}}
    refField = f0
    fieldList = [None]*numFields
    # Define k1th field
    for k1 in range(numFields):
        fieldList[k1] = refField
        # Compute field parameters
        vin = minSize*grwRatioS**(k1)
        vout = maxSize
        xmin = func_xmin(x0, minWidth, grwRatioX, k1)
        ymin = func_ymin(y0, minHeight, grwRatioY, k1)
        xmax = func_xmax(x0, minWidth, grwRatioX, k1)
        ymax = func_ymax(y0, minHeight, grwRatioY, k1)
        thickness = minSize*grwRatioS
        # Make the field
        mesh.field.add("Box", refField)
        mesh.field.setNumber(refField, "VIn"      , vin)
        mesh.field.setNumber(refField, "VOut"     , vout)
        mesh.field.setNumber(refField, "XMin"     , xmin)
        mesh.field.setNumber(refField, "XMax"     , xmax)
        mesh.field.setNumber(refField, "YMin"     , ymin)
        mesh.field.setNumber(refField, "YMax"     , ymax)
        mesh.field.setNumber(refField, "Thickness", thickness)
        refField += 1
    return fieldList
# }}}

# Mesh volume 1D {{{
def MeshVolume1D(nodes, elements):
    # Get initial information {{{
    numEles = len(elements)
    nnod = len(elements[0])
    if not nnod == 2:
        raise("Error: for the moment only linear elements are available.")
    # }}}
    # Compute the length of each element and sum it up
    volume = 0.0
    for ele in elements:
        n1 = nodes[ele[0] - 1]
        n2 = nodes[ele[1] - 1]
        x_c = n2[0] - n1[0]
        y_c = n2[1] - n1[1]
        z_c = n2[2] - n1[2]
        volume += math.sqrt(x_c*x_c + y_c*y_c + z_c*z_c)
    return volume
# }}}

# Define fields around an arbitrary contact line {{{
def MeshFieldArbitraryContact(mesh, f0, curveList, xmin0, xmax0, ymin0,
        ymax0, sizemin, grw_xmin, grw_xmax, grw_ymin, grw_ymax,
        grw_size, sizemax, lenratio, numFields = 25, sampling = 10000,
        VOut = 1.0):
    # Initialisation {{{
    fieldList = [None]*numFields
    ff = f0
    # }}}
    # Make base distance field {{{
    di_fi = ff
    mesh.field.add("Distance", ff)
    mesh.field.setNumbers(ff, "CurvesList", curveList)
    mesh.field.setNumber(ff, "Sampling", sampling)
    # }}}
    # Make boxes and thresholds {{{
    xmin = xmin0
    xmax = xmax0
    ymin = ymin0
    ymax = ymax0
    size = sizemin
    for k1 in range(numFields):
        # Set box
        ff += 1
        mesh.field.add("Box", ff)
        mesh.field.setNumber(ff, "VIn"      , 0.0)
        mesh.field.setNumber(ff, "VOut"     , VOut)
        mesh.field.setNumber(ff, "XMin"     , xmin)
        mesh.field.setNumber(ff, "XMax"     , xmax)
        mesh.field.setNumber(ff, "YMin"     , ymin)
        mesh.field.setNumber(ff, "YMax"     , ymax)
        mesh.field.setNumber(ff, "Thickness", 0.0)
        # Set threshold base field
        ff += 1
        mesh.field.add("MathEval", ff)
        mesh.field.setString(ff, "F", "F%d + F%d" % (ff - 1, di_fi))
        # Set mesh size field as threshold
        ff += 1
        fieldList[k1] = ff
        mesh.field.add("Threshold", ff)
        mesh.field.setNumber(ff, "InField", ff - 1)
        mesh.field.setNumber(ff, "SizeMin", size)
        mesh.field.setNumber(ff, "SizeMax", sizemax)
        mesh.field.setNumber(ff, "DistMin", size*lenratio)
        mesh.field.setNumber(ff, "DistMax", size*lenratio*2.0)
        # Update sizes
        xmin *= grw_xmin
        xmax *= grw_xmax
        ymin *= grw_ymin
        ymax *= grw_ymax
        size *= grw_size
    # }}}
    # Set thresholds as background mesh {{{
    bacFieldNum = ff + 1
    mesh.field.add("Min", bacFieldNum)
    mesh.field.setNumbers(bacFieldNum, "FieldsList", fieldList)
    mesh.field.setAsBackgroundMesh(bacFieldNum)
    fieldList.append(bacFieldNum)
    # }}}
    return fieldList
# }}}

# Mesh field distance--threshold to curve {{{
def MeshFieldDistanceCurve(mesh, f0, curveList, sizeMin, distMin,
        sizeGrwRate, distGrwRate, sizeMax, numFields,
        growthStyle = 'linear', sampling = 10000, at = "CurvesList"):
    # Initialisation {{{
    fieldList = [None]*numFields
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
        fieldList[k1] = ff
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
        else:
            raise("Error: growth style not available")
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
# }}}

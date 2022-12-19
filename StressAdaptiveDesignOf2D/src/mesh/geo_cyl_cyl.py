# Libraries {{{
import os
import sys
import gmsh
import getopt
import json
import math
import numpy as np

# In-house modules
sys.path.append(os.path.dirname(os.path.abspath(__file__))[:-4] + 'PythonUtilities/')
import ReadGmsh as rgmsh
# }}}

# Get output file name and popup call {{{
try:
    opts, args = getopt.getopt(sys.argv[1:], 'n:pi:frd')
    if len(args)>0:
        raise getopt.GetoptError('Too many arguments!')
except getopt.GetoptError:
    print('Wrong call! The excecution can have the following arguments:')
    print('To indicate the export name: -n "fileName"')
    print('To indicate popup True use: -p')
    print('To indicate the input parameters: -i "fileName" (.json)')
    print('To indicate full geometry -f')
    print('To indicate remeshing -r')
    print('To indicate remeshing of both the cylinder and the plane -d')
    raise

inParams = False
export = False
popup = False
fullGeo = False
remesh = False
double = False
for opt, arg in opts:
    # File name
    if opt in ['-n']:
        export = True
        fileName = arg
    elif opt in ['-p']:
        popup = True
    elif opt in ['-f']:
        fullGeo = True
    elif opt in ['-r']:
        remesh = True
    elif opt in ['-d']:
        remesh = True
        double = True
    elif opt in ['-i']:
        parFile = arg
        inParams = True
if not inParams:
    print("Wrong call! The input file with the parameters is not " +
            "present. Use -i 'fileName' (.json)")
    raise
# }}}

# Parameters {{{
# Default {{{
baseDir = os.getcwd()
surfName        = "_surf"
recoFolder      = "/RECONSTRUCTION"
repoFile        =  "/report.tab"
if export:
    fileName, fileNameExt = os.path.splitext(fileName)
base_y0_Name    = "bas_y0"
base_left_Name  = "bs_left"
base_cont_Name  = "bas_cont"
baseName        = "base"
wp_x0_Name      = "wp_x0"
wp_load_Name    = "wp_load"
wp_cont_Name    = "wp_cont"
wpName          = "wp"
wpGroups = [wp_x0_Name,
            wp_load_Name,
            wp_cont_Name,
            wpName]
bsGroups = [base_y0_Name,
            base_left_Name,
            base_cont_Name,
            baseName]
wpx0ids = [0]
wpldids = [1]
wpcntid = 2
wpnmid  = 3
wpNoGy_int_ids = [0, 2]
bsy0ids = [0]
bsx0ids = [1]
bscntid =  2
bsnmid  =  3
bsNoGy_int_ids = [1, 2]
# }}}
# External {{{
with open(parFile, 'r') as inFile:
    params = json.load(inFile)
meshParams = params["mesh"]
fileNames  = params["fileNames"]

# Geometry type and parameters {{{
geoParams = meshParams["geoParams"]
workpiece_w = geoParams[0]
# }}}

lcMin       = meshParams["lcMin"]
lcMax       = meshParams["lcMax"]
fine_y      = meshParams["fine_y"]
growthRatio = meshParams["growthRatio"]
lengthRatio = meshParams["lengthRatio"]
try:
    numFields_c = meshParams["numFields_c"]
except KeyError:
    numFields_c = 25
try:
    numFields_p = meshParams["numFields_p"]
except KeyError:
    numFields_p = 25
baseWidth   = 10.0*workpiece_w

lcMaxwp = lcMax
lcMaxbs = lcMax

asterFolder = fileNames["asterFolder"]
meshGroupFileName = fileNames["meshGroupFileName"]
workDir = baseDir + asterFolder
# }}}
# }}}

# Write mesh groups or read contact area {{{
if not remesh:
    meshGroups = {}
    meshGroups["wpGroups"] = wpGroups
    meshGroups["bsGroups"] = bsGroups
    meshGroups["wpx0ids"] = wpx0ids
    meshGroups["wpldids"] = wpldids
    meshGroups["wpcntid"] = wpcntid
    meshGroups["wpnmid"] = wpnmid
    meshGroups["wpNoGy_int_ids"] = wpNoGy_int_ids
    meshGroups["bsy0ids"] = bsy0ids
    meshGroups["bsx0ids"] = bsx0ids
    meshGroups["bscntid"] = bscntid
    meshGroups["bsnmid"] = bsnmid
    meshGroups["bsNoGy_int_ids"] = bsNoGy_int_ids
    with open(workDir + meshGroupFileName, 'w') as ofile:
        json.dump(meshGroups, ofile, indent = 4)
else:
    with open(workDir + repoFile, 'r') as ifile:
        contArea = ifile.readlines()
    contArea = contArea[-1]
    contArea = contArea.split(',')
    contArea = eval(contArea[4])
# }}}

# Organisation of the parameters {{{
# Be careful! All the following items must have consistent Ids. That
# means that $$$FileNams[i], $$$GroNams[i] and $$$Lc[i] together
# describe the ith contour of the $$$ domain.
wpNumGrps = len(wpGroups) - 1
wpFileNams = [None]*wpNumGrps
for k1 in range(wpNumGrps):
    wpFileNams[k1] = workDir + recoFolder + '/' + wpGroups[k1] + '.msh'
wpGrpNams = wpGroups[:-1]
wpLc =  [lcMaxwp]*wpNumGrps

bsNumGrps = len(bsGroups) - 1
bsFileNams = [None]*bsNumGrps
for k1 in range(bsNumGrps):
    bsFileNams[k1] = workDir + recoFolder + '/' + bsGroups[k1] + '.msh'
bsGrpNams = bsGroups[:-1]
bsLc = [lcMaxbs]*bsNumGrps
# }}}

# GMSH {{{
if remesh:
    # Remesh workpiece geometry {{{
    # Initialisation {{{
    gmsh.initialize()
    gmsh.clear()
    gmsh.model.add("geo")
    model = gmsh.model
    geo = model.geo
    mesh = model.mesh
    # }}}

    # Mesh reconstruction {{{
    resu = rgmsh.GmshOutput(wpFileNams[2])
    nds, eles = resu.views[0].GetMeshFromLines('scalar', 0)
    nds, eles, _ = rgmsh.OrganiseNodesIn1DMesh(nds, eles)
    # Create points
    points = []
    pointIds = []
    xx = nds[0][0]
    yy = nds[0][1]
    for node in nds:
        xi = node[0]
        yi = node[1]
        if yi < workpiece_w and xi < 2.0*workpiece_w:
            xa = xx
            ya = yy
            xx = xi
            yy = yi
            points.append([xi, yi, 0.0])
            pointIds.append(geo.addPoint(xi, yi, 0.0, lcMax))
    dx = xx - xa
    xi = xx + dx
    yi = workpiece_w
    points.append([xi, yi, 0.0])
    pointIds.append(geo.addPoint(xi, yi, 0.0, lcMax))
    points = np.array(points)
    # Make line
    wp_l_cont = geo.addBSpline(pointIds)
    # Top line
    yi = points[-1, 1]
    wp_p_tr = pointIds[-1]
    wp_p_tl = geo.addPoint(0., yi, 0.0, lcMax)
    wp_l_t = geo.addLine(wp_p_tr, wp_p_tl)
    # Left line
    wp_p_bl = pointIds[0]
    wp_l_l = geo.addLine(wp_p_tl, wp_p_bl)
    wp_loop = geo.addCurveLoop([wp_l_cont, wp_l_t, wp_l_l])
    wp_face = geo.addPlaneSurface([wp_loop])
    model.geo.synchronize()
    # }}}

    # Create mesh groups {{{
    wp_fixed_x = model.addPhysicalGroup(1, [wp_l_l])
    model.setPhysicalName(1, wp_fixed_x, wp_x0_Name)
    wp_load = model.addPhysicalGroup(1, [wp_l_t])
    model.setPhysicalName(1, wp_load, wp_load_Name)
    wp_cont = model.addPhysicalGroup(1, [wp_l_cont])
    model.setPhysicalName(1, wp_cont, wp_cont_Name)

    wp = model.addPhysicalGroup(2, [wp_face])
    model.setPhysicalName(2, wp, wpName)
    # }}}

    # Find lowest x and y in wp contact line {{{
    nds = points
    wpx0 = nds[:, 0].min()
    wpy0id = np.argmin(nds[:, 1])
    wpy0 = nds[:, 1][wpy0id]
    wpy0x = nds[:, 0][wpy0id]
    #geo.translate([(2,1)], 0.0, -wpy0, 0.0)
    # Find fine_x (fine mesh from (wpx0, wpy0) to (fine_x, fine_y))
    min_dist = 1000000.0
    fine_y_mv = fine_y + wpy0
    for node in nds:
        xx = node[0]
        yy = node[1]
        if xx >= wpy0x:
            if min_dist > abs(yy - fine_y_mv):
                min_dist = abs(yy - fine_y_mv)
                fine_x = node[0]
                nnn = node
    # }}}

    # Create fields for mesh size {{{
    length_x = nds[:, 0].max()
    length_y = nds[:, 1].max()
    lcMaxwp = lcMax*max(length_x, length_y)/workpiece_w
    fieldList = rgmsh.MeshFieldArbitraryContact(mesh, 1, [wp_l_cont],
            -contArea, 1.25*contArea, -1000, 1000, lcMin, 1.05, 1.05, 1.0, 1.0,
            1.11, lcMaxwp, lengthRatio, numFields = numFields_c,
            VOut = max(length_x, length_y))
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    # }}}

    # Generate mesh {{{
    gmsh.option.setNumber("Mesh.Algorithm", 8)
    mesh.generate(2)
    mesh.recombine()
    gmsh.option.setNumber("Mesh.Smoothing", 10)
    # }}}

    # Save mesh {{{
    if export:
        gmsh.write(fileName + 'wp' + fileNameExt)
    # }}}

    # Pop up {{{
    if popup:
        gmsh.fltk.run()
    # }}}

    # End
    gmsh.finalize()
    # }}}
else:
    # Make workpiece geometry {{{
    # Initialisation {{{
    gmsh.initialize()
    gmsh.clear()
    gmsh.model.add("geo")
    model = gmsh.model
    geo = model.geo
    mesh = model.mesh
    # }}}

    # Create geometry {{{
    if fullGeo:
        # Points
        wp_tl = geo.addPoint(-workpiece_w, workpiece_w, 0.0, 10.0)
        wp_tr = geo.addPoint(workpiece_w, workpiece_w, 0.0, 10.0)
        wp_tc = geo.addPoint(0.0, workpiece_w, 0.0, 10.0)
        # Lines
        wp_c = geo.addCircleArc(wp_tl, wp_tc, wp_tr)
        wp_t = geo.addLine(wp_tr, wp_tl)
        wp_clist = [wp_c]
    else:
        # Points
        wp_bl = geo.addPoint(0.0, 0.0, 0.0, 10.0)
        wp_tr = geo.addPoint(workpiece_w, workpiece_w, 0.0, 10.0)
        wp_tl = geo.addPoint(0.0, workpiece_w, 0.0, 10.0)
        # Lines
        wp_c = geo.addCircleArc(wp_bl, wp_tl, wp_tr)
        wp_t = geo.addLine(wp_tr, wp_tl)
        wp_l = geo.addLine(wp_tl, wp_bl)
        wp_clist = [wp_c]
    # Workpiece loop and face
    if fullGeo:
        wp_loop = geo.addCurveLoop([wp_t] + wp_clist)
    else:
        wp_loop = geo.addCurveLoop([wp_t, wp_l] + wp_clist)
    wp_face = geo.addPlaneSurface([wp_loop])
    # Synchronize
    geo.synchronize()
    # }}}

    # Create fields for mesh size {{{
    fieldList = rgmsh.MeshFieldArbitraryContact(mesh, 1, wp_clist, -baseWidth,
            baseWidth, 0.0, fine_y, lcMin, 1.0, 1.0, 1.0, 1.1, 1.1, lcMax,
            lengthRatio, numFields = numFields_c)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    # }}}

    # Create mesh groups {{{
    wp_fixed_x = model.addPhysicalGroup(1, [wp_l])
    model.setPhysicalName(1, wp_fixed_x, wp_x0_Name)
    wp_load = model.addPhysicalGroup(1, [wp_t])
    model.setPhysicalName(1, wp_load, wp_load_Name)
    wp_cont = model.addPhysicalGroup(1, wp_clist)
    model.setPhysicalName(1, wp_cont, wp_cont_Name)

    wp = model.addPhysicalGroup(2, [wp_face])
    model.setPhysicalName(2, wp, wpName)
    # }}}

    # Generate mesh {{{
    gmsh.option.setNumber("Mesh.Algorithm", 8)
    mesh.generate(2)
    mesh.recombine()
    gmsh.option.setNumber("Mesh.Smoothing", 10)
    # }}}

    # Save mesh {{{
    if export:
        gmsh.write(fileName + 'wp' + fileNameExt)
    # }}}

    # Pop up {{{
    if popup:
        gmsh.fltk.run()
    # }}}

    # End
    gmsh.finalize()
    # }}}
if double:
    # Remesh base geometry {{{
    # Initialisation {{{
    gmsh.initialize()
    gmsh.clear()
    gmsh.model.add("geo")
    model = gmsh.model
    geo = model.geo
    mesh = model.mesh
    # }}}

    # Mesh reconstruction {{{
    model = rgmsh.Remesh([bsGrpNams], [bsFileNams], [bsLc], [baseName], model,
            recombine = False)
    # Synchronize
    model.geo.synchronize()
    # }}}

    # Find the coordinates of the left top point {{{
    resu = rgmsh.GmshOutput(bsFileNams[2])
    lcoord, _ = resu.views[0].GetMeshFromLines('scalar', 0)
    lcoord = np.array([list(nod) for nod in lcoord])
    maxlefty = lcoord[:, 1].max()
    geo.translate([(2,1)], 0.0, -maxlefty, 0.0)
    # }}}

    # Create fields for mesh size {{{
    resu = rgmsh.GmshOutput(bsFileNams[2])
    nds, _ = resu.views[0].GetMeshFromLines('scalar', 0)
    nds = np.array([list(nod) for nod in nds])
    length_x = nds[:, 0].max()
    length_y = abs(nds[:, 1].min())
    lcMaxbs = lcMax*max(length_x, length_y)/workpiece_w
    fieldList = rgmsh.MeshFieldArbitraryContact(mesh, 1, [bscntid + 1],
            -contArea, contArea, -1000, 1000, lcMin, 1.05, 1.05, 1.0,
            1.0, 1.11, lcMax, lengthRatio, numFields = numFields_p,
            VOut = max(length_x, length_y))
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    # }}}

    # Generate mesh {{{
    gmsh.option.setNumber("Mesh.Algorithm", 8)
    mesh.generate(2)
    mesh.recombine()
    gmsh.option.setNumber("Mesh.Smoothing", 10)
    # }}}

    # Save mesh {{{
    if export:
        gmsh.write(fileName + 'bs' + fileNameExt)
    # }}}

    # Pop up {{{
    if popup:
        gmsh.fltk.run()
    # }}}

    # End
    gmsh.finalize()
    # }}}
else:
    # Make base geometry {{{
    # Initialisation {{{
    gmsh.initialize()
    gmsh.clear()
    gmsh.model.add("geo")
    model = gmsh.model
    geo = model.geo
    mesh = model.mesh
    # }}}

    # Create geometry {{{
    if fullGeo:
        # Points
        bs_bl = geo.addPoint(-workpiece_w, -workpiece_w, 0.0, 10.0)
        bs_br = geo.addPoint(workpiece_w, -workpiece_w, 0.0, 10.0)
        bs_bc = geo.addPoint(0.0, -workpiece_w, 0.0, 10.0)
        # Lines
        bs_c = geo.addCircleArc(bs_br, bs_bc, bs_bl)
        bs_b = geo.addLine(bs_bl, bs_br)
        bs_clist = [bs_c]
    else:
        # Points
        bs_tl = geo.addPoint(0.0, 0.0, 0.0, 10.0)
        bs_br = geo.addPoint(workpiece_w, -workpiece_w, 0.0, 10.0)
        bs_bl = geo.addPoint(0.0, -workpiece_w, 0.0, 10.0)
        # Lines
        bs_c = geo.addCircleArc(bs_tl, bs_bl, bs_br)
        bs_b = geo.addLine(bs_br, bs_bl)
        bs_l = geo.addLine(bs_bl, bs_tl)
        bs_clist = [bs_c]
    # Workpiece loop and face
    if fullGeo:
        bs_loop = geo.addCurveLoop([bs_b] + bs_clist)
    else:
        bs_loop = geo.addCurveLoop([bs_b, bs_l] + bs_clist)
    bs_face = geo.addPlaneSurface([bs_loop])
    # Synchronize
    geo.synchronize()
    # }}}

    # Create fields for mesh size {{{
    if remesh:
        fieldList = rgmsh.MeshFieldArbitraryContact(mesh, 1, bs_clist,
                -contArea, 1.25*contArea, -1000, 1000, lcMin, 1.05, 1.05, 1.0, 1.0,
                1.11, lcMax, lengthRatio, numFields = numFields_p,
                VOut = max(length_x, length_y))
    else:
        fieldList = rgmsh.MeshFieldArbitraryContact(mesh, 1, bs_clist,
                -baseWidth, baseWidth, -fine_y, 0.0, lcMin, 1.0, 1.0, 1.1, 1.0,
                1.1, lcMax, lengthRatio, numFields = numFields_p)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    # }}}

    # Create mesh groups {{{
    bs_fixed_y = model.addPhysicalGroup(1, [bs_b])
    model.setPhysicalName(1, bs_fixed_y, base_y0_Name)
    base_cont = model.addPhysicalGroup(1, bs_clist)
    model.setPhysicalName(1, base_cont, base_cont_Name)
    base_left = model.addPhysicalGroup(1, [bs_l])
    model.setPhysicalName(1, base_left, base_left_Name)

    semi = model.addPhysicalGroup(2, [bs_face])
    model.setPhysicalName(2, semi, baseName)
    # }}}

    # Generate mesh {{{
    gmsh.option.setNumber("Mesh.Algorithm", 8)
    mesh.generate(2)
    mesh.recombine()
    gmsh.option.setNumber("Mesh.Smoothing", 10)
    # }}}

    # Save mesh {{{
    if export:
        gmsh.write(fileName + 'bs' + fileNameExt)
    # }}}

    # Pop up {{{
    if popup:
        gmsh.fltk.run()
    # }}}

    # End
    gmsh.finalize()
    # }}}
# }}}

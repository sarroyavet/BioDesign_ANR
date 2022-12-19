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
if export:
    fileName, fileNameExt = os.path.splitext(fileName)
base_y0_Name    = "bas_y0"
base_left_Name  = "bs_left"
base_right_Name = "bs_right"
base_cont_Name  = "bas_cont"
baseName        = "base"
wp_x0_Name      = "wp_x0"
bs_x0_Name      = "bs_x0"
wp_load_Name    = "wp_load"
wp_cont_Name    = "wp_cont"
wpName          = "wp"
if fullGeo:
    wpGroups = [wp_load_Name,
                wp_cont_Name,
                wpName]
else:
    wpGroups = [wp_x0_Name,
                wp_load_Name,
                wp_cont_Name,
                wpName]
bsGroups = [base_y0_Name,
            base_left_Name,
            base_right_Name,
            base_cont_Name,
            baseName]
if fullGeo:
    wpx0ids = []
    wpldids = [0]
    wpcntid = 1
    wpnmid  = 2
    wpNoGy_int_ids = []
else:
    wpx0ids = [0]
    wpldids = [1]
    wpcntid = 2
    wpnmid  = 3
    wpNoGy_int_ids = [0, 2]
bsy0ids = [0]
bsx0ids = [1, 2]
bscntid =  3
bsnmid  =  4
if fullGeo:
    bsNoGy_int_ids = [1, 2]
else:
    bsNoGy_int_ids = [1, 3]
# }}}
# External {{{
with open(parFile, 'r') as inFile:
    params = json.load(inFile)
meshParams = params["mesh"]
fileNames  = params["fileNames"]

# Geometry type and parameters {{{
geoParams = meshParams["geoParams"]
workpiece_w = geoParams[0]
baseWidth   = geoParams[1]
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
lcMaxwp = lcMax
lcMaxbs = lcMax*baseWidth/workpiece_w
fine_x = math.sqrt(2.0*fine_y*workpiece_w - fine_y**2.0)

asterFolder = fileNames["asterFolder"]
meshGroupFileName = fileNames["meshGroupFileName"]
workDir = baseDir + asterFolder
# }}}
# }}}

# Write mesh groups {{{
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
    os.makedirs(workDir, exist_ok = True)
    with open(workDir + meshGroupFileName, 'w') as ofile:
        json.dump(meshGroups, ofile, indent = 4)
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
    model = rgmsh.Remesh([wpGrpNams], [wpFileNams], [wpLc], [wpName], model,
            recombine = False)
    # Synchronize
    model.geo.synchronize()
    # }}}

    # Find lowest x and y in wp contact line {{{
    resu = rgmsh.GmshOutput(wpFileNams[2])
    nds, _ = resu.views[0].GetMeshFromLines('scalar', 0)
    nds = np.array([list(nod) for nod in nds])
    wpx0 = nds[:, 0].min()
    wpy0id = np.argmin(nds[:, 1])
    wpy0 = nds[:, 1][wpy0id]
    wpy0x = nds[:, 0][wpy0id]
    geo.translate([(2,1)], 0.0, -wpy0, 0.0)
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
    fieldList = rgmsh.MeshFieldArbitraryContact(mesh, 1, [wpcntid + 1], -baseWidth,
            baseWidth, 0.0, fine_y, lcMin, 1.0, 1.0, 1.0, 1.1, 1.1, lcMax,
            lengthRatio, numFields = numFields_c)
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
    lcMax = lcMaxwp
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
        wp_bc = geo.addPoint(0.0, 0.0, 0.0, 10.0)
        # Lines
        wp_c_l = geo.addCircleArc(wp_tl, wp_tc, wp_bc)
        wp_c_r = geo.addCircleArc(wp_bc, wp_tc, wp_tr)
        wp_t = geo.addLine(wp_tr, wp_tl)
        wp_clist = [wp_c_l, wp_c_r]
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
    model.geo.synchronize()
    if fullGeo:
        #wp_fixed_x = model.addPhysicalGroup(0, [wp_bc])
        #model.setPhysicalName(0, wp_fixed_x, wp_x0_Name)
        pass
    else:
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
        #gmsh.option.setNumber("Mesh.SaveGroupsOfNodes", -1)
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
    resu = rgmsh.GmshOutput(bsFileNams[1])
    lcoord, _ = resu.views[0].GetMeshFromLines('scalar', 0)
    lcoord = np.array([list(nod) for nod in lcoord])
    maxlefty = lcoord[:, 1].max()
    geo.translate([(2,1)], 0.0, -maxlefty, 0.0)
    # }}}

    # Create fields for mesh size {{{
    fieldList = rgmsh.MeshFieldArbitraryContact(mesh, 1, [bscntid + 1],
            -baseWidth, fine_x, -fine_y, fine_y, lcMin, 1.0, 1.1, 1.2, 1.2,
            1.25, lcMaxbs, lengthRatio, numFields = numFields_p,
            VOut = baseWidth)
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
    lcMax = lcMaxbs
    lcMin = lcMin
    # Initialisation {{{
    gmsh.initialize()
    gmsh.clear()
    gmsh.model.add("geo")
    model = gmsh.model
    geo = model.geo
    mesh = model.mesh
    # }}}

    # Create geometry {{{
    base_w = baseWidth
    base_h = 4.0*workpiece_w
    if fullGeo:
        base_bl = geo.addPoint(-base_w, -base_h, 0.0, lcMax)
        base_cl = geo.addPoint(-base_w,  0.0, 0.0, lcMax)
        base_cc = geo.addPoint(0.0,  0.0, 0.0, lcMax)
    else:
        base_bl = geo.addPoint(0.0, -base_h, 0.0, lcMax)
        base_cl = geo.addPoint(0.0,  0.0, 0.0, lcMax)
    base_br = geo.addPoint(base_w, -base_h, 0.0, lcMax)
    base_cr = geo.addPoint(base_w,  0.0, 0.0, lcMax)

    base_b  = geo.addLine(base_bl, base_br)
    base_r  = geo.addLine(base_br, base_cr)
    if fullGeo:
        base_c_r = geo.addLine(base_cr, base_cc)
        base_c_l = geo.addLine(base_cc, base_cl)
        base_clist = [base_c_r, base_cl]
    else:
        base_c = geo.addLine(base_cr, base_cl)
        base_clist = [base_c]
    base_l  = geo.addLine(base_cl, base_bl)

    base_loop = geo.addCurveLoop([base_b, base_r] + base_clist + [base_l])

    base_face = geo.addPlaneSurface([base_loop])
    # Synchronize
    geo.synchronize()
    # }}}

    # Create fields for mesh size {{{
    # Make box fields for the base geometry
    if fullGeo:
        bsFields = rgmsh.BoxFieldsIncreasingSize(mesh, 1, 0.0,
                0.0, fine_x, lcMin*lengthRatio, lcMin, lcMax, 1.05,
                1.2, 1.25, numFields_p, 'cc')
    else:
        bsFields = rgmsh.BoxFieldsIncreasingSize(mesh, 1, 0.0,
                0.0, fine_x, lcMin*lengthRatio, lcMin, lcMax, 1.05,
                1.2, 1.25, numFields_p, 'tl')
    fieldList = bsFields
    # Threshold as background mesh
    bacFieldNum = bsFields[-1] + 1
    mesh.field.add("Min", bacFieldNum)
    mesh.field.setNumbers(bacFieldNum, "FieldsList", fieldList)
    mesh.field.setAsBackgroundMesh(bacFieldNum)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    # }}}

    # Create mesh groups {{{
    if fullGeo:
        #bs_fixed_x = model.addPhysicalGroup(0, [base_cc])
        #model.setPhysicalName(0, bs_fixed_x, bs_x0_Name)
        pass
    base_left = model.addPhysicalGroup(1, [base_l])
    model.setPhysicalName(1, base_left, base_left_Name)
    base_right = model.addPhysicalGroup(1, [base_r])
    model.setPhysicalName(1, base_right, base_right_Name)
    base_fixed_y = model.addPhysicalGroup(1, [base_b])
    model.setPhysicalName(1, base_fixed_y, base_y0_Name)
    base_cont = model.addPhysicalGroup(1, base_clist)
    model.setPhysicalName(1, base_cont, base_cont_Name)

    semi = model.addPhysicalGroup(2, [base_face])
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
#        gmsh.option.setNumber("Mesh.SaveGroupsOfNodes", -1)
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

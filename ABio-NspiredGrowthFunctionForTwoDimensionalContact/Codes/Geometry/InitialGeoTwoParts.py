# Libraries {{{
import os
import sys
import gmsh
import getopt
import json
import math
import numpy as np

# In-house modules
sys.path.append('PythonUtilities/')
import ReadGmsh as rgmsh
# }}}

# Get output file name and popup call {{{
try:
    opts, args = getopt.getopt(sys.argv[1:], 'n:pi:f')
    if len(args)>0:
        raise getopt.GetoptError('Too many arguments!')
except getopt.GetoptError:
    print('Wrong call! The excecution can have the following arguments:')
    print('To indicate the export name: -n "fileName"')
    print('To indicate popup True use: -p')
    print('To indicate the input parameters: -i "fileName" (.json)')
    print('To indicate full geometry -f')
    raise

inParams = False
export = False
popup = False
fullGeo = False
for opt, arg in opts:
    # File name
    if opt in ['-n']:
        export = True
        fileName = arg
    elif opt in ['-p']:
        popup = True
    elif opt in ['-f']:
        fullGeo = True
    elif opt in ['-i']:
        parFile = arg
        inParams = True
if not inParams:
    print("Wrong call! The input file with the parameters is not " +
            "present. Use -i 'fileName' (.json)")
    raise
# }}}

# Parameters {{{
if export:
    fileName, fileNameExt = os.path.splitext(fileName)
with open(parFile, 'r') as inFile:
    params = json.load(inFile)
meshParams = params["mesh"]

# Geometry type and parameters {{{
geoType   = meshParams["geoType"]
geoParams = meshParams["geoParams"]
if geoType == 'circ':
    workpiece_w = geoParams[0]
    if len(geoParams) == 2:
        numExtras = geoParams[1]
    else:
        numExtras = 10
elif geoType == 'lund':
    workpiece_w = geoParams[0]
    workpiece_f = geoParams[1]
elif geoType == 'para':
    workpiece_w = geoParams[0]
    workpiece_o = geoParams[1]
elif geoType == 'hype':
    workpiece_w = geoParams[0]
    workpiece_o = geoParams[1]
    workpiece_p = geoParams[2]
    if workpiece_p > 1.0:
        raise("Error: For the 'hype' the third parameter must belong to (0, 1.0]")
elif geoType == 'rec-arc':
    workpiece_w =geoParams[0]
    workpiece_r =geoParams[1]
    workpiece_h =geoParams[2]
else:
    raise("Error: unknown geoType, use either: 'hype', 'circ', 'para', 'rec-arc' or 'lund'.")
# }}}

lcMin       = meshParams["lcMin"]
lcMax       = meshParams["lcMax"]
fine_y      = meshParams["fine_y"]
growthRatio = meshParams["growthRatio"]
lengthRatio = meshParams["lengthRatio"]
baseWidth = meshParams["baseWidth"]

baseGroups = meshParams["baseGroups"]
base_x0_Name   = baseGroups[0]
base_y0_Name   = baseGroups[1]
base_cont_Name = baseGroups[2]
baseName       = baseGroups[3]

wpGroups = meshParams["wpGroups"]
wp_x0_Name     = wpGroups[0]
wp_load_Name   = wpGroups[1]
wp_cont_Name   = wpGroups[2]
wpName         = wpGroups[3]

surfName   = meshParams["surfName"]
# }}}

# GMSH {{{
# Workpiece geometry {{{
# Initialisation {{{
gmsh.initialize()
gmsh.clear()
gmsh.model.add("geo")
model = gmsh.model
geo = model.geo
mesh = model.mesh
# }}}

# Create geometry {{{
# Contact line
def Hype(x, a, b, n):
    term = 1.0 - abs(x/a)**n
    return -b*(term**(1.0/n) - 1.0)
def Poly(x, a, n):
    return (a*x)**n
def Lund_f_of_x(x, w, c = math.pi/2.0):
    return c*math.log(1.0/(1.0 - (x/(w))**2.0))
def Lund_f_of_y(y, w, c = math.pi/2.0):
    try:
        term = math.exp(y/c)
    except OverflowError:
        term = float('inf')
    return w*math.sqrt((1.0 - 1.0/term))
def Circ(x, r):
    return r - math.sqrt(r**2.0 - x**2.0)
numPoints = int(2*workpiece_w/lcMin)
if geoType == 'circ':
    if fullGeo:
        # Points
        wp_tl = geo.addPoint(-workpiece_w, workpiece_w, 0.0, 10.0)
        wp_tr = geo.addPoint(workpiece_w, workpiece_w, 0.0, 10.0)
        wp_tc = geo.addPoint(0.0, workpiece_w, 0.0, 10.0)
        # Lines
        wp_c = geo.addCircleArc(wp_tl, wp_tc, wp_tr)
        wp_t = geo.addLine(wp_tr, wp_tl)
        wp_clist = [wp_c]
        # Get fine_x
        fine_x = math.sqrt(2.0*fine_y*workpiece_w - fine_y**2.0)
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
        # Get fine_x
        fine_x = math.sqrt(2.0*fine_y*workpiece_w - fine_y**2.0)
elif geoType == 'rec-arc':
    # Points
    wp_tl = geo.addPoint(0.0, workpiece_h, 0.0, 10.0)
    wp_bl = geo.addPoint(0.0, 0.0, 0.0, 10.0)
    wp_tr = geo.addPoint(workpiece_w, workpiece_h, 0.0, 10.0)
    wp_cb = geo.addPoint(workpiece_w - workpiece_r, 0.0, 0.0, 10.0)
    wp_ct = geo.addPoint(workpiece_w, workpiece_r, 0.0,
            10.0)
    wp_cr = geo.addPoint(workpiece_w - workpiece_r, workpiece_r, 0.0, 10.0)
    # Lines
    wp_t = geo.addLine(wp_tr, wp_tl)
    wp_l = geo.addLine(wp_tl, wp_bl)
    wp_c1 = geo.addLine(wp_bl, wp_cb)
    wp_c2 = geo.addCircleArc(wp_cb, wp_cr, wp_ct)
    wp_c3 = geo.addLine(wp_ct, wp_tr)
    wp_clist = [wp_c1, wp_c2, wp_c3]
    # Get fine_x
    fine_x = math.sqrt(2.0*fine_y*workpiece_w - fine_y**2.0) + \
            workpiece_w
else:
    if geoType == 'lund':
        fine_x = 1.1*workpiece_w
        wp_tl = geo.addPoint(0.0, workpiece_w, 0.0, 10.0)
        wp_tr = geo.addPoint(workpiece_w, workpiece_w, 0.0, 10.0)
        y = np.linspace(0.0, 1.0, numPoints*10)
        x = np.array([Lund_f_of_y(yy, workpiece_w, workpiece_f) for yy in y])
    else:
        x = np.linspace(0.0, workpiece_w, numPoints)
        if geoType == 'para':
            wp_tl = geo.addPoint(0.0, 1.0, 0.0, 10.0)
            y = np.array([Poly(xx, 1.0/workpiece_w, workpiece_o) for xx in x])
        elif geoType == 'hype':
            a = workpiece_w
            b = workpiece_p
            n = workpiece_o
            wp_tl = geo.addPoint(0.0, b, 0.0, 10.0)
            y = np.array([Hype(xx, a, b, n) for xx in x])
    wp_spl_pois = []
    for k1 in range(numPoints):
        wp_spl_pois.append(geo.addPoint(x[k1], y[k1], 0.0, lcMax))
    if geoType == 'lund':
        wp_r = geo.addLine(wp_spl_pois[-1], wp_tr)
    wp_c = geo.addSpline(wp_spl_pois)
    # Top and left lines
    if geoType == 'lund':
        wp_t = geo.addLine(wp_tr, wp_tl)
    else:
        wp_t = geo.addLine(wp_spl_pois[-1], wp_tl)
    wp_l = geo.addLine(wp_tl, wp_spl_pois[0])
    wp_clist = [wp_c]
# Workpiece loop and face
if fullGeo:
    wp_loop = geo.addCurveLoop([wp_t] + wp_clist)
else:
    if geoType == 'lund':
        wp_loop = geo.addCurveLoop([wp_r,wp_t, wp_l] + wp_clist)
    else:
        wp_loop = geo.addCurveLoop([wp_t, wp_l] + wp_clist)
wp_face = geo.addPlaneSurface([wp_loop])
# Synchronize
geo.synchronize()
# }}}

# Create fields for mesh size {{{
# Fields from border and contact box {{{
numFields = 25
# Base distance field
ff = 1
di_fi = ff
mesh.field.add("Distance", ff)
mesh.field.setNumbers(ff, "CurvesList", wp_clist)
mesh.field.setNumber(ff, "Sampling", 10000)
fieldList = [None]*numFields
sizey = fine_y
minH = lcMin
grw_sizey = 1.1
grw_lcMin = 1.1
for k1 in range(numFields):
    # Set box
    ff += 1
    mesh.field.add("Box", ff)
    mesh.field.setNumber(ff, "VIn"      , 0.0)
    mesh.field.setNumber(ff, "VOut"     , 1.0)
    mesh.field.setNumber(ff, "XMin"     , 0.0)
    if fullGeo:
        mesh.field.setNumber(ff, "XMin"     , -baseWidth)
    else:
        mesh.field.setNumber(ff, "YMin"     , 0.0)
    mesh.field.setNumber(ff, "XMax"     , baseWidth)
    mesh.field.setNumber(ff, "YMax"     , sizey)
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
    mesh.field.setNumber(ff, "SizeMin", minH)
    mesh.field.setNumber(ff, "SizeMax", lcMax)
    mesh.field.setNumber(ff, "DistMin", minH*lengthRatio)
    mesh.field.setNumber(ff, "DistMax", minH*lengthRatio*2.0)
    # Update size
    sizey *= grw_sizey
    minH *= grw_lcMin
# }}}
# Threshold as background mesh
bacFieldNum = ff + 1
mesh.field.add("Min", bacFieldNum)
mesh.field.setNumbers(bacFieldNum, "FieldsList", fieldList)
mesh.field.setAsBackgroundMesh(bacFieldNum)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
# }}}

# Create mesh groups {{{
if fullGeo:
    wp_fixed_x = model.addPhysicalGroup(0, [wp_tr])
    model.setPhysicalName(0, wp_fixed_x, wp_x0_Name)
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
#gmsh.option.setNumber("Mesh.ElementOrder", 2)
#gmsh.option.setNumber("Mesh.SecondOrderLinear", 0)
gmsh.option.setNumber("Mesh.Algorithm", 8)
#mesh.setRecombine(2, surfName)
mesh.generate(2)
mesh.recombine()
#mesh.setOrder(2)
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

# Base geometry {{{
lcMax = 10.0*lcMax
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
else:
    base_bl = geo.addPoint(0.0, -base_h, 0.0, lcMax)
    base_cl = geo.addPoint(0.0,  0.0, 0.0, lcMax)
base_br = geo.addPoint(base_w, -base_h, 0.0, lcMax)
base_cr = geo.addPoint(base_w,  0.0, 0.0, lcMax)

base_b  = geo.addLine(base_bl, base_br)
base_r  = geo.addLine(base_br, base_cr)
base_c = geo.addLine(base_cr, base_cl)
base_l  = geo.addLine(base_cl, base_bl)

base_loop = geo.addCurveLoop([base_b, base_r, base_c, base_l])

base_face = geo.addPlaneSurface([base_loop])
# Synchronize
geo.synchronize()
# }}}

# Create fields for mesh size {{{
# Make box fields for the base geometry
if fullGeo:
    bsFields = rgmsh.BoxFieldsIncreasingSize(mesh, 1, 0.0,
            0.0, fine_x, lcMin*lengthRatio, lcMin, lcMax, 1.05,
            1.2, 1.25, 40, 'cc')
else:
    bsFields = rgmsh.BoxFieldsIncreasingSize(mesh, 1, 0.0,
            0.0, fine_x, lcMin*lengthRatio, lcMin, lcMax, 1.05,
            1.2, 1.25, 40, 'tl')
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
base_fixed_x = model.addPhysicalGroup(1, [base_l, base_r])
model.setPhysicalName(1, base_fixed_x, base_x0_Name)
base_fixed_y = model.addPhysicalGroup(1, [base_b])
model.setPhysicalName(1, base_fixed_y, base_y0_Name)
base_cont = model.addPhysicalGroup(1, [base_c])
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

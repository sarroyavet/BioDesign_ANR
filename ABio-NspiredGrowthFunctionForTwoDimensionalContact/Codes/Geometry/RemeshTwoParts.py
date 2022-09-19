# Libraries {{{
import os
import sys
import json
import getopt
import gmsh
import numpy as np

# In-house modules
sys.path.append('PythonUtilities/')
import AsterStudyUtilities as asus
import ReadGmsh as rgmsh
# }}}

# Get output file name and popup call {{{
try:
    opts, args = getopt.getopt(sys.argv[1:], 'n:pi:')
    if len(args)>0:
        raise getopt.GetoptError('Too many arguments!')
except getopt.GetoptError:
    print('Wrong call! The excecution can have the following arguments:')
    print('To indicate the name use: -n "fileName"')
    print('To indicate popup True use: -p')
    print('To indicate the input parameters: -i "fileName" (.json)')
    raise
inParams = False
export = False
popup = False
for opt, arg in opts:
    # File name
    if opt in ['-n']:
        export = True
        fileName = arg
    elif opt in ['-p']:
        popup = True
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
fileNamesParams = params["fileNames"]
geoParams = meshParams["geoParams"]
workpiece_w = geoParams[0]

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

asterFolder = fileNamesParams["asterFolder"]
recoFolder  = fileNamesParams["recoFolder"]
# }}}

# Organisation of the parameters {{{
# Be careful! All the following items must have consistent Ids. That
# means that $$$FileNams[i], $$$GroNams[i] and $$$Lc[i] together
# describe the ith contour of the $$$ domain.
wpNumGrps = len(wpGroups) - 1
workDir = os.getcwd() + asterFolder + recoFolder
wpFileNams = [None]*wpNumGrps
for k1 in range(wpNumGrps):
    wpFileNams[k1] = workDir + '/' + wpGroups[k1] + '.msh'
wpGrpNams = wpGroups[:-1]
wpLc =  [lcMax]*wpNumGrps
# }}}

# Remesh {{{
# Workpiece geometry {{{
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
# Fields from border and contact box {{{
numFields = 25
# Base distance field
ff = 1
di_fi = ff
mesh.field.add("Distance", ff)
mesh.field.setNumbers(ff, "CurvesList", [3])
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
    mesh.field.setNumber(ff, "XMin"     , -1.0)
    mesh.field.setNumber(ff, "YMin"     , -1.0)
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

# Base geometry {{{
lcMax = 10.0*lcMax
lcMin = 2.0*lcMin
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
base_bl = geo.addPoint(0.0, -base_h, 0.0, lcMax)
base_br = geo.addPoint(base_w, -base_h, 0.0, lcMax)
base_cr = geo.addPoint(base_w,  0.0, 0.0, lcMax)
base_cl = geo.addPoint(0.0,  0.0, 0.0, lcMax)

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

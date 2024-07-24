#!/usr/bin/env python3
# Libraries {{{
import os
import sys
import json
import getopt
import time
import psutil
import gmsh
import numpy as np
from scipy.interpolate import UnivariateSpline
from matplotlib import pyplot as plt
from shapely.geometry import LineString
from scipy import integrate

# In-house modules
srcFolder = os.path.join(os.getcwd(), '../../src/')
sys.path.append(os.path.join(srcFolder, 'PythonUtilities'))
import mesh
from datfile import datfile
# }}}

# Get the file with the parameters {{{
try:
    opts, args = getopt.getopt(sys.argv[1:], 'i:r')
    if len(args)>0:
        raise getopt.GetoptError('Too many arguments!')
except getopt.GetoptError:
    print('Wrong call! The excecution can have the following arguments:')
    print('To indicate the input parameters: -i "fileName" (.json)')
    print('To remesh: -r')
    raise
# Initialisation
inParams = False
newGeo = True
for opt, arg in opts:
    # File name
    if opt in ['-i']:
        parFile = arg
        inParams = True
    # New geometry
    if opt in ['-r']:
        newGeo = False
if not inParams:
    message = "Wrong call! The input file with the parameters is not "
    message += "present. Use -i 'fileName' (.json)"
    raise ValueError(message)
# }}}

# Read parameters {{{
# Default
with open(os.path.join(srcFolder, 'config.json')) as fle:
    defaultParameters = json.load(fle)
locals().update(defaultParameters)
# Specific
with open(parFile) as fle:
    params = json.load(fle)

initialGeometry = params["initialGeometry"]
meshParams = params["mesh"]
fileNames.update(params.get("fileNames", {}))

locals().update(fileNames)
workDir = os.path.dirname(parFile)
# }}}

# Set mesh parameters {{{
sizeMin = meshParams["lcMin"]
sizeGrwRate = meshParams["sizeGrwRate"]
distMin = meshParams.get("distMin", sizeMin*5.0)
distGrwRate = meshParams.get("distGrwRate", sizeGrwRate*5.0)
sizeMax = meshParams.get("lcMax", 10.0*sizeMin)
fine_x_rate = meshParams.get("fine_x_rate", 2.0)
# }}}

if newGeo:
    # Make geometries {{{
    # Make meshes
    for name, geo in initialGeometry.items():
        mainParams = geo["main"]
        auxParams  = geo["aux"]
        main = mesh.Geometry(name, **mainParams)
        aux = mesh.Geometry("aux", **auxParams)
        main.MakeGeometry()
        aux.MakeGeometry(makeFace = False, clearAll = False)
        curveList = [aux.lines[0]["id"]]
        mesh.MeshFieldDistanceCurve(0, curveList, sizeMin, distMin,
                sizeGrwRate, distGrwRate, sizeMax)
        main.MakeMesh()
        main.Export(folder = workDir, fmt = geomFMT)
    # }}}
else:
    # Remesh {{{
    # Get xlims
    meshReport = datfile(workDir + meshSizeFile).datablocks['0'].variables
    xlims = [meshReport["xlim0"][-1], meshReport["xlim1"][-1]]
    xmin = min(xlims)
    xmax = max(xlims)
    xdist = xmax - xmin
    xmin -= xdist*fine_x_rate
    xmax += xdist*fine_x_rate
    # Set old mesh file name
    oldMeshFile = workDir + recoFile
    # Reconstruct and remesh geometry
    for name, geo in initialGeometry.items():
        # Read initial geometry parameters
        mainParams = geo["main"]
        max_x = mainParams.get("max_x", 10.0e+9)
        flip_boundaries = mainParams.get("flip_boundaries", [])
        auxParams  = geo["aux"]
        oldpoints = mainParams["points"]
        oldlines = mainParams["lines"]
        refLine = auxParams["refLine"]
        # Set reference line list
        if isinstance(refLine, str):
            refLine = [refLine]
        # Reconstruct geometries
        main = mesh.Reconstruct(name, oldpoints, oldlines, oldMeshFile)
        cutGeom_points, cutGeom_lines = main.CutGeometry(xmax = max_x,
                flip_boundaries = flip_boundaries)
        #main = mesh.Geometry(name, cutGeom_points, cutGeom_lines)
        # Create cut reference lines
        cutLines = []
        for refLine_i in refLine:
            cutLines.append(main.CutLine(refLine_i, xmin = xmin,
            xmax = xmax))
        # Create auxiliary lines
        auxLines = []
        k_i = 0
        for cutLine_i in cutLines:
            if len(cutLine_i[0]) == 0:
                pass
            else:
                auxLines.append(mesh.Geometry("aux{}".format(k_i), *cutLine_i))
        # Make geometries
        main.MakeGeometry()
        curveList = []
        for aux in auxLines:
            aux.MakeGeometry(makeFace = False, clearAll = False)
            curveList.append(aux.lines[0]["id"])
        # Make mesh
        mesh.MeshFieldDistanceCurve(0, curveList, sizeMin, distMin,
                sizeGrwRate, distGrwRate, sizeMax)
        main.MakeMesh()
        main.Export(folder = workDir, fmt = geomFMT)
    # }}}

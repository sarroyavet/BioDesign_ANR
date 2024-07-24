#!/usr/bin/env python3
# Libraries {{{
import sys
import os
import numpy as np
import copy
import json
import gmsh

# In-house modules
srcFolder = os.path.join(os.getcwd(), '../../src/')
sys.path.append(os.path.join(srcFolder, 'PythonUtilities'))
import mesh
# }}}

# Definition of parameters {{{
RadialForceBaseDir = lambda be_i, st_i: os.getcwd() + \
        "/Tests/radial_force/Be_{}_St_{}/params.json".format(be_i, st_i)
        # "/Tests/Bearing/RadialForce/DOE_Fixed_Base_GrSt_0.625_fine_mesh/Be_{}_St_{}/params.json".format(be_i, st_i)
meshNames = ["esc", "mai"]
flipReocPoints = {"esc" : False,
                  "mai" : True}
height_sign = {"esc" :  1.0,
               "mai" : -1.0}
# Load GROUP_MA
load_group_ma = ["esc_dy"]
# Period if any
load_period = 32.0
# load_period = 16.0
# Time (without cycles)
# cyclic_deltaT = 0.25
fixed_deltaT = 1.0
t_f = 640.0
prints = 64
# Reference St
St_ref = 1.0
# Width and height
width_ref = 0.5
height_ref = 1.0
# Growth St factor
Gr_St_factor = 1.25 # Such that Gr = St*Gr_St_factor
# Gr_St_factor = 3.125 # Such that Gr = St*Gr_St_factor
# Gr_St_factor = 6.25 # Such that Gr = St*Gr_St_factor
# }}}

# Design matrix {{{
designMatrix = np.array([
    # Be, St, Mf, cycle, width reference, maxDis, deltaT
    [0.5, 0.1, 0.5, -1.0, 0.6, 0.025, 1.0],
    [0.5, 1.0, 0.5, -1.0, 0.6, 0.025, 1.0],
    [0.5, 10.0, 0.5,-1.0, 0.6, 0.025, 1.0],
    [0.5, 0.1, 0.5, 1.0, 0.6,  2.5, 1.0],
    [0.5, 1.0, 0.5, 1.0, 0.6,  2.5, 1.0],
    [0.5, 10.0, 0.5, 1.0, 0.6, 2.5, 1.0],
    ])
# }}}

# Intermediary functions {{{
# Get initial geometry {{{
def GetInitialGeometry(parFile, H, W):
    # Get working directory (of the parameter file)
    workDir = os.path.dirname(parFile)
    # Read parameters {{{
    # Default
    with open(os.path.join(srcFolder, 'config.json')) as fle:
        defaultParameters = json.load(fle)
    fileNames = defaultParameters["fileNames"]
    # Specific
    with open(parFile) as fle:
        params = json.load(fle)
    initialGeometry = params["initialGeometry"]
    fileNames.update(params.get("fileNames", {}))
    # Old mesh file
    recoFile = fileNames["recoFile"]
    oldMeshFile = workDir + recoFile
    # }}}
    # Reconstruct each geometry {{{
    geomDict = dict()
    for name in meshNames:
        # Get previous initial geometry {{{
        geo = initialGeometry[name]
        # Read initial geometry parameters
        mainParams = geo["main"]
        auxParams  = geo["aux"]
        oldpoints = mainParams["points"]
        oldlines = mainParams["lines"]
        refLine = auxParams["refLine"]
        # }}}
        # Make and cut reference contact zone {{{
        # Reconstruct full geometry
        main = mesh.Reconstruct(name, oldpoints, oldlines, oldMeshFile)
        # Cut and make reference contact zone
        aux = mesh.Geometry("aux", *main.CutLine(refLine[0], xmax = W))
        # Get reconstructed contact points
        recoPoints = aux.points
        if flipReocPoints[name]:
            recoPoints = np.flip(recoPoints, 0)
        # Set x0 = 0 and y0 =0
        recoPoints[:, 0] -= recoPoints[0, 0]
        recoPoints[:, 1] -= recoPoints[0, 1]
        # }}}
        # Make new contact points {{{
        newPoints = []
        order_cont_l = []
        order_aux = []
        id_i = 0
        # Make the left contour
        for k1 in range(len(recoPoints) - 1):
            xi = recoPoints[-(1 + k1), 0]
            # Scale y
            yi = recoPoints[-(1 + k1), 1]*St_ref
            # Add new point and assign it an id
            newPoints.append([-xi, yi])
            order_cont_l.append(id_i)
            order_aux.append(id_i)
            id_i += 1
        order_cont_l.append(id_i)
        order_cont_r = []
        # Make the right contour
        for k1 in range(len(recoPoints)):
            xi = recoPoints[k1, 0]
            # Scale y
            yi = recoPoints[k1, 1]*St_ref
            # Add new point and assign it an id
            newPoints.append([xi, yi])
            order_cont_r.append(id_i)
            order_aux.append(id_i)
            id_i += 1
        # }}}
        # Make other lines {{{
        H_i = H*height_sign[name]
        # Make right-hand-side line
        order_dxR = [id_i - 1, id_i]
        newPoints.append([W, H_i])
        id_i += 1
        # Make top/bottom line
        order_dy = [id_i - 1, id_i]
        newPoints.append([-W, H_i])
        # Make left-hand-side line
        order_dxL = [id_i, 0]
        # }}}
        # Line setup {{{
        newLines = [{"physical_name" : name + "_conL",
                     "type" : "bspline",
                     "ordered_points" : order_cont_l},
                    {"physical_name" : name + "_conR",
                     "type" : "bspline",
                     "ordered_points" : order_cont_r},
                    {"physical_name" : name + "_dxR",
                     "type" : "straight",
                     "ordered_points" : order_dxR},
                    {"physical_name" : name + "_dy",
                     "type" : "straight",
                     "ordered_points" : order_dy},
                    {"physical_name" : name + "_dxL",
                     "type" : "straight",
                     "ordered_points" : order_dxL}]
        auxLines = [{"physical_name" : name + "_aux",
                     "type" : "bspline",
                     "ordered_points" : order_aux}]
        # }}}
        # Make dictionary {{{
        geomDict[name] = {"main" : {"points" : newPoints,
                                    "lines" : newLines},
                          "aux" : {"points" : newPoints,
                                   "lines" : auxLines,
                                   "refLine" : [name + "_conL",
                                                name + "_conR"]}}
        # }}}
    # }}}
    return geomDict
# }}}
# Contact load {{{
def ContactLoad(St, Mf, period):
    load = {"FY" : -St,
            "M" : St*Mf/2.0,
            "GROUP_MA" : load_group_ma}
    if period == 1:
        # load["FACTOR_FUNCTION"] = "lambda t, PERIOD: np.sin(2.0*np.pi*t/PERIOD)"
        load["FACTOR_FUNCTION"] = "lambda t, PERIOD : np.concatenate((np.linspace(-1.0, 1.0, int(PERIOD/2)), np.linspace(1.0, -1.0, int(PERIOD/2))))[int(np.floor(PERIOD*(t/PERIOD - np.floor(t/PERIOD))))]"
        # load["FACTOR_FUNCTION"] = "lambda t, PERIOD: np.linspace(-1.0, 1.0, PERIOD)[int(np.floor(PERIOD*(t/PERIOD - 0.5 - np.floor(t/PERIOD - 0.5))))]"
        load["FACTOR_PARAMETERS"] = {"PERIOD" : load_period}
    elif period == -1:
        pass
    else:
        message  = "Invalid period input. Use 1 to add the cyclic function"
        message += " and -1 to not to add."
        raise ValueError(message)
    return [load]
# }}}
# Time parameters {{{
def TimeParameters(period, **kwargs):
    timeParams = {"final" : t_f,
                  "prints" : prints}
    if period == 1:
        timeParams["deltaT"] = kwargs["cyclic_deltaT"]
    elif period == -1:
        timeParams["deltaT"] = fixed_deltaT
    else:
        message  = "Invalid period input. Use 1 to add the smaller deltaT"
        message += " and -1 to not to add."
        raise ValueError(message)
    return timeParams
# }}}
# }}}

# Main function {{{
def ParameterCases(baseParameters, unit0, asterStudyBase):
    # Get number of cases
    numCases = designMatrix.shape[0]
    # Base parameter shallow to copy
    shallow_params = copy.deepcopy(baseParameters)
    # Make parameter dictionary for each case {{{
    cases = []
    u = unit0
    for k1 in range(numCases):
        # Get case parameters
        Be_i, St_i, Mf_i, cycle_i, width_ref_i, maxDis_i, deltaT_i = designMatrix[k1]
        # Initialisation
        params_ij = copy.deepcopy(baseParameters)
        initialGeometry_ij = GetInitialGeometry(RadialForceBaseDir(Be_i, St_i),
                height_ref, width_ref_i)
        contactLoad_ij = ContactLoad(St_i, Mf_i, cycle_i)
        timeParams_ij = TimeParameters(cycle_i, cyclic_deltaT = deltaT_i)
        # Add initial geometry
        params_ij["initialGeometry"] = initialGeometry_ij
        # Aster folder
        asterFolder_ij = os.path.join(asterStudyBase,
                "Be_{}_St_{}_Mf_{}_Cy_{}_widRef_{}_maxDis_{}_deltaT_{}".format(Be_i,
                    St_i, Mf_i, cycle_i, width_ref_i, maxDis_i, deltaT_i))
        params_ij["fileNames"]["asterFolder"] = asterFolder_ij
        # Add time parameters
        params_ij["model"]["timeParams"] = timeParams_ij
        # Add contact load
        params_ij["boundary_condition"]["contactLoad"] = contactLoad_ij
        # Add Gr
        params_ij["model"]["growthParams"]["Gr"] = St_i*Gr_St_factor
        # Add maxDis
        params_ij["model"]["growthParams"]["maxDis"] = maxDis_i
        # Add read/write unit
        params_ij["c_a_unit"] = u
        #
        cases.append(params_ij)
        u += 1
    # }}}
    return cases
# }}}

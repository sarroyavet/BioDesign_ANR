# Libraries {{{
import sys
import os
import numpy as np
import copy
# }}}

# Definition of parameters {{{
# Design variables min and max range
# Ge = w/lf
Ge_default = 1.0
# Be = l0/lf array
Be_arr = np.array([0.2, 0.8])
# Sg = lg/lf
Sg_arr = 10.0
# height (circle) and width (rectangle) factor
h_factor = 2.0
w_factor = 3.0
Gr_St_factor = 0.625 # Such that Gr = St*Gr_St_factor
numContEles = 50.0
refLcMin = Be_arr[0]/numContEles # 0.1 = min(Be/1.0), Ge = 1.0
refLcMax = 40.0*refLcMin
maxDepl_factor = 50.0
a_dot_ref = refLcMin/3.0
# }}}

# Design matrix {{{
# Order of variables: Be, St
designMatrix = np.array([
        [0.5, 0.1],
        [0.5, 1.0],
        [0.5, 10.0],
        ])
# }}}

# Intermediary functions {{{
# Radius from St and Be {{{
def Radius(St, Be, Ge, nu):
    return (1000.0*np.pi/32.0)*(Be**2.0)/(St*Ge*(1.0 - nu**2.0))
# }}}
# Circle bottom y {{{
def CircleBottomY(x, r, cx, cy):
    return cy - np.sqrt(r**2.0 - (x - cx)**2.0)
# }}}
# Make points {{{
def MakePoints(r, Ge, Be, h_fact, w_fact):
    # Compute length
    lfstar = 1.0/Ge
    l0star = Be/Ge
    # Get circle point based on its size and on the length
    if r >= lfstar:
        circ_x = lfstar
        circ_y = CircleBottomY(circ_x, r, 0.0, r)
    else:
        circ_x = r
        circ_y = r
    # Make points
    cir_points = [[0.0, 0.0],
                  [circ_x, circ_y],
                  [circ_x, h_fact*lfstar],
                  [0.0, h_fact*lfstar],
                  [0.0, r]]
    rec_points = [[0.0, 0.0],
                  [w_fact*lfstar, 0.0],
                  [w_fact*lfstar, -w_fact*lfstar],
                  [0.0, -w_fact*lfstar]]
    aux_points = [[0.0, 0.0],
                  [1.1*l0star, 0.0]]
    dict_points = {"esc" : cir_points,
                   "mai" : rec_points,
                   "aux" : aux_points}
    return dict_points
# }}}
# }}}

# Main function {{{
def ParameterCases(baseParameters, unit0, asterStudyBase):
    # Get number of cases
    numCases = designMatrix.shape[0]
    # Get nu
    nu = baseParameters["model"]["materials"]["nu"]
    # Set common entries of initial geometry {{{
    initialGeometry = {
            "esc" : {
                "main" : {
                    "lines" : [{"physical_name" : "esc_con",
                                "type" : "circlearc",
                                "ordered_points" : [0, 4, 1]},
                               {"physical_name" : "esc_dxR",
                                "type" : "straight",
                                "ordered_points" : [1, 2]},
                               {"physical_name" : "esc_dy",
                                "type" : "straight",
                                "ordered_points" : [2, 3]},
                               {"physical_name" : "esc_dxM",
                                "type" : "straight",
                                "ordered_points" : [3, 0]}]},
                "aux" : {
                    "lines" : [{"physical_name" : "esc_aux",
                                "type" : "straight",
                                "ordered_points" : [0, 1]}],
                    "refLine" : ["esc_con"]}},
            "mai" : {
                "flip_boundaries" : ["mai_con"],
                "main" : {
                    "lines" : [{"physical_name" : "mai_con",
                                "type" : "straight",
                                "ordered_points" : [0, 1]},
                               {"physical_name" : "mai_dxR",
                                "type" : "straight",
                                "ordered_points" : [1, 2]},
                               {"physical_name" : "mai_dy",
                                "type" : "straight",
                                "ordered_points" : [2, 3]},
                               {"physical_name" : "mai_dxM",
                                "type" : "straight",
                                "ordered_points" : [3, 0]}]},
                "aux" : {
                    "lines" : [{"physical_name" : "mai_aux",
                                "type" : "straight",
                                "ordered_points" : [0, 1]}],
                    "refLine" : ["mai_con"]}}
                }
    contactLoad = [{"M" : 0.0,
                    "GROUP_MA" : ["esc_dy"]}]
    # }}}
    # Base parameter shallow to copy
    shallow_params = copy.deepcopy(baseParameters)
    shallow_initial = copy.deepcopy(initialGeometry)
    # Make parameter dictionary for each case {{{
    cases = []
    u = unit0
    for k1 in range(numCases):
        Ge_i = Ge_default
        Be_i, St_i = designMatrix[k1]
        l0star = Be_i/Ge_i
        lfstar = 1.0/Ge_i
        # Initialisation
        params_ij = copy.deepcopy(baseParameters)
        initialGeometry_ij = copy.deepcopy(initialGeometry)
        contactLoad_ij = copy.deepcopy(contactLoad)
        # Get points
        r_ij = Radius(St_i, Be_i, Ge_i, nu)
        dict_points_ij = MakePoints(r_ij, Ge_i, Be_i, h_factor,
                w_factor)
        # Add points to initial geometry
        initialGeometry_ij["esc"]["main"]["points"] = dict_points_ij["esc"]
        initialGeometry_ij["esc"]["aux"]["points"] = dict_points_ij["aux"]
        initialGeometry_ij["esc"]["max_x"] = lfstar*1.1
        initialGeometry_ij["mai"]["main"]["points"] = dict_points_ij["mai"]
        initialGeometry_ij["mai"]["aux"]["points"] = dict_points_ij["aux"]
        initialGeometry_ij["mai"]["max_x"] = w_factor*lfstar*1.1
        # Add initial geometry
        params_ij["initialGeometry"] = initialGeometry_ij
        # Aster folder
        asterFolder_ij = os.path.join(asterStudyBase,
                "Be_{}_St_{}".format(Be_i, St_i))
        params_ij["fileNames"]["asterFolder"] = asterFolder_ij
        # Add contact load
        contactLoad_ij[0]["FY"] = -(St_i/2.0)/Ge_i # Divided 2.0 because of the symmetry
        params_ij["boundary_condition"]["contactLoad"] = contactLoad_ij
        # Add lcMin
        params_ij["mesh"]["lcMax"] = refLcMax/Ge_i
        params_ij["mesh"]["lcMin"] = refLcMin/Ge_i
        # Add Gr
        params_ij["model"]["growthParams"]["Gr"] = St_i*Gr_St_factor
        # Add maxDis
        maxDis = Sg_arr/Ge_i
        params_ij["model"]["growthParams"]["maxDis"] = maxDis
        # Add a_f which refers to the final length/area
        params_ij["model"]["growthParams"]["a_f"] = 0.5/Ge_i
        # Add a_f maxDepl
        params_ij["model"]["algoParams"]["maxDepl"] = refLcMin*St_i/(Ge_i*maxDepl_factor)
        params_ij["model"]["algoParams"]["a_dot_ref"] = a_dot_ref
        # Add read/write unit
        params_ij["c_a_unit"] = u
        #
        cases.append(params_ij)
        u += 1
    # }}}
    return cases
# }}}

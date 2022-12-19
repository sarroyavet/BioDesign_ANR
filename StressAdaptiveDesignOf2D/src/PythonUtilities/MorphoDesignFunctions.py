# Libraries {{{
import math
import numpy as np
from NormAnalysis import LpNorm_1D
# }}}

# Functions {{{
# Invariants of a matrix {{{
def firstInv(mat):
    # Test square matrix
    if not((mat.ndim == 2) and (len(mat)**2 == mat.size)):
        raise('Error: the given matrix is not square!')
    # Initialisation
    inv1 = 0.0
    # Loop for computing (s_{ii})
    dim = int(math.sqrt(mat.size))
    for i in range(dim):
        inv1 += mat[i][i]
    return inv1
def secondInv(mat):
    # Test square matrix
    if not((mat.ndim == 2) and (len(mat)**2 == mat.size)):
        raise('Error: the given matrix is not square!')
    # Initialisation
    inv2 = 0.0
    # Loop for computing (s_{ii}*s_{jj} - s_{ij}*s_{ji})/2.0
    dim = int(math.sqrt(mat.size))
    for i in range(dim):
        for j in range(dim):
            inv2 += (mat[i][i]*mat[j][j] -
                    mat[i][j]*mat[j][i])
    inv2 *= 0.5
    return inv2
# }}}

# Octahedral stresses in plane strain {{{
def Shyd_pln_strain(SIXX, SIYY, SIZZ):
    return (SIXX + SIYY + SIZZ)/3.0
def Soct_pln_strain(SIXX, SIYY, SIZZ, SIXY):
    mat = np.array([[SIXX, SIXY, 0.0], [SIXY, SIYY, 0.0], [0.0, 0.0,
        SIZZ]])
    I1 = firstInv(mat)
    I2 = secondInv(mat)
    J2 = (I1*I1/3.0 - I2)
    return math.sqrt(2.0*J2/3.0)
def Shyd_pln_stress(SIXX, SIYY):
    return (SIXX + SIYY + 0.0)/3.0
def Soct_pln_stress(SIXX, SIYY, SIXY):
    mat = np.array([[SIXX, SIXY, 0.0], [SIXY, SIYY, 0.0], [0.0, 0.0,
        0.0]])
    I1 = firstInv(mat)
    I2 = secondInv(mat)
    J2 = (I1*I1/3.0 - I2)
    return math.sqrt(2.0*J2/3.0)
# }}}

# Von Mises {{{
def Svmis_pln_strain(SIXX, SIYY, SIZZ, SIXY):
    mat = np.array([[SIXX, SIXY, 0.0], [SIXY, SIYY, 0.0], [0.0, 0.0,
        SIZZ]])
    I1 = firstInv(mat)
    I2 = secondInv(mat)
    J2 = (I1*I1/3.0 - I2)
    return math.sqrt(3.0*J2)
def Svmis_pln_stress(SIXX, SIYY, SIXY):
    mat = np.array([[SIXX, SIXY, 0.0], [SIXY, SIYY, 0.0], [0.0, 0.0,
        0.0]])
    I1 = firstInv(mat)
    I2 = secondInv(mat)
    J2 = (I1*I1/3.0 - I2)
    return math.sqrt(3.0*J2)
# }}}

# Murray terms {{{
def fn1(x, a, b):
    return a/(x + b)
def fn2(x, a, b):
    return a*x/(x*x + b*b)
# }}}

# Growth stress and soft field {{{
def sig1(sig, tau, c = 0.0):
    sic = -sig if sig < 0.0 else 0.0
    term = sic - c*tau
    if term < 0.0:
        term = 0.0
    return term
def sig2(sig, tau, n = 1.0):
    return abs(tau*sig)
def f_nli(s, b, n = 1.0):
    return s**n/(s**n + b**n)
def Sg_nli(sig, tau, a, b1, b2, c = 0, n = 1.0):
    s1 = sig1(sig, tau, c)
    s2 = sig2(sig, tau)
    return f_nli(s1, b1, n) - a*f_nli(s2, b2, n)
#    return 1.0 - a*f_nli(s2, b2)
def Sg_lin(sig, tau, r, c = 0.0, n = 1.0):
    s1 = sig1(sig, tau, c)
    s2 = sig2(sig, tau, n)
    return s1 - r*s2
def Sg_Ss_Sa(sig, tau, maxTau, alpha, beta, gamma, n, tlim, Rp, vel = 1.0):
    # Compressive hydrostatic stress
    sc = sig1(sig, tau)
    # Shear stress (Von Mises)
    sv = sig2(sig, tau, n)
    # Area term (look for sv = tlim)
    sa = (2.0*Sigmoid(maxTau, vel = vel, cen = tlim) - 1.0)
    # Specialisation (homogenise contact pressure)
    ss = (1.0 - Rp)*(sc - gamma*sv)
    return alpha*(ss + beta*sa)
def Sg_Ss_Sa_t1(sig, tau, maxTau, alpha, beta, tlim, Rp, vel = 1.0):
    # Compressive hydrostatic stress
    sc = sig1(sig, 0.0)
    # Area term (look for sv = tlim)
    sa = (2.0*Sigmoid(maxTau, vel = vel, cen = tlim) - 1.0)
    # Specialisation (homogenise contact pressure)
    ss = abs((1.0 - Rp)*sc)
    return alpha*(-ss + beta*sa)
def Sgrowth(Shyd, Sshr, alpha, shrlim = 0, hydlim = 0.0,
        fun = 'sigmoid', vel = 20.0):
    if fun == 'heaviside':
        fun = Heaviside
    elif fun == 'sigmoid':
        fun = lambda x: Sigmoid(x, vel = vel)
    else:
        raise("Error: fun must be either heaviside or sigmoid.")
    sig = -Shyd if Shyd < 0.0 else 0.0
    tau = Sshr
    sa = alpha*fun(shrlim - tau)*fun(hydlim - sig)*sig
    return sa
def SgrowthDesign(Shyd, Sshr, alpha, Rp, Rv, eta, shrlim = 0, hydlim = 0.0,
        fun = 'sigmoid', vel = 20.0):
    sa = Sgrowth(Shyd, Sshr, 1.0, shrlim = shrlim, hydlim = hydlim,
            fun = fun, vel = vel)
    sg = alpha*(sa*(1.0 - Rp) + eta*(1.0 - Rv))
    return sg
def Sgrowth_flag(Shyd, Sshr, shrlim = 0, hydlim = 0.0,
        fun = 'heaviside', vel = 1.0):
    if fun == 'heaviside':
        fun = Heaviside
    elif fun == 'sigmoid':
        fun = lambda x: Sigmoid(x, vel = vel)
    else:
        raise("Error: fun must be either heaviside or sigmoid.")
    sa = fun(shrlim - Sshr)*fun(hydlim - Shyd)
    return sa
# }}}

# Sigmoid and step function # {{{
def Sigmoid(x, vel = 1.0, cen = 0.0):
    try:
        term = math.exp(-vel*(x - cen))
    except OverflowError:
        term = float('inf')
    return 1.0/(1.0 + term)
def Heaviside(x):
    if x > 0:
        return 1.0
    else:
        return 0.0
def SigmoidThresholds(x, thr1, thr2, per = 0.01):
    """A Sigmoid function defined by two points: (thr1, 1 -per ) and
    (thr2, per). f(x) = 1/(1 + exp(v(x - c)))"""
    # Tests inputs {{{
    if thr1 == thr2:
        raise("Error: thr1 and thr2 must be different.")
    if per <= 0.0 or per >= 0.5:
        raise("Error: 0.0 < per < 0.5.")
    if thr1 > thr2:
        thr1, thr2 = thr2, thr1
    # }}}
    # Set points
    bt = 1.0 - per
    bb = per
    v, c = FindSigThrPars(bt, bb, thr1, thr2)
    # Compute f(x)
    val = v*(x - c)
    if val > 709.0:
        return 0.0
    else:
        return 1.0/(1.0 + math.exp(val))
def FindSigThrPars(bt, bb, thr1, thr2):
    # Find c
    termt = math.log(1.0/bt - 1.0)
    termb = math.log(1.0/bb - 1.0)
    c = (termt*thr2 - termb*thr1)/(termt - termb)
    # Find v
    v = termt/(thr1 - c)
    return v, c
# }}}

# One-dimensional linear element methods {{{
# 1D linear element size
def Size1DLinearElement(nod1, nod2):
    x1 = nod1[0]
    y1 = nod1[1]
    x2 = nod2[0]
    y2 = nod2[1]
    return math.sqrt((x2 - x1)**2.0 + (y2 - y1)**2.0)
# Maximum 1D length in a group of elements
# mail_py is a mesh object of the type mail_py and groupname is the
# name of the group of one-dimensional elements whose maximun size
# will be calculated.
def Maximum1DLengthInElementGroup(mail_py, groupName):
    # Get mesh info
    grpEles = mail_py.gma[groupName]
    eles = mail_py.co
    nodes = mail_py.cn
    # Get maximum size
    numEles = len(grpEles)
    maxSize = 0.0
    for eleId in grpEles:
        eleNods = eles[eleId]
        if not(len(eleNods) == 2):
            raise("Error: the element must have exactly two nodes. " +
                    "It has " + str(len(eleNods)))
        nod1 = [nodes[eleNods[0]][0], nodes[eleNods[0]][1]]
        nod2 = [nodes[eleNods[1]][0], nodes[eleNods[1]][1]]
        eleSize = Size1DLinearElement(nod1, nod2)
        if eleSize > maxSize:
            maxSize = eleSize
    return maxSize
# Length of a 1D mesh group
# mail_py is a mesh object of the type mail_py and groupname is the
# name of the group of one-dimensional elements whose maximun size
# will be calculated.
def Length1DElementGroup(mail_py, groupName):
    # Get mesh info
    grpEles = mail_py.gma[groupName]
    eles = mail_py.co
    nodes = mail_py.cn
    # Get maximum size
    numEles = len(grpEles)
    length = 0.0
    for eleId in grpEles:
        eleNods = eles[eleId]
        if not(len(eleNods) == 2):
            raise("Error: the element must have exactly two nodes. " +
                    "It has " + str(len(eleNods)))
        nod1 = [nodes[eleNods[0]][0], nodes[eleNods[0]][1]]
        nod2 = [nodes[eleNods[1]][0], nodes[eleNods[1]][1]]
        eleSize = Size1DLinearElement(nod1, nod2)
        length += eleSize
    return length
# }}}

# Contact area from ordered data{{{
def ContactArea(absc_curv, lags_c, zero = 1e-6):
    """Sum area as long as lags_c > zero in the element."""
    # Test input variables {{{
    if not (isinstance(absc_curv, np.ndarray)):
        raise("Error: absc_curv must be a numpy array.")
    if not (absc_curv.ndim == 1):
        raise("Error: absc_curv.ndim must be 1.")
    if not (isinstance(lags_c, np.ndarray)):
        raise("Error: lags_c must be a numpy array.")
    if not (lags_c.ndim == 1):
        raise("Error: lags_c.ndim must be 1.")
    if not (absc_curv.size == lags_c.size):
        raise("Error: the size of absc_curv and lags_c must be equal.")
    # }}}
    numNods = absc_curv.size
    area = 0.0
    for k1 in range(numNods - 1):
        p1 = lags_c[k1]
        p2 = lags_c[k1 + 1]
        p = 0.5*(p1 + p2)
        if p > zero:
            x1 = absc_curv[k1]
            x2 = absc_curv[k1 + 1]
            area += abs(x2 - x1)
    return area
def IntegralToxFromOrderedData(x, f):
    # Test numpy arrays
    if not (isinstance(x, np.ndarray) and
            isinstance(f, np.ndarray)):
        raise("Error: x and f must be numpy arrays.")
    if not (x.ndim == 1 and f.ndim == 1):
        raise("Error: x and f must be 1D numpy arrays.")
    if not (x.shape[0] == f.shape[0]):
        raise("Error: x and f must have the same length.")
    length = x.shape[0]
    g = np.zeros(length)
    for k1 in range(length - 1):
        x1 = x[k1]
        x2 = x[k1 + 1]
        f1 = f[k1]
        f2 = f[k1 + 1]
        g[k1 + 1] = 0.5*(f1 + f2)*(x2 - x1)
    return g
# }}}

# Maximum length of in-contact elements {{{
def MaxLengthInContactElements(absc_curv, lags_c, zero = 1e-6):
    """Take into account only the elements with lags_c > zero."""
    # Test input variables {{{
    if not (isinstance(absc_curv, np.ndarray)):
        raise("Error: absc_curv must be a numpy array.")
    if not (absc_curv.ndim == 1):
        raise("Error: absc_curv.ndim must be 1.")
    if not (isinstance(lags_c, np.ndarray)):
        raise("Error: lags_c must be a numpy array.")
    if not (lags_c.ndim == 1):
        raise("Error: lags_c.ndim must be 1.")
    if not (absc_curv.size == lags_c.size):
        raise("Error: the size of absc_curv and lags_c must be equal.")
    # }}}
    numNods = absc_curv.size
    maxLength = 0.0
    for k1 in range(numNods - 1):
        p1 = lags_c[k1]
        p2 = lags_c[k1 + 1]
        p = 0.5*(p1 + p2)
        if p > zero:
            x1 = absc_curv[k1]
            x2 = absc_curv[k1 + 1]
            length = abs(x2 - x1)
            if length > maxLength:
                maxLength = length
    return maxLength
# }}}

# Pressure weighted variance {{{
def PressureWeightedVariance(absc_curv, lags_c,
        zero = 1e-6):
    """Weighted variance of nonzero pressure. The weight of each value
    is given by the size of the elements"""
    # Test input variables {{{
    if not (isinstance(absc_curv, np.ndarray)):
        raise("Error: absc_curv must be a numpy array.")
    if not (absc_curv.ndim == 1):
        raise("Error: absc_curv.ndim must be 1.")
    if not (isinstance(lags_c, np.ndarray)):
        raise("Error: lags_c must be a numpy array.")
    if not (lags_c.ndim == 1):
        raise("Error: lags_c.ndim must be 1.")
    if not (absc_curv.size == lags_c.size):
        raise("Error: the size of absc_curv and lags_c must be equal.")
    # }}}
    # Isolate the nonzero elements and their pressure {{{
    lags_c = abs(lags_c)
    numNods = absc_curv.size
    weight = []
    press = []
    for k1 in range(numNods - 1):
        p1 = lags_c[k1]
        p2 = lags_c[k1 + 1]
        p = 0.5*(p1 + p2)
        if p > zero:
            x1 = absc_curv[k1]
            x2 = absc_curv[k1 + 1]
            press.append(p)
            weight.append(abs(x2 - x1))
    weight = np.array(weight)
    press = np.array(press)
    # }}}
    # Normalise the weight {{{
    size = weight.sum()
    weight = np.array([wei/size for wei in weight])
    # }}}
    # Compute weighted variance {{{
    average = np.average(press, weights = weight)
    variance = np.average((press - average)**2.0, weights = weight)
    # }}}
    return (variance, average)
# }}}

# Lundberg profile {{{
def LundProfile(x, b, St, nu = 0.3):
    fac = 0.1*(1.0 - nu*nu)/(math.pi*b*St)
    try:
        term = math.log(1.0 - (x/b)**2.0)
    except ValueError:
        print(x, b)
        print((x/b)**2.0)
        print(1.0 - (x/b)**2.0)
        raise
    return -fac*term
# }}}

# Contact area from CONT_ELEM 1D {{{
def ContactAreaCONT_ELEM1D(mail_py, cont_elem):
    # Set up input variables
    cont = cont_elem.valeurs
    maille = cont_elem.maille
    nodes = mail_py.cn[:, :2]
    eles = mail_py.co
    contSet = np.array([mm - 1 for mm in maille])
    # Compute area
    area = 0.0
    for k1 in range(contSet.size):
        cc = cont[k1]
        if cc > 0.0:
            ele = eles[contSet[k1]]
            node1 = nodes[ele[0]]
            node2 = nodes[ele[1]]
            dx = node2[0] - node1[0]
            dy = node2[1] - node1[1]
            area += math.sqrt(dx*dx + dy*dy)
    return area
# }}}

# Active nodes from CONT_ELEM 1D {{{
def ActiveNodesCONT_ELEM1D(mail_py, cont_elem):
    # Set up input variables
    cont = cont_elem.valeurs
    maille = cont_elem.maille
    nodes = mail_py.cn[:, :2]
    eles = mail_py.co
    contSet = np.array([mm - 1 for mm in maille])
    # Active nodes
    activeNodeIds = []
    for k1 in range(contSet.size):
        cc = cont[k1]
        if cc > 0.0:
            ele = eles[contSet[k1]]
            activeNodeIds.append(ele[0])
            activeNodeIds.append(ele[1])
    activeNodeIds = list(set(activeNodeIds))
    activeNodes = [None]*len(activeNodeIds)
    for k1 in range(len(activeNodeIds)):
        ID = activeNodeIds[k1]
        xx = nodes[ID, 0]
        yy = nodes[ID, 1]
        activeNodes[k1] = [ID, xx, yy]
    return np.array(activeNodes)
# }}}

# Active nodes from JEU {{{
def ActiveNodesJEU(mail_py, jeu, maxJeu = 0.001):
    # Set up input variables
    cont = jeu.valeurs
    maille = jeu.maille
    nodes = mail_py.cn[:, :2]
    eles = mail_py.co
    contSet = np.array([mm - 1 for mm in maille])
    # Active nodes
    activeNodeIds = []
    for k1 in range(contSet.size):
        cc = cont[k1]
        if abs(cc) < maxJeu:
            ele = eles[contSet[k1]]
            activeNodeIds.append(ele[0])
            activeNodeIds.append(ele[1])
    activeNodeIds = list(set(activeNodeIds))
    activeNodes = [None]*len(activeNodeIds)
    for k1 in range(len(activeNodeIds)):
        ID = activeNodeIds[k1]
        xx = nodes[ID, 0]
        yy = nodes[ID, 1]
        activeNodes[k1] = [ID, xx, yy]
    return np.array(activeNodes)
# }}}

# Contact area from CONT_ELEM 1D {{{
def MaxLengthCONT_ELEM1D(mail_py, cont_elem):
    # Set up input variables
    cont = cont_elem.valeurs
    maille = cont_elem.maille
    nodes = mail_py.cn[:, :2]
    eles = mail_py.co
    contSet = np.array([mm - 1 for mm in maille])
    # Compute area
    length = 0.0
    for k1 in range(contSet.size):
        cc = cont[k1]
        if cc > 0.0:
            ele = eles[contSet[k1]]
            node1 = nodes[ele[0]]
            node2 = nodes[ele[1]]
            dx = node2[0] - node1[0]
            dy = node2[1] - node1[1]
            length_d = math.sqrt(dx*dx + dy*dy)
            if length < length_d:
                length = length_d
    return length
# }}}

# Compute Rp from ordered data{{{
def Rp_OrderedData(absc_curv, lags_c, area, load, atEnd = False):
    # Set optimum profile
    lenA = absc_curv.size
    p_o = load/area
    popt = np.zeros(lenA)
    min_curv = absc_curv.min()
    max_curv = absc_curv.max()
    length = 0.0
    if atEnd:
        x0 = absc_curv[-1]
        popt[-1] = p_o
        for k1 in range(lenA - 1):
            x1 = absc_curv[lenA - k1 - 2]
            length += abs(x1 - x0)
            if length > area:
                break
            else:
                popt[lenA - k1 - 2] = p_o
                x0 = x1
    else:
        x0 = absc_curv[0]
        popt[0] = p_o
        for k1 in range(lenA - 1):
            x1 = absc_curv[k1 + 1]
            length += abs(x1 - x0)
            if length > area:
                break
            else:
                popt[k1 + 1] = p_o
                x0 = x1
    # Compute norms
    difNorm = LpNorm_1D(absc_curv, popt - lags_c)
    optNorm = LpNorm_1D(absc_curv, popt)
    appNorm = LpNorm_1D(absc_curv, lags_c)
    # Return Rp
    return 1.0 - difNorm/(optNorm + appNorm)

# }}}
# }}}

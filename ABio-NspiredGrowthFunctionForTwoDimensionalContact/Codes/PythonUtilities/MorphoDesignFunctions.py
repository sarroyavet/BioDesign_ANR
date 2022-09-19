import math
import numpy as np

# Functions
# Lp norm in one--dimensional arrays {{{
def LpNorm_1D(x, f, p = 2.0):
    # Check numpy arrays
    if isinstance(x, np.ndarray) and isinstance(f, np.ndarray):
        pass
    else:
        raise Exception("x and f must be ndarrays!")
    # Check for size consistency
    if x.size == f.size:
        numPoints = x.size
        pass
    else:
        raise Exception("x and f must have the same dimension!")
    # Check x to be sorted
    if np.all(np.diff(x) >= 0.0):
        pass
    else:
        raise Exception("x must be sorted!")
    # Copute integral: int_\Omega (f)^p dx; with \Omega = [min(x),
    # max(x)]
    norm = 0.0
    for k1 in range(numPoints - 1):
        x1 = x[k1]
        x2 = x[k1 + 1]
        f1 = f[k1]
        f2 = f[k1 + 1]
        norm += (((f2 + f1)/2.0)**p)*(x2 - x1)
    # Compute norm
    norm = norm**(1.0/p)
    return norm
# }}}

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

# Hydrostatic stress in plane strain {{{
def Shyd_pln_strain(SIXX, SIYY, SIZZ):
    return (SIXX + SIYY + SIZZ)/3.0
# }}}

# Von Mises {{{
def Svmis_pln_strain(SIXX, SIYY, SIZZ, SIXY):
    mat = np.array([[SIXX, SIXY, 0.0], [SIXY, SIYY, 0.0], [0.0, 0.0,
        SIZZ]])
    I1 = firstInv(mat)
    I2 = secondInv(mat)
    J2 = (I1*I1/3.0 - I2)
    return math.sqrt(3.0*J2)
# }}}

# Growth stress {{{
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
def SgrowthDesign(Shyd, Sshr, alpha, Qp, Rv, eta, shrlim = 0, hydlim = 0.0,
        fun = 'sigmoid', vel = 20.0):
    sa = Sgrowth(Shyd, Sshr, 1.0, shrlim = shrlim, hydlim = hydlim,
            fun = fun, vel = vel)
    sg = alpha*(sa*(1.0 - Qp) + eta*(1.0 - Rv))
    return sg
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

# Compute Qp from ordered data{{{
def Qp_OrderedData(absc_curv, lags_c, area, load, atEnd = False):
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
    # Return Qp
    return 1.0 - difNorm/(optNorm + appNorm)
# }}}

# Libraries {{{
from paraview.simple import *
import json
# }}}

# Initialisation {{{
results = FindSource('resu*')
inp = results
# }}}

# Apply displacements {{{
appDisp = WarpByVector(registrationName='appDisp', Input=inp)
appDisp.Vectors = ['POINTS', 'disp']
inp = appDisp
# }}}

# Normalised stresses {{{
# Shear stress
norTau = PythonCalculator(registrationName='norTau', Input=inp)
inp = norTau
norTau.Expression = 'tauMis/max(tauMis)'
norTau.ArrayName = 'norTau'
# Hydrostatic stress
absSig = PythonCalculator(registrationName='absSig', Input=inp)
inp = absSig
absSig.Expression = 'abs(sigHyd)'
absSig.ArrayName = 'absSig'
norSig = PythonCalculator(registrationName='norSig', Input=inp)
inp = norSig
norSig.Expression = 'sigHyd/abs(min(sigHyd))'
norSig.ArrayName = 'norSig'
# Growth stress
norGrw = PythonCalculator(registrationName='norGrw', Input=inp)
inp = norGrw
norGrw.Expression = '(growthData - min(growthData))/(max(growthData) - min(growthData))'
norGrw.ArrayName = 'norGrw'
# }}}

# Growth displacement {{{
grwDis = PythonCalculator(registrationName='grwDis', Input=inp)
inp = grwDis
grwDis.Expression = 'sqrt(growth_DX**2.0 + growth_DY**2.0)'
grwDis.ArrayName = 'grwDis'
# }}}

# Libraries {{{
from paraview.simple import *
import json
# }}}

# Initialisation {{{
results = FindSource('results*')
inp = results
# }}}

# Apply displacements {{{
appDisp = WarpByVector(registrationName='appDisp', Input=inp)
appDisp.Vectors = ['POINTS', 'disp']
inp = appDisp
# }}}

# Normalised stresses {{{
# Hydrostatic stress
norTau = PythonCalculator(registrationName='norTau', Input=inp)
inp = norTau
norTau.Expression = 'tauMis/max(tauMis)'
norTau.ArrayName = 'norTau'
# Shear stress
norSig = PythonCalculator(registrationName='norSig', Input=inp)
inp = norSig
norSig.Expression = 'sigHyd/abs(min(sigHyd))'
norSig.ArrayName = 'norSig'
# Growth stress
norGrw = PythonCalculator(registrationName='norGrw', Input=inp)
inp = norGrw
norGrw.Expression = 'growthData/abs(max(growthData))'
norGrw.ArrayName = 'norGrw'
# }}}

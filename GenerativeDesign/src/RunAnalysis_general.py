#!/usr/bin/env python3
# Libraries {{{
import os
import sys
import json
import getopt
import time
import psutil
import numpy as np
import multiprocessing

# In-house modules
abspath = os.path.dirname(os.path.abspath(__file__))
sys.path.append(abspath + '/PythonUtilities/')
import AsterStudyUtilities as asus
from datfile import datfile
# }}}

# Start time
startTime = time.time()
# Get the file with the parameters {{{
try:
    opts, args = getopt.getopt(sys.argv[1:], 'i:kcpad')
    if len(args)>0:
        raise getopt.GetoptError('Too many arguments!')
except getopt.GetoptError:
    print('Wrong call! The excecution can have the following arguments:')
    print('To indicate the input parameters: -i "fileName" (.json)')
    print('To indicate to keep the existing aster study folder: -k')
    print('To indicate to continue an aster study folder: -c')
    print('To indicate to run the post-processing: -p')
    print('To indicate to code_aster with the alternative root: -a')
    print('To indicate to use an existing mesh (useful when debugging): -d')
    raise
# Initialisation
inParams = False
remove = True
continuee = False
onlyPost = False
defRoot = True
makeMesh = True
for opt, arg in opts:
    # File name
    if opt in ['-i']:
        parFile = arg
        inParams = True
    if opt in ['-k']:
        remove = False
    if opt in ['-c']:
        continuee = True
    if opt in ['-p']:
        onlyPost = True
    if opt in ['-a']:
        defRoot = False
    if opt in ['-d']:
        makeMesh = False
        remove = False
if not inParams:
    print("Wrong call! The input file with the parameters is not " +
            "present. Use -i 'fileName' (.json)")
    raise
# }}}

# Parameters {{{
# Load default parameter file
srcFolder = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(srcFolder, 'config.json')) as fle:
    defaultParameters = json.load(fle)
locals().update(defaultParameters)
memory_limit_default = psutil.virtual_memory().available*0.1/1000000.0
# Load specific parameter file
with open(parFile) as fle:
    params = json.load(fle)
fileNames.update(params.get("fileNames", {}))
codeParams = params.get("code_aster", {})
modelParams     = params["model"]

locals().update(fileNames)
# }}}

# Set up iteration parameters {{{
timeParams = modelParams["timeParams"]
final   = timeParams["final"]
deltaT  = timeParams["deltaT"]
prints  = timeParams["prints"]
numItes = int(abs(final/deltaT))
if numItes < prints:
    raise ValueError("Error: timeParams are not consistent. There are too many prints.")
if not numItes%prints == 0:
    raise ValueError("The number of prints (" + str(prints) + ") is not a multiple of the number of iterations (" + str(numItes) + ").")
# }}}

# Compute necessary paths {{{
baseDir = os.getcwd()
if asterFolder[:5] == "/home":
    workDir = asterFolder
else:
    workDir = baseDir + asterFolder
if parFile[:5] == '/home':
    paramsPath = parFile
else:
    if continuee or not remove:
        paramsPath = os.path.join(baseDir, parFile)
    else:
        paramsPath = os.path.join(workDir, parFile)

geoGenerator = baseDir + geoGen
# }}}

# Aster study {{{
commFile = dummFile[:-5] + '.comm'
dict_P = {"memory_limit" : memory_limit_default}
dict_P_ext = codeParams.get("dict_P", {})
dict_P.update(dict_P_ext)
dict_F = {"comm" : commFile}
dict_A = {"memjeveux" : 60}
study = asus.AsterStudy(dict_F = dict_F, dict_P = dict_P, dict_A = dict_A,
        defaultRoot = defRoot)
study._workDir = workDir
# }}}

# Initialisation of the growth process {{{
if not continuee:
    # Set up aster study
    study.CreateStudy(deletePrevious = remove, workDir = workDir)
    with open(paramsPath, 'w') as fle:
        json.dump(params, fle, indent = 4)
    os.makedirs(workDir + resuFolder, exist_ok = True)

    dummies = [('#workDir', workDir),
            ('#parFile', paramsPath),
            ('#srcDir', srcFolder)]
    asus.ReplaceDummyFile(srcFolder + '/codeasters' + dummFile,
            workDir + commFile, dummies)
    # Create mesh
    if makeMesh:
        fail = os.system("python3 {} -i {}".format(geoGenerator, paramsPath))
        if fail != 0:
            raise ValueError("Unsuccessful meshing.")
    makeMesh = True
    # Run aster study
    fail = study.RunStudy(outSalome = True)
# }}}

# Continue simulating if necessary {{{
# Get the last time simulated
meshReport = datfile(workDir + meshSizeFile).datablocks['0'].variables
ite = meshReport["ite"][-1]
run = True if ite < numItes else False
# Set prevIte to avoid error loops
prevIte = -1
while run:
    # Rewrite export
    study.CreateStudy(deletePrevious = False, workDir = workDir)
    # Remesh if necessary
    if makeMesh:
        fail = os.system("python3 {} -i {} -r".format(geoGenerator, paramsPath))
        if fail != 0:
            raise ValueError("Unsuccessful remeshing.")
    makeMesh = True
    # Run aster study
    fail = study.RunStudy(outSalome = True)
    if fail:
        if prevIte == ite:
            raise ValueError("Error: unsuccesful aster study, even after remesh.")
        else:
            prevIte = ite
    # Get the last time simulated
    meshReport = datfile(workDir + meshSizeFile).datablocks['0'].variables
    ite = meshReport["ite"][-1]
    run = True if ite < numItes else False
# }}}


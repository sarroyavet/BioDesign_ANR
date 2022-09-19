# Libraries {{{
import os
import sys
import json
import getopt
import time

# In-house modules
sys.path.append('PythonUtilities/')
import AsterStudyUtilities as asus
# }}}

# Start time
startTime = time.time()
# Get the file with the parameters {{{
try:
    opts, args = getopt.getopt(sys.argv[1:], 'i:kc')
    if len(args)>0:
        raise getopt.GetoptError('Too many arguments!')
except getopt.GetoptError:
    print('Wrong call! The excecution can have the following arguments:')
    print('To indicate the input parameters: -i "fileName" (.json)')
    print('To indicate to keep the existing aster study folder: -k')
    print('To indicate to continue an aster study folder: -c')
    raise
# Initialisation
inParams = False
remove = True
continuee = False
for opt, arg in opts:
    # File name
    if opt in ['-i']:
        parFile = arg
        inParams = True
    if opt in ['-k']:
        remove = False
    if opt in ['-c']:
        continuee = True
if not inParams:
    print("Wrong call! The input file with the parameters is not " +
            "present. Use -i 'fileName' (.json)")
    raise
# }}}

# Read parameters {{{
with open(parFile, 'r') as inFile:
    params = json.load(inFile)
fileNamesParams = params["fileNames"]
modelParams = params["model"]

commFile        = fileNamesParams["commFile"]
dummCommFile    = fileNamesParams["dummCommFile"]
asterFolder     = fileNamesParams["asterFolder"]
resuFolder      = fileNamesParams["resuFolder"]
genGeomFile     = fileNamesParams["genGeomFile"]
geomFile        = fileNamesParams["geomFile"]
resuFile        = fileNamesParams["resuFile"]
recoFolder        = fileNamesParams["recoFolder"]
recGeomFile     = fileNamesParams["recGeomFile"]
meshSizeFile = fileNamesParams["meshSizeFile"]

final  = modelParams["final"]
alpha  = modelParams["alpha"]
prints = modelParams["prints"]
numItes = int(abs(final/alpha))
if not numItes%prints == 0:
    numItes += (prints - numItes%prints)
# }}}

# Compute necessary paths {{{
baseDir = os.getcwd()
workDir = baseDir + asterFolder
paramsPath = parFile if parFile[0] == '/' else baseDir + '/' + parFile
resuName, resuExt = os.path.splitext(workDir + resuFolder + resuFile)
# }}}

# Aster study {{{
study = asus.AsterStudy()
study.workDir = workDir
study.comm = commFile
#study.memory_limit = 12000.0
# }}}

# Initialisation of the growth process {{{
if not continuee:
    study.CreateStudy(deletePrevious = remove)
    os.makedirs(workDir + resuFolder, exist_ok = True)
    os.makedirs(workDir + recoFolder, exist_ok = True)

    dummies = [('#workDir', workDir),
            ('#parFile', paramsPath)]
    asus.ReplaceDummyFile(baseDir + dummCommFile, workDir + commFile,
            dummies)
    ##############################
    # First run with original mesh
    ##############################
    # Create mesh
    fail = os.system('python3 ' + baseDir + genGeomFile + ' -n "' +
            baseDir + asterFolder + geomFile + '" -i ' + paramsPath)
    if fail:
        raise("Error: unsuccesful mesh generation")
    # Run aster study
    fail = study.RunStudy(outSalome = True)
    if fail:
        raise("Error: unsuccesful aster study")
# See the number of iterations
with open(workDir + meshSizeFile, 'r') as inFile:
    line = inFile.readlines()
    line = line[-1]
    line = line.split(',')
    line = line[0]
    ite = eval(line)
run = True if ite < numItes else False
# }}}

# Continuation of the growth process if necessary {{{
#####################################################################
# If the system has not reach the total number of iterations, run the
# analysis with the remesh algorithm
#####################################################################
while run:
    # Create mesh
    fail = os.system('python3 ' + baseDir + recGeomFile + ' -n "' +
             baseDir + asterFolder + geomFile +
             '" -i ' + paramsPath)
    if fail:
        raise("Error: unsuccesful remeshing")
    # Run aster study
    fail = study.RunStudy(outSalome = True)
    if fail:
        raise("Error: unsuccesful aster study")
    # Update run flag
    with open(workDir + meshSizeFile, 'r') as inFile:
        line = inFile.readlines()
        line = line[-1]
        line = line.split(',')
        line = line[0]
        ite = eval(line)
    run = True if ite < numItes else False
os.system('cp {} {}'.format(parFile, asterFolder))
# }}}
# End time
endTime = time.time()
print("Execution time: %g" % (endTime - startTime))

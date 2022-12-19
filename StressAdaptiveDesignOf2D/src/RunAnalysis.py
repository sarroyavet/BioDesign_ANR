# Libraries {{{
import os
import sys
import json
import getopt
import time

# In-house modules
sys.path.append('../../src/PythonUtilities/')
import AsterStudyUtilities as asus
from PostProcessing import MorDesResults as mdr
# }}}

# Start time
startTime = time.time()
# Get the file with the parameters {{{
try:
    opts, args = getopt.getopt(sys.argv[1:], 'i:kcp')
    if len(args)>0:
        raise getopt.GetoptError('Too many arguments!')
except getopt.GetoptError:
    print('Wrong call! The excecution can have the following arguments:')
    print('To indicate the input parameters: -i "fileName" (.json)')
    print('To indicate to keep the existing aster study folder: -k')
    print('To indicate to continue an aster study folder: -c')
    print('To indicate to run the post-processing: -p')
    raise
# Initialisation
inParams = False
remove = True
continuee = False
onlyPost = False
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
if not inParams:
    print("Wrong call! The input file with the parameters is not " +
            "present. Use -i 'fileName' (.json)")
    raise
# }}}

# Parameters {{{
# Default {{{
resuFolder      = "/RESULTS"
resuFile        = "/results.vtk"
recoFolder      = "/RECONSTRUCTION"
meshSizeFile =  "/meshSize.tab"
# }}}
# Read parameters {{{
with open(parFile, 'r') as inFile:
    params = json.load(inFile)
fileNamesParams = params["fileNames"]
modelParams = params["model"]
memory_limit = params["memory_limit"]

asterFolder     = fileNamesParams["asterFolder"]
geomFile        = fileNamesParams["geomFile"]
srcFolder       = fileNamesParams["srcFolder"]
genGeomFile     = fileNamesParams["genGeomFile"]
recGeomFile     = fileNamesParams["recGeomFile"]
dummFile        = fileNamesParams["dummFile"]
commFile        = dummFile[:-5] + '.comm'

final  = modelParams["final"]
alpha  = modelParams["alpha"]
prints = modelParams["prints"]
numItes = int(abs(final/alpha))
if not numItes%prints == 0:
    numItes += (prints - numItes%prints)
# }}}
# }}}

# Compute necessary paths {{{
baseDir = os.getcwd()
workDir = baseDir + asterFolder
paramsPath = parFile if parFile[0] == '/' else baseDir + '/' + parFile
resuName, resuExt = os.path.splitext(workDir + resuFolder + resuFile)
# }}}
if not onlyPost:
    # Aster study {{{
    study = asus.AsterStudy()
    study.workDir = workDir
    study.comm = commFile
    study.memory_limit = memory_limit
    # }}}

    # Initialisation of the growth process {{{
    if not continuee:
        study.CreateStudy(deletePrevious = remove)
        os.system('cp {} {}'.format(parFile, workDir))
        os.makedirs(workDir + resuFolder, exist_ok = True)
        os.makedirs(workDir + recoFolder, exist_ok = True)

        dummies = [('#workDir', workDir),
                ('#parFile', paramsPath),
                ('#srcDir', srcFolder)]
        asus.ReplaceDummyFile(srcFolder + '/codeasters' + dummFile,
                workDir + commFile, dummies)
        ##############################
        # First run with original mesh
        ##############################
        # Create mesh
        fail = os.system('python3 ' + srcFolder + '/mesh' + genGeomFile + ' -n "' +
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
    # alanysis with the remesh algorithm
    #####################################################################
    while run:
        # Create mesh
        fail = os.system('python3 ' + srcFolder + '/mesh' + recGeomFile + ' -n "' +
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
    # }}}
    # End time
    endTime = time.time()
    print("Execution time: %g" % (endTime - startTime))

# Post-process {{{
study = mdr(workDir, masterName = "/bas_cont.msh", slaveName = "/wp_cont.msh")
study.add_relativeContour()
cm = 1.0/2.54
study.Compute_Curvature()
study.PlotContContours(figsize = [8*cm, 8*0.75*cm])
study.PlotBestPressure()
study.PlotMeshSize()
study.PlotReport()
def figConf(ax):
    ax.set_xlim(45.0, 55.0)
study.PlotPressure()
# }}}

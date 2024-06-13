#!/usr/bin/env python3
# Libraries {{{
import os
import sys
import json
import getopt
import time
import psutil

# In-house modules
abspath = os.path.dirname(os.path.abspath(__file__))
sys.path.append(abspath + '/PythonUtilities/')
import AsterStudyUtilities as asus
# }}}

# Start time
startTime = time.time()
# Get the file with the parameters {{{
try:
    opts, args = getopt.getopt(sys.argv[1:], 'i:kcpa')
    if len(args)>0:
        raise getopt.GetoptError('Too many arguments!')
except getopt.GetoptError:
    print('Wrong call! The excecution can have the following arguments:')
    print('To indicate the input parameters: -i "fileName" (.json)')
    print('To indicate to keep the existing aster study folder: -k')
    print('To indicate to continue an aster study folder: -c')
    print('To indicate to run the post-processing: -p')
    print('To indicate to code_aster with the alternative root: -a')
    raise
# Initialisation
inParams = False
remove = True
continuee = False
onlyPost = False
defRoot = True
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
if not inParams:
    print("Wrong call! The input file with the parameters is not " +
            "present. Use -i 'fileName' (.json)")
    raise
# }}}

# Parameters {{{
# Read parameters {{{
with open(parFile, 'r') as inFile:
    params = json.load(inFile)
fileNamesParams   = params["fileNames"]
codeParams = params.get("code_aster", {})

asterFolder       = fileNamesParams["asterFolder"]
dummFile          = fileNamesParams["dummFile"]
geomFile = fileNamesParams["geomFile"]

memory_limit_default = psutil.virtual_memory().available*0.1/1000000.0
# }}}
# }}}

# Compute necessary paths {{{
baseDir = os.getcwd()
workDir = baseDir + asterFolder
if parFile[:5] == '/home':
    paramsPath = parFile
else:
    paramsPath = workDir + '/' + parFile
# }}}

# Set up study {{{
commFile          = dummFile[:-5] + '.comm'
dict_P = {"memory_limit" : memory_limit_default}
dict_P_ext = codeParams.get("dict_P", {})
dict_P.update(dict_P_ext)
dict_F = {"comm" : commFile}
dict_A = {"memjeveux" : 60}
study = asus.AsterStudy(dict_F = dict_F, dict_P = dict_P, dict_A = dict_A,
        defaultRoot = defRoot)
# Set work directory
study._workDir = workDir
study.CreateStudy(deletePrevious = remove, workDir = workDir)
# Change dummies
srcFolder = os.path.dirname(os.path.abspath(__file__))
dummies = [('#workDir', workDir),
        ('#parFile', paramsPath),
        ('#srcDir', srcFolder)]
asus.ReplaceDummyFile(srcFolder + '/codeasters' + dummFile,
        workDir + commFile, dummies)
# Copy meshes
for key in geomFile:
    meshPathIn = baseDir + geomFile[key]
    print(os.path.basename(meshPathIn))
    os.system('cp {} {}'.format(meshPathIn, workDir))
# Copy parameters
os.system('cp {} {}'.format(parFile, workDir))
# }}}

fail = study.RunStudy(outSalome = True)

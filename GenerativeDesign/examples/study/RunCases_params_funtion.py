#!/usr/bin/env python3
# Libraries {{{
import sys
import os
import json
import getopt
import multiprocessing
import time
import numpy as np
import shutil

# In-house modules
abspath = os.path.dirname(os.path.abspath(__file__))
sys.path.append(abspath + '/PythonUtilities/')
# }}}

# Get external arguments and read parameters {{{
try:
    opts, args = getopt.getopt(sys.argv[1:], 'i:')
    if len(args)>0:
        raise getopt.GetoptError('Too many arguments!')
except getopt.GetoptError:
    print('Wrong call! The excecution can have the following arguments:')
    print('To indicate the input parameters: -i "fileName" (.json)')
    raise
# Initialisation
parFile = None
for opt, arg in opts:
    # File name
    if opt in ['-i']:
        parFile = arg
if not parFile:
    print("Wrong call! The input file with the parameters is not " +
            "present. Use -i 'fileName' (.json)")
    raise
# }}}

# Parameters {{{
# Read parameters
with open(parFile, 'r') as inFile:
    runParams = json.load(inFile)
parFile     = runParams["parFile"]
studyFolder = runParams["studyFolder"]
makeParameterModule = runParams["makeParameterModule"]
numOfProcessors = runParams["numOfProcessors"]
runSimultations = runParams["runSimultations"]
try:
    alternativeRoot = runParams["alternativeRoot"]
except KeyError:
    alternativeRoot = False
RunAnalysis = runParams["RunAnalysis"]
U = runParams.get("unit0", 10)
continuePreviews = runParams.get("continuePreviews", False)
# }}}

# Functions {{{
# Function to run analyses.
# It can be used to run in parallel.
def RunCase(paramFileName):
    exe = 'python3 ' + abspath + '{} -k'.format(RunAnalysis)
    # Add options
    if alternativeRoot:
        exe += ' -a'
    if continuePreviews:
        exe += ' -c'
    # Add path
    exe += ' -i ' + paramFileName
    # Execute
    os.system(exe)
    return
# }}}

# Create params.json for each case {{{
masterDir = os.getcwd()
# Add function to create parameter files
exec("from {} import ParameterCases".format(makeParameterModule))
# Read base parameters
with open(masterDir + parFile, 'r') as File:
    params = json.load(File)
# Make list with the parameters of each case
paramsList = ParameterCases(params, U, masterDir + studyFolder)
# Write parameters
paramsPaths = []
for params_i in paramsList:
    # Get aster folder
    asterFolder_i = params_i["fileNames"]["asterFolder"]
    if not continuePreviews:
        if os.path.exists(asterFolder_i):
            shutil.rmtree(asterFolder_i)
        os.makedirs(asterFolder_i)
    path_i = os.path.join(asterFolder_i, "params.json")
    paramsPaths.append(path_i)
    with open(path_i, 'w') as fle:
        json.dump(params_i, fle, indent = 4)
# }}}

# Run analyses {{{
if runSimultations:
    paramsPaths.reverse()
    # Run in parallel
    with multiprocessing.Pool(processes = numOfProcessors) as pool:
        for case in pool.imap(RunCase, paramsPaths):
            pass
## }}}

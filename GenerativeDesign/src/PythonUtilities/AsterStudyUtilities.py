############################################################
# Module with utilities to create code aster studies and run
# them
############################################################
# Libraries {{{
import os
import sys
# }}}

# Class for Code_aster {{{
class AsterStudy:
    # Properties {{{
    @property
    def dict_F(self):
        return self._dict_F
    @property
    def dict_P(self):
        return self._dict_P
    @property
    def dict_A(self):
        return self._dict_A
    @property
    def workDir(self):
        return self._workDir
    @property
    def defaultRoot(self):
        return self._defaultRoot
    # }}}
    # Initialiser with default values {{{
    def __init__(self, dict_F, dict_P = {}, dict_A = {}, defaultRoot = True):
        # Find aster_root
        # Search for parallel version
        if os.path.isfile('/opt/aster146p/bin/as_run'):
            aster_root = '/opt/aster146p/bin/as_run'
        else:
            # Search for stable version
            if os.path.isfile('/opt/aster146/bin/as_run'):
                aster_root = '/opt/aster146/bin/as_run'
            else:
                if os.path.isfile('/opt/aster/bin/as_run'):
                    aster_root = '/opt/aster/bin/as_run'
                else:
                    raise SystemError("Code_aster version not found.")
        # Set version
        if ("mpi_nbcpu" in dict_F.keys()) or ("mpi_nbnoeud" in dict_F.keys()):
            version = "PAR14.6MUPT"
        else:
            version = "stable"
        # Set default settings {{{
        self._dict_P = {
                "actions" : "make_etude",
                "aster_root" : aster_root,
                "version" : version,
                "time_limit" : 720000
                }
        self._dict_A = {
                }
        self._defaultRoot = defaultRoot
        if defaultRoot == False:
            self._dict_P["aster_root"] = (os.getenv('HOME') + '/salome_meca' +
                                  '/appli_V2019.0.3_universal' +
                                  '/salome shell -- as_run')
        # }}}
        self._dict_F = dict_F
        # Update P and A dictionaries
        self._dict_P.update(dict_P)
        self._dict_A.update(dict_A)
        return
    # }}}
    # Write export file and directory {{{
    def CreateStudy(self, deletePrevious = True, workDir = None):
        # Test for existence of required names
        if (workDir == None):
            print('Warning: workDir has not been defined. '
                    + 'os.getcwd() will be used instead')
            workDir = os.getcwd()
        self._workDir = workDir
        comm = self.dict_F["comm"]
        # Create the main directory
        if deletePrevious:
            os.system('rm -r {}'.format(workDir))
        os.makedirs(self.workDir, exist_ok = True)
        # Create mess and export directories
        os.makedirs(workDir + '/MESS', exist_ok = True)
        os.makedirs(workDir + '/EXPORT', exist_ok = True)
        # Write export file
        with open(workDir + '/EXPORT/Study.export', 'w') as fle:
            # Write P parameters
            for key, value in self.dict_P.items():
                fle.write('P {} {}\n'.format(key, value))
            # Write A parameters
            for key, value in self.dict_A.items():
                fle.write('A {} {}\n'.format(key, value))
            # Write comm and mess
            fle.write('F comm {} D 1\n'.format(workDir + comm))
            fle.write('F mess {} R 6\n'.format(workDir + '/MESS/mess.mess'))
    # }}}
    # Run Study {{{
    # delete$$$s allows the deletion of previous results
    # outSalome to indicate if the command Salome needs to be used. If False,
    # the program assumes that the salome shell is open and if True,
    # it assumes that it is not open.
    def RunStudy(self, deleteRMEDs = False,
            deleteEPSs = False,
            deleteTABs = False,
            deleteGMSHs = False,
            deleteMSHs = False,
            outSalome = False):
        # Test for the existence of the export file
        fullExportPath = self.workDir + '/EXPORT/Study.export'
        isFile = os.path.isfile(fullExportPath)
        if isFile:
            # Delete previous *.rmed files if demanded
            if deleteRMEDs:
                os.system(f'rm {self.workDir}/*.rmed')
            if deleteEPSs:
                os.system(f'rm {self.workDir}/*.eps')
            if deleteTABs:
                os.system(f'rm {self.workDir}/*.tab')
            if deleteGMSHs:
                os.system(f'rm {self.workDir}/*.gmsh')
            if deleteMSHs:
                os.system(f'rm {self.workDir}/*.msh')
            # Run aster study
            if outSalome:
                aster_root = self.dict_P["aster_root"]
                suc = os.system(aster_root + ' ' + fullExportPath)
            else:
                suc = os.system('as_run ' + fullExportPath)
        else:
            sys.exit('Error: it is not possible to run the'
                    + ' study. The export file has not been'
                    + ' created yet.')
        return suc
    # }}}
# }}}

# Create file from dummy and add work directory {{{
def ReplaceDummyFile(dummFileName, outFileName, strs):
    dummFile = open(dummFileName, 'r')
    outFile  = open(outFileName, 'w')
    lines = dummFile.readlines()
    for line in lines:
        writeLine = True
        for dummStr, outStr in strs:
            if (dummStr in line):
                outFile.write(line.replace(dummStr, outStr))
                writeLine = False
        if writeLine:
            outFile.write(line)
    dummFile.close()
    outFile.close()
# }}}

# Run a study in a list of studies {{{
#  Useful to run in parallel
def RunStudyInList(study, deleteRMEDs = False):
    study.RunStudy(deleteRMEDs = deleteRMEDs)
# }}}

# Run salome at an specific port and then kill the port {{{
#  Default port: 2910
def RunSalomeAtPort(fileName, port = 2910):
    runSalome = f'runSalome.py --port {port} -t {fileName}'
    killSalome = f'killSalomeWithPort.py {port}'
    os.system(runSalome)
    os.system(killSalome)
# }}}

# Delete tmp {{{
def DeleteTmpFiles(exportFile, fileNames):
    # Get tmp directory {{{
    with open(exportFile) as fle:
        readlines = fle.readlines()
        # Search rep_trav
        for line in readlines:
            if "rep_trav" in line:
                rep_trav_line = line
                break
        # Get rep_trav input
        rep_trav_line = rep_trav_line.split()
        rep_trav_dir = rep_trav_line[-1]
    # }}}
    # Delete each file in fileNames {{{
    for filename in fileNames:
        filePath = os.path.join(rep_trav_dir, filename)
        try:
            os.system("rm " + filePath)
        except FileNotFoundError:
            pass
    # }}}
    return
# }}}

# Test {{{
if __name__ == '__main__':
    study = AsterStudy()
    study.workDir = os.getcwd() + '/TryModule'
    study.comm = '/study.comm'
    study.CreateStudy()
# }}}

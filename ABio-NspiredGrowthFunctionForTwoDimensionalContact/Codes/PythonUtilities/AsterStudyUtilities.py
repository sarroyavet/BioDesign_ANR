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
    # Initialiser with default values {{{
    def __init__(self):
        self.workDir        = None
        self.comm           = None
        self.aster_root     = (os.getenv('HOME') + '/salome_meca' +
                              '/appli_V2019.0.3_universal' +
                              '/salome')
        self.export         = '/EXPORT/Study.export'
        self.mess           = '/MESS/mess.mess'
        self.actions        = 'make_etude'
        self.version        = 'stable'
        self.consbtc        = 'oui'
        self.corefilesize   = 'unlimited'
        self.cpresok        = 'RESNOOK'
        self.debug          = 'nodebug'
        self.follow_output  = 'yes'
        self.facmtps        = 1
        self.lang           = 'en'
        self.mpi_nbcpu      = 1
        self.mpi_nbnoeud    = 1
        self.mode           = 'interactif'
        self.memjob         = 5120000
        self.memory_limit   = 20000.0
        self.time_limit     = 72000.0
        self.tpsjob         = 1201
        self.memjeveux      = 20000.0
        self.tpmax          = 72000.0
    # }}}
    # Write export file and directory {{{
    def CreateStudy(self, deletePrevious = True):
        workDir     = self.workDir
        comm        = self.comm
        aster_root  = self.aster_root
        export      = self.export
        mess        = self.mess
        # Test for existence of required names
        if (workDir == None):
            print('Warning: workDir has not been defined. '
                    + 'os.getcwd() will be used instead')
            workDir = os.getcwd()
            self.workDir = workDir
        if (comm == None):
            sys.exit('Error: the .comm file has not been ' +
                    'defined')
        # Create the main directory
        if deletePrevious:
            os.system(f'rm -r {self.workDir}')
        os.makedirs(self.workDir, exist_ok = True)
        # Create mess and export directories
        messDir = os.path.dirname(workDir + mess)
        exportDir = os.path.dirname(workDir + export)
        os.makedirs(exportDir, exist_ok = True)
        os.makedirs(messDir, exist_ok = True)
        # Write export file
        outFile = open(self.workDir + self.export, 'w')
        outFile.write('P actions %s\n' % (self.actions))
        outFile.write('P aster_root %s shell -- as_run\n' % (aster_root))
        outFile.write(('P consbtc %s\n' % (self.consbtc)))
        outFile.write(('P corefilesize %s\n' % (self.corefilesize)))
        outFile.write(('P cpresok %s\n' % (self.cpresok)))
        outFile.write(('P debug %s\n' % (self.debug)))
        outFile.write(('P facmtps %g\n' % (self.facmtps)))
        outFile.write(('P follow_output %s\n' % (self.follow_output)))
        outFile.write(('P lang %s\n' % (self.lang)))
        outFile.write(('P memjob %g\n' % (self.memjob)))
        outFile.write(('P memory_limit %g\n' % (self.memory_limit)))
        outFile.write(('P mode %s\n' % (self.mode)))
        outFile.write(('P mpi_nbcpu %g\n' % (self.mpi_nbcpu)))
        outFile.write(('P mpi_nbnoeud %g\n' % (self.mpi_nbnoeud)))
        outFile.write(('P time_limit %g\n' % (self.time_limit)))
        outFile.write(('P tpsjob %g\n' % (self.tpsjob)))
        outFile.write(('P version %s\n' % (self.version)))
        outFile.write(('A memjeveux %g\n' % (self.memjeveux)))
        outFile.write(('A tpmax %g\n' % (self.tpmax)))
        outFile.write(('F comm %s D  1\n'
            % (self.workDir + self.comm)))
        outFile.write(('F mess %s R  6\n'
            % (self.workDir + self.mess)))
        outFile.close()
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
        fullExportPath = self.workDir + self.export
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
                suc = os.system(os.getenv('HOME') + '/salome_meca/appli_V2019.0.3_universal/salome shell -- as_run ' + fullExportPath)
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

# Test {{{
if __name__ == '__main__':
    study = AsterStudy()
    study.workDir = os.getcwd() + '/TryModule'
    study.comm = '/study.comm'
    study.CreateStudy()
# }}}

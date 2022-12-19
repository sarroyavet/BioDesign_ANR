# Libraries {{{
import sys
import os
import json
import re
import copy
import multiprocessing
import time
import glob
import numpy as np
from scipy.signal import savgol_filter
from attrdict import AttrMap
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"],
    "font.size" : 9})
import matplotlib.animation as animation
from cycler import cycler

# In-house modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import ReadGmsh as rgmsh
from Geoanalysis import StraightLine as st
# }}}

# Classes {{{
# MorphoCases class {{{
class MorphoCases(object):

    # Properties {{{
    @property
    def maxVal(self):
        return self._maxVal
    @property
    def valRatio(self):
        return self._valRatio
    @property
    def numCases(self):
        return self._numCases
    @property
    def valName(self):
        return self._valName
    @property
    def cases(self):
        return self._cases
    # }}}

    # Initialiser {{{
    def __init__(self, maxVal, valRatio, numCases, valName,
            params, folderNames = [], paramsName = 'params.json',
            root = None, summ = False):
        self._maxVal = maxVal
        self._valRatio = valRatio
        self._numCases = numCases
        self._valName = valName
        # Create cases
        if summ:
            self._cases = [CaseWithParameters(maxVal + valRatio*k1,
                valName, copy.deepcopy(params)) for k1 in
                range(numCases)]
        else:
            self._cases = [CaseWithParameters(maxVal/(valRatio**k1),
                valName, copy.deepcopy(params)) for k1 in
                range(numCases)]
        # Add folders
        if (len(folderNames) == 0) and root:
            folNames = self.FolderNamesFromCases(root)
        elif (not len(folderNames) == 0) and not root:
            folNames = folderNames
        else:
            raise("Error: one and only one of the following should be defined: folderNames and root")
        for k1 in range(numCases):
            self.cases[k1].folderName = folNames[k1]
        return
    # }}}

    # Folder names from cases {{{
    def FolderNamesFromCases(self, root):
        if not(isinstance(root, str)):
            raise('Error: root must be a string.')
        # Create file names as root + cases[i]
        folderNames = [None]*self.numCases
        for k1 in range(self.numCases):
            folderNames[k1] = root + str(self.cases[k1].val)
        self.folderNames = folderNames
        return folderNames
    # }}}

    # Smooth data {{{
    @staticmethod
    def SmoothData(data, factor = 1, order = 5):
        """This function applies savgol_filter to a numpy array of data.
         considers each column of the array as a data source
        separately.  returns a numpy array of same shape with the
        smoothed data."""
        if not isinstance(data, np.ndarray):
            raise("Error: data must be a numpy array")
        rows, cols = data.shape
        smoothData = [None]*cols
        # Set filter window width
        length = int(rows/factor)
        length = length - 1 if length%2 == 0 else length
        # Smooth each column
        for col in range(cols):
            dataCol = data[:, col]
            smoothData[col] = savgol_filter(dataCol, length, order)
        smoothData = np.column_stack(tuple(smoothData))
        return smoothData
    # }}}

    # Post-process data and smooth data{{{
    def PostProcessDataAndSmoothedData(self, reportName, xlabels,
            ylabels, figExt = '.eps'):
        """This function executes the post-processing of a study of
        cases. Mainly, it provides graphs of contact pressure,
        maximum octahedral shear stress and minimum octahedral
        hydrostatic stress along with their smoothed version. In
        addition, it provides a graph with the standard deviation
        between the raw and smoothed data for each case."""
        if not isinstance(reportName, str):
            raise("Error: reportName must be a string.")
        if not reportName[0] == '/':
            reportName = '/' + reportName
        # Load and smooth data
        for k1 in range(self.numCases):
            case = self.cases[k1]
            data = case.add_dataInTable(reportName,
                    'rawReport')
            smData = MorphoCases.SmoothData(data)
            case.add_data(smData, 'smoothedReport')
        # Plot data
        repName, repExt = os.path.splitext(reportName)
        figPath = repName + figExt
        for k1 in range(self.numCases):
            case = self.cases[k1]
            case.MorphoPlotToFile_SmoothedData(figPath,
                    ('rawReport', 'smoothedReport'), xlabels, ylabels)
        return
    # }}}
# }}}

# Case class {{{
class Case(object):

    # Properties {{{
    @property
    def data(self):
        return self._data
    @property
    def folderName(self):
        if self._folderName:
            if self._folderName[:5] == '/home':
                return self._folderName
            else:
                return self._folderName[1:]
        else:
            return self._folderName
    @folderName.setter
    def folderName(self, folderName):
        if folderName == None or isinstance(folderName, str):
            self._folderName = folderName
        else:
            raise("Error: folderName can only be either None or a string")
        return
    # }}}

    # Initialiser {{{
    def __init__(self, folderName = None):
        self._data = {}
        self.folderName = folderName
        return
    # }}}

    # Data-related methods {{{
    def add_data(self, data, key):
        if isinstance(data, np.ndarray) and isinstance(key, str):
            self._data[key] = data
        else:
            raise("Error: data must be a numpy array and key a string.")
        return

    def delete_from_data(self, key, obj, axis):
        newData = np.delete(self.data[key], obj, axis)
        self.add_data(newData, key)
        return

    def add_dataInTable(self, fileName, key, numLinesOmit = 0,
            flipCols = [], flipColMax = []):
        if not isinstance(flipColMax, list):
            raise("Error: flipColMax must be a list.")
        if not isinstance(flipCols, list):
            raise("Error: flipCols must be a list.")
        if not isinstance(fileName, str):
            raise("Error: fileName must be a string.")
        if not fileName[0] == '/':
            fileName = '/' + fileName
        if self.folderName:
            path = self.folderName + fileName
        else:
            print("Warning: no folder name has been indicated.")
            path = fileName[1:]
        # Read file
        data = []
        with open(path, 'r') as File:
            for _ in range(numLinesOmit):
                File.readline()
            reading = True
            while reading:
                line = File.readline()
                if not line:
                    reading = False
                else:
                    # Convert non-commented line in a list of float
                    if not line[0] == '#':
                        line = re.split(',|\s', line)
                        line = list(filter(None, line))
                        line = [eval(val) for val in line]
                        if not len(line) == 0:
                            data.append(line)
        try:
            data = np.array(data)
        except:
            raise('Error: not possible to convert data list into a numpy array. Check that: all rows have the same amount of columns and that the non numeric rows have a # at the beginning')
        # Flip columns {{{
        if not len(flipCols) == 0:
            if data.ndim == 1:
                numCols = data.size
            else:
                _, numCols = data.shape
            cols = list(range(numCols))
            dataCols = [None]*numCols
            for col in cols:
                if col in flipCols:
                    dataCols[col] = np.flipud(data[:, col])
                else:
                    dataCols[col] = data[:, col]
            data = np.column_stack(tuple(dataCols))
        # }}}
        # Flip column maximum {{{
        if not len(flipColMax) == 0:
            if data.ndim == 1:
                numCols = data.size
            else:
                _, numCols = data.shape
            cols = list(range(numCols))
            dataCols = [None]*numCols
            for col in cols:
                if col in flipColMax:
                    maxCol = data[:, col].max()
                    dataCols[col] = maxCol - data[:, col]
                    #dataCols[col] = -data[:, col]
                else:
                    dataCols[col] = data[:, col]
            data = np.column_stack(tuple(dataCols))
        # }}}
        self.add_data(data, key)
        return data
    # }}}

    # Post-process data block{{{
    def PostProcessData(self, key, figName, xlabel, ylabel,
            legends = None, lineStyles = None, figExt = '.eps',
            figsize = [6, 4.5], figConf = lambda ax: None):
        if not isinstance(figExt, str):
            raise("Error: figExt must be a string.")
        if not isinstance(figName, str):
            raise("Error: figName must be a string.")
        if not isinstance(xlabel, str):
            raise("Error: xlabel must be a string.")
        if not isinstance(ylabel, str):
            raise("Error: ylabel must be a string.")
        if not isinstance(key, str):
            raise("Error: key must be a string.")
        data = self.data[key]
        rows, cols = data.shape
        if cols == 1:
            raise("Error: the numpy array in data[key] must have at least 2 columns.")
        if lineStyles:
            if not len(lineStyles) == cols - 1:
                raise("Error: lineStyles must be a list or tuple whose length must be the number of columns of data minus one.")
        else:
            lineStyles = ['-']*(cols - 1)
        if legends:
            if not len(legends) == cols - 1:
                raise("Error: legends must be a list or tuple whose length must be the number of columns of data minus one.")
        x = data[:, 0]
        data = data[:, 1:]
        fig = plt.figure(figsize = figsize)
        ax = plt.subplot(111)
        for k1 in range(cols - 1):
            ax.plot(x, data[:, k1], lineStyles[k1])
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if legends:
            ax.legend(legends)
        figConf(ax)
        if not figName[0] == '/':
            figName = '/' + figName
        fname, fext = os.path.splitext(figName)
        fig.tight_layout()
        if self.folderName:
            fig.savefig(self.folderName + fname + figExt)
        else:
            print("Warning: no folder name has been indicated.")
            fig.savefig(fileName[1:])
        plt.close()
        return
    # }}}

    # Post-process data block with fig configuration {{{
    def PostProcessDataFigConf(self, key, figName, figExt = '.eps',
            figsize = [6, 4.5], figConf = None, markers = None):
        data = self.data[key]
        rows, cols = data.shape
        # Test input {{{
        if not isinstance(figExt, str):
            raise("Error: figExt must be a string.")
        if not isinstance(figName, str):
            raise("Error: figName must be a string.")
        if not isinstance(key, str):
            raise("Error: key must be a string.")
        if not isinstance(figsize, list):
            raise("Error: figsize must be a list.")
        if not len(figsize) == 2:
            raise("Error: the length of figsize must be 2.")
        if cols == 1:
            raise("Error: the numpy array in data[key] must have at least 2 columns.")
        if markers == None:
            markers = [","]*(cols - 1)
        else:
            if not isinstance(markers, list):
                raise("Error: markers must be a list.")
            if not len(markers) == cols - 1:
                raise("Error: the length of markers must be equal to the number of columns of data minus one (the x-axes). If you want to let a curve without markers use: ','.")
        # }}}
        # Plot data {{{
        x = data[:, 0]
        data = data[:, 1:]
        fig = plt.figure(figsize = figsize)
        ax = plt.subplot(111)
        for k1 in range(cols - 1):
            ax.plot(x, data[:, k1], marker = markers[k1])
        if not figConf == None:
            figConf(ax)
        # }}}
        # Save file {{{
        if not figName[0] == '/':
            figName = '/' + figName
        fname, fext = os.path.splitext(figName)
        if self.folderName:
            fig.savefig(self.folderName + fname + figExt)
        else:
            print("Warning: no folder name has been indicated.")
            fig.savefig(fileName[1:])
        plt.close()
        # }}}
        return
    # }}}

    # PostProcessing2Ys {{{
    def PostProcess2Ys(self, key, figName, figExt = '.eps',
            figsize = [6, 4.5], figConf1 = None, figConf2 = None,
            markers = None, colors = None, x = 0, y1 = [1], y2 = [2]):
        data = self.data[key]
        rows, cols = data.shape
        # Test input {{{
        if not isinstance(figExt, str):
            raise("Error: figExt must be a string.")
        if not isinstance(figName, str):
            raise("Error: figName must be a string.")
        if not isinstance(key, str):
            raise("Error: key must be a string.")
        if not isinstance(figsize, list):
            raise("Error: figsize must be a list.")
        if not len(figsize) == 2:
            raise("Error: the length of figsize must be 2.")
        if cols == 1:
            raise("Error: the numpy array in data[key] must have at least 2 columns.")
        if markers == None:
            markers = [","]*(cols - 1)
        else:
            if not isinstance(markers, list):
                raise("Error: markers must be a list.")
            if not len(markers) >= cols - 1:
                raise("Error: the length of markers must be equal to the number of columns of data minus one (the x-axes). If you want to let a curve without markers use: ','.")
        if colors == None:
            colors = ["tab:blue",
                      "tab:orange",
                      "tab:green",
                      "tab:red",
                      "tab:purple",
                      "tab:brown",
                      "tab:pink",
                      "tab:gray",
                      "tab:olive",
                      "tab:cyan"]
        else:
            if not isinstance(colors, list):
                raise("Error: colors must be a list.")
            if not len(colors) >= cols - 1:
                raise("Error: the length of markers must be equal to the number of columns of data minus one (the x-axes). If you want to let a curve without markers use: ','.")
        # }}}
        # Plot data {{{
        fig, ax1 = plt.subplots(figsize = figsize)
        ax2 = ax1.twinx()
        k1 = 0
        for yy in y1:
            ax1.plot(data[:, x], data[:, yy], marker = markers[k1], color = colors[k1])
            k1 += 1
        for yy in y2:
            ax2.plot(data[:, x], data[:, yy], marker = markers[k1], color = colors[k1])
            k1 += 1
        if not figConf1 == None:
            figConf1(ax1)
        if not figConf2 == None:
            figConf2(ax2)
        # }}}
        # Save file {{{
        if not figName[0] == '/':
            figName = '/' + figName
        fname, fext = os.path.splitext(figName)
        if self.folderName:
            fig.savefig(self.folderName + fname + figExt)
        else:
            print("Warning: no folder name has been indicated.")
            fig.savefig(fileName[1:])
        plt.close()
        # }}}
    # }}}

    # Morpho data and smoothed data plot to file {{{
    def MorphoPlotToFile_SmoothedData(self, figName, keys, xlabels,
            ylabels, figFormat = '.eps'):
        if not isinstance(figName, str):
            raise('Error: figName must be a string')
        fileName, fileExt = os.path.splitext(figName)
        fileName += figFormat
        if not fileName[0] == '/':
            fileName = '/' + fileName
        x      = self.data[keys[0]][:, 0]
        data   = self.data[keys[0]][:, 1:]
        smData = self.data[keys[1]][:, 1:]
        rows, cols = data.shape
        fig = plt.figure(figsize = [7, 9.5])
        axes = [None]*cols
        if not len(xlabels) == cols:
            raise("Error: there are not enough xlables. If you do not want a label use None within the list of labels")
        if not len(ylabels) == cols:
            raise("Error: there are not enough ylables. If you do not want a label use None within the list of labels")
        for k3 in range(cols):
            axes[k3] = plt.subplot(cols, 1, k3 + 1)
            axes[k3].plot(x, data[:, k3], '.')
            axes[k3].plot(x, smData[:, k3])
            if xlabels[k3]:
                axes[k3].set_xlabel(xlabels[k3])
            if ylabels[k3]:
                axes[k3].set_ylabel(ylabels[k3])
        if self.folderName:
            fig.savefig(self.folderName + fileName)
        else:
            print("Warning: no folder name has been indicated.")
            fig.savefig(fileName[1:])
        plt.close()
        return
    # }}}

    # Animate data {{{
    def AnimateData(self, figName, keys = None, savekwargs = None,
            figConf = lambda ax: None, figExt = '.mp4',
            shareAbscissa = True, numLines = 1, markers = None,
            lineStyles = '-'):
        """Animation of 2D arrays (x, y)"""
        # Test inputs {{{
        if not isinstance(figExt, str):
            raise("Error: figExt must be a string.")
        if not isinstance(figName, str):
            raise("Error: figName must be a string.")
        if keys == None:
            keys = self.data.keys()
        else:
            if not isinstance(keys, list):
                raise("Error: keys must be a list or None (in which case all keys from self.data are used).")
        if not callable(figConf):
            raise("Error: figConf must be a callable function.")
        if not (isinstance(savekwargs, dict) or savekwargs == None):
            raise("Error: savekwargs must be either a dictionary or None.")
        if markers == None:
            markers = [None]*numLines
        if not (isinstance(markers, list)):
            raise("Error: markers must be either None or a list.")
        if not(len(markers) == numLines):
            raise("Error: the length of markers must be equal to numLines.")
        if lineStyles == '-':
            lineStyles = ['-']*numLines
        if not (isinstance(lineStyles, list)):
            raise("Error: lineStyles must be either '-' or a list.")
        if not(len(lineStyles) == numLines):
            raise("Error: the length of lineStyles must be equal to numLines.")
        # }}}
        # Initialisation
        fig, ax = plt.subplots()
        lines = [None]*numLines
        for k1 in range(numLines):
            lines[k1], = plt.plot([], [], marker = markers[k1],
                    linestyle = lineStyles[k1])
        # Setup of the initial configuration
        def init():
            figConf(ax)
            return lines
        # Add the x-y data
        def animate(key):
            data = self.data[key]
            rows, cols = data.shape
            if shareAbscissa:
                numYs = cols - 1
                if not numYs == numLines:
                    raise("Error: the data cannot generate the indicated number of lines")
                x = data[:, 0]
                for k1 in range(numYs):
                    y = data[:, k1 + 1]
                    lines[k1].set_data(x, y)
            else:
                if not cols%2 == 0:
                    raise("Error: if shareAbscissa is false, the number of columns of each data block must be even.")
                numYs = int(cols/2)
                if not numYs == numLines:
                    raise("Error: the data cannot generate the indicated number of lines")
                for k1 in range(numYs):
                    x = data[:, k1*2]
                    y = data[:, k1*2 + 1]
                    lines[k1].set_data(x, y)
            return lines
        # Create animation
        ani = animation.FuncAnimation(fig, animate, frames = keys,
                blit = True, init_func = init)
        # Save animation
        if not figName[0] == '/':
            figName = '/' + figName
        fname, fext = os.path.splitext(figName)
        if self.folderName:
            figFullName = self.folderName + fname + figExt
        else:
            print("Warning: no folder name has been indicated.")
            figFullName = fileName[1:]
        ani.save(figFullName)
        plt.close()
        return
    # }}}

    # Normalised data {{{
    def NormaliseData(self, key, cols = []):
        """Normalise the data by column. Key indicates the data set to
        be normalised and cols indicates the columns of the data block
        to be normalised. If cols == [] then all the columns will be
        normalised."""
        # Tests input {{{
        if not isinstance(key, str):
            raise("Error: key must be a string.")
        if not key in self.data.keys():
            raise("Error: there is no data for the indicated key: " + key)
        data = self.data[key]
        if data.ndim == 1:
            numDatas = 1
        else:
            numDatas = data.shape[1]
        allowedDatas = list(range(numDatas))
        if not isinstance(cols, list):
            raise("Error: cols must be a list.")
        if len(cols) == 0:
            cols = allowedDatas
        else:
            for col in range(1, numDatas):
                allowedDatas.append(-col)
            for col in cols:
                if not col in allowedDatas:
                    raise("Error: at least one of the elements of cols refers to a column out of bounds.")

        # }}}
        # Normalisation of each column {{{
        normCols = []
        for col in range(numDatas):
            if col in cols:
                normCols.append(Case.NormaliseVector(data[:, col]))
            else:
                normCols.append(data[:, col])
        normCols = np.column_stack(tuple(normCols))
        # }}}
        self.add_data(normCols, key)
        return
    @staticmethod
    def NormaliseVector(vec):
        minVec = vec.min()
        maxVec = vec.max()
        dif = maxVec - minVec
        m = 1.0/dif
        b = -minVec/dif
        size = vec.shape[0]
        normVec = [None]*size
        for k1 in range(size):
            normVec[k1] = m*vec[k1] + b
        normVec = np.array(normVec)
        return normVec
    # }}}
# }}}

# Class case with parameters {{{
class CaseWithParameters(Case):
    # Properties {{{
    @property
    def val(self):
        return self._val
    @property
    def valName(self):
        return self._valName
    @property
    def params(self):
        return self._params
    @params.setter
    def params(self, params):
        # Check params is a dictionary
        if not(isinstance(params, dict)):
            raise('Error: params must be a dictionary')
        parameters = AttrMap(params.copy())
        exec("parameters.%s = %g" %(self.valName, self.val))
        self._params = dict(parameters)
        return
    @property
    def folderName(self):
        return super().folderName
    @folderName.setter
    def folderName(self, folderName):
        super(CaseWithParameters,
                self.__class__).folderName.fset(self, folderName)
        if hasattr(self, '_params'):
            self._params["fileNames"]["asterFolder"] = folderName
        return
    # }}}

    # Initialiser {{{
    def __init__(self, val, valName, params, folderName = None):
        super().__init__(folderName)
        self._val = val
        self._valName = valName
        self.params = params
        self.folderName = folderName
        return
    # }}}

    # Write parameters {{{
    def WriteParameters(self, fileName, ensureNewFolder = False,
            unit = 10):
        if not isinstance(fileName, str):
            raise('Error: fileName must be a string')
        if not fileName[0] == '/':
            fileName = '/' + fileName
        if ensureNewFolder:
            os.system('rm -r ' + self.folderName)
        if self.folderName:
            os.makedirs(self.folderName, exist_ok = True)
            path = self.folderName + fileName
        else:
            print("Warning: no folder name has been indicated.")
            path = fileName[1:]
        # Modify unit
        self.params["c_a_unit"] = unit
        with open(path, 'w') as File:
            json.dump(self.params, File, indent = 4)
        return path
    # }}}
# }}}

# MorDesResults {{{
class MorDesResults(object):
    # Properties {{{
    @property
    def asterFolder(self):
        return self._asterFolder
    @property
    def genResults(self):
        return self._genResults
    @property
    def pressure(self):
        return self._pressure
    @property
    def contContours(self):
        return self._contContours
    # }}}

    # Initialiser {{{
    def __init__(self, asterFolder,
            bestPressureName = '/bestPressure.tab',
            meshSizeName = '/meshSize.tab', reportName = '/report.tab',
            pressureName = '/RESULTS/pressure.tab',
            contourSubDir = '/RECONSTRUCTION', masterName = None,
            slaveName  = None, contAid = 4):
        # Add aster folder
        self._asterFolder = asterFolder
        # Add general results {{{
        self._genResults = Case(asterFolder)
        # Add one-file results
        self.genResults.add_dataInTable(bestPressureName, "bestPressure",
                numLinesOmit = 7, flipColMax = [0])
        self.genResults.add_dataInTable(meshSizeName, "meshSize")
        self.genResults.add_dataInTable(reportName, "report")
        # }}}
        # Add pressures {{{
        # Add pressure at each time
        self._pressure = Case(asterFolder)
        presName, presExt = os.path.splitext(pressureName)
        lenExt  = len(presExt)
        if asterFolder[:6] == '/home/':
            presNames = glob.glob(asterFolder + presName + '*')
            lenName = len(asterFolder + presName)
        else:
            presNames = glob.glob('.' + asterFolder + presName + '*')
            lenName = len('.' + asterFolder + presName)
        sortedNames = []
        for name in presNames:
            num = name[lenName: -lenExt]
            sortedNames.append((eval(num), num))
        sortedNames.sort()
        for name in sortedNames:
            fileName = presName + name[1] + presExt
            self.pressure.add_dataInTable(fileName, name[1],
                    numLinesOmit = 7, flipColMax = [0])
        # }}}
        # Add contact contours {{{
        contArea = self.genResults.data["report"][-1, contAid]
        #  Master {{{
        if masterName == None:
            master_nds = None
        else:
            resu = rgmsh.GmshOutput(asterFolder + contourSubDir + masterName)
            nds, _ = resu.views[0].GetMeshFromLines('scalar', 0)
            nds = np.array([list(nod) for nod in nds])
            arcLine = 0.0
            dxm = -nds[0, 0]
            dym = -nds[0, 1]
            master_nds = [[nds[0, 0] + dxm, nds[0, 1] + dym]]
            pnode = nds[0]
            for node in nds[1:]:
                if arcLine > contArea*1.25:
                    break
                dx = node[0] - pnode[0]
                dy = node[1] - pnode[1]
                ds = np.sqrt(dx**2.0 + dy**2.0)
                arcLine += ds
                master_nds.append([node[0] + dxm, node[1] + dym])
                pnode = node
            master_nds = np.array(master_nds)
        # }}}
        #  Slave {{{
        if slaveName == None:
            slave_nds = None
        else:
            resu = rgmsh.GmshOutput(asterFolder + contourSubDir + slaveName)
            nds, _ = resu.views[0].GetMeshFromLines('scalar', 0)
            nds = np.array([list(nod) for nod in nds])
            arcLine = 0.0
            dxs = -nds[0, 0]
            dys = -nds[0, 1] + 0.
            slave_nds = [[nds[0, 0] + dxs, nds[0, 1] + dys]]
            pnode = nds[0]
            for node in nds[1:]:
                if arcLine > contArea*1.25:
                    break
                dx = node[0] - pnode[0]
                dy = node[1] - pnode[1]
                ds = np.sqrt(dx**2.0 + dy**2.0)
                arcLine += ds
                slave_nds.append([node[0] + dxs, node[1] + dys])
                pnode = node
            slave_nds = np.array(slave_nds)
        # }}}
        #
        self._contContours = {"master" : master_nds,
                "slave" : slave_nds,
                "contArea" : contArea}
        # }}}
        return
    # }}}

    # Add relative contour {{{
    def add_relativeContour(self):
        master = self.contContours["master"]
        slave = self.contContours["slave"]
        # Test master and slave arrays
        if not (isinstance(master, np.ndarray) and
                isinstance(slave, np.ndarray)):
            raise("Error: add_relativeContour requires both master and slave be defined.")
        # Initialisation
        lenx, _ = master.shape
        x = np.zeros(lenx)
        y = np.zeros(lenx)
        # Define curved x as the arc length of the master curve {{{
        for k1 in range(lenx - 1):
            node1 = master[k1]
            node2 = master[k1 + 1]
            dx = node2[0] - node1[0]
            dy = node2[1] - node1[1]
            ds = np.sqrt(dx**2.0 + dy**2.0)
            x[k1 + 1] = x[k1] + ds
        # }}}
        # Define y as perpendicular straight line to x {{{
        #  Normal lines {{{
        normal_lines = [None]*lenx
        normal = np.zeros((lenx, 2))
        for k1 in range(lenx):
            if k1 == 0:
                line = st(p1 = master[0], p2 = master[1])
                normal_lines[k1] = st(m = -1.0/line.m, p1 = master[k1])
            elif k1 == lenx - 1:
                line = st(p1 = master[-2], p2 = master[-1])
                normal_lines[k1] = st(m = -1.0/line.m, p1 = master[k1])
            else:
                p1 = master[k1 - 1]
                p2 = master[k1]
                p3 = master[k1 + 1]
                xx = [p1[0], p2[0], p3[0]]
                yy = [p1[1], p2[1], p3[1]]
                mm, _ = np.polyfit(xx, yy, 1)
                normal_lines[k1] = st(m = -1.0/mm, p1 = master[k1])
        #  }}}
        #  Distance from curved x to the slave surface {{{
        for k1 in range(lenx):
            # Line perpendicular to the point
            line1 = normal_lines[k1]
            # Closest two points of the slave curve to line1
            p1, _, ID = line1.ClosestPointFromSetOfPoints(slave)
            slave2 = np.delete(slave, ID, 0)
            p2, _, _ = line1.ClosestPointFromSetOfPoints(slave2)
            # line2 is a line that passes through p1  and p2
            line2 = st(p1 = p1, p2 = p2)
            # Point that corresponds to the image of the master point
            pc = st.Intersection(line1, line2)
            node = master[k1]
            dx = pc[0] - node[0]
            dy = pc[1] - node[1]
            y[k1] = np.sqrt(dx**2.0 + dy**2.0)
        #  }}}
        relative_nds = np.zeros((lenx, 2))
        for k1 in range(lenx):
            relative_nds[k1, 0] = x[k1]
            relative_nds[k1, 1] = y[k1]
        self._contContours["relative"] = np.array(relative_nds)
        # }}}
        return
    # }}}

    # Compute curvature {{{
    def Compute_Curvature(self):
        data = self.contContours
        contArea = data["contArea"]
        nds = data["relative"]
        x = nds[:, 0]
        y = nds[:, 1]
        x_t = np.gradient(x, edge_order = 2)
        y_t = np.gradient(y, edge_order = 2)
        xx_t = np.gradient(x_t, edge_order = 2)
        yy_t = np.gradient(y_t, edge_order = 2)
        curv = np.abs(xx_t*y_t - x_t*yy_t)/(x_t*x_t + y_t*y_t)**1.5
        self._contContours["curvature"] = curv
        return
    # }}}

    # Plot best pressure {{{
    def PlotBestPressure(self):
        contArea = self.genResults.data["report"][:, 4]
        def figConf(ax):
            ax.set_xlabel(r"$l\mathrm{[mm]}$")
            ax.set_ylabel(r"$p\mathrm{[GPa]}$")
            ax.set_xlim(0.0, contArea.max()*1.1)
            return ax
        self.genResults.PostProcessDataFigConf("bestPressure", "bestPressure",
                figExt = ".pdf", figConf = figConf)
        return
    # }}}

    # Plot mesh size {{{
    def PlotMeshSize(self):
        def figConf1(ax):
            ax.set_xlabel("Iteration")
            ax.legend([r"$\min(h_{\Gamma^c})$"], loc = "upper left")
            return ax
        def figConf2(ax):
            ax.legend([r"CPU-time"], loc = "upper right")
            return ax
        self.genResults.PostProcess2Ys("meshSize", "meshSize",
                figExt = '.pdf', y1 = [1], y2 = [4], figConf1 = figConf1,
                figConf2 = figConf2)
        return
    # }}}

    # Plot report {{{
    def PlotReport(self):
        def figConf1(ax):
            ax.set_xlabel("$\hat{t}$")
            ax.legend([r"$p_\mathrm{max}[\mathrm{GPa}]$",
                       r"${\sigma_v}^\mathrm{max}[\mathrm{GPa}]$",
                       r"$\sigma_\mathrm{hyc}^\mathrm{max}[\mathrm{GPa}]$"], loc = "upper left")
            return ax
        def figConf2(ax):
            ax.legend([r"$a_c[\mathrm{mm}]$"], loc = "upper right")
            return ax
        self.genResults.PostProcess2Ys("report", "report-1",
                figExt = '.pdf', y1 = [1, 2, 3], y2 = [4], figConf1 = figConf1,
                figConf2 = figConf2)
        dummCase = copy.deepcopy(self.genResults)
        dummCase.delete_from_data("report", 1, 1)
        dummCase.delete_from_data("report", 1, 1)
        dummCase.delete_from_data("report", 1, 1)
        dummCase.delete_from_data("report", 1, 1)
        def figConf3(ax):
            def forward(x):
                return x**2.0
            def inverse(x):
                return x**0.5
            ax.set_ylim(bottom = 0.0)
            ax.set_yscale('function', functions = (forward, inverse))
            ax.grid()
            return ax
        dummCase.PostProcessData("report", "report-2", r"$\hat{t}$",
                "", legends = [r"$Q_v$", r"$Q_p$"], figExt = ".pdf",
                figConf = figConf3)
        return
    # }}}

    # Plot pressure {{{
    def PlotPressure(self, numLines = 0, figConf = lambda x: x):
        data = self.pressure.data
        if numLines == 0:
            numLines = len(data)
        lkeys = list(data.keys())
        delta = int(len(lkeys)/numLines)
        if delta == 0:
            numLines = len(data)
            lkeys = list(data.keys())
            delta = int(len(lkeys)/numLines)
        graphs = [None]*numLines
        time = [None]*numLines
        for k1 in range(numLines - 1):
            key = lkeys[k1*delta]
            graphs[k1] = data[key]
            time[k1] = eval(key)
        key = lkeys[-1]
        graphs[-1] = data[key]
        time[-1] = eval(key)
        cm = 1.0/2.54
        fig, ax = plt.subplots(1, 2, figsize = [8*cm, 6*cm],
                gridspec_kw={'width_ratios': [10, 1]})
        colors = plt.cm.coolwarm(np.linspace(0, 1, numLines))
        for k1 in range(numLines):
            ax[0].plot(graphs[k1][:, 0], graphs[k1][:, 1], '.-', color = colors[k1])
        cmap = mpl.cm.coolwarm
        norm = mpl.colors.Normalize(vmin = 0.0, vmax = time[-1])
        cb   = mpl.colorbar.ColorbarBase(ax[1], cmap = cmap, norm = norm,
                boundaries = time)
        contArea = self.genResults.data["report"][:, 4]
        ax[0].set_xlim(0.0, contArea.max()*1.1)
        ax[0].set_xlabel(r'$l[\mathrm{mm}]$')
        ax[0].set_ylabel(r'$p[\mathrm{GPa}]$')
        ax[1].set_title(r'$\hat{t}_m$')
        figConf(ax[0])
        fig.tight_layout()
        plt.savefig(self.pressure.folderName + '/pressure.pdf')
        return
    # }}}

    # Plot contact contours {{{
    def PlotContContours(self, figsize = [4, 3], figConfms = lambda x: x,
            figConfr = lambda x: x, figExt = ".pdf"):
        data = self.contContours
        fig, ax = plt.subplots(figsize = figsize)
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(r"$y$")
        if data["master"] is None:
            if data["slave"] is None:
                print("Warning: no contact contours given.")
            else:
                ax.plot(data["slave"][:, 0], data["slave"][:, 1], 'k')
                ax.legend(["Slave"])
        else:
#            ax.plot(data["master"][:, 0], data["master"][:, 1],
#                    linewidth = 0.5, color = 'k')
            ax.fill_between(x = data["master"][:, 0],
                    y1 = data["master"][:, 1],
                    y2 = -100.0,
                    color = '#e3e3cdff',
                    alpha = 1.0)
            if data["slave"] is None:
                ax.legend(["Master"])
            else:
#                ax.plot(data["slave"][:, 0], data["slave"][:, 1],
#                    linewidth = 0.5, color = 'k')
                ax.fill_between(x = data["slave"][:, 0],
                        y1 = data["slave"][:, 1],
                        y2 = 1.0,
                        color = '#9ebeffff',
                        alpha = 1.0)
                ax.legend(["Master", "Slave"])
        ax.set_xlim(0.0, 1.05*data["contArea"])
        ax.set_ylim(-1.05*data["contArea"], 0.0)
        figConfms(ax)
        fig.tight_layout()
        fig.savefig(self.asterFolder + "/ContactGeometry" + figExt)
        if "relative" in data:
            contArea = data["contArea"]
            xx = data["relative"][:, 0]
            yy = data["relative"][:, 1]
            diff_x_cont = np.asarray(xx)
            idx = (np.abs(xx - contArea)).argmin()
            maxY = yy[idx]*2.0
            fig, ax = plt.subplots(figsize = figsize)
            ax.set_xlabel(r"$l[\mathrm{mm}]$")
            ax.set_ylabel(r"$n[\mathrm{mm}]$")
            ax.plot(data["relative"][:, 0], data["relative"][:, 1],
                    color = '#b40426ff')
            ax.plot([data["contArea"], data["contArea"]], [0.0, 1.0],
                    'k:')
            ax.set_xlim(0.0, 1.1*data["contArea"])
            xticks = ax.get_xticks()
            xticklabels = ax.get_xticklabels()
            newxticks = []
            newxticklabels = []
            for k1 in range(len(xticks) - 1):
                if xticks[k1] < contArea:
                    if xticks[k1 + 1] < contArea:
                        newxticks.append(xticks[k1])
                        newxticklabels.append(r'$%g$' % (xticks[k1]))
            newxticks.append(contArea)
            newxticklabels.append(r'$a_c$')
            ax.set_xticks(newxticks)
            ax.set_xticklabels(newxticklabels)
            ax.set_ylim(0.0, maxY)
            #ax.set_yscale('log')
            figConfr(ax)
            fig.tight_layout()
            fig.savefig(self.asterFolder + "/RelativeContactGeometry" + figExt)
        return
    # }}}

# }}}

# Functions
# Lp norm of a 1D function {{{
def LpNorm1D(x, y, p = 2.0):
    # Test input {{{
    if not isinstance(x, np.ndarray):
        raise("Error: x must be a numpy array.")
    if not isinstance(y, np.ndarray):
        raise("Error: y must be a numpy array.")
    if not x.ndim == 1:
        raise("Error: x must be a 1D array.")
    if not y.ndim == 1:
        raise("Error: y must be a 1D array.")
    if not x.size == y.size:
        raise("Error: x and y must have the same size.")
    if not (isinstance(p, float) or isinstance(p, int)):
        raise("Error: p must be a float or an integer.")
    if p <= 0:
        raise("Error: p must be greater than zero. But if you need L_0 or L_inf you can add it to this code.")
    # }}}
    norm = 0.0
    for k1 in range(x.size - 1):
        xa = x[k1]
        xb = x[k1+ 1]
        ya = y[k1]
        yb = y[k1+ 1]
        dx = abs(xb - xa)
        ym = 0.5*(ya + yb)
        norm = ym**p*dx
    norm = norm**(1.0/p)
    return norm
# }}}
# }}}

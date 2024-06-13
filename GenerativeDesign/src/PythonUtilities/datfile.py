# Libraries {{{
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"],
    "font.size" : 9})
# }}}

# Class datfile {{{
class datfile(object):
    # Properties {{{
    @property
    def filename(self):
        return self._filename
    @property
    def title(self):
        return self._title
    @property
    def datablocks(self):
        return self._datablocks
    # }}}

    # Initialiser {{{
    def __init__(self, filename, fmt = '5.6f'):
        # Read lines
        fle = open(filename, 'r')
        # Add title
        titleline = fle.readline()
        try:
            exec(titleline.replace(" ",""), globals(), globals())
            self._title = TITLE
        except:
            print("Given line: " + titleline)
            raise("Error: the first line must contain TITLE = '...'.")
        # Add blocks
        keepreading = True
        self._datablocks = {}
        timeline = fle.readline().split(',')
        while keepreading:
            # Get time
            try:
                exec(timeline[0].replace(" ",""), globals(), globals())
                time = TIME
            except:
                print("Given line: " + timeline[0])
                raise("Error: TIME was not defined.")
            # Get data and variable name lines
            if time == None:
                exec(fle.readline().replace(" ",""), globals(), globals())
                datalines = fle.readlines()
                time = 0
                fmt = '1d'
            else:
                try:
                    exec(timeline[1].replace(" ",""), globals(), globals())
                    numPoints = N
                except:
                    print("Given line: " + timeline[1])
                    raise("Error: TIME is different to None and N is not defined.")
                exec(fle.readline().replace(" ",""), globals(), globals())
                datalines = []
                for k1 in range(numPoints):
                    datalines.append(fle.readline())
            # Add datablock
            timestring = ('{:' + fmt + '}').format(time)
            self._datablocks[timestring] = datfile.datablock(VARIABLES, datalines)
            # Read next time line
            timeline = fle.readline().split(',')
            if timeline == ['']:
                keepreading = False
    # }}}

    # get_data_at_time {{{
    def get_data_at_time(self, time, tol = 0.000001):
        for datakey in self.datablocks:
            ti = eval(datakey)
            if abs(ti - time) < tol:
                return self.datablocks[datakey]
        print("Error: there is no data for the indicated time with tol = " + str(tol))
        raise
        return
    # }}}

    # Data block class {{{
    class datablock(object):
        # Properties {{{
        @property
        def numVars(self):
            return self._numVars
        @property
        def variables(self):
            return self._variables
        # }}}

        # add_data {{{
        def add_data(self, data, key):
            if isinstance(data, np.ndarray) and isinstance(key, str):
                self._variables[key] = data
            else:
                raise("Error: data must be a numpy array and key a string.")
            return
        # }}}

        # Initialiser {{{
        def __init__(self, variables, rawLines):
            # Get the number of variables in variable list/tuple
            numVars = len(variables)
            self._numVars = numVars
            # Get the data in rawLines as a numpy array
            data = []
            for k1 in range(numVars):
                data.append([])
            for line in rawLines:
                strings = line.split(',')
                if strings[-1] == "\n":
                    strings.remove("\n")
                row = [eval(string) for string in strings]
                for k1 in range(numVars):
                    data[k1].append(row[k1])
            # Save each column data with the variable name in a dictionary
            self._variables = {}
            for k1 in range(numVars):
                self.add_data(np.array(data[k1]), variables[k1])
        # }}}

        # Plot {{{
        def Plot(self, xkey, ykeys, **kwargs):
            variables = self.variables
            keys = kwargs.keys()
            # Get variables to plot {{{
            xkey = xkey.replace(" ", "")
            x = variables[xkey]
            if isinstance(ykeys, str):
                ykeys = ykeys.replace(" ", "")
                ys = [variables[ykeys]]
            else:
                ys = []
                for ykey in ykeys:
                    ys.append(variables[ykey.replace(" ", "")])
            numLines = len(ys)
            # }}}
            # Set up fig {{{
            if 'figsize' in keys:
                figsize = kwargs['figsize']
            else:
                cm = 1.0/2.54
                figsize = [8*cm, 6*cm]
            fig, ax = plt.subplots(figsize = figsize, layout = "constrained")
            # }}}
            # Set up colors {{{
            if 'colormap' in keys:
                colormap = kwargs['colormap']
            else:
                colormap = plt.cm.coolwarm
            colors = colormap(np.linspace(0, 1, numLines))
            # }}}
            # Set up markers {{{
            if 'markers' in kwargs:
                markers = kwargs['markers']
                if isinstance(markers, str) or isinstance(markers, int):
                    markers = [markers]*numLines
            else:
                markers = ['-']*numLines
            # }}}
            # Make plot {{{
            def_kwargs = [{}]*numLines
            plot_kwargs = kwargs.get("plot_kwargs", def_kwargs)
            if isinstance(plot_kwargs, dict):
                plot_kwargs = [plot_kwargs]
            for k1 in range(numLines):
                ax.plot(x, ys[k1], markers[k1], color = colors[k1],
                        **plot_kwargs[k1])
            # }}}
            # Run figConf if defined {{{
            if 'figConf' in keys:
                figConf = kwargs["figConf"]
                figConf(ax)
            # }}}
            # Save figure if outname is defined, else show() {{{
            # fig.tight_layout()
            if 'outname' in keys:
                fig.savefig(kwargs['outname'])
                plt.close(fig)
            else:
                plt.show()
            # }}}
            return
        # }}}
    # }}}
# }}}

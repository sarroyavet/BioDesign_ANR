#!/usr/bin/env python3
# Libraries {{{
import os
import sys
import json
import getopt
import time
import psutil
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"],
    "font.size" : 9})
import matplotlib as mpl
from pdb import set_trace
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.linear_model import LinearRegression
from scipy.signal import savgol_filter

# In-house modules
srcFolder = os.path.join(os.getcwd(), '../../src/')
sys.path.append(os.path.join(srcFolder, 'PythonUtilities'))
from datfile import datfile
# }}}

# Get the file with the parameters {{{
try:
    opts, args = getopt.getopt(sys.argv[1:], 'i:x:y:mbcs:z')
    if len(args)>0:
        raise getopt.GetoptError('Too many arguments!')
except getopt.GetoptError:
    print('Wrong call! The excecution can have the following arguments:')
    print('To indicate the input parameters: -i "fileName" (.json)')
    print('To indicate the x limits: -x "[xmin, xmax]"')
    print('To indicate the y limits: -y "[ymin, ymax]"')
    print('To mirror a symmetric field: -m')
    print('Not to plot each pressure: -b')
    print('To correct Q_p: -c')
    print('To smooth report: -s')
    print('To center the x at zero: -z')
    raise
# Initialisation
inParams = False
in_xlim = False
in_ylim = False
mirror = False
pPlot = True
cPlot = False
QpCorrection = False
smooth = False
for opt, arg in opts:
    # File name
    if opt in ['-i']:
        parFile = arg
        inParams = True
    # x lim
    if opt in ['-x']:
        in_xlim = True
        xlim = eval(arg)
    # y lim
    if opt in ['-y']:
        in_ylim = True
        ylim = eval(arg)
    # Mirror
    if opt in ['-m']:
        mirror = True
    # Plot pressure
    if opt in ['-b']:
        pPlot = False
    # Centre at zero
    if opt in ['-z']:
        cPlot = True
    # Q_p correction
    if opt in ['-c']:
        QpCorrection = True
    # Smooth
    if opt in ['-s']:
        smooth = True
        smooth_par = eval(arg)
if not inParams:
    print("Wrong call! The input file with the parameters is not " +
            "present. Use -i 'fileName' (.json)")
    raise
# }}}

# Functions {{{
# Smooth {{{
# Taken from https://stackoverflow.com/questions/20618804/how-to-smooth-a-curve-for-a-dataset
def Smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth
# }}}
# }}}

# Parameters {{{
# Default
with open(os.path.join(srcFolder, 'config.json')) as fle:
    defaultParameters = json.load(fle)
locals().update(defaultParameters)
# Specific
with open(parFile) as fle:
    params = json.load(fle)
with open(parFile, 'r') as inFile:
    params = json.load(inFile)
fileNames.update(params.get("fileNames", {}))

locals().update(fileNames)

if asterFolder[:5] == "/home":
    workDir = asterFolder
else:
    workDir = os.getcwd() + asterFolder
userhome = os.path.expanduser('~')
if not userhome in workDir:
    workDir = os.path.join(os.getcwd(), os.path.dirname(parFile))
# }}}

# Report {{{
report = datfile(workDir + repoFile)
report = report.datablocks['0']

# Set up variables {{{
contArea = report.variables["cont_area"]
p_c_max = report.variables["p_c_max"]
a_c0 = contArea[0]
p_c_max0 = p_c_max[0]
normalised_area = contArea/a_c0
normalised_p_c_max = p_c_max/p_c_max0
report.add_data(normalised_area, "norm_cont_area")

mea_kesc = report.variables["mea_kesc"]
mea_kmai = report.variables["mea_kmai"]
mea_kdif = report.variables["mea_kdif"]
mea_krelDif = 1.0 - (mea_kdif)/(mea_kesc + mea_kmai)
report.add_data(mea_krelDif, "mea_krelDif")

# Correction of Q_p
if QpCorrection:
    old_Q_p = report.variables["Q_p"]
    report.variables["Q_p"] = 2.0*old_Q_p - 1.0

if smooth:
    time = report.variables["Time"]
    rough_Q_p = report.variables["Q_p"]
    rough_normalised_area = report.variables["norm_cont_area"]
    report.variables["Q_p"] = savgol_filter(rough_Q_p, smooth_par, 3)
    report.variables["norm_cont_area"] = savgol_filter(rough_normalised_area,
                                                       smooth_par, 3)
    smooth_Q_p = report.variables["Q_p"]
    print("Q_p:\trough(min, max): {}, {}".format(rough_Q_p.min(), rough_Q_p.max()))
    print("\tsmooth(min, max): {}, {}".format(smooth_Q_p.min(), smooth_Q_p.max()))
# }}}

# Figure size
cm = 1.0/2.54
fig_w = 5*cm
fig_h = fig_w*3/4
figsize = [fig_w, fig_h]

# Time vs Q_p {{{
def figConf_time_Qp(ax):
    ax.set_xlabel(r"$\mathrm{Time}$")
    ax.set_ylabel(r"$Q_p$")
    #ax.set_ylim([0.85, 0.97])
    return ax

plotKey = "Q_p"
outname = workDir + '/report_' + plotKey + '.pdf'
report.Plot("Time", plotKey, outname = outname,
        figConf = figConf_time_Qp)
# }}}

# Time vs maximum contact pressure {{{
def figConf_time_Pcmax(ax):
    ax.set_xlabel(r"$\mathrm{Time}$")
    ax.set_ylabel(r"$\frac{p_c^\mathrm{max}}{p_c^\mathrm{max}(0)}$")
    #ax.set_ylim([0.96, 2.5])
    return ax

plotKey = "p_c_max"
outname = workDir + '/report_' + plotKey + '.pdf'
report.Plot("Time", plotKey, outname = outname,
        figConf = figConf_time_Pcmax)
# }}}

# Time vs normalised contact area {{{
def figConf_time_Na(ax):
    ax.set_xlabel(r"$\mathrm{Time}$")
    ax.set_ylabel(r"$\frac{l_c}{l_{c_{0}} }$")
    ax.set_ylabel(r"${l_c}/{l_{c_{0}} }$")
    #ax.set_ylim([0.96, 2.5])
    return ax

plotKey = "norm_cont_area"
outname = workDir + '/report_' + plotKey + '.pdf'
report.Plot("Time", plotKey, outname = outname,
        figConf = figConf_time_Na, colormap = plt.cm.gray,
        figsize = figsize)
# }}}

# Normalised area vs Q_p {{{
def figConf_Na_Qp(ax):
    ax.set_xlabel(r"$\frac{l_c}{l_{c_{0}} }$")
    ax.set_xlabel(r"${l_c}/{l_{c_{0}} }$")
    ax.set_ylabel(r"$Q_p$")
    #ax.set_xlim([0.96, 2.5])
    #ax.set_ylim([0.85, 0.97])
    return ax

plotKey = "Q_p"
outname = workDir + '/report_Na_' + plotKey + '.pdf'
report.Plot("norm_cont_area", plotKey, outname = outname,
        figConf = figConf_Na_Qp, colormap = plt.cm.gray,
        figsize = figsize),
# }}}

# Curvatures "mea_kesc", "mea_kmai", "mea_kdif" {{{
def figConf_time_kappa(ax):
    ax.set_xlabel(r"$\mathrm{Time}$")
    ax.legend([r"$||\kappa_\mathrm{esc}||$",
        r"$||\kappa_\mathrm{mai}||$",
        r"$||\kappa_\mathrm{esc} - \kappa_\mathrm{mai}||$"])
    return ax

plotKey = ["mea_kesc", "mea_kmai", "mea_kdif"]
outname = workDir + '/report_time_' + "kappa" + '.pdf'
report.Plot("Time", plotKey, outname = outname,
        figConf = figConf_time_kappa)
# }}}

# Curvatures "mea_krelDif" {{{
def figConf_time_kappa_relDif(ax):
    ax.set_xlabel(r"$\mathrm{Time}$")
    ax.set_ylabel(r"$1 - ||(\kappa_\mathrm{esc} - \kappa_\mathrm{mai})||/(||\kappa_\mathrm{esc}|| + ||\kappa_\mathrm{mai}||)$")
    return ax

plotKey = "mea_krelDif"
outname = workDir + '/report_time_' + plotKey + '.pdf'
report.Plot("Time", plotKey, outname = outname,
        figConf = figConf_time_kappa_relDif)
# }}}

# Curvatures "alpha_i" {{{
def figConf_time_alpha(ax):
    ax.set_xlabel(r"$\mathrm{Time}$")
    ax.set_ylabel(r"$\alpha$")
    ax.set_yscale("log")
    return ax

plotKey = "alpha_i"
outname = workDir + '/report_time_' + plotKey + '.pdf'
report.Plot("Time", plotKey, outname = outname,
        figConf = figConf_time_alpha)
# }}}
# }}}

# Pressure  {{{
# Figure configuration
def figConf(ax):
    ax.set_xlim([0.0, 1.1*max(normalised_area)])
    ax.set_ylim([0.0, 1.1*max(normalised_p_c_max)])
    if in_xlim:
        ax.set_xlim(xlim)
    if in_ylim:
        ax.set_ylim(ylim)
    ax.legend([r"$p_c/p_{c_0}^\mathrm{max}$",
               r"$p_c^\mathrm{ref}/p_{c_0}^\mathrm{max}$"])
    ax.set_xlabel(r"$l_c/b_0$")
    return ax
# Load data
pressure = datfile(workDir + presResuFile)
# Plot data
for key, p_i in pressure.datablocks.items():
    p_i.add_data(p_i.variables["l_c[mm]"]/a_c0, "norm(l_c)")
    p_i.add_data(p_i.variables["p_c[GPa]"]/p_c_max0, "norm(p_c)")
    p_i.add_data(p_i.variables["p_r[GPa]"]/p_c_max0, "norm(p_r)")
    if pPlot:
        p_i.Plot("norm(l_c)", ["norm(p_c)", "norm(p_r)"],
            outname = workDir + resuFolder + '/PressureTime' + key + '.png',
            figsize = figsize, figConf = figConf)
    if cPlot:
        x_i = p_i.variables["l_c[mm]"]
        x_i_xmin = x_i.min()
        x_i_xmax = x_i.max()
        x_i_xmean = (x_i_xmin + x_i_xmax)/2.0
        p_i.add_data(x_i - x_i_xmean, "cl_c")
# }}}

# Pressure all in one {{{
# Set size
cm = 1.0/2.54
fig_w = 5*cm
fig_h = fig_w*3/4
figsize = [fig_w, fig_h]
# Make figure
fig, ax = plt.subplots(1, 2, figsize = figsize, layout = "constrained",
        gridspec_kw={'width_ratios': [20, 1]})
# Colors
colormap = plt.cm.YlOrRd
# Add lines
try:
    with open(workDir + "/plot_time_list.json", "r") as fle:
        time_list = json.load(fle)["time_list"]
except:
    time_list = pressure.datablocks.keys()
numLines = len(time_list)
colorRange = np.linspace(0.1, 1, numLines)
colors = colormap(colorRange)
k1 = 0
for key in time_list:
    p_i = pressure.datablocks[key]
    if cPlot:
        x = p_i.variables["cl_c"]
    else:
        x = p_i.variables["l_c[mm]"]
    y = p_i.variables["norm(p_c)"]
    if mirror:
        x = np.concatenate((-np.flip(x), x))
        y = np.concatenate((np.flip(y), y))
    ax[0].plot(x, y, color = colors[k1])
    k1 += 1
ax[0].set_xlabel(r"$l_{c}$")
ax[0].set_ylabel(r"$p_c/p_{c_0}^\mathrm{max}$")
if in_xlim:
    ax[0].set_xlim(xlim)
if in_ylim:
    ax[0].set_ylim(ylim)
# Add colour bar
cb   = mpl.colorbar.ColorbarBase(ax[1], cmap = colormap,
        # norm = norm,
        orientation = "vertical")
ax[1].set_ylim([colorRange.min(), 1.0])
ax[1].set_yticks([colorRange.min(), 1.0])
ax[1].set_yticklabels([r'$t_0$', r'$t_f$'])
# Save
# fig.tight_layout()
fig.savefig(os.path.join(workDir, "pressure.png"), dpi = 360)
fig.savefig(os.path.join(workDir, "pressure.pdf"))
fig.savefig(os.path.join(workDir, "pressure.svg"))
# }}}

#!/usr/bin/env /proj/sot/ska/bin/python

#########################################################################################################
#                                                                                                       #
#   create_count_rate_plots.py: create count rate and background rate plots                             #
#                                                                                                       #
#           author: t. isobe (tisobe@cfa.harvard.edu)                                                   #
#                                                                                                       #
#           last update: Jun 01, 2017                                                                   #
#                                                                                                       #
#########################################################################################################

import os
import sys
import re
import string
import time
import numpy

import matplotlib as mpl

if __name__ == '__main__':

    mpl.use('Agg')

from pylab import *
import matplotlib.pyplot       as plt
import matplotlib.font_manager as font_manager
import matplotlib.lines        as lines

import mpld3
from mpld3 import plugins, utils

#
#--- reading directory list
#
path = '/data/aschrc6/wilton/isobe/Project8/ArLac/Scripts3/house_keeping/dir_list'

f    = open(path, 'r')
data = [line.strip() for line in f.readlines()]
f.close()

for ent in data:
    atemp = re.split(':', ent)
    var  = atemp[1].strip()
    line = atemp[0].strip()
    exec "%s = %s" %(var, line)
html_top = html_top.replace('#', ':')
#
#--- append a path to a private folder to python directory
#
sys.path.append(bin_dir)
sys.path.append(mta_dir)
#
#--- converTimeFormat contains MTA time conversion routines
#
import convertTimeFormat    as tcnv
import mta_common_functions as mcf
import find_moving_average  as fmv
import robust_linear        as robust

hrci_offset = [[0,0], [2,0], [0,2], [-2,0], [0,-2], [2,2], [-2,2], [-2,-2], [2,-2], [4,0], [0,4],[-4,0], [0,-4], [6,0], [0,6], [-6,0], [0,-6], [10,10], [-10,10], [-10,-10], [10,-10]]

hrcs_offset = [[0,0], [2,0], [0,2], [-2,0], [0,-2], [2,2], [-2,2], [-2,-2], [2,-2], [4,0], [-4,0], [4,2], [-4,2], [-4, -2], [4,-2], [6,0], [-6,0], [10,2], [-10,2], [-10,-2],[10,-2]]

mday = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

#---------------------------------------------------------------------------------------------------
#-- create_count_rate_plots: creates create count rate and background rate plots                  --
#---------------------------------------------------------------------------------------------------

def create_count_rate_plots(chk):
    """
    creates create count rate and background rate plots
    input: none but read data from <data_dir> 
    output: plots in <web_dir>/Plots3/Count Rate/   
    """
#
#--- check whether new data came in
#
    ifile = house_keeping + 'hrci_new'
    sfile = house_keeping + 'hrcs_new'
    if chk == "":
        if (not os.path.isfile(ifile)) and (not os.path.isfile(sfile)):
            exit(1)

    xmin  = 1999
    xmax  = time.strftime('%Y', time.localtime())
    xmax  = int(float(xmax)) + 1
    xname = 'Time (year)'
    ymin  = 0
    ymax  = 20
    yname = 'Count/Sec/Pix * 10^6'
    yname2= 'Count/Sec'

    for col in ['pi', 'samp']:
        for inst in ['s', 'i']:
            infile  = data_dir + col  + '_list_' + inst
            [time_list, cnt_list, err_list, bkg_list, berr_list, obsid_list] = read_file(infile, inst, col)

            for pos in range(0, 21):
                outname = web_dir + 'Plots3/Count_rates/' + col +'_' +   inst + '_' + str(pos) + '.png'
#
#--- if there are not enough data, show a blank page
#
                if len(time_list) < 3:
                    cmd = 'cp ' + house_keeping + 'no_data.png ' + outname
                    os.system(cmd)
                else:
#
#--- main plot
#
                    if time_list[pos] > 3:
                        plot_panel(xmin, xmax, ymin, ymax, time_list[pos], cnt_list[pos], err_list[pos], bkg_list[pos], berr_list[pos], xname, yname)
                        cmd    = 'mv out_plot.png ' + outname
                    else:
                        cmd    = 'cp ' +house_keeping + 'no_data.png ' + outname

                    os.system(cmd)

#---------------------------------------------------------------------------------------------------
#-- read_file: read arlac data file                                                               --
#---------------------------------------------------------------------------------------------------

def read_file(infile, inst, col):
    """
    read arlac data file
    input:  infile      --- data file name
    output: time_list   --- a list of time data
            cnt_list    --- a list of the count rates
            bkg_list    --- a list of the background rates
    """
    if inst.lower() == 'i':
        offset_list = hrci_offset
    else:
        offset_list = hrcs_offset
    
    data  = read_data(infile)
    if data == []:
        return False
    
    time_list  = []
    cnt_list   = []
    err_list   = []
    bkg_list   = []
    berr_list  = []
    obsid_list = []
    for k in range(0, 21):
        time_list.append([])
        cnt_list.append([])
        err_list.append([])
        bkg_list.append([])
        berr_list.append([])
        obsid_list.append([])
    
    for ent in data:
        atemp = re.split('\s+', ent)
        time  = get_frac_year(atemp[2])
        if time < 2000:
            continue
#
#--- find where to put the data; 0 to 20 data slots
#
        yoffset = round(float(atemp[10]))
        zoffset = round(float(atemp[11]))
        for k in range(0, 21):
            if offset_list[k][0] == yoffset and  offset_list[k][1] == zoffset:
                pos = k
                break
    
        avg    = float(atemp[12])
        err    = float(atemp[13])
        crate  = float(atemp[15])
        brate  = float(atemp[16])
        carea  = float(atemp[26])
        barea  = float(atemp[27])

        expos  = float(atemp[3])
#
#--- create count rate errors 
#--- the count rate is per seconds; so to compute the count error, bring back the total counts
#--- by multiplying the exposure time before computing the count error.
#
        err    = math.sqrt(crate * expos) / expos 
        crate -= brate * carea / barea
        berr   = math.sqrt(brate * expos) / expos / barea
        berr  *= 10**6
        brate  = brate / barea
        brate *= 10**6

        if avg == -999 or err == -999:
            continue
        else:
            time_list[pos].append(time)
            cnt_list[pos].append(crate)
            err_list[pos].append(err)
            bkg_list[pos].append(brate)
            berr_list[pos].append(berr)
            obsid_list[pos].append(atemp[1])
    
    return [time_list, cnt_list, err_list,  bkg_list, berr_list,  obsid_list]

#---------------------------------------------------------------------------------------------------
#-- get_frac_year: convert date into a fractional year format                                     --
#---------------------------------------------------------------------------------------------------

def get_frac_year(date):
    """
    convert date into a fractional year format
    input: date    --- date in the format of <yyyy>-<mm>-<dd>T<hh>:<mm>:<ss>
    output: fyear  --- date in fractional year. less than a day will be ignored
    """
    atemp = re.split('T',  date)
    btemp = re.split('\-', atemp[0])
    fyear = float(btemp[0])
    mon   = int(float(btemp[1]))
    day   = int(float(btemp[2]))

    if mon == 2 and tcnv.isLeapYear(fyear) == 1:
        lday = mday[mon-1] + 1
    else:
        lday = mday[mon-1]

    mon   += day / lday
    fyear += mon / 12.0

    return fyear

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

def read_data(infile, remove=0):

    try:
        f    = open(infile, 'r')
        data = [line.strip() for line in f.readlines()]
        f.close()

        if remove == 1:
            mcf.rm_file(infile)
    except:
        data = []

    return data

#---------------------------------------------------------------------------------------------------
#-- plot_single_panel: plot a single data set on a single panel                                  ---
#---------------------------------------------------------------------------------------------------

def plot_panel(xmin, xmax, ymin, ymax, xdata, ydata, yerr, ydata2, yerr2, xname, yname):

    """
    plot a single data set on a single panel
    Input:  xmin    --- min x
            xmax    --- max x
            ymin    --- min y
            ymax    --- max y
            xdata   --- independent variable
            ydata   --- dependent variable
            yerr    --- error of the ydata
            ydata2  --- dependent variable 2
            xname   --- x axis label
            ynane   --- y axis label
    Output: png plot named <outname>
    """

    colorList = ('blue', 'green', 'red', 'aqua', 'lime', 'fuchsia', 'maroon', 'black', 'yellow', 'olive')
    fsize = 12 
    lsize = 2
    resolution = 100.0
    connect = 0
#
#--- fit line --- use robust method
#
    xx = []
    for m in range(0, len(xdata)):
        xx.append(xdata[m] - 1999)

    [sint, slope, sa, serr] =  line_fit(xx, ydata, yerr)
    lint  =  '%2.3f' % (round(sint,  3))

    if slope < 0:
        sign = -1
    else:
        sign = 1

    lslope = '%2.3f' % (round(abs(slope), 3))
    lerr   = '%2.3f' % (round(serr,  3))
#
#--- second fit
#
    [sint2, slope2, sa, serr2] =  line_fit(xx, ydata2, yerr)
    if sint2 < 0.0001:
        lint2  =  '%2.4f' % (round(sint2,  4))
    else:
        lint2  =  '%2.3f' % (round(sint2,  3))

    if slope2 < 0:
        sign2 = -1
    else:
        sign2 = 1

    lslope2 = '%2.3f' % (round(abs(slope2*1e3), 3))
    lerr2   = '%2.3f' % (round(serr2*1e3,  3))
#
#--- close everything opened before
#
    plt.close('all')
#
#--- set font size
#
    mpl.rcParams['font.size'] = fsize
    props = font_manager.FontProperties(size=9)
#
#---- the first panel
#
#--- set plotting range
#
    ychk = ydata + ydata2
    [ymin, ymax] = set_y_range(ychk)
    ax1  = plt.subplot(211)
    ax1.set_autoscale_on(False)
    ax1.set_xbound(xmin,xmax)
    ax1.set_xlim(xmin=xmin, xmax=xmax, auto=False)
    ax1.set_ylim(ymin=ymin, ymax=ymax, auto=False)
#
#--- plot data
#
    p, = plt.plot(xdata, ydata, color=colorList[3], marker='o', markersize=4, lw = connect)
    plt.errorbar(xdata,  ydata, yerr=yerr, lw = 0, elinewidth=1)
#
#--- plot fitted line
#
    start = sint + slope * (xmin - 1999)
    stop  = sint + slope * (xmax - 1999)
#
#--- label axes
#
    yname2 = 'Counts/Sec'
    plt.ylabel(yname2, size=fsize)
#
#--- add what is plotted on this plot
#
    xdiff = xmax - xmin
    xpos  = xmin + 0.1 * xdiff
    ydiff = ymax - ymin
    ypos  = ymax - 0.09 * ydiff

    label = 'Net Rate'
    plt.text(xpos, ypos, label, size=fsize)
#
#---- second panel
#
#
#--- set plotting range
#
    [ymin, ymax] = set_y_range(ydata2)
    ax2  = plt.subplot(212)
    ax2.set_autoscale_on(False)
    ax2.set_xbound(xmin,xmax)
    ax2.set_xlim(xmin=xmin, xmax=xmax, auto=False)
    ax2.set_ylim(ymin=ymin, ymax=ymax, auto=False)
#
#--- plot data
#
    p, = plt.plot(xdata, ydata2, color=colorList[4], marker='o', markersize=4, lw = connect)
    plt.errorbar(xdata,  ydata2, yerr=yerr2, lw = 0, elinewidth=1)
#
#--- plot fitted line
#
    start = sint2 + slope2 * (xmin - 1999)
    stop  = sint2 + slope2 * (xmax - 1999)
#
#--- label axes
#
    plt.xlabel(xname, size=fsize)
    plt.ylabel(yname, size=fsize)
#
#--- add what is plotted on this plot
#
    xdiff = xmax - xmin
    xpos  = xmin + 0.1 * xdiff
    ydiff = ymax - ymin
    ypos  = ymax - 0.08 * ydiff

    label = 'Background Rate'
    plt.text(xpos, ypos, label, size=fsize)
#
#--- set the size of the plotting area in inch
#
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(10.0, 6.0)
    outname = 'zplot.png'
    plt.savefig(outname, format='png', dpi=resolution)

    cmd = 'convert zplot.png -trim out_plot.png'
    os.system(cmd)
    mcf.rm_file('zplot.png')

    plt.close('all')


#---------------------------------------------------------------------------------------------------
#-- set_y_range: set min and max range for the plotting                                           --
#---------------------------------------------------------------------------------------------------

def set_y_range(data):
    """
    set min and max range for the plotting
    input:  data    --- data list
    output: dmin    --- min of the plotting range
            dmax    --- max of the plotting range
    """

    if len(data) == 0:
        return[0, 1]

    dmin = min(data)
    dmax = max(data)
    diff = dmax - dmin

    dmax = dmax + 0.3 * diff
    dmin = dmin - 0.1 * diff 
    if diff > 1:
        dmax = int(10 * dmax) + 5
        dmin = int(10 * dmin) - 5
        dmax = 0.1 * dmax
        dmin = 0.1 * dmin

    if dmin < 0:
        dmin = 0

    return [dmin, dmax]


#---------------------------------------------------------------------------------------------------
#-- line_fit: fit a weighted linear line fit                                                      --
#---------------------------------------------------------------------------------------------------

def line_fit(x, y, e):
    """
    fit a weighted linear line fit
    input:  x       --- independent data
            y       --- dependent data
            e       --- y error
    output: a       --- intercept
            b       --- slope
            siga    --- error on the intercept
            sigb    --- error on the slope
    """

    suma  = 0
    sumx  = 0
    sumy  = 0
    sumx2 = 0
    sumy2 = 0
    sumxy = 0

    dlen = len(x)
    if dlen < 3:
        return [0, 0, 0, 0]

    for k in range(0, dlen):
        try:
            weight = 1.0 / e[k]**2
        except:
            weight = 1.0
        suma  += weight
        sumx  += weight * x[k]
        sumy  += weight * y[k]
        sumx2 += weight * x[k] * x[k]
        sumy2 += weight * y[k] * y[k]
        sumxy += weight * x[k] * y[k]

    delta = suma * sumx2 - sumx* sumx
    a     = (sumx2 * sumy - sumx * sumxy) / delta
    b     = (sumxy * suma - sumx * sumy ) / delta
    if dlen <= 2:
        siga = 0
        sigb = 0
    else:
        var   = (sumy2 + a * a * suma + b * b * sumx2 - 2.0 *(a * sumy + b * sumxy - a * b * sumx)) / (len(x) -2)
        siga  = math.sqrt(var * sumx2 / delta)
        sigb  = math.sqrt(var * suma  / delta)

    return [a, b, siga, sigb]


#---------------------------------------------------------------------------------------------------

if __name__ == "__main__":

    if len(sys.argv) == 2:
        chk = sys.argv[1]
    else:
        chk = ''

    create_count_rate_plots(chk)


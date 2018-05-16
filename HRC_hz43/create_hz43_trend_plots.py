#!/usr/bin/env /proj/sot/ska/bin/python

#########################################################################################################
#                                                                                                       #
#       create_hz43_trend_plots.py: creates binned trend plots                                          #
#                                                                                                       #
#           author: t. isobe (tisobe@cfa.harvard.edu)                                                   #
#                                                                                                       #
#           last update: May 26, 2017                                                                   #
#                                                                                                       #
#########################################################################################################

import os
import sys
import re
import string
import time

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
path = '/data/aschrc6/wilton/isobe/Project8/HZ43/Scripts/house_keeping/dir_list'

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

mday = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

year_div = 2012

exclude_list= ['1009']
    
#---------------------------------------------------------------------------------------------------
#-- create_plots: creates binned trend plots                                                      --
#---------------------------------------------------------------------------------------------------

def create_plots():
    """
    creates binned trend plots 
    input: none but read data from <data_dir> 
    output: plots in <web_dir>/Plots/   such as samp_n_s_14.png or pi_n_s_14.png
    """

    xmin  = 1999
    xmax  = time.strftime('%Y', time.localtime())
    xmax  = int(float(xmax)) + 1
    xname = 'Time (year)'

    for col in ['pi_', 'samp_']:
        if col == 'pi_':
            yname = 'PI'
            ymin  = 40
            ymax  = 190
        else:
            yname = 'Scaled SumAmp'
            ymin  = 0
            ymax  = 150

        for inst in ['s', 'i']:
#
#---  wavelength part plots
#
            for sign in ['n_', 'p_']:
                oname = data_dir + col + inst + '_' + sign + 'fitting_results'
                fo    = open(oname, 'w')

                for step in range(0, 20):
                    infile       = data_dir + col + sign + 'list_' + str(step) + '_' +  inst
                    outname      = web_dir + 'Plots/' + col + sign + inst + '_' + str(step) + '.html'
                    profile_page = web_dir + 'Plots/Profile_page/' + col +  sign + inst + '_' + str(step) + '.html'
#
#--- read data
#
                    try:
                        [time_list, data_list, err_list, obsid_list] = read_file(infile)
                    except:
                        cmd = 'cp ' + house_keeping + 'no_data.html ' + outname
                        os.system(cmd)
                        continue
#
#--- if there are not enough data, show a blank page
#
                    if len(time_list) < 3:
                        cmd = 'cp ' + house_keeping + 'no_data.html ' + outname
                        os.system(cmd)
                    else:
                        [info_list, hlink_list]  = create_info_link(obsid_list, col, sign, step)
#
#--- main plot
#
                        [fig, lint, lslope, lerr, lint2, lslope2, lerr2] \
                                    = plot_single_panel(xmin, xmax, ymin, ymax, time_list, data_list, err_list, \
                                                        info_list, xname, yname, label='')
#
#--- make html page for it
#
                        make_html_page(fig, outname, col, inst, sign, step, profile_page)

                        oline = str(step) + '\t' + lint  + '\t' + lslope  + '\t' + lerr  + '\t'
                        oline = oline     + '\t' + lint2 + '\t' + lslope2 + '\t' + lerr2 + '\n'
                        fo.write(oline)
#
#--- thumb nail plot
#
                        create_thumb_nail_plots(xmin, xmax, ymin, ymax, time_list, data_list, err_list, \
                                                    info_list, xname, yname, label='')

                        p_name = web_dir + 'Plots/Thumb_plots/' + col + sign + inst + '_' + str(step) + '_thumb_plot.png'
                        cmd    = 'mv out_plot.png ' + p_name
                        os.system(cmd)

                        if inst =='i':
                            title = 'HRC I ' +  yname + ': Location ' + str(step)
                        else:
                            title = 'HRC S ' +  yname + ': Location ' + str(step)

                        make_profile_page(time_list, hlink_list, title, profile_page, outname)
                fo.close()
#
#--- center part plots
#
            infile       = data_dir + col + 'center_list_' + inst
            outname      = web_dir + 'Plots/'  + col + 'center_list_' + inst + '.html'
            profile_page = web_dir + 'Plots/Profile_page/' + col +  'center_list_' + inst  + '.html'
            try:
                [time_list, data_list, err_list, obsid_list] = read_file(infile)
            except:
                continue 

            [info_list, hlink_list] = create_info_link(obsid_list, col, sign, 'center')
            [fig, lint, lslope, lerr, lint2, lslope2, lerr2] \
                                = plot_single_panel(xmin, xmax, ymin, ymax, time_list, data_list, err_list, \
                                                    info_list, xname, yname, label='')

            make_html_page(fig, outname, col, inst, 'center', 'center', profile_page)

            create_thumb_nail_plots(xmin, xmax, ymin, ymax, time_list, data_list, err_list, \
                                        info_list, xname, yname, label='')
            p_name = web_dir + 'Plots/Thumb_plots/' + col + inst + '_' + 'center_thumb_plot.png'
            cmd    = 'mv out_plot.png ' + p_name
            os.system(cmd)

            if inst == 'i':
                title = 'HRC I' + yname + 'Center'
            else:
                title = 'HRC S' + yname + 'Center'

            make_profile_page(time_list, hlink_list, title, profile_page, outname)

            oname = data_dir + col + inst + '_center_fitting_results'
            fo    = open(oname, 'w')
            oline = lint  + '\t' + lslope + '\t' + lerr + '\t'
            oline = oline + '\t' + lint2  + '\t' + lslope2 + '\t' + lerr2 + '\n'
            fo.write(oline)
            fo.close()


#---------------------------------------------------------------------------------------------------
#-- read_file: read hz 43 data file                                                               --
#---------------------------------------------------------------------------------------------------

def read_file(infile):
    """
    read hz 43 data file
    input:  infile      --- data file name
    output: time_list   --- a list of time data
            data_list   --- a list of value data either sumamp or pi
            err_list    --- a list of the error of the value
            obsid_list  --- a list of obsid of the observations
    """

    data  = read_data(infile)
    if data == []:
        return False

    time_list  = []
    data_list  = []
    err_list   = []
    obsid_list = []

    for ent in data:
        atemp = re.split('\s+', ent)
        time  = get_frac_year(atemp[2])
        if time < 2000:
            continue

        obsid = str(atemp[1].strip())
        if  obsid in exclude_list:
            continue

        avg   = float(atemp[4])
        err   = float(atemp[5])
        med   = float(atemp[6])
        if avg == -999 or err == -999 or pi == -999:
            continue
        else:
            time_list.append(time)
            data_list.append(avg)
            err_list.append(err)
            obsid_list.append(atemp[1])

    zmin = min(data_list)
    zmax = max(data_list)

    return [time_list, data_list, err_list, obsid_list]

#---------------------------------------------------------------------------------------------------
#-- create_info_link: create a list of links used by an interactive plot                          --
#---------------------------------------------------------------------------------------------------

def create_info_link(obsid_list, col, sign, step):
    """
    create a list of links used by an interactive plot
    input:  obsid_list  --- a list of obsids
            col         --- pi or samp
            sign        --- n or p (negative or positive sides of the arms)
            step        --- position in the oarm
    output: outlist     --- a list of html address to a distribution plot
    """

    html_plot = html_top + '/Plots/Indivisual_Plots/'
    outlist = []
    hlist   = []
    for obsid in obsid_list:
        if step == 'center':
            ofile =  html_plot + str(obsid)   + '/' + col +  'center_list_vfit.png' 
            olink = '<p><img src="' + ofile + '" width= 600px></p>'
        else:
            ofile =  html_plot +  str(obsid)   + '/' + col +  sign + 'list_' +  str(step) + '_vfit.png' 
            olink = '<p><img src="' + ofile + '" width= 600px></p>'

        hlist.append(ofile)
        outlist.append(olink)

    return [outlist, hlist]

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

def plot_single_panel(xmin, xmax, ymin, ymax, xdata, ydata, yerror, info_list,  xname, yname, \
                      label, fsize = 9, psize = 60.0, marker = 's', pcolor =0, lcolor=0,\
                      lsize=1, resolution=100, linefit=1, connect=0):

    """
    plot a single data set on a single panel
    Input:  xmin    --- min x
            xmax    --- max x
            ymin    --- min y
            ymax    --- max y
            xdata   --- independent variable
            ydata   --- dependent variable
            yerror  --- error in y axis; if it is '0', no error bar
            info_list   --- a list of information to be display on the each data point
            xname   --- x axis label
            ynane   --- y axis label
            label   --- a text to indecate what is plotted
            fsize   ---  font size, default = 9
            psize   ---  plotting point size, default = 2.0
            marker  ---  marker shape, defalut = 'o'
            pcolor  ---  color of the point, default= 7 ---> 'black'
            lcolor  ---  color of the fitted line, default = 7 ---> 'black'
                colorList = ('blue', 'red', 'green', 'aqua', 'fuchsia','lime', 'maroon', 'black', 'olive', 'yellow')
            lsize:      fitted line width, defalut = 1
            resolution-- the resolution of the plot in dpi
            linefit  --- if it is 1, fit a line estimated by robust method
            connect  --- if it is > 0, lsize data point with a line, the larger the value thinker the line
    Output: png plot named <outname>
    """
    colorList = ('blue', 'green', 'red', 'aqua', 'lime', 'fuchsia', 'maroon', 'black', 'yellow', 'olive')
#
#--- this css is used for the pop up page
#
    css = """
        body{
            width:600px;
            height:300px;
        }
        p{
            text-align:center;
        }
    """
#
#--- fit line --- use robust method
#
    if linefit == 1:
        xx = []
        for m in range(0, len(xdata)):
            xx.append(xdata[m] - 1999)

        [sint, slope, sa, serr, sint2, slope2, sa2, serr2] =  line_fit_prep(xx, ydata, yerror)
        lint  =  '%2.3f' % (round(sint,  3))
        lint2 =  '%2.3f' % (round(sint2, 3))

        if slope < 0:
            sign = -1
        else:
            sign = 1

        if slope2 < 0:
            sign2 = -1
        else:
            sign2 = 1

        lslope  = '%2.3f' % (round(abs(slope), 3))
        lerr    = '%2.3f' % (round(serr,  3))

        lslope2 = '%2.3f' % (round(abs(slope2), 3))
        lerr2   = '%2.3f' % (round(serr2,  3))
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
#--- set plotting range
#
    fig, ax = plt.subplots()
    ax.set_xlim(xmin=xmin, xmax=xmax, auto=False)
    ax.set_ylim(ymin=ymin, ymax=ymax, auto=False)
    ax.yaxis.set_label_coords(-0.05, 0.5)
    ax.xaxis.set_label_coords(0.5, -0.05)
    ax.set_xlabel(xname)
    ax.set_ylabel(yname)
#
#--- plot data
#
    pv = ax.plot(xdata, ydata, color=colorList[pcolor], marker=marker, markersize=10, lw = connect)
    plugins.connect(fig, mpld3.plugins.PointHTMLTooltip(pv[0], info_list, css=css, hoffset=-350))
#
#--- plot error bar
#
    if yerror != 0:
        plt.errorbar(xdata, ydata, yerr=yerror, lw = 0, elinewidth=1)
#
#--- plot fitted line
#
    if linefit == 1:
        start = sint + slope * (xmin - 1999)
        if slope2 == 0:
            stop  = sint + slope * (xmax - 1999)
        else:
            stop  = sint + slope * (year_div - 1999)

        plt.plot([xmin, year_div], [start, stop], color=colorList[lcolor], lw = lsize)
#
        if slope2 != 0:
            start = sint2 + slope2 * (year_div - 1999)
            stop  = sint2 + slope2 * (xmax - 1999)
            plt.plot([year_div, xmax], [start, stop], color=colorList[lcolor], lw = lsize)
#
#--- add what is plotted on this plot
#
    xdiff = xmax - xmin
    xpos  = xmin + 0.1 * xdiff
    ydiff = ymax - ymin
    ypos  = ymax - 0.08 * ydiff
    ypos2 = ymax - 0.12 * ydiff

    if linefit == 1:
        if slope2 == 0:
            if sign >  0:
                atext = 'Slope: '  + str(lslope) 
                atext = atext + '+/-' + lerr
            else:
                atext = 'Slope: -'  + str(lslope) 
                atext = atext + '+/-' + lerr

            plt.text(xpos, ypos,  atext,  size=fsize)
        else:
            if sign >  0:
                atext = 'Slope (before ' + str(year_div) + '): '  + str(lslope) 
                atext = atext + '+/-' + lerr
            else:
                atext = 'Slope (before ' + str(year_div) + '):  -'  + str(lslope) 
                atext = atext + '+/-' + lerr 
    
            if sign2 > 0:
                atext2 = 'Slope (after ' + str(year_div) + '): '  + str(lslope2) 
                atext2 = atext2 + '+/-' + lerr2
            else:
                atext2 = 'Slope (after ' + str(year_div) + '): -'  + str(lslope2) 
                atext2 = atext2 + '+/-' + lerr2

            plt.text(xpos, ypos,  atext,  size=fsize)
            plt.text(xpos, ypos2, atext2, size=fsize)
#
#--- set the size of the plotting area in inch
#
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(10.5, 6.0)
    fig.tight_layout()

    plt.close('all')

    if sign < 0:
        lslope = '-' + lslope
    if sign2 < 0:
        lslope2 = '-' + lslope2

    return [fig, lint, lslope, lerr, lint2, lslope2, lerr2]



#---------------------------------------------------------------------------------------------------
#-- create_thumb_nail_plots: plot a single data set on a single panel for thumbmail plot          --
#---------------------------------------------------------------------------------------------------

def create_thumb_nail_plots(xmin, xmax, ymin, ymax, xdata, ydata, yerror, info_list,  xname, yname, \
                      label, fsize = 0, psize = 30.0, marker = 's', pcolor =0, lcolor=0,\
                      lsize=1, resolution=100, linefit=1, connect=0):

    """
    plot a single data set on a single panel for thumbmail plot
    Input:  xmin    --- min x
            xmax    --- max x
            ymin    --- min y
            ymax    --- max y
            xdata   --- independent variable
            ydata   --- dependent variable
            yerror  --- error in y axis; if it is '0', no error bar
            info_list   --- a list of information to be display on the each data point
            xname   --- x axis label
            ynane   --- y axis label
            label   --- a text to indecate what is plotted
            fsize   ---  font size, default = 9
            psize   ---  plotting point size, default = 2.0
            marker  ---  marker shape, defalut = 'o'
            pcolor  ---  color of the point, default= 7 ---> 'black'
            lcolor  ---  color of the fitted line, default = 7 ---> 'black'
                colorList = ('blue', 'red', 'green', 'aqua', 'fuchsia','lime', 'maroon', 'black', 'olive', 'yellow')
            lsize:      fitted line width, defalut = 1
            resolution-- the resolution of the plot in dpi
            linefit  --- if it is 1, fit a line estimated by robust method
            connect  --- if it is > 0, lsize data point with a line, the larger the value thinker the line
    Output: png plot named <outname>
    """
    colorList = ('blue', 'green', 'red', 'aqua', 'lime', 'fuchsia', 'maroon', 'black', 'yellow', 'olive')
#
#--- fit line --- use robust method
#
    if linefit == 1:
        xx = []
        for m in range(0, len(xdata)):
            xx.append(xdata[m] - 1999)

        [sint, slope, sa, serr, sint2, slope2, sa2, serr2] =  line_fit_prep(xx, ydata, yerror)
        lint  =  '%2.3f' % (round(sint,  3))

        if slope < 0:
            sign = -1
        else:
            sign = 1

        lslope = '%2.3f' % (round(abs(slope), 3))
        lerr   = '%2.3f' % (round(serr,  3))
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
#--- set plotting range
#
    fig, ax = plt.subplots(1)
    ax.set_xbound(xmin,xmax)
    ax.set_xlim(xmin=xmin, xmax=xmax, auto=False)
    ax.set_ylim(ymin=ymin, ymax=ymax, auto=False)
#
#--- plot data
#
    ax.plot(xdata, ydata, color=colorList[pcolor], marker=marker, markersize=3, lw = connect)
#
#--- plot error bar
#
    if yerror != 0:
        plt.errorbar(xdata, ydata, yerr=yerror, lw = 0, elinewidth=1)
#
#--- plot fitted line
#
    if linefit == 1:
        if slope2 == 0:
            start = sint + slope * (xmin - 1999)
            stop  = sint + slope * (xmax - 1999)
            plt.plot([xmin, xmax], [start, stop], color=colorList[lcolor], lw = lsize)
        else:
            start = sint + slope * (xmin - 1999)
            stop  = sint + slope * (year_div - 1999)
            plt.plot([xmin, year_div], [start, stop], color=colorList[lcolor], lw = lsize)
#
            start = sint2 + slope2 * (year_div - 1999)
            stop  = sint2 + slope2 * (xmax - 1999)
            plt.plot([year_div, xmax], [start, stop], color=colorList[lcolor], lw = lsize)
#
#--- remove tickers
#
        line = ax.get_xticklabels()
        for label in line:
            label.set_visible(False)

        line = ax.get_yticklabels()
        for label in line:
            label.set_visible(False)
#
#--- set the size of the plotting area in inch
#
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(1.0, 0.5)
    plt.savefig('zout.png', format='png', dpi=200)

    plt.close('all')

    cmd = 'convert zout.png -trim out_plot.png'
    os.system(cmd)
    mcf.rm_file('zout.png')


#---------------------------------------------------------------------------------------------------
#-- set_min_max: set plotting range                                                              ---
#---------------------------------------------------------------------------------------------------

def set_min_max(xdata, ydata, xtime = 0, ybot = -999):

    """
    set plotting range
    Input:  xdata   ---- xdata
            ydata   ---- ydata
            xtime   ---- if it is >0, it set the plotting range from 1999 to the current in year
            ybot    ---- if it is == 0, the ymin will be 0, if the ymin computed is smaller than 0
    Output: [xmin, xmax, ymin, ymax]
    """

    xmin  = min(xdata)
    xmax  = max(xdata)
    xdiff = xmax - xmin
    xmin -= 0.1 * xdiff
    xmax += 0.2 * xdiff

    if xtime > 0:
        xmin  = 1999
        tlist = tcnv.currentTime()
        xmax  = tlist[0] + 1

    ymin  = min(ydata)
    ymax  = max(ydata)
    ydiff = ymax - ymin
    ymin -= 0.1 * ydiff
    ymax += 0.2 * ydiff

    if ybot == 0:
        if ymin < 0:
            ymin = 0

    return [xmin, xmax, ymin, ymax]

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

def line_fit_prep(x, y, e):

    cut = year_div - 1999
    if len(x) < 4:
        [a1, b1, siga1, sigb1] = line_fit(x, y, e)
        [a2, b2, siga2, sigb2] = [ 0, 0, 0, 0]
    else:
        x1  = []
        y1  = []
        e1  = []
        x2  = []
        y2  = []
        e2  = []
        for k in range(0, len(x)):
            if x[k] < cut:
                x1.append(x[k])
                y1.append(y[k])
                e1.append(e[k])
            else:
                x2.append(x[k])
                y2.append(y[k])
                e2.append(e[k])

        [a1, b1, siga1, sigb1] = line_fit(x1, y1, e1)
        [a2, b2, siga2, sigb2] = line_fit(x2, y2, e2)

    return [a1, b1, siga1, sigb1, a2, b2, siga2, sigb2]


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
#-- make_html_page: create trend html page with an interactive plot                               --
#---------------------------------------------------------------------------------------------------

def make_html_page(fig, outname, col, inst, sign, step, profile_page):
    """
    create trend html page with an interactive plot
    input:  fig     --- a plot code
            outname --- output name
            col     --- samp or pi
            inst    --- instrument s or i
            sign    --- negative or positive side (n/p)
            step    --- position at arm 
    output: <web_dir>/Plots/outname
    """
#
#--- if the step is empty, don't make any page
#
    if step == "":
        return False
#
#--- set step width
#
    if step != 'center':
        if inst == 'i':
            swidth = 20
            start  = step * swidth + 10
        else:
            swidth = 10
            start  = step * swidth
#
#--- header part
#
    out = '<!DOCTYPE html>\n'
    out = out + '<html>\n'
    out = out + '<head>\n'
    out = out + '<title> HRC Energy Trending: Radial Distribution</title>\n'

    out = out + '<style>\n'
    out = out + '\ta{color:#F5F5DC;}\n'

    out = out + '.fraction {\n'
    out = out + '    display: inline-block;\n'
    out = out + '    vertical-align: middle; \n'
    out = out + '    margin: 0 0.2em 0.4ex;\n'
    out = out + '    text-align: center;\n'
    out = out + '}\n'
    out = out + '.fraction > span {\n'
    out = out + '    display: block;\n'
    out = out + '    padding-top: 0.15em;\n'
    out = out + '}\n'
    out = out + '.fraction span.fdn {border-top: thin solid black;}\n'
    out = out + '.fraction span.bar {display: none;}\n'

    out = out + '</style>\n'

    out = out + '</head>\n'
    out = out + '<body style="background-color:#F5F5DC;width:95%;margin-left:10px; margin-right;10px">\n'
#
#--- creating the title
#
    if inst == 'i':
        device = ' HRC I: '
    else:
        device = ' HRC S: '

    if col  == 'samp_':
        title =  device + 'Scaled Sum Amp '
    else:
        title =  device + 'PI '

    if sign == 'center':
        title  = title + " Center"
        stitle = title

    else:
        if sign == 'p_':
            title = title + 'Positive ' 
        else:
            title = title + 'Negative ' 
    
        stop   = start  + swidth
        stitle = title  + 'Range: ' +  str(start) + ' - ' + str(stop) 
        title  = stitle + '&#8491;\n'

    out = out + '<h2 style="background-color:blue; color:#F5F5DC;">' + title + '</h2>\n'
    out = out + '<div style="text-align:right;">\n'
#
#--- create links to previous and next wavelength pages
#
    if step != 'center' and inst != 'i':
        paddress = html_top + 'Plots/' + col + sign + inst + '_' +  str(step -1) + '.html'
        pwave    = str(start-swidth)  + ' - ' + str(start) + '&#8491;'
    
        naddress = html_top + 'Plots/' + col + sign + inst + '_' +  str(step +1) + '.html'
        nwave    = str(stop)  + ' - ' + str(stop+swidth) + '&#8491;'

        out = out + '<span style="background-color:#006400; color:#F5F5DC;">'
        if step == 5:
            out = out + '<a href="' + naddress + '"><b>' + nwave + '&gt;&gt;</b></a>\n'
    
        elif step == 16 and sign == 'n_':
            out = out + '<a href="' + paddress + '"><b>&lt;&lt;' + pwave + '</b></a>\n'
    
        elif step == 17:
            out = out + '<a href="' + paddress + '"><b>&lt;&lt;' + pwave + '</b></a>\n'
    
        else:
            out = out + '<a href="' + paddress + '"><b>&lt;&lt;' + pwave + '</b></a>\n'
    
            out = out + '</span><span style="color:#F5F5DC">&nbsp;&nbsp;</span>\n' 
    
            out = out + '<span style="background-color:#006400; color:#F5F5DC;">'
            out = out + '<a href="' + naddress + '"><b>' + nwave + '&gt;&gt;</b></a>\n'
        out = out + '</span>'
    
        out = out + '<br/>&nbsp;&nbsp;<br />'
#
#--- back to the top page
#
    out = out + '<span style="background-color:#006400; color:#F5F5DC;">'
    out = out + '<a href=' + html_top + 'hz43.html'
    out = out + ' style="text-align:right;"><b>Back to Top Page</b></a>\n'
    out = out + '</span>\n'
    out = out + '</div>\n'


    out = out + '<p style="text-align:left"><b>\n'
    out = out + 'Hover the mouse over  the data point on the plot below to see the information about the data point.</b><br />\n'
    out = out + '(Note: it may take a while to load the interactive plot.)  </p>\n'
    out = out + '<p style="text-align:left">\n'
    out = out + 'The fitted line on the popup window is a Gamma distribtuion:'

    out = out + '<div style="padding-top:5px; padding-bottom:5px;padding-left:80px;">'

    out = out + 'g(x; k, &theta;) = '
    out = out + '<div class="fraction">'
    out = out + '<span class="fup"><i>&theta;<sup>k</sup></i></span> <span class="bar"></span>'
    out = out + '<span class="fdn"><i>&Gamma;(k)</i></span>'
    out = out + '</div>'
    out = out + 'x<sup>k-1</sup> exp(-&theta;* x)'
    out = out + '</div>'

    out = out + 'where <em>k</em> is the shape parameter and <em>&theta;</em> is the scale parameter.'

    out = out + '</p>'
    #out = out + '<p style="text-align:left">\n'
    #out = out + 'The intercept of the slope fitted below starts at Year 1999 not Year 0.</p>\n'
#
#--- convert the figure code into html 
#
    out = out + '<h3>' + stitle + '</h3>\n'
    out = out + '<div style="padding-left:80px;">\n'
    out = out + mpld3.fig_to_html(fig)
    out = out + '</div>'
#
#--- open a profile page
#
    hprofile = profile_page.replace("/proj/web-cxc-dmz/htdocs", "http://cxc.cfa.harvard.edu")
    out = out + '<h3><a href="' + hprofile + '" style=\'color:blue;\'>Open Fitted Distribution Page</a></h3>\n\n'

#
#--- count plots
#
    out = out + '<h3>Count Rate: ' + stitle + '</h3>'

    #out = out + '<p style="text-align:left"> Note: The count errors are generally larger than the plotting range '
    #out = out + ' and they are not shown on the plots.</p>\n'

    out = out + '<div style="margin-left:70px;">\n'
    out = out +'<img src="' + html_top + 'Plots/Count_rates/'
    if sign == 'center':
        out = out + col + 'center_'+ inst + '.png">\n'
    else:
        out = out + col + sign + inst + '_' + str(step) + '.png">\n'
    out = out + "</div>\n"
#
#--- back to the top page
#
    out = out + '<div style="text-align:right;padding-bottom:40px;">\n'
    out = out + '<span style="background-color:#006400; color:#F5F5DC;">'
    out = out + '<a href=' + html_top + 'hz43.html'
    out = out + ' style="text-align:right;"><b>Back to Top Page</b></a>\n'
    out = out + '</span>\n'
    out = out + '</div>\n'
#
#--- the tail part
#
    out = out + '<hr />'
    out = out + '<p style="text-align:left;padding-top:10px; padding-bottom:20px"> \n'
    out = out + '<i>If you have any questions about this page, please contact \n'
    out = out + '<a href="mailto:tisobe@cfa.harvard.edu" style="background-color:#006400;">'
    out = out + 'tisobe@cfa.harvard.edu</a>.\n'
    out = out + '</i></p>\n'

    out = out + '</body>\n'
    out = out + '</html>\n'

    fo  = open(outname, 'w')
    fo.write(out)
    fo.close()

#---------------------------------------------------------------------------------------------------
#-- make_profile_page: create a distribution profile plot page                                   ---
#---------------------------------------------------------------------------------------------------

def make_profile_page(time_list, info_list, title, outname, backpage):
    """
    create a distribution profile plot page
    imput:  time_list   --- a list of time in year
            info_list   --- a list of information line 
            title       --- a title of the ouput
            outname     --- output file name
            backpage    --- a link to the page this page is linked from
    output: outpname
    """
#
#--- if the step is empty, don't make any page
#
    if step == "":
        return False
#
#--- header part
#
    out = '<!DOCTYPE html>\n'
    out = out + '<html>\n'
    out = out + '<head>\n'
    out = out + '<title>  HZ43 HRC Energy Trending: Radial Distribution</title>\n'

    out = out + '<style>\n'
    out = out + '\ta{color:#F5F5DC;}\n'

    out = out + '.fraction {\n'
    out = out + '    display: inline-block;\n'
    out = out + '    vertical-align: middle; \n'
    out = out + '    margin: 0 0.2em 0.4ex;\n'
    out = out + '    text-align: center;\n'
    out = out + '}\n'
    out = out + '.fraction > span {\n'
    out = out + '    display: block;\n'
    out = out + '    padding-top: 0.15em;\n'
    out = out + '}\n'
    out = out + '.fraction span.fdn {border-top: thin solid black;}\n'
    out = out + '.fraction span.bar {display: none;}\n'

    out = out + '</style>\n'

    out = out + '<script language="JavaScript">   \n'
    out = out + '    function WindowOpener(imgname) {   \n'
    out = out + '    msgWindow = open("","displayname","toolbar=no,directories=no,menubar=no,location=no,scrollbars=no,status=no,,width=1000,height=500,resize=no");   \n'
    out = out + '    msgWindow.document.clear();   \n'
    out = out + '    msgWindow.document.write("<html><title>Trend plot:   \'+imgname+\'</TITLE>");   \n'
    out = out + '    msgWindow.document.write("<body style=\'background-color:#F5F5DC;\'>");   \n'
    out = out + '    msgWindow.document.write("<img src=\'../Indivisual_Plots/"+imgname+"\' border =0 ><p></p></body></html>")   \n'
    out = out + '    msgWindow.document.close();   \n'
    out = out + '    msgWindow.focus();   \n'
    out = out + '    }   \n'
    out = out + '</script>   \n'

    out = out + '</head>\n'
    out = out + '<body style="background-color:#F5F5DC;width:95%;margin-left:10px; margin-right;10px">\n'

    out = out + '<h2>Distribution Page: ' + title + '</h2>\n'

    hprofile = backpage.replace("/proj/web-cxc-dmz/htdocs", "http://cxc.cfa.harvard.edu")

    out = out + '<div style="text-align:right;padding-bottom:5px;">\n'
    out = out + '<h3><a href="' + hprofile + '" style=\'color:blue\'>Back to Previous Page</a></h3>\n'
    out = out + '<h3><a href="' + html_top + 'hz43.html"  style=\'color:blue\'>Back to Top Page</a></h3>\n'

    out = out + '</div>\n'

    out = out + '<p>Click a figure to enlarge</p>\n'
    
    out = out + '<table border=1 cellpadding=2>\n'

    ddict = {}
    temp  = []
    for k in range(0, len(time_list)):
        ctime  = float(round(time_list[k], 2))
        ddict[ctime] = info_list[k]
        temp.append(ctime)
    time_sorted = sorted(temp)

    for k in range(0, len(time_sorted)):

        if k % 5 == 0:
            if k != 0:
                out = out + '</tr>\n'
            out = out + '<tr>\n'

        out = out + '<td>\n'
        out = out + '<b>Year: ' + str(time_sorted[k]) + '</b><br />' 
        xtemp = re.split('\/', ddict[time_sorted[k]])
        pname = xtemp[-2] + '/' + xtemp[-1]
        out = out + '<a href="javascript:WindowOpener(\'' + pname + '\')"><img src="' + ddict[time_sorted[k]] + '" style="width:300px;"></a></td>\n'

    
    out = out + '</tr>\n'
    out = out + '</table>\n'

    out = out + '<div style="text-align:right;padding-bottom:20px;">\n'
    out = out + '<h3><a href="' + hprofile + '" style=\'color:blue\'>Back to Previous Page</a></h3>\n'
    out = out + '<h3><a href="' + html_top + 'hz43.html"  style=\'color:blue\'>Back to Top Page</a></h3>\n'
    out = out + '</div>\n'
    out = out + '</body>\n</html>\n'

    fo  = open(outname, 'w')
    fo.write(out)
    fo.close()

#---------------------------------------------------------------------------------------------------

if __name__ == "__main__":

    create_plots()

#!/usr/bin/env /proj/sot/ska/bin/python

#########################################################################################################
#                                                                                                       #
#       create_arlac_trend_plots.py: creates binned trend plots                                         #
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

mday     = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
year_div = 2012
    
#---------------------------------------------------------------------------------------------------
#-- create_plots: creates binned trend plots                                                      --
#---------------------------------------------------------------------------------------------------

def create_plots(chk=''):
    """
    creates binned trend plots 
    input: none but read data from <data_dir> 
    output: plots in <web_dir>/Plots3/   such as samp_n_s_14.png or pi_n_s_14.png
    """
#
#--- check whether new data came in
#
    ifile = house_keeping + 'hrci_new'
    sfile = house_keeping + 'hrcs_new'
    if chk == '':
        if (not os.path.isfile(ifile)) and (not os.path.isfile(sfile)):
            exit(1)


    xmin  = 1999
    xmax  = time.strftime('%Y', time.localtime())
    xmax  = int(float(xmax)) + 1
    xname = 'Time (year)'

    for col in ['pi', 'samp']:
        if col == 'pi':
            yname = 'PI'
            ymin  = 40
            ymax  = 190
        else:
            yname = 'Scaled SumAmp'
            ymin  = 0
            ymax  = 150

        for inst in ['s', 'i']:

            infile  = data_dir + col + '_list_' + inst
            [time_list, data_list, err_list, cnt_list, bkg_list, obsid_list] = read_file(infile, inst)

            oname = data_dir + col + '_' + inst + '_fitting_results'
            fo    = open(oname, 'w')
            for pos in range(0, 21):
                outname      = web_dir + 'Plots3/'  + col + '_' + inst + '_' + str(pos) +  '.html'
                profile_page = web_dir + 'Plots3/Profile_page/' + col + '_' + inst + '_' + str(pos) + '.html'

                [info_list, hlink_list]  = create_info_link(obsid_list[pos], col)

                [fig, lint, lslope, lerr, lint2, lslope2, lerr2] \
                            = plot_single_panel(xmin, xmax, ymin, ymax, time_list[pos], data_list[pos], err_list[pos], \
                                                info_list, xname, yname, label='')
    
                make_html_page(fig, outname, col, inst, pos, profile_page)
    
                create_thumb_nail_plots(xmin, xmax, ymin, ymax, time_list[pos], data_list[pos], err_list[pos], \
                                            info_list, xname, yname, label='')

                p_name = web_dir + 'Plots3/Thumb_plots/' + col + '_' +  inst + '_' + str(pos) +  '_thumb_plot.png'
                cmd    = 'mv out_plot.png ' + p_name
                os.system(cmd)
    
                oline = lint  + '\t' + lslope + '\t' + lerr + '\t'
                oline = oline + '\t' + lint2  + '\t' + lslope2 + '\t' + lerr2 + '\n'
                fo.write(oline)


                if inst == 'i':
                    osets =  hrci_offset
                    title = 'HRC I ' +  yname + ': Location ' + str(osets[pos])
                else:
                    osets =  hrcs_offset
                    title = 'HRC S ' +  yname + ': Location ' + str(osets[pos])
    
                make_profile_page(time_list[pos], hlink_list, title, profile_page, outname)

            fo.close()

#---------------------------------------------------------------------------------------------------
#-- read_file: read arlac data file                                                               --
#---------------------------------------------------------------------------------------------------

def read_file(infile, inst):
    """
    read arlac data file
    input:  infile      --- data file name
    output: time_list   --- a list of time data
            data_list   --- a list of value data either sumamp or pi
            err_list    --- a list of the error of the value
            obsid_list  --- a list of obsid of the observations
    """
    if inst.lower() == 'i':
        offset_list = hrci_offset
    else:
        offset_list = hrcs_offset

    data  = read_data(infile)
    if data == []:
        return False

    time_list  = []
    data_list  = []
    err_list   = []
    cnt_list   = []
    bkg_list   = []
    obsid_list = []
    for k in range(0, 21):
        time_list.append([])
        data_list.append([])
        err_list.append([])
        cnt_list.append([])
        bkg_list.append([])
        obsid_list.append([])

    for ent in data:
        atemp = re.split('\s+', ent)
        time  = get_frac_year(atemp[2])
        if time < 2000:
            continue
#
#--- find where to put the data; 0 to 20 data slots
#
        yoffset = int(round(float(atemp[10])))
        zoffset = int(round(float(atemp[11])))
        pos = -1
        for k in range(0, 21):
            if (offset_list[k][0] == yoffset) and  (offset_list[k][1] == zoffset):
                pos = k
                break
        if pos < 0:
            continue

        avg   = float(atemp[12])
        err   = float(atemp[13])
        crate = float(atemp[15])
        brate = float(atemp[16])
        if avg == -999 or err == -999:
            continue
        else:
            time_list[pos].append(time)
            data_list[pos].append(avg)
            err_list[pos].append(err)
            cnt_list[pos].append(crate)
            bkg_list[pos].append(brate)
            obsid_list[pos].append(atemp[1])

    return [time_list, data_list, err_list, cnt_list, bkg_list, obsid_list]

#---------------------------------------------------------------------------------------------------
#-- create_info_link: create a list of links used by an interactive plot                          --
#---------------------------------------------------------------------------------------------------

def create_info_link(obsid_list, col):
    """
    create a list of links used by an interactive plot
    input:  obsid_list  --- a list of obsids
            col         --- pi or samp
    output: outlist     --- a list of html address to a distribution plot
    """

    html_plot = html_top + '/Plots3/Indivisual_Plots/'
    outlist  = []
    hlist    = []
    for obsid in obsid_list:
        ofile = html_plot + str(obsid)   + '_' + col +  '_list_vfit.png' 
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

def make_html_page(fig, outname, col, inst, step, profile_page):
    """
    create trend html page with an interactive plot
    input:  fig     --- a plot code
            outname --- output name
            col     --- samp or pi
            inst    --- instrument s or i
            sign    --- negative or positive side (n/p)
            step    --- position at arm 
    output: <web_dir>/Plots3/outname
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
    out = out + '<title> AR LacHRC Energy Trending: Radial Distribution</title>\n'

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
        olist  = hrci_offset
    else:
        device = ' HRC S: '
        olist  = hrcs_offset
    
    if col  == 'samp':
        title =  device + 'Scaled Sum Amp '
    else:
        title =  device + 'PI '

    title = title + 'Offset Position: ' +  str(olist[step]) + '\n'

    out = out + '<h2 style="background-color:blue; color:#F5F5DC;">' + title + '</h2>\n'
    out = out + '<div style="text-align:right;">\n'
#
#--- create links to previous and next wavelength pages
#
    paddress = html_top + 'Plots3/' + col + '_' + inst + '_' +  str(step -1) + '.html'
    pwave    = 'Location: ' + str(olist[step-1])

    naddress = html_top + 'Plots3/' + col + '_' + inst + '_' +  str(step +1) + '.html'
    try:
        nwave    = 'Location: ' + str(olist[step+1])
    except:
        nwave    = ''

    out = out + '<span style="background-color:#006400; color:#F5F5DC;">'
    if step == 0:
        out = out + '<a href="' + naddress + '"><b>' + nwave + '&gt;&gt;</b></a>\n'

#    elif step == 16 and sign == 'n_':
#        out = out + '<a href="' + paddress + '"><b>&lt;&lt;' + pwave + '</b></a>\n'

    elif step == 20:
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
    out = out + '<a href=' + html_top + 'arlac_trends.html'
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

    out = out + 'where <em>k</em> is the shape parameter and <em>&theta;</em> is the scale parameter. '
    out = out + 'Shape parameter is defined by: (mean / sigma)<sup>2</sup>, and '
    out = out + 'scale parameter id defined by: mean / sigma<sup>2</sup>. </p>\n'

    out = out + '</p>'
#
#--- convert the figure code into html 
#
    out = out + '<h3>' + title + '</h3>'
    out = out + '<div style="padding-left:80px;">\n'
    out = out + mpld3.fig_to_html(fig)
    out = out + '</div>\n\n'
#
#--- open a profile page
#
    hprofile = profile_page.replace("/proj/web-cxc-dmz/htdocs", "http://cxc.cfa.harvard.edu")
    out = out + '<h3><a href="' + hprofile + '" style=\'color:blue;\'>Open Fitted Distribution Page</a></h3>\n\n'
#
#--- count plots
#
    out = out + '<h3>' + title + ':  Count Rate</h3>'

    out = out + '<div style="margin-left:70px;">\n'
    out = out +'<img src="' + html_top + 'Plots3/Count_rates/'
    if sign == 'center':
        out = out + col + 'center_'+ inst + '.png">\n'
    else:
        out = out + col + '_' + inst + '_' + str(step) + '.png">\n'
    out = out + "</div>\n"
#
#--- back to the top page
#
    out = out + '<div style="text-align:right;padding-bottom:40px;">\n'
    out = out + '<span style="background-color:#006400; color:#F5F5DC;">'
    out = out + '<a href=' + html_top + 'arlac_trends.html'
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
    out = out + '<title> AR LacHRC Energy Trending: Radial Distribution</title>\n'

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
    out = out + '<h3><a href="' + html_top + 'arlac_trends.html"  style=\'color:blue\'>Back to Top Page</a></h3>\n'

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
        pname = xtemp[-1]
        out = out + '<a href="javascript:WindowOpener(\'' + pname + '\')"><img src="' + ddict[time_sorted[k]] + '" style="width:300px;"></a></td>\n'

    
    out = out + '</tr>\n'
    out = out + '</table>\n'

    out = out + '<div style="text-align:right;padding-bottom:20px;">\n'
    out = out + '<h3><a href="' + hprofile + '" style=\'color:blue\'>Back to Previous Page</a></h3>\n'
    out = out + '<h3><a href="' + html_top + 'arlac_trends.html"  style=\'color:blue\'>Back to Top Page</a></h3>\n'
    out = out + '</div>\n'
    out = out + '</body>\n</html>\n'

    fo  = open(outname, 'w')
    fo.write(out)
    fo.close()

#---------------------------------------------------------------------------------------------------

if __name__ == "__main__":

    if len(sys.argv) == 2:
        chk = sys.argv[1]
    else:
        chk = ''

    create_plots(chk)

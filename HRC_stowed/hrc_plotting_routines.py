#!/usr/bin/env /proj/sot/ska/bin/python

#############################################################################################
#                                                                                           #
#       hrc_plotting_routines.py: create hrc plots                                          #
#                                                                                           #
#               author: t. isobe (tisobe@cfa.harvard.edu)                                   #
#                                                                                           #
#               last update: Apr 14, 2016                                                   #
#                                                                                           #
#############################################################################################

import os
import sys
import re
import string
import time
import math
import numpy
import astropy.io.fits  as pyfits

import matplotlib as mpl

if __name__ == '__main__':
    mpl.use('Agg')

from pylab import *
import matplotlib.pyplot       as plt
import matplotlib.font_manager as font_manager
import matplotlib.lines        as lines
#
#--- reading directory list
#
path = '/data/aschrc6/wilton/isobe/Project1/Scripts/house_keeping/dir_list'

f= open(path, 'r')
data = [line.strip() for line in f.readlines()]
f.close()

for ent in data:
    atemp = re.split(':', ent)
    var   = atemp[1].strip()
    line  = atemp[0].strip()
    exec "%s = %s" %(var, line)

sys.path.append(mta_dir)
sys.path.append(bin_dir)
#
#--- import several functions
#
import convertTimeFormat          as tcnv       #---- contains MTA time conversion routines
import mta_common_functions       as mcf        #---- contains other functions commonly used in MTA scripts
import robust_linear              as robust     #---- robust linear fit routine
from kapteyn import kmpfit                      #---- least sq fitting routine 
import hrc_stowed_common_function as hcf
#
#--- temp writing file name
#
import random
rtail   = int(10000 * random.random())       #---- put a romdom # tail so that it won't mix up with other scripts space
rtail   = int(time.time())
zspace  = '/tmp/zspace' + str(rtail) + 'hrc'
#
#--- set a few other things for plottings
#
color = ('blue', 'red', 'green', 'aqua', 'fuchsia','lime', 'maroon', 'black', 'olive', 'yellow')
markerList = ('o',    '*',     '+',   '^',    's',    'D',       '1',      '2',     '3',      '4')
labels     = ['s2hvst=5', 's2hvst=6', 's2hvst=7', 's2hvst=8']

fsize= 9

#-----------------------------------------------------------------------------------------------
#-- plot_hrc_trend_data: create all trending plots                                           ---
#-----------------------------------------------------------------------------------------------

def plot_hrc_trend_data():
    """
    create all trending plots
    input:  none but read from <data_dir>
    output: hrc_i_115_cnt_trend.png    
            hrc_i_115_ephin_trend.png  
            hrc_i_115_pha_trend.png        
            hrc_i_115_shield_valid.png  
            hrc_i_115_ephin_total.png  
            hrc_i_115_ephin_valid.png  
            hrc_i_115_s2hvst_vt_ratio.png  etc
    """
    for sfile in ['hrc_i_115_stat_results',      'hrc_s_125_hi_stat_results',   \
                  'hrc_s_125_sec1_stat_results', 'hrc_s_125_sec2_stat_results', \
                  'hrc_s_125_sec3_stat_results']:

        infile = data_dir + 'Stats/' + sfile
        mc     = re.search('hrc_i', sfile)
        mc2    = re.search('hrc_s_125_hi', sfile)

        if mc is not None:
            inst = 1
        elif mc2 is not None:
            inst = 2
        else:
            inst = 3
#
#--- read data
#
        indata = read_stat_data(infile)
#
#--- plot trend data
#
        plot_count_trend(sfile,      indata)
        plot_ephin_trend(sfile,      indata, inst)
        plot_pha_trend(sfile,        indata, inst)
        plot_ephin_comparison(sfile, indata, inst)
        plot_s2hvst_vt_ratio(sfile,  indata)

        plot_dead_time_trend()

#-----------------------------------------------------------------------------------------------
#-- plot_count_trend: plot valid, total, and shield rate time trend data                      --
#-----------------------------------------------------------------------------------------------

def plot_count_trend(sfile, indata):
    """
    plot valid, total, and shield rate time trend data
    input:  sfile   --- a data file name
            indata  --- a list of data lists
    output: <head_name>_cnt_trend.png, e.g. hrc_i_115_cnt_trend.png
    """
#
#--- setting up the plot surface
#
    prep_plot()
#
#--- set plotting range
#
    xmin = 1999
    xmax = int(max(indata[0])) + 1
    ymin = 0
#
#--- valid rate trend
#
    a1     = plt.subplot(311)
    valid  = separate_data_s2hvst(indata[0], indata[12], indata[14], indata[21]) 
    ymax   = max(indata[12]) * 1.2
    title  = 'Mean of Valid Count Rate'
    plot_panel(valid, xmin, xmax, ymin, ymax, labels, title)
#
#--- total rate trend
#
    a2     = plt.subplot(312)
    total  = separate_data_s2hvst(indata[0], indata[9],  indata[11], indata[21]) 
    ymax   = max(indata[9])  * 1.2
    title  = 'Mean of Total Count Rate'
    plot_panel(total, xmin, xmax, ymin, ymax, labels, title)
#
#--- sheild rate trend
#
    a3     = plt.subplot(313)
    shield = separate_data_s2hvst(indata[0], indata[15], indata[17], indata[21]) 
    ymax   = 10
    title  = 'Mean of Shield Count Rate (x 10^-3)'
    plot_panel(shield, xmin, xmax, ymin, ymax, labels, title)
#
#--- add x ticks label only on the panels of the last row
#
    plt.setp(a1.get_xticklabels(), visible=False)
    plt.setp(a2.get_xticklabels(), visible=False)
#
#--- add the legend and x axis name
#
    add_legend(a1)

    a1.set_ylabel("Counts per Sec")
    a2.set_ylabel("Counts per Sec")
    a3.set_ylabel("Counts per Sec")
    a3.set_xlabel("Time (Year)", size=fsize)
#
#--- save the plot in file
#
    save_plot_in_file(sfile, 'cnt_trend.png')

#-----------------------------------------------------------------------------------------------
#-- add_legend: add legend depending on ncol                                                 ---
#-----------------------------------------------------------------------------------------------

def add_legend(a1):
    """
    add legend depending on ncol
    """

    try:
        a1.legend(loc='upper left', bbox_to_anchor=(0.0, 1.20), fancybox=False, shadow=False,ncol=4)
    except:
        try:
            a1.legend(loc='upper left', bbox_to_anchor=(0.0, 1.20), fancybox=False, shadow=False,ncol=3)
        except:
            try:
                a1.legend(loc='upper left', bbox_to_anchor=(0.0, 1.20), fancybox=False, shadow=False,ncol=2)
            except:
                pass

#-----------------------------------------------------------------------------------------------
#-- prep_plot: setup the matplot plotting panel                                               --
#-----------------------------------------------------------------------------------------------

def prep_plot():
    """
    setup the matplot plotting panel
    input:  none
    output: none
    """
#
#--- clean up the plotting device
#
    plt.close('all')
#
#---- set a few parameters
#
    mpl.rcParams['font.size'] = fsize
    props = font_manager.FontProperties(size=10)
    plt.subplots_adjust(hspace=0.04)

#-----------------------------------------------------------------------------------------------
#-- save_plot_in_file: save the plot in file                                                  --
#-----------------------------------------------------------------------------------------------

def save_plot_in_file(sfile, tail):
    """
    save the plot in file
    input:  sfile   --- input data name
            tail    --- tail to finish the plot file
    output: <head_name>_<tail>, e.g. hrc_i_115_cnt_trend.png
    """
    outname = sfile.replace('stat_results', tail)
    outname = plot_dir + outname
    create_plot_file(outname)

#-----------------------------------------------------------------------------------------------
#-- plot_ephin_trend: plot ephin data and valid/total count rate ratios                       --
#-----------------------------------------------------------------------------------------------

def plot_ephin_trend(sfile, indata, inst):
    """
    plot ephin data and valid/total count rate ratios
    input:  sfile   --- a name of the data file
            indata  --- a list of data lists
    output: a plot file (e.g. hrc_i_115_ephin_trend.png)
    """
#
#--- setting up a plotting surface
#
    prep_plot()
#
#--- set plotting range
#
    xmin = 1999
    xmax = int(max(indata[0])) + 1
    ymin = 0
#
#--- ephin integral trend ----------------------
#
    a1     = plt.subplot(211)
    ephin  = separate_data_s2hvst(indata[0], indata[23], indata[24], indata[21]) 
    ymax   = max(indata[12]) * 1.2
    title  = ''

    plot_panel(ephin, xmin, xmax, ymin, ymax, labels, title)
#
#--- special treatment for semi-log plot
#
    plt.ylim(0.001,10)
    plt.yscale('log')
    title  = 'EPHIN Integral Channel Flux'
    tfsize = fsize + 1
    plt.text(2000, 5, title, size=tfsize, weight='bold')
#
#--- valid / total ratio trend ----------------
#
#
#--- compute ratio: valid/total
#
    [ratio, std] = compute_ratio(indata[12], indata[9], indata[14], indata[11])
#
#--- plotting the data
#
    a2     = plt.subplot(212)
    rate   = separate_data_s2hvst(indata[0], ratio, std, indata[21]) 
    if inst == 1:
        ymin = 0.1
        ymax = 0.7
        deg  = 2
    elif inst == 2:
        ymin = 0.7
        ymax = 1.3
        deg  = 1
    elif inst == 3:
        ymin = 0.2
        ymax = 0.8
        deg  = 1

    title  = 'Mean of Valid Counts/Total Count'

    plot_panel(rate, xmin, xmax, ymin, ymax, labels, title)
#
#--- polynomial fit on the data with degree of 2
#
    [coeff, xsave, ysave] = poly_fit(indata[0], ratio, xmin, xmax, deg=deg,  shift=-2000)
    plt.plot(xsave, ysave, color='lime', lw =2)
#
#--- position the display equation of the polynomial fit
#
    xpos = xmin + 0.05 * (xmax - xmin)
    ypos = ymax - 0.2 * (ymax-ymin)
    if deg == 2:
        line = set_poly_eqation_disp(coeff, shift=-2000)
    else:
        line = 'Slope: ' + str(round(coeff[1],3))

    plt.text(xpos, ypos , line, size=fsize, weight='bold')
#
#--- add x ticks label only on the panels of the last row
#
    plt.setp(a1.get_xticklabels(), visible=False)
#
#--- add the legend and x axis name
#
    #a1.legend(loc='upper left', bbox_to_anchor=(0.0, 1.20), fancybox=False, shadow=False,ncol=5)
    a1.set_ylabel("Counts per Sec (in Log)")
    a2.set_ylabel("Vaild / Total")
    a2.set_xlabel("Time (Year)", size=fsize)
#
#--- save the plot in file
#
    save_plot_in_file(sfile, 'ephin_trend.png')


#-----------------------------------------------------------------------------------------------
#-- plot_pha_trend: plot median peak postion and peak width trends                           ---
#-----------------------------------------------------------------------------------------------

def plot_pha_trend(sfile, indata, inst=1):
    """
    plot median peak postion and peak width trends
    input:  sfile   --- a data file name
            indata  --- a list of data lists
    output: <head_name>_pha_trend.png, e.g.: hrc_i_115_pha_trend.png
    """
#
#--- setting a plotting surface
#
    prep_plot()
#
#--- set plotting range
#
    xmin = 1999
    xmax = int(max(indata[0])) + 1
#
#--- since no error bar, create 0 error bar list
#
    err_bar = []
    for i in range(0, len(indata[0])):
        err_bar.append(0)
#
#--- median pha peak position trend ---------------------------------
#
    a1 = plt.subplot(211)
    peak  = separate_data_s2hvst(indata[0], indata[7], err_bar, indata[21]) 
    if inst == 1:
        ymin   = 40
        ymax   = 130
        bpoint = 2010
    elif inst == 2:
        ymin   = 100
        ymax   = 220
        bpoint = 2012
    else:
        ymin   = 80
        ymax   = 190
        bpoint = 2012


    title  = 'Median PHA Peak Position'
    plot_panel(peak, xmin, xmax, ymin, ymax, labels, title,psize=4, ebar=0, tpos = 2)
#
#--- two section line fits
#
    bpoints = [bpoint]
    [int_list, slope_list, serr_list, xa_list, ya_list]  \
                    = multi_section_linear_fit(indata[0], indata[7], xmin, xmax, bpoints, shift=-2000)
    for k in range(0, 2):
        plt.plot(xa_list[k], ya_list[k], color='lime', lw =2)
    
    add_slope_description(xmin, xmax, ymin, ymax,  slope_list, str(bpoint))
#
#--- peak width trend ---------------------------------
#
    a2 = plt.subplot(212)
    width  = separate_data_s2hvst(indata[0], indata[8],  err_bar, indata[21]) 
    ymin   = 60
    ymax   = 90
    title  = 'Width of PHA Peak'
    plot_panel(width, xmin, xmax, ymin, ymax, labels, title, psize=4, ebar=0, tpos = 2)
#
#--- two section line fits
#
    bpoints = [bpoint]
    [int_list, slope_list, serr_list, xa_list, ya_list]  \
                    = multi_section_linear_fit(indata[0], indata[8], xmin, xmax, bpoints, shift=-2000)
    for k in range(0, 2):
        plt.plot(xa_list[k], ya_list[k], color='lime', lw =2)

    add_slope_description(xmin, xmax, ymin, ymax,  slope_list, str(bpoint))
#
#-----------------------------------------
# 

#
#--- add x ticks label only on the panels of the last row
#
    plt.setp(a1.get_xticklabels(), visible=False)
#
#--- add the legend and x axis name
#
    a1.legend(loc='upper left', bbox_to_anchor=(0.0, 1.20), fancybox=False, shadow=False,ncol=5)
    a1.set_ylabel('PHA')
    a2.set_ylabel('PHA')
    a2.set_xlabel("Time (Year)", size=fsize)
#
#--- save the plot in file
#
    save_plot_in_file(sfile, 'pha_trend.png')

#-----------------------------------------------------------------------------------------------
#-- add_slope_description: add two slope values on the plot                                   --
#-----------------------------------------------------------------------------------------------

def add_slope_description(xmin, xmax, ymin, ymax, slope_list, breakp, tfsize=9):
    """
    add two slope values on the plot
    input:  xmin        --- min of x
            xmax        --- max of x
            ymin        --- min of y
            ymax        --- max of y
            slope_list  --- a list of slope vlaues (2 entries)
            breakp      --- an x value of breaking point
            tfsize      --- a font size
    output: add slope values to a plot
    """

    xdiff = xmax - xmin
    ydiff = ymax - ymin
    xpos  = xmin + 0.1 * xdiff

    ypos  = ymin + 0.4 * ydiff
    title = 'Slope: %3.3f' % slope_list[0]
    title = title + ' (before ' + str(breakp) + ')' 
    plt.text(xpos, ypos, title, size=tfsize, weight='bold')

    ypos  = ymin + 0.3 * ydiff
    title = 'Slope: %3.3f ' % slope_list[1]
    title = title + ' (after ' + str(breakp) + ')' 
    plt.text(xpos, ypos, title, size=tfsize, weight='bold')

#-----------------------------------------------------------------------------------------------
#-- plot_ephin_comparison: plot count rates against ephin count and shield rate               --
#-----------------------------------------------------------------------------------------------

def plot_ephin_comparison(sfile, indata, inst=1):
    """
    plot count rates against ephin count and shield rate
    input:  sfile   --- a data file name
            indata  --- a list of data lists
    output: <head_name>_ephin_valid.png, <head_name>_shield_valid.png
    """
#
#--- setting a plotting surface
#
    prep_plot()
#
#---- against ephin rate --------------------------------------
#
    xmin   = 0.005
    xmax   = 6
#
#--- indata[12] --- valid data
#--- indata[9]  --- total data
#
    [ymin1, ymax1, ymin2, ymax2] = find_plotting_range(indata[12], indata[9])

    psize = 3
    tfsize = fsize -1

    labels= ['Ephin < 1.0', 'Ephin >1.0']
    a1 = plt.subplot(211)
    valid  = separate_data_time(indata[23], indata[12], indata[14], indata[0]) 
    title  = 'Mean Valid Count Rate'

    a1.set_xscale("log", nonposx='clip')
    plot_panel(valid, xmin, xmax, ymin1, ymax1, labels, title, step=2, txpos=1)

    a2 = plt.subplot(212)
    valid  = separate_data_time(indata[23], indata[9], indata[11], indata[0]) 
    title  = 'Mean Total Count Rate'

    a2.set_xscale("log", nonposx='clip')
    plot_panel(valid, xmin, xmax, ymin2, ymax2, labels, title, step=2, txpos=1)
#
#--- add the legend and x axis name
#
    a2.set_xlabel("EPIN Integral Channel Flux", size=fsize)

    a1.set_ylabel("Valid Count Rate", size=fsize)
    a2.set_ylabel("Total Count Rate", size=fsize)
#
#--- add x ticks label only on the panels of the last row
#
    plt.setp(a1.get_xticklabels(), visible=False)
#
#--- save the plot in file
#
    save_plot_in_file(sfile,  'ephin_valid.png')
#
#---- against shield rate --------------------
#
#
#--- setting a plotting surface
#
    prep_plot()

    xmin   = 2
    xmax   = 8
    psize  = 3
    tfsize = fsize -1

    a1     = plt.subplot(211)
    valid  = separate_data_ephin(indata[15], indata[12], indata[14], indata[23]) 
    title  = 'Mean Valid Count Rate'
    plot_panel(valid, xmin, xmax, ymin1, ymax1, labels, title, step=2)

    a2     = plt.subplot(212)
    total  = separate_data_ephin(indata[15], indata[9], indata[11], indata[23]) 
    title  = 'Mean Total Count Rate'
    labels= ['Ephin < 1.0', 'Ephin >1.0']
    plot_panel(total, xmin, xmax, ymin2, ymax2, labels, title, step=2)
# 
#
#--- add the legend and x axis name
#
    a2.set_xlabel("Mean Shield Count Rate (x10^-3)", size=fsize)

    a1.set_ylabel("Valid Count Rate", size=fsize)
    a2.set_ylabel("Total Count Rate", size=fsize)

    add_legend(a1)
#
#--- add x ticks label only on the panels of the last row
#
    plt.setp(a1.get_xticklabels(), visible=False)
#
#--- save the plot in file
#
    save_plot_in_file(sfile,  'shield_valid.png')

#-----------------------------------------------------------------------------------------------
#-- find_plotting_range: set min and max range on two different data sets with the same size interval
#-----------------------------------------------------------------------------------------------

def find_plotting_range(d1, d2, step=10):
    """
    set min and max range on two different data sets with the same size interval
    input:  d1  --- data 1
            d2  --- data 2
            step    --- min inteval size. if step = 10, min and max are 10 * i i = 0, 1, 2, ...
    ouput:  min1, max1, min2, max2
    """

    min1  = min(d1)
    max1  = max(d1)
    diff1 = max1 - min1

    min2  = min(d2)
    max2  = max(d2)
    diff2 = max2 - min2

    if diff1 >= diff2:

        min1 -= 0.15 * diff1
        val = int(min1 / step)
        min1  = step * (val-1)
        if min1 < 0:
            min1 = 0

        max1 += 0.15 * diff1
        val = int(max1 / step)
        max1  = step * (val+1)

        diff  = max1 - min1
        da    = numpy.array(d2)
        dm    = step * int(da.mean()/step)
        min2  = dm - 0.5 * diff
        max2  = dm + 0.5 * diff
        if min2 < 0:
            min2 = 0
            max2 = diff
    else:

        min2 -= 0.15 * diff2
        val = int(min2 / step)
        min2  = step * (val-1)
        if min2 < 0:
            min2 = 0

        max2 += 0.15 * diff2
        val = int(max2 / step)
        max2  = step * (val+1)

        diff  = max2 - min2
        da    = numpy.array(d1)
        dm    = step * int(da.mean()/step)
        min1  = dm - 0.5 * diff
        max1  = dm + 0.5 * diff
        if min1 < 0:
            min1 = 0
            max1 = diff

    return [min1, max1, min2, max2]


#-----------------------------------------------------------------------------------------------
#-- plot_s2hvst_vt_ratio: create s2hvst - valid rate/total rate plot                         ---
#-----------------------------------------------------------------------------------------------

def plot_s2hvst_vt_ratio(sfile, indata):
    """
    create s2hvst - valid rate/total rate plot
    input:  sfile   --- input file name
            indata  --- a list of lists of data
    output: <head_name>_s2hvst_vt_ratio.png
    """
#
#--- create ratio data from vald and total count rate
#
    [rdata, rerr] =  compute_ratio(indata[12], indata[9], indata[14], indata[11])
#
#--- setting plotting surface
#
    mpl.rcParams['font.size'] = 14
    props = font_manager.FontProperties(size=12)

    xmin   = 4
    xmax   = 9
    psize  = 9
    fsize  = 14
    tfsize = fsize +2

    ymin   = 0.20
    ymax   = 0.50
    a1     = plt.subplot(111)
    plt.axis([xmin, xmax, ymin, ymax])
    plt.plot(indata[21], rdata, lw=0, marker='+',  markersize=psize)
# 
#--- add fitting line
#
    (sint, slope, serror) = robust.robust_fit(indata[21], rdata, iter=100)
    start = sint + slope * xmin
    stop  = sint + slope * xmax
    plt.plot([xmin, xmax],[start,stop], lw=2.0, marker='+',  markersize=psize)

    xdiff = xmax - xmin
    ydiff = ymax - ymin
    xpos  = xmin + 0.1 * xdiff
    ypos  = ymax - 0.1 * ydiff

    text  = 'Slope: ' + str(round(slope, 4) ) + '+/-' + str(round(serror, 4))
    plt.text(xpos, ypos, text, size=tfsize, weight='bold')

#
#--- add the legend and x axis name
#
    a1.set_xlabel("S2HVST", size=fsize)
    a1.set_ylabel("Valid Rate / Total Rate", size=fsize)
#
#--- save the plot in file
#
    save_plot_in_file(sfile, 's2hvst_vt_ratio.png')

#-----------------------------------------------------------------------------------------------
#-- compute_ratio: compute ratio of two entries and the std of the ratio                     ---
#-----------------------------------------------------------------------------------------------

def compute_ratio(vlist1, vlist2, slist1, slist2):
    """
    compute ratio of two entries and the std of the ratio
    input:  vlist1  --- a list of data which holds values of numerator
            vlist2  --- a list of data which holds values of denominator
            slist1  --- a list of std of the data in vlist1
            slist2  --- a list of std of the data in vlist2
    output: ratio   --- a list of ratio data
            std     --- a list of std of the ratio data
    """

    ratio = []
    std   = []
    for k in range(0, len(vlist1)):
        try:
#
#--- compute ratio
#
            rval = vlist1[k] / vlist2[k]
            ratio.append(rval)
#
#--- compute std
#
            sval = sqrt(((slist1[k]/vlist1[k])**2 + (slist2[k]/vlist2[k])**2) * rval**2)
            std.append(sval)
        except:
            ratio.append(0)
            std.append(0)

    return [ratio, std]

#-----------------------------------------------------------------------------------------------
#-- plot_dead_time_trend: update hrc_i_115_dtf_plot.png plot file                             --
#-----------------------------------------------------------------------------------------------

def plot_dead_time_trend():
    """ 
    update hrc_i_115_dtf_plot.png plot file
    input:  none, but read from hrc_i_115_stat_results
    output: <plot dir>/hrc_i_115_dtf_plot.png
    """ 
#
#--- read data
#
    dfile = data_dir + 'Stats/hrc_i_115_stat_results'
    f     = open(dfile, 'r')
    data  = [line.strip() for line in f.readlines()]
    f.close()

    time     = []
    cnt_rate = []
    v_rate   = []
    dtf      = []
#
#--- computer dead time correciton fuctor
#
    for ent in data:
        atemp     = re.split('\s+', ent)
        date      = tcnv.sectoFracYear(float(atemp[1]))
        a_rate    =  float(atemp[5]) / float(atemp[13])

        time.append(date)
        cnt_rate.append(float(atemp[5]))
        v_rate.append(float(atemp[13]))
        dtf.append(a_rate)

#
#--- setting a plotting surface
#
    prep_plot()
    mpl.rcParams['font.size'] = 11
    props = font_manager.FontProperties(size=11)
    psize  = 8
    tfsize =12 
#
#--- start plotting
#
    xmin   = 1999 
    xmax   = int(max(time)) + 1
#
#--- Count Rate per Sec
#
    a1 = plt.subplot(311)
    ymin   = 3.45
    ymax   = 3.55
    title  = 'Count Rate per Sec'
    plt.axis([xmin, xmax, ymin, ymax])
    plt.plot(time, cnt_rate, color='blue', lw=0, marker='+',  markersize=psize)
    plt.text(2000, 3.538, title, size=tfsize, weight='bold')
#
#--- Valid Count Rate per Sec
#
    a2 = plt.subplot(312)
    ymin   = 0
    ymax   = 200 
    title  = 'Valid Count Rate per Sec'
    plt.axis([xmin, xmax, ymin, ymax])
    plt.plot(time, v_rate, color='blue', lw=0, marker='+',  markersize=psize)
    plt.text(2000, 160, title, size=tfsize, weight='bold')
#
#--- Dead Time Correction Fator
#
    a3 = plt.subplot(313)
    ymin   = 0
    ymax   = 0.1 
    title  = 'Dead Time Correctin Factor'
    plt.axis([xmin, xmax, ymin, ymax])
    plt.plot(time, dtf, color='blue', lw=0, marker='+',  markersize=psize)
    plt.text(2000, 0.02, title, size=tfsize, weight='bold')
#
#--- add x ticks label only on the panels of the last row
#
    plt.setp(a1.get_xticklabels(), visible=False)
    plt.setp(a2.get_xticklabels(), visible=False)

    a3.set_xlabel("Time (Year)", size=11)
#
#--- save the plot in file
#
    outname = 'hrc_i_115_dtf_plot.png'
    outname = plot_dir + outname
    create_plot_file(outname)


#-----------------------------------------------------------------------------------------------
#-- poly_fit: polynomial fit routine                                                          --
#-----------------------------------------------------------------------------------------------

def poly_fit(xin, yin, xmin, xmax, deg=2, snum=200, shift=0):
    """
    polynomial fit routine
    input:  xin     --- a list of independent data
            yin     --- a list of dependent data
            xmin    --- min of x axis (could be different from min(xin))
            xmax    --- max of x axis
            deg     --- degree of polynomial    default: 2
            snum    --- numbers of steps        defalut: 200
            shift   --- data shift to be added  defalut: 0
    """
#
#--- convert lists into numpy arrays
#
    x     = numpy.array(xin) + shift
    y     = numpy.array(yin) 
#
#--- compute poli fit coefficients
#
    coeff = list(numpy.polyfit(x, y, deg))
#
#--- create x and y lists for plotting based on coeff
#
    fsnum = float(snum)
    step  = (xmax - xmin) / fsnum
    xsave = []
    ysave = []
    for i in range(0, int(snum)):
        xp = step * i
#
#--- estimate y value
#
        yp = 0
        for k in range(0, deg+1):
            pos = deg -k
            yp = yp + coeff[pos] * xp**k
#
#--- shift back to the original position
#
        xp -= shift

        xsave.append(xp)
        ysave.append(yp)

    return [coeff, xsave, ysave]

#-----------------------------------------------------------------------------------------------
#-- set_poly_eqation_disp: create line fit equeation for display                              --
#-----------------------------------------------------------------------------------------------

def set_poly_eqation_disp(coeff, xd = 'T', iround = 4, shift= 0):
    """
    create line fit equeation for display
    input:  coeff   --- a list of coefficient. highest order first
            xd      --- independent val key     default: "T"
            iround  --- rounding position       default: 4
            shift   --- data shifted quantity   default: 0
    output: line    --- created line, e.g.: Line Fit: 0.100 + 0.100*T + 0.100*T^2
    """
#
#--- round the coefficients for display 
#
    out   = []
    for ent in coeff:
        val = round(ent, iround)
        out.append(val)

    dim = len(coeff)

    line = 'Line Fit: ' 

    for k in range(1, dim+1):
#
#--- first order
#
        if k == 1:
            line = line + str(out[dim-k] - shift)
        else:
#
#--- if the value is negative, change it to positive and set the sign to ' - '
#
            val = out[dim-k]
            if val < 0.0:
                val *= -1
                sind = ' - ' 
            else:
                sind = ' + ' 
#
#--- second order
#
            if k == 2:
                td = '*T'
#
#--- order higher than second
#
            else:
                td = '*T^' + str(k-1)
            line = line + sind + str(val) + td

    return line

#-----------------------------------------------------------------------------------------------
#-- plot_panel: plotting each panel                                                          ---
#-----------------------------------------------------------------------------------------------

def plot_panel(data, xmin, xmax, ymin, ymax, entLabels, title, lsize=1, psize=2, fsize=9, ebar=1, tpos = 1, step=4, txpos=''):
    """
    plotting each panel
    input:  data    --- list of data lists. there are <step> sets of data 
                        for example, if the data set depending of s2hvst value
                            structure: [time1, time2.., time4, val1,..., val4, sigma1,..., sigma4]
                            each entry is a list
            xmin    --- x axis min
            xmax    --- x axis max
            ymin    --- y axis min
            ymax    --- y axis max
            entLabels --- a list of label of each data
            title   --- the title of the plot
            lsize   --- line size   default: 1
            psize   --- point size  default: 2
            fsize   --- font size   default: 9
            ebar    --- indicator to show whether add error bar default: 1 --- yes
            tpos    --- position of title 1: left top, 2: right top, 3: left bottom, 4:right bottom
            step    --- how many data set in the same cateogry (see data above for the example)
            txpos   --- indicator of xpos setting for title
    output: return panel for the farther operation
    """
#
#--- set the plotting surface
#
    plt.axis([xmin, xmax, ymin, ymax])

    for k in range (0, step):
#
#--- set data positions
#
        pos1 = k                    #--- position of the time data
        pos2 = pos1 + step          #--- position of the main data
        pos3 = pos2 + step          #--- position of the sigma
#
#--- plot data with error bars; fmt='o' indicates no connecting lines between data points
#--- if you want to connect, use '-o'
#
        if ebar == 1:
            plt.errorbar(data[pos1], data[pos2], yerr=data[pos3], color=color[k*2], fmt='o', markersize=psize, label=entLabels[k])
        else:
            plt.plot(data[pos1], data[pos2], color=color[k*2], lw=0, marker='+',  markersize=psize, label=entLabels[k])
#
#--- add the title to the plot
#
    if txpos == '':
        xpos = xmin + 0.05 * (xmax - xmin)
    else:
        xpos = xmin * 1.05
    ypos = ymax - 0.10 * (ymax - ymin)
    if tpos == 2:
        xpos = xmax - 0.3 * (xmax - xmin)
    if tpos == 3:
        ypos = ymin + 0.3 * (ymax - ymin)
    if tpos == 4:
        xpos = xmax - 0.3 * (xmax - xmin)
        ypos = ymin + 0.3 * (ymax - ymin)
        
    tfsize = fsize + 1
    plt.text(xpos, ypos, title, size=tfsize, weight='bold')


#-----------------------------------------------------------------------------------------------
#-- read_stat_data: read a stat data file and return a list of lists of data                  --
#-----------------------------------------------------------------------------------------------

def read_stat_data(infile):
    """
    read a stat data file and return a list of lists of data
    input:  infile    --- input file name (with a full path)
    output:  0: time            --- time in year date
             1: tstart          --- start time in seconds from 1998.1.1
             2: tstop           --- stop time in seconds from 1998.1.1
             3: duration        --- duration in seconds
             4: total_count     --- total counts
             5: count_per_sec   --- counts per seonds
             6: pha_mean        --- pha mean
             7: pha_median      --- pha median
             8: pha_sigma       --- pha sigma
             9: t_mean          --- total count rate mean
            10: t_median        --- total count rate median
            11: t_sigma         --- total count rate sigma
            12: v_mean          --- valid count rate mean
            13: v_median        --- valid count rate median
            14: v_sigma         --- valid count rate sigma
            15: s_mean          --- shield count rate mean
            16: s_median        --- shield count rate median
            17: s_sigma         --- shield count rate sigma
            18: anti_co_mean    --- anti conicidnece rate mean
            19: anti_co_median  --- anti conicidnece rate median
            20: anti_co_sigma   --- anti conicidnece rate sigma
            21: s2hvst          --- s2hvst value
            22: s2hvlv          --- s2jvlv valie
            23: scint           --- scint 
            24: scint_std       --- scint sigma
    """

    data = hcf.read_file_data(infile)

    date           = []
    tstart         = []
    tstop          = []
    duration       = []
    total_count    = []
    count_per_sec  = []
    pha_mean       = []
    pha_median     = []
    pha_sigma      = []
    t_mean         = []
    t_median       = []
    t_sigma        = []
    v_mean         = []
    v_median       = []
    v_sigma        = []
    s_mean         = []
    s_median       = []
    s_sigma        = []
    anti_co_mean   = []
    anti_co_median = []
    anti_co_sigma  = []
    s2hvst         = []
    s2hvlv         = []
    scint          = []
    scint_std      = []
    time           = []

    for ent in data:
        aval  = re.split('\s+', ent)
        aval3 = float(aval[3])
#
#--- if the observation interval is shorter than 900 sec, drop the data set.
#--- these data set are not accurate enough.
#
        if aval3 < 900:
            continue
        else:

            try:
                date.append(aval[0])
                start = float(aval[1])
                tstart.append(start)
                tstop.append(float(aval[2]))
                duration.append(aval3)
                total_count.append(float(aval[4]))
                count_per_sec.append(float(aval[5]))
                pha_mean.append(float(aval[6]))
                pha_median.append(float(aval[7]))
                pha_sigma.append(float(aval[8]))
                t_mean.append(float(aval[9]))
                t_median.append(float(aval[10]))
                t_sigma.append(float(aval[11]))
                v_mean.append(float(aval[12]))
                v_median.append(float(aval[13]))
                v_sigma.append(float(aval[14]))
#
#--- changing a plotting range for an easy view
#
                s_mean.append(0.001   * float(aval[15]))
                s_median.append(0.001 * float(aval[16]))
                s_sigma.append(0.001  * float(aval[17]))
         
                anti_co_mean.append(float(aval[18]))
                anti_co_median.append(float(aval[19]))
                anti_co_sigma.append(float(aval[20]))
         
                s2hvst.append(int(float(aval[21])))
                s2hvlv.append(int(float(aval[22]) + 0.5))
         
                scint.append(float(aval[23]))
                scint_std.append(float(aval[24]))
    
                tval = tcnv.sectoFracYear(start)
                time.append(tval)
            except:
                continue

    return [time, tstart, tstop, duration, total_count, count_per_sec, \
            pha_mean,  pha_median, pha_sigma,  \
            t_mean, t_median, t_sigma, \
            v_mean, v_median, v_sigma, \
            s_mean, s_median, s_sigma, \
            anti_co_mean, anti_co_median, anti_co_sigma, \
            s2hvst, s2hvlv, scint, scint_std]

#-----------------------------------------------------------------------------------------------
#-- separate_data_s2hvst: separate data depending on s2hvst value                             --
#-----------------------------------------------------------------------------------------------

def separate_data_s2hvst(fdata, sdata, sigma, s2hvst):
    """
    separate data depending on s2hvst value
    input:  fdata   --- first data
            sdata   --- second data
            sigma   --- sigma of sdata
            s2hvst  --- s2hvst
    output: [fdata5, fdata6, fdata7, fdata8, sdata5, sdata6, sdata7, sdata8, sigma5, sigma6, sigma7, sigma8]
    """
#
#--- initialize
#
    for k in range(5, 9):
        exec "fdata%s = []" % (k)
        exec "sdata%s = []" % (k)
        exec "sigma%s = []" % (k)
#
#--- separate data by s2hvst value; ignore s2hvst other than 5, 6, 7, and 8
#
    for k in range(0, len(fdata)):
        val = str(s2hvst[k])
        try:
            exec "fdata%s.append(%f)" % (val, fdata[k])
            exec "sdata%s.append(%f)" % (val, sdata[k])
            exec "sigma%s.append(%f)" % (val, sigma[k])
        except:
            pass

    return [fdata5, fdata6, fdata7, fdata8, sdata5, sdata6, sdata7, sdata8, sigma5, sigma6, sigma7, sigma8]

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

def separate_data_time(fdata, sdata, sigma, time):
#
#--- initialize
#
    for k in range(0, 2):
        exec "fdata%s = []" % (k)
        exec "sdata%s = []" % (k)
        exec "sigma%s = []" % (k)
#
#--- separate data by  time (2006 - 2009)
#
    for k in range(0, len(time)):
        #if time[k] >= 252460799 and time[k] <= 378691197:
        if time[k] >= 2006 and time[k] <= 2010:
            fdata1.append(fdata[k])
            sdata1.append(sdata[k])
            sigma1.append(sigma[k])
        else:
            fdata0.append(fdata[k])
            sdata0.append(sdata[k])
            sigma0.append(sigma[k])

    return [fdata0, fdata1, sdata0, sdata1, sigma0, sigma1]



#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

def separate_data_ephin(fdata, sdata, sigma, ephin):
#
#--- initialize
#
    for k in range(0, 2):
        exec "fdata%s = []" % (k)
        exec "sdata%s = []" % (k)
        exec "sigma%s = []" % (k)
#
#--- separate data by  ephin (>1)
#
    for k in range(0, len(ephin)):
        if ephin[k] >= 1.0:
            fdata1.append(fdata[k])
            sdata1.append(sdata[k])
            sigma1.append(sigma[k])
        else:
            fdata0.append(fdata[k])
            sdata0.append(sdata[k])
            sigma0.append(sigma[k])

    return [fdata0, fdata1, sdata0, sdata1, sigma0, sigma1]


#-----------------------------------------------------------------------------------------------
#-- create_plot_file: save the plot in a file                                                 --
#-----------------------------------------------------------------------------------------------

def create_plot_file(outname, xsize=10.0, ysize=9.0, resolution=300, fmt='png'):
    """
    save the plot in a file
    input:  outname     --- output file name
            xsize       --- x axis size in inch     default: 10 in
            ysize       --- y axis size in inch     defalut:  8 in
            resolution  --- resolution of the plot  defalut: 300 ppi
            fmg         --- format of the plot      defalut: png
    output: outname     --- plot in <fmt>
    """
#
#--- set the size of the plotting area in inch
#
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(xsize, ysize)
#
#--- save the plot in png format
#
    plt.savefig(outname, format=fmt, dpi=resolution)


#-----------------------------------------------------------------------------------------------
#-- multi_section_linear_fit: compute two line fit on a data and return lists of fitted line data 
#-----------------------------------------------------------------------------------------------

def multi_section_linear_fit(x, y, xmin, xmax, bpoints, snum=100, shift=0):
    """
    compute two line fit on a data and return lists of fitted line data
    input:  x           --- independent data
            y           --- dependent data
            xmin        --- min of x
            xmax        --- max of x
            bpoints     --- a list of breaking point e.g, [2000]
            snum        --- numbers of data points to return
            shift       --- whether shift data. useful when x is far from 0
    output: int_list    --- a list of intercept values
            slope_list  --- a list of slope values
            serr_list   --- a list of errors of slope
            xa_list     --- a list of x points for the fitted line
            ya_list     --- a list of y points for the fitted line
    """
#
#--- add xmax as the last "breaking" point
#
    bpoints.append(xmax)
    bpnum = len(bpoints) 
#
#--- separate data into sections
#
    for k in range(0, bpnum):
        sk = str(k)
        exec "x%s = []" % (str(k))
        exec "y%s = []" % (str(k))
#
#--- when the x is far from 0, robust method does not work well;
#--- so add "shift" (possibly negative) to make the data near 0
#
    for i in range(0, len(x)):
        for k in range(0, bpnum):
            if x[i] <= bpoints[k]:
                exec "x%s.append(%f)" % (str(k), x[i] + shift)
                exec "y%s.append(%f)" % (str(k), y[i])
                break
#
#--- compute intercept and slope for each section
#
    int_list   = []
    slope_list = []
    serr_list  = []
    xa_list    = []
    ya_list    = []
    for k in range(0, bpnum):
        exec "xs = x%s" % (str(k))
        exec "ys = y%s" % (str(k))

        (sint, slope, serror) = robust.robust_fit(xs, ys)

        int_list.append(sint)
        slope_list.append(slope)
        serr_list.append(serror)
#
#--- create plotting data sets
#
        smin = min(xs)
        smax = max(xs)
        step = (smax - smin) / float(snum)
        if k > 0:
            add = bpoints[k-1] + shift
        else:
            add = 0

        tx_list = []
        ty_list = []

        for j in range(0, snum):
            tx  = step * j + add
            ty  = sint + slope * tx
#
#--- shift back x value
#
            tx -= shift
            tx_list.append(tx)
            ty_list.append(ty)

        xa_list.append(tx_list)
        ya_list.append(ty_list)

    return [int_list, slope_list, serr_list, xa_list, ya_list]


#-----------------------------------------------------------------------------------------------

if __name__ == "__main__":

    plot_hrc_trend_data()

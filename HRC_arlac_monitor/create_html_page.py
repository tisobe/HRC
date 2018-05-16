#!/usr/bin/env /proj/sot/ska/bin/python

#############################################################################################################
#                                                                                                           #
#           create_html_page.py: update arlac monitoring html page                                           #
#                                                                                                           #
#           author: t. isobe (tisobe@cfa.harvard.edu)                                                       #
#                                                                                                           #
#           Last Update: Mar 02, 2017                                                                       #
#                                                                                                           #
#############################################################################################################

import sys
import os
import string
import re
import math
import unittest
import time
import numpy
import astropy.io.fits  as pyfits
from datetime import datetime
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
path = '/data/aschrc6/wilton/isobe/Project8/ArLac/Scripts/house_keeping/dir_list'

f    = open(path, 'r')
data = [line.strip() for line in f.readlines()]
f.close()

for ent in data:
    atemp = re.split(':', ent)
    var  = atemp[1].strip()
    line = atemp[0].strip()
    exec "%s = %s" %(var, line)
#
#--- append path to a private folders
#
sys.path.append(mta_dir)
sys.path.append(bin_dir)

import mta_common_functions     as mcf
import convertTimeFormat        as tcnv
import find_coordinate          as fc

#
#--- temp writing file name
#
rtail  = int(time.time())
zspace = '/tmp/zspace' + str(rtail)

#-----------------------------------------------------------------------------------------
#-- create_html_and_plot: update arlac monitoring html page                              --
#-----------------------------------------------------------------------------------------

def create_html_and_plot():
    """
    update arlac monitoring html page
    input:  none:
    output: <html_page>/arlac_vis_montior.html
            <html_dir>/Plots/*.png
    """

#
#--- read html page template
#
    template = house_keeping + '/arlac_template'
    f        = open(template, 'r')
    text     = f.read()
    f.close()
#
#--- go though different settings; set input data, plot file name
#
    for pos in ('i_0', 'i_10', 'i_25', 'i_m10', 'i_m25', 's_0', 's_10', 's_25', 's_m10', 's_m25'):
        dfile = data_dir + 'hrc_' + pos + '_results'
        out   = html_dir + 'Plots/hrc_' + pos + '_rate.png'
        try:
            data  = read_data(dfile)
        except:
            continue
#
#--- plot data
#
        if len(data) >=3:
            plot_data(data, out, pos)
#
#--- substitute the updated data table into the html page
#
            line   = create_html_table(data)
            marker = '#TABLE_' + pos.upper() + '#'
            text   = text.replace(marker, line)
#
#--- change the updated date
#
    ctime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    text = text.replace('#UPDATE#', ctime)
#
#--- print out the html page
#
    hfile = html_dir + 'arlac_vis_montior.html'
    fo    = open(hfile, 'w')
    fo.write(text)
    fo.close()


#-----------------------------------------------------------------------------------------
#-- plot_data: plotting data                                                            --
#-----------------------------------------------------------------------------------------

def plot_data(data, outname, pos):
    """
    plotting data
    input:  data    --- data
            outname --- output file name
            pos     --- indicator of position of aiming
    output  outname --- png plot of the data
    """

    date = []
    cnt  = []
    err  = []
    for ent in data:
        atemp = re.split('\s+', ent)
        #try:
        ytime = convert_to_ytime(atemp[2])
        exp   = float(atemp[3])
        val   = float(atemp[4])
        if val < 10:
            continue
        val  /= exp
        sig   = float(atemp[5])/exp

#        if (pos == 's_25' or pos == 's_m25') and (val <= 10.0):
#            continue

        date.append(ytime)
        cnt.append(val)
        err.append(sig)
        #except:
        #    continue
#
#--- fit a linear line
#
    [a, b, sa, sb] =  line_fit(date, cnt, err)
#
#--- set min max
#
    xmin = 1999.0
    ta   = time.localtime()
    xmax = ta.tm_year + 1
    ymin = 0.0
    ymax = 10
#
#--- start plotting
#
    plt.close('all')
    mpl.rcParams['font.size'] = 9
    props = font_manager.FontProperties(size=9)

    ax  = plt.subplot(111)
    ax.set_autoscale_on(False)
    ax.set_xbound(xmin,xmax)
    ax.set_xlim(xmin=xmin, xmax=xmax, auto=False)
    ax.set_ylim(ymin=ymin, ymax=ymax, auto=False)

    plt.errorbar(date, cnt, yerr=err, fmt='o', lw=1)
    plt.xlabel('Time (year)')
    plt.ylabel('Source Rate (cnt/s)')
#
#--- plot fitting line
#
    start = a + b * xmin
    stop  = a + b * xmax
    plt.plot([xmin, xmax], [start, stop], lw =1, color='blue')
#
#--- write the fitting equation
#
    xdiff = xmax - xmin
    ydiff = ymax - ymin
    xpos  = xmin + 0.1  * xdiff
    ypos  = ymax - 0.1  * ydiff
    ac    = '%1.2f' % (round(a,  2))
    if abs(b) < 0.01:
        bc    = '%1.4f' % (round(abs(b),  4))
        ec    = '%1.4f' % (round(sb, 4))
    else:
        bc    = '%1.2f' % (round(abs(b),  2))
        ec    = '%1.2f' % (round(sb, 2))
    if b > 0 :
        text  = '(Source Rate) = ' + ac  + ' + (' + bc + '+/-' + ec + ') * Time'
    else:
        text  = '(Source Rate) = ' + ac  + ' - (' + bc + '+/-' + ec + ') * Time'
    plt.text(xpos, ypos, text)
#
#--- set the size of the plotting area in inch (width: 10.0in, height 5 in)
#
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(10.0, 5.0)

#
#--- save the plot in png format
#
    plt.savefig(outname, format='png', dpi=100)


#-----------------------------------------------------------------------------------------
#-- convert_to_ytime: convert <yyyy>-<mm>-<dd>T<hh>:<mm>:<ss> into fractional year format 
#-----------------------------------------------------------------------------------------

def convert_to_ytime(sdate):
    """
    convert <yyyy>-<mm>-<dd>T<hh>:<mm>:<ss> into fractional year format. 
    hours, minuts, and seconds are ignored
    input:  sdate   date in yyyy>-<mm>-<dd>T<hh>:<mm>:<ss> format
    output: ytime   date in fractional year format
    """

    m_add = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]

    atemp = re.split('T', sdate)
    btemp = re.split('\-', atemp[0])
    year  = float(btemp[0])
    mon   = float(btemp[1])
    imon  = int(mon)
    day   = float(btemp[2])

    ydate = day + m_add[imon -1]

    base  = 365
    chk   = 4 * int(0.25 * year)
    if(chk == year):
        base = 366
        if mon > 2:
            ydate = ydate + 1

    ytime = year + ydate / base

    return ytime

#-----------------------------------------------------------------------------------------
#-- line_fit: fit a weighted linear line fit                                            --
#-----------------------------------------------------------------------------------------

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

#-----------------------------------------------------------------------------------------
#-- create_html_table: create a html table based on the given data                      --
#-----------------------------------------------------------------------------------------

def create_html_table(data):
    """
    create a html table based on the given data
    input:  data    --- a table data
    output: line    --- a html table 
    """

    line = '<table border=1 cellpadding=2 cellspacing=2>\n'
    line = line + '<tr>\n'
    line = line + '<th>ObsID</th>\n'
    line = line + '<th>Filename</th>\n'
    line = line + '<th>Date</th>\n'
    line = line + '<th>Exposure</th>\n'
    line = line + '<th>Net Counts</th>\n'
    line = line + '<th>Error</th>\n'
    line = line + '<th>DeadTime Correction</th>\n'
    line = line + '</tr>\n'
     
    out  = sorted_by_time(data)

    for ent in out:
        atemp = re.split('\s+', ent)
        line = line + '<tr>\n'
        line = line + '<td style="text-align:right">'  + str(atemp[0]) + '</td>\n'
        line = line + '<td style="text-align:center">' + str(atemp[1]) + '</td>\n'
        line = line + '<td style="text-align:center">' + str(atemp[2]) + '</td>\n'
        line = line + '<td style="text-align:center">' + str(atemp[3]) + '</td>\n'
        line = line + '<td style="text-align:center">' + str(atemp[4]) + '</td>\n'
        line = line + '<td style="text-align:center">' + str(atemp[5]) + '</td>\n'
        line = line + '<td style="text-align:center">' + str(atemp[6]) + '</td>\n'
        line = line + '</tr>\n'

    line = line + '</table>\n\n'

    return line

#-----------------------------------------------------------------------------------------
#-- sorted_by_time: rearrange the data in time order                                    --
#-----------------------------------------------------------------------------------------

def sorted_by_time(data):
    """
    rearrange the data in time order
    input:  data    --- assume that the third entry is time in the form of <yyyy>-<mm>-<dd>T<hh>:<mm>:<ss>
    ouput:  out     --- time sorted data
    """

    tlist = []
    tdict = {}
    for ent in data:
        atemp = re.split('\s+', ent)
        tm = atemp[2]
        tm = tm.replace('-', '')
        tm = tm.replace('T', '')
        tm = tm.replace(':', '')
        tm = int(float(tm))
        tlist.append(tm)
        tdict[tm] = ent

    tlist.sort()

    out = []
    for ent in tlist:
        out.append(tdict[ent])

    return out


#-----------------------------------------------------------------------------------------
#-- read_data: read data file                                                           --
#-----------------------------------------------------------------------------------------

def read_data(infile, remove=0):

    f    = open(infile, 'r')
    data = [line.strip() for line in f.readlines()]
    f.close()

    if remove == 1:
        mcf.rm_file(infile)

    return data


#-----------------------------------------------------------------------------------------

if __name__ == '__main__':

    create_html_and_plot()

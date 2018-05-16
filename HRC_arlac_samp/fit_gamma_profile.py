#!/usr/bin/env /proj/sot/ska/bin/python

#####################################################################################################
#                                                                                                   #
#           fit_gamma_profile.py: fitting gamma function on a given data                            #
#                                                                                                   #
#           author: t. isobe(tisobe@cfa.harvard.edu)                                                #
#                                                                                                   #
#           Last Update:    May 24, 2017                                                            #
#                                                                                                   #
#####################################################################################################

import os
import sys
import re
import string
import operator
import math
import numpy
import astropy.io.fits  as pyfits
import scipy
import scipy.stats      as stats 
import time

from pylab import *
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import matplotlib.lines as lines

import matplotlib as mpl

if __name__ == '__main__':
        mpl.use('Agg')
#
#--- reading directory list
#
path = '/data/aschrc6/wilton/isobe/Project8/HZ43/Scripts/house_keeping/dir_list'

f= open(path, 'r')
data = [line.strip() for line in f.readlines()]
f.close()

for ent in data:
    atemp = re.split(':', ent)
    var  = atemp[1].strip()
    line = atemp[0].strip()
    exec "%s = %s" %(var, line)
#
#--- append  pathes to private folders to a python directory
#
sys.path.append(bin_dir)
sys.path.append(mta_dir)
#
#--- import several functions
#
import convertTimeFormat          as tcnv       #---- contains MTA time conversion routines
import mta_common_functions       as mcf        #---- contains other functions commonly used in MTA scripts

from kapteyn import kmpfit

#
#--- temp writing file name
#
rtail  = int(time.time())
zspace = '/tmp/zspace' + str(rtail)

#---------------------------------------------------------------------------------------------------
#--- fit_gamma_profile: fitting Gamman function to the given data                                ---
#---------------------------------------------------------------------------------------------------

def  fit_gamma_profile(x, y, avg, std, plot_title='', xname='',yname=''):
    """
    fitting Gamman function to the given data
    input:  x           --- x data
            y           --- y data
            avg         --- average of the distribution
            std         --- std of the distribution
            plot_title  --- title of the plots; if it is not empty, a plot is created
            xname       --- name of the x axis
            yname       --- name of the y axis
    output: [a, b, h]   --- shape, sale parameters and normlization factor
            if plot_title is not "", it also create a png plot ('out.png')
    """
#
#--- make the initial guess of shape and scale parameters
#
    a0    = (avg/std)**2
    b0    = avg/(std*std)
    m0    = max(y)

    [a, b, h] = gamma_fit(x, y, a0, b0, m0)
#
#--- plot the data
#
    if plot_title != '':
            plot_gamma(x, y, a, b, h, plot_title, xname, yname)


    return  [a, b, h]

#---------------------------------------------------------------------------------------------------
#-- plot_gamma: plotting a fitted result on the data                                             ---
#---------------------------------------------------------------------------------------------------

def plot_gamma(x, y, a,b,h, name='',xname='', yname='', outfile='out.png'):

    """
    plotting a fitted result on the data
    Input:      x   ---- a list of bin value
                y   ---- a list of phas
                a
                b
                amp ---- amp of the normal distribution
                center ---- center of the normal distribution
                width  ---- width of the normal distribution
                name   ---- title of the plot; default ""
                outfile --- output file name; defalut "out.png"
    Output:     outfile in <plot_dir>
    """
#
#--- set plotting range
#
    xmin = int(min(x))
    xmax = int(max(x)) + 1
    ymin = 0
    ymax = int(1.1 * max(y))
#
#--- create Gaussian distribution
#
    p    = (a, b)
    est_v = []
    est_list = []
    for v in range(0, xmax):
        est = gfunction(v, a, b, h)
        est_v.append(v)
        est_list.append(est)
#
#--- clean up the plotting device
#
    plt.close('all')
#
#---- set a few parameters
#
    mpl.rcParams['font.size'] = 12
    props = font_manager.FontProperties(size=12)
    
    plt.axis([xmin, xmax, ymin, ymax])
    plt.xlabel(xname, size=12)
    plt.ylabel(yname, size=12)
    plt.title(name, size=12)
#
#--- set information text postion and content
#
    tx  = 0.5 * xmax
    ty  = 0.9 * ymax
    avg = a/b
    sig = math.sqrt(a/(b*b))
    line = 'Mean: ' + str(round(avg,2)) + '   Sigma: ' + str(round(sig,2)) + '   Norm: ' + str(round(h,2))


    plt.text(tx, ty, line)
#
#--- actual plotting
#
    p, = plt.plot(x, y, marker='.', markersize=4.0, lw = 0)
    p, = plt.plot(est_v, est_list, marker= '', lw=3)
#
#--- set the size of the plotting area in inch (width: 10.0in, height 2.08in x number of panels)
#
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(10.0, 5.0)
#
#--- save the plot in png format
#
    plt.savefig(outfile, format='png', dpi=100)

#---------------------------------------------------------------------------------------------------
#-- gfunction: gamman function                                                                    --
#---------------------------------------------------------------------------------------------------

def gfunction(x, a, b, h):
    """
    gamman function 
    input:  x   --- x value
            a   --- shape parameter
            b   --- scale parameter
            h   --- normalization
    output: g   --- value returned
    """
    gam = scipy.special.gamma(a)
    g = h * x**(a-1) * exp(-x*b) * (b**a/gam)

    return g

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

def residualsG(p, data):
    """
    Return weighted residuals of Gamma function
    input:   p       --- parameter list (a, b, h)
            data    --- data set (x, y , err); err is ignored
    output: residual
    """

    x, y, err = data
    a, b, h   = p

    return (y- gfunction(x, a, b, h)) / err

#---------------------------------------------------------------------------------------------------
#-- gamma_fit: fitting gamma profile to the data                                                 ---
#---------------------------------------------------------------------------------------------------

def gamma_fit(x, y, a0, b0, m0):

    """
    fitting gammaf profile to the data
    Input:      x   --- independent var
                y   --- dependent var
    Output: [a, b]
    """

    N = len(y)
    x   = numpy.array(x)
    y   = numpy.array(y)
    err = numpy.ones(N)
#
#--- initial guess
#
    
    p0 = [a0, b0, m0]

#
#--- fit the model
#
    try:
        fitter = kmpfit.Fitter(residuals=residualsG, data=(x,y,err))
        #fitter.parinfo = [{}, {}, {'fixed':True}]  # Take zero level fixed in fit
        fitter.fit(params0=p0)
    
        a, b, h = fitter.params
    except:
        a, b, h = p0

    return [a, b, h]

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

def read_data(dfile):

    f     = open(dfile, 'r')
    data  = [line.strip() for line in f.readlines()]
    f.close()
    fdata = []
    for ent in data:
        fdata.append(int(float(ent)))

    avg = numpy.mean(fdata)
    std = numpy.std(fdata)

    x = []
    y = []
    for k in range(1, 257):
        x.append(k)
        y.append(0)

    for ent in fdata:
        val = int(float(ent))
        y[val] += 1
#

    return [x, y, avg, std]

#--------------------------------------------------------------------
if __name__ == "__main__":


    if len(sys.argv) == 2:
        dfile = sys.argv[1]

        [x, y] = read_data(dfile)

        fit_gamma_profile(x, y, avg, std,  plot_title='test', xname='Channel',yname='Count Rate Rate')

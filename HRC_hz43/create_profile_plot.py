#!/usr/bin/env /proj/sot/ska/bin/python

#####################################################################################################
#                                                                                                   #
#           create_profile_plot.py: fit voigt profile on the data and creat a plot                  #
#                                                                                                   #
#           author: t. isobe(tisobe@cfa.harvard.edu)                                                #
#                                                                                                   #
#           Last Update:    May 02, 2017                                                            #
#                                                                                                   #
#####################################################################################################

import os
import sys
import re
import string
import random
import operator
import math
import numpy
import time
import astropy.io.fits  as pyfits
import unittest

#
#--- from ska
#
from Ska.Shell import getenv, bash

ascdsenv = getenv('source /home/ascds/.ascrc -r release', shell='tcsh')
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
html_top = html_top.replace('#', ':')
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
import fit_voigt_profile          as voigt
import gamma_function             as gamma

#
#--- temp writing file name
#
rtail  = int(time.time())
zspace = '/tmp/zspace' + str(rtail)

#---------------------------------------------------------------------------------------------------
#-- run_for_all_profile_plot: read which obsids are available for profile plotting and plot them   -
#---------------------------------------------------------------------------------------------------

def run_for_all_profile_plot():
    """
    read which obsids are available for profile plotting and plot them
    input: none
    output: <web_dir>/Plots/Indivisual_Plots/<obsid>/pi_p_list_<#>_vifits.png etc
    """
    sdir  = data_dir + 'Fittings/'
    cmd   = 'ls -d ' + sdir + '/*  > ' + zspace
    os.system(cmd)
    f     = open(zspace, 'r')
    olist = [line.strip() for line in f.readlines()]
    f.close()
    mcf.rm_file(zspace)

    for ent in olist:
        atemp = re.split('\/', ent)
        obsid = atemp[-1]

        print 'ObsID: ' + str(obsid)
        create_profile_plot(obsid)

#---------------------------------------------------------------------------------------------------
#-- create_profile_plot: fit voigt profile on the data and creat a plot                          ---
#---------------------------------------------------------------------------------------------------

def  create_profile_plot(obsid):
    """
    fit voigt profile on the data and creat a plot
    input:  obsid   --- obsid. the data is read from <data_dir>
    output: <web_dir>/Plots/Indivisual_Plots/<obsid>/<head>_<col #>_vfits.png
    """
    head_list = ['samp_center_list', 'pi_center_list', 'samp_p_list_', 'samp_n_list_', 'pi_p_list_', 'pi_n_list_']

    sdir = data_dir + 'Fittings/' + str(obsid) + '/'

    for m in range(0, len(head_list)):
        out  = sdir + head_list[m] +  'fit_results'
        fo   = open(out, 'w')
        if m < 2:
            fname = sdir + head_list[m]
            line  = create_plot_and_table_entry(head_list[m], '', obsid,  fname)
            if line:
                fo.write(line)
            else:
                continue


        else:
            for k in range(0, 20):
                fname = sdir + head_list[m] + str(k)
                line  = create_plot_and_table_entry(head_list[m], k, obsid,  fname)
                if line:
                    fo.write(line)
                else:
                    continue

        fo.close()

            
#---------------------------------------------------------------------------------------------------
#-- create_plot_and_table_entry: fit a voigt profile on the data, create a plot, and return the fitting results 
#---------------------------------------------------------------------------------------------------

def create_plot_and_table_entry(head, pos, obsid, fname):
    """
    fit a voigt profile on the data, create a plot, and return the fitting results
    input:  head    --- head part of the file
            pos     --- position of the data. either center of cell <#>
            obsid   --- obsid
            fname   --- the data file name
    output: <web_dir>/Plots/Indivisual_Plots/<obsid>/<head>_<col #>_vfits.png
            line    --- voight fitted results in a table column form
    """

    if not  os.path.isfile(fname):
        return False
    try:
        out = binning_data(fname)
        if out == False:
            return False
        else:
            [abin, data, avg, std] = out
    except:
        return False
#
#--- fit a voigt distribution on the data
#
    title = 'ObsID: ' + str(obsid)
    [a, b, h]  = gamma.fit_gamma_profile(abin, data, avg, std, plot_title=title)
    #[center, width, amp, alphaD, alphaL, I, a_back, b_back]  \
        #= voigt.fit_voigt_profile(abin, data, type='voigt', plot_title=title)
#
#--- rename a plotting file
#            
    outdir  = web_dir + 'Plots/Indivisual_Plots/' + str(obsid) + '/' 

    if not os.path.isdir(outdir):
        cmd = 'mkdir ' + outdir
        os.system(cmd)

    if not os.path.isdir(outdir):
        cmd = 'mkdir ' + outdir
        os.system(cmd)

    outfile = outdir + head + str(pos) + '_vfit.png'
    cmd = 'mv out.png ' + outfile
    os.system(cmd)
#
#--- keep the fitting results
#
    line = str(pos) + '\t'
    #line = line + str(center) + '\t'
    #line = line + str(width)  + '\t'
    #line = line + str(amp)    + '\t'
    #line = line + str(alphaD) + '\t'
    #line = line + str(alphaL) + '\t'
    #line = line + str(I)      + '\t'
    #line = line + str(a_back) + '\t'
    #line = line + str(b_back) + '\n'
    line = line + str(a) + '\t'
    line = line + str(b) + '\t'
    line = line + str(h) + '\n'

    return line

#---------------------------------------------------------------------------------------------------
#-- binning_data: binning the data into 256 bins                                                 ---
#---------------------------------------------------------------------------------------------------

def binning_data(fname):
    """
    binning the data into 256 bins
    input:  fname   --- data file name
    output: abin    --- a list of bin # 0 - 255
            out     --- a list of the counts of each bin
    """

    f     = open(fname, 'r')
    data  = [line.strip() for line in f.readlines()]
    f.close()

    out = [0] * 256

    fdata = []
    for ent in data:
        fdata.append(int(float(ent)))

    avg = numpy.mean(fdata)
    std = numpy.std(fdata)

    for val in fdata:
        if val >= 256:
            continue
        out[val]  += 1

    abin = []
    for k in range(0, 256):
        abin.append(k)

    if str(avg) == 'nan':
        return False

    else:
        return [abin, out, avg, std]



#--------------------------------------------------------------------

if __name__ == '__main__':

    if len(sys.argv) > 1:
        obsid = sys.argv[1]
        obsid = int(float(obsid))

        create_profile_plot(obsid)
    else:
        run_for_all_profile_plot()



#!/usr/bin/env /proj/sot/ska/bin/python

#############################################################################################################
#                                                                                                           #
#       save_proccessed_data.py: save cerated data in subdirectries                                         #
#                                                                                                           #
#           author: t. isobe (tisobe@cfa.harvard.edu)                                                       #
#                                                                                                           #
#           Last Update: Jan 27, 2017                                                                       #
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
#
#--- from ska
#
from Ska.Shell import getenv, bash

ascdsenv = getenv('source /home/ascds/.ascrc -r release; source /home/mta/bin/reset_param ', shell='tcsh')
#
#--- reading directory list
#
path = '/data/mta/Script/Python_script2.7/house_keeping/dir_list'

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

import mta_common_functions as mcf
import convertTimeFormat    as tcnv
#
#--- temp writing file name
#
rtail  = int(time.time())
zspace = '/tmp/zspace' + str(rtail)

NULL   = 'NULL'

top1 = '/data/aschrc1/GENHRC/RAW/HRC_ENG/'
top2 = '/data/aschrc1/GENHRC/RAW/HRC_HK0/'

#--------------------------------------------------------------------------------
#-- move_files: create a monthly combined rate data file/move the data into a yearly sub directory
#--------------------------------------------------------------------------------

def move_files():
    """
    create a monthly combined rate data file and also move the data into a yearly
    sub directory at the beginning of the year. 
    input:  none
    output: rates_<yyyy>_<mm>.fits.gz
            HRC_ENG/<yyyy>/Eng0/*.fits.gz
            HRC_ENG/<yyyy>/Rates/*.fits
            HRC_HK0/<yyyy>/*.fits
    """
#
#--- find today's date
#
    tlist = time.localtime()
    tyear = int(float(tlist[0]))
    tmon  = int(float(tlist[0]))
    tday  = int(float(tlist[0]))

#
#--- make a monthly combined rates file at every month around day 5 or 6
#

    if tday >= 5 and tday <= 7:
        year = tyear
        mon  = tmon -1
        if mon < 1:
            mon = 12
            year -= 1
        
        input_head  = top1 + 'hrc_rates_'
        output_head = top1 + 'rates_'
        combine_rates_data(year, mon, input_head, output_head)
        
        input_head  = top1 + 'hrc_4_eng0_'
        output_head = top1 + 'eng0_4_'
        combine_rates_data(year, mon, input_head, output_head)
        
        input_head  = top1 + 'hrc_5_eng0_'
        output_head = top1 + 'eng0_5_'
        combine_rates_data(year, mon, input_head, output_head)
        
        input_head  = top2 + 'hrc_hk0_'
        output_head = top2 + 'hk0_'
        combine_rates_data(year, mon, input_head, output_head)


#------------------------------------------

        lmon = str(mon)
        if mon < 10:
            lmon = '0' + lmon
#
#--- ENG0 dir; create new save directories
#
        ddir = top1 + str(year) 
        if not os.path.isdir(ddir):
            cmd = 'mkdir ' + ddir
            os.system(cmd)

            cmd = 'mkdir ' + ddir + '/Rates/'
            os.system(cmd)

            cmd = 'mkdir ' + ddir + '/Eng0/'
            os.system(cmd)
#
#--- rates files
#
        cmd = 'ls ' + top1 + 'hrc_rates_' + str(year) + lmon + '*.fits* > ' + zspace
        os.system(cmd)
        f   = open(zspace, 'r')
        data = [line.strip() for line in f.readlines()]
        f.close()
        mcf.rm_file(zspace)

        for ent in data:
            cmd = 'mv ' + str(ent) + ' ' + ddir + '/Rates/. '
            os.system(cmd)
#
#--- eng files
#
        cmd = 'ls ' + top1 + 'hrc_*_eng0_' + str(year) + lmon + '*.fits* > ' + zspace
        os.system(cmd)
        f   = open(zspace, 'r')
        data = [line.strip() for line in f.readlines()]
        f.close()
        mcf.rm_file(zspace)

        for ent in data:
            cmd = 'mv ' + str(ent) + ' ' + ddir + '/Eng0/. '
            os.system(cmd)
#
#--- HK0 dir; create a new directory
#
        ddir = top2 + str(year) 
        if not os.path.isdir(ddir):
            cmd = 'mkdir ' + ddir
            os.system(cmd)

        cmd = 'ls ' + top2 + 'hrc_hk0_' + str(year) + lmon + '*.fits* > ' + zspace
        os.system(cmd)
        f   = open(zspace, 'r')
        data = [line.strip() for line in f.readlines()]
        f.close()
        mcf.rm_file(zspace)

        for ent in data:
            cmd = 'mv ' + str(ent) + ' ' + ddir + '/. '
            os.system(cmd)

#------------------------------------------

#
#--- reate monthly combined files; it will be moved only when at a new year
#
        if mon == 1:
            cmd = 'ls ' + top1 + 'rates_' + str(year) + '_*.fits* > ' + zspace
            os.system(cmd)
            f   = open(zspace, 'r')
            data = [line.strip() for line in f.readlines()]
            f.close()
            mcf.rm_file(zspace)
    
            for ent in data:
                cmd = 'mv ' + str(ent) + ' ' + ddir + '/. '
                os.system(cmd)

#--------------------------------------------------------------------------------
#-- combine_rates_data: combine hrc_reates data into month long data set      ---
#--------------------------------------------------------------------------------

def combine_rates_data(year, mon, input_head, output_head):
    """
    combine specified daily data set into a month long data set
    input:  year    --- year
            mon     --- mon
    output: rates_<yyyy>_<mm>.fits.gz in ..../HRC_ENG/ directory
    """

    lmon = str(mon)
    if mon < 10:
        lmon = '0' + lmon

    cmd  = 'ls ' + input_head + str(year) + lmon  + '*fits.gz > ' + zspace
    os.system(cmd)
    f    = open(zspace, 'r')
    data = [line.strip() for line in f.readlines()]
    f.close()
    mcf.rm_file(zspace)

    if len(data) > 0:

        cmd = 'cp ' + data[0] + ' temp.fits.gz'
        os.system(cmd)
        cmd = 'gzip -d temp.fits.gz'
        os.system(cmd)
    
        cmd1 = '/usr/bin/env PERL5LIB="" ' 
        for k in range(1, len(data)):
            cmd2 = ' dmmerge "temp.fits,' + data[k] + '"  out.fits'
            cmd  = cmd1 + cmd2
            try:
                bash(cmd, env=ascdsenv)
            except:
                try:
                    bash(cmd, env=ascdsenv)
                except:
                    pass
    
            cmd = 'mv out.fits temp.fits'
            try:
                os.system(cmd)
            except:
                pass
    
        out = output_head + str(year) + '_' + lmon + '.fits'
    
        cmd = 'mv temp.fits ' + out
        os.system(cmd)
        cmd =  'gzip ' + out
        os.system(cmd)
        cmd = 'chmod 755 ' + out

#--------------------------------------------------------------------------------

if __name__ == "__main__":

    #move_files()
    year = 2016
    mon  = 12
    for year in range(2012, 2017):
        cmd = 'mkdir  /data/aschrc1/GENHRC/RAW/HRC_HK0/' + str(year)  + '/HK0/'
        os.system(cmd)
        cmd = 'mv /data/aschrc1/GENHRC/RAW/HRC_HK0/' + str(year) + '/hrc_hk0_*.fits.gz '
        cmd = cmd + ' /data/aschrc1/GENHRC/RAW/HRC_HK0/' + str(year)  + '/HK0/.' 
        os.system(cmd)
        for mon in range(1, 13):
            print "Year: " + str(year) + '   Month: ' + str(mon)


            input_head  = '/data/aschrc1/GENHRC/RAW/HRC_HK0/' + str(year) + '/HK0/hrc_hk0_'
            output_head = '/data/aschrc1/GENHRC/RAW/HRC_HK0/' + str(year) + '/hrc_hk0_'
            combine_rates_data(year, mon, input_head, output_head)




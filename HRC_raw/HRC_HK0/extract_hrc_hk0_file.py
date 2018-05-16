#!/usr/bin/env /proj/sot/ska/bin/python

#############################################################################################
#                                                                                           #
#       extract_hrc_hk0_file.py: extract hrc hk0 fits files                                 #
#                                                                                           #
#           author: t. isobe (tisobe@cfa.harvard.edu)                                       #
#                                                                                           #
#           last update: Jan 13, 2017                                                       #
#                                                                                           #
#############################################################################################

import os
import sys
import re
import string
import random
import operator
import time
import datetime
import math
import numpy
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
path = '/data/aschrc1/GENHRC/TOOLS/HRC_HK0/house_keeping/dir_list_py'

f    = open(path, 'r')
data = [line.strip() for line in f.readlines()]
f.close()

for ent in data:
    atemp = re.split(':', ent)
    var  = atemp[1].strip()
    line = atemp[0].strip()
    exec "%s = %s" %(var, line)
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
import fits_operation       as mfits

#
#--- temp writing file name
#
rtail  = int(10000 * random.random())       #---- put a romdom # tail so that it won't mix up with other scripts space
zspace = '/tmp/zspace' + str(rtail)
#outdir = '/data/aschrc1/GENHRC/TOOLS/HRC_HK0/Outdir/'

#-----------------------------------------------------------------------------
#-- run_hrc_scripts:  extract hrc hk0 fits files                           ---
#-----------------------------------------------------------------------------

def run_hrc_scripts():
    """
     extract hrc hk0 fits files
    input: none but extract from the database
    output: hrc_4_eng_<yyyy><mm><dd>.fits
            hrc_rates_<yyyy><mm><dd>.fits
            hrc_5_eng_<yyyy><mm><dd>.fits
    """

#
#--- find today's time in seconds from epoch
#
    today = time.mktime(time.localtime())
#
#--- set the last data collection day to two days ago
#
    tlim  = today - 3 * 86400.0 
#
#--- find the last data file created
#
    [year, mon, day] = find_the_last_data_created_date()
#
#--- convert the last data created date into seconds from epoch
#
    dt    = datetime.datetime(year, mon, day, 0, 0)
    ltime = time.mktime(dt.timetuple())

#
#--- extract data and update tables
#
    while(ltime < tlim):
        ltime += 86400
#
#--- 883630800.0 is the seconds difference between python epoch and 1998.1.1
#
        atime  = ltime  - 883630800.0
        ftime  = time.localtime(ltime)
        year   = int(float(time.strftime("%Y", ftime)))
        mon    = int(float(time.strftime("%m", ftime)))
        day    = int(float(time.strftime("%d", ftime)))

        print str(year) + '<-->' + str(mon) + '<-->' + str(day)

        run_get_hrc_hk0(atime, year, mon, day)

#-----------------------------------------------------------------------------
#-- find_the_last_data_created_date: find the data of the last created fits data file
#-----------------------------------------------------------------------------

def find_the_last_data_created_date():
    """
    find the data of the last created fits data file.
    input:  none
    output: year    --- year
            mon     --- momth
            day     --- day of the month
    """

    cmd   = 'ls /data/aschrc1/GENHRC/RAW/HRC_HK0/hrc_hk0*fits* >' + zspace
    os.system(cmd)

    f     = open(zspace, 'r')
    data  = [line.strip() for line in f.readlines()]
    f.close()
    mcf.rm_file(zspace)

    date  = []
    for ent in data:
        atemp = re.split('hrc_hk0_', ent)
        btemp = re.split('\.fits', atemp[1])
        date.append(btemp[0])

    date.sort()
    lent  = str(date[-1])
    year  = lent[0] + lent[1] + lent[2] + lent[3]
    year  = int(float(year))
    mon   = lent[4] + lent[5]
    mon   = int(float(mon))
    day   = lent[6] + lent[7]
    day   = int(float(day))

    return [year, mon, day]

#-----------------------------------------------------------------------------
#-- run_get_hrc_hk0: extract hrc hk0 fits file and creates combined hk0 file -
#-----------------------------------------------------------------------------

def run_get_hrc_hk0(tstart, year, mon, day):
    """
    extract hrc hk0 fits file and creates combined hk0 file
    input:  tstart  --- time in seconds from 1998.1.1
            year    --- year
            mon     --- month
            day     --- day of the month
    output: hrc_hk0_<yyy><mm><dd>.fits.gz
            bad_2smtratm.rdb    --- updated
    """

    cyear = str(year)
    cmon  = str(mon)
    cday  = str(day)

    if mon < 10:
        cmon = '0' + cmon
    if day < 10:
        cday = '0' + cday
#
#--- set the data extract interval to a day
#
    tstop =  tstart + 86400.0
#
#--- extract  hrc4eng data with arc5gl
#
    run_arc5gl(tstart, tstop)
#
#--- merge all hrc4eng data into one fits file
#
    cmd  = 'ls hrcf*_hk0.fits.gz >dat.lis'
    os.system(cmd)

    cmd1 = "/usr/bin/env PERL5LIB="
    cmd2 = ' dmmerge infile=@dat.lis"[time=' + str(tstart) + ':' + str(tstop) + ',scidpren=0000xxxx000xxxxx]'
    cmd2 = cmd2 + '[cols \!ce00atm,\!ce01atm,\!quality]"'
    cmd2 = cmd2 + ' outfile=hrc_hk0_'+ cyear + cmon + cday +'.fits mode=h'
    cmd  = cmd1 + cmd2
    bash(cmd,  env=ascdsenv)
#
#--- append the results to rdb data tables
#
    cmd2 = ' dmlist hrc_hk0_'+ cyear + cmon + cday +'.fits'
    cmd2 = cmd2 + '"[smtratm=55:57]" blocks | grep HSKPNG >' + zspace
    cmd  = cmd1 + cmd2
    bash(cmd,  env=ascdsenv)
    
    f    = open(zspace, 'r')
    data = f.read()
    f.close()
    mcf.rm_file(zspace)
    data.strip()
    atemp = re.split('\s+', data)
    val   = atemp[7]
    ofile = outdir + 'bad_2smtratm.rdb'
    f     = open(ofile, 'a')
    oline = 'hrc_hk0_'+ cyear + cmon + cday +'.fits.gz\t' + str(val) + '\n'
    f.write(oline)
    f.close()
#
#--- move the fits files to the saving direcotry
#
    cmd   = 'gzip hrc_hk0*fits'
    os.system(cmd)

    cmd   = 'mv hrc_hk0*fits.gz ' + outdir
    os.system(cmd)
#
#--- clean up the original files
#
    cmd  = 'rm `cat dat.lis`'
    os.system(cmd)
    mcf.rm_file('./dat.lis')


#-----------------------------------------------------------------------------
#-- run_arc5gl: run arc5gl                                                  --
#-----------------------------------------------------------------------------

def run_arc5gl(tstart, tstop):
    """
    run arc5gl
    input:  start   --- start time
            stop    --- stop time
    output: extracted fits files
    """

    fo    = open(zspace, 'w')

    line  = 'operation=retrieve\n'
    line  = line + 'tstart=' + str(tstart) + '\n'
    line  = line + 'tstop='  + str(tstop)  + '\n'
    line  = line + 'dataset=flight\n'
    line  = line + 'level=0\n'
    line  = line + 'detector=hrc\n'
    line  = line + 'filetype=hrchk\n'
    line  = line + 'go\n'

    fo.write(line)
    fo.close()

    cmd1 = "/usr/bin/env PERL5LIB="
    cmd2 = ' /proj/axaf/simul/bin/arc5gl -user isobe -script ' + zspace
    cmd  = cmd1 + cmd2
    bash(cmd,  env=ascdsenv)

#-----------------------------------------------------------------------------

if __name__ == '__main__':

    run_hrc_scripts()

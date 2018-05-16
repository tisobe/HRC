#!/usr/bin/env /proj/sot/ska/bin/python

#############################################################################################
#                                                                                           #
#       extract_hrc_eng_files.py: extract hrc eng4 and eng5 data and make summary           #
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
path = '/data/aschrc1/GENHRC/TOOLS/HRC_ENG/house_keeping/dir_list_py'

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
#outdir = '/data/aschrc1/GENHRC/RAW/HRC_ENG/'

#-----------------------------------------------------------------------------
#-- run_hrc_scripts: extract hrc eng4 and eng5 data and make summary       ---
#-----------------------------------------------------------------------------

def run_hrc_scripts():
    """
    extract hrc eng4 and eng5 data and make summary
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

        run_get_hrc_4_eng(atime, year, mon, day)
        run_get_hrc_5_eng(atime, year, mon, day)


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

    cmd   = 'ls /data/aschrc1/GENHRC/RAW/HRC_ENG/hrc_rates*fits* >' + zspace
    os.system(cmd)

    f     = open(zspace, 'r')
    data  = [line.strip() for line in f.readlines()]
    f.close()
    mcf.rm_file(zspace)

    date  = []
    for ent in data:
        atemp = re.split('hrc_rates_', ent)
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
#-- run_get_hrc_4_eng: extract hrc eng4 fits file and creates combined eng4 file and rate file
#-----------------------------------------------------------------------------

def run_get_hrc_4_eng(tstart, year, mon, day):
    """
    extract hrc eng4 fits file and creates combined eng4 file and rate file
    input:  tstart  --- time in seconds from 1998.1.1
            year    --- year
            mon     --- month
            day     --- day of the month
    output: hrc_4_eng0_<yyy><mm><dd>.fits
            hrc_rates_<yyy><mm><dd>.fits
            mcptot_a_stats.rdb  --- updated
            shield_a_stats.rdb  --- updated
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
    run_arc5gl(tstart, tstop, 'hrc4eng')
#
#--- merge all hrc4eng data into one fits file
#
    cmd  = 'ls hrcf*_4_eng0.fits.gz >dat.lis'
    os.system(cmd)

    cmd1 = "/usr/bin/env PERL5LIB="
    cmd2 = ' dmmerge infile=@dat.lis"[time=' + str(tstart) + ':' + str(tstop) + ',quality=0000000000000000000,mnf=0]" '
    cmd2 = cmd2 + 'outfile=hrc_4_eng0_'+ cyear + cmon + cday +'.fits mode=h'
    cmd  = cmd1 + cmd2
    bash(cmd,  env=ascdsenv)
#
#--- calculate rate stat
#
    cmd  = 'rm `cat dat.lis`'
    os.system(cmd)
    mcf.rm_file('./dat.lis')

    cmd2 = ' dmtcalc infile=hrc_4_eng0_' + cyear + cmon + cday + '.fits"[cols '
    cmd2 = cmd2 + '\!QUALITY,\!2CEAHVPT,\!2N15VAVL,\!2N15VBVL,\!2P05VAVL,\!2P05VBVL,'
    cmd2 = cmd2 + '\!2P15VAVL,\!2P15VBVL,\!2P24VAVL,\!2P24VBVL,\!2UVLSPXT]" '
    cmd2 = cmd2 + 'outfile=hrc_rates_' + cyear + cmon + cday + '.fits expression=@col_calc.lis'
    cmd  = cmd1 + cmd2
    bash(cmd,  env=ascdsenv)
#
#--- append the results to rdb data tables
#
    ofile = outdir + 'shield_a_stats.rdb'
    append_to_rdb('shield_a', ofile, cyear, cmon, cday)

    ofile = outdir + 'mcptot_a_stats.rdb'
    append_to_rdb('mcptot_a', ofile, cyear, cmon, cday)
#
#--- move the fits files to the saving direcotry
#
    cmd   = 'mv hrc_4_eng0*fits hrc_rate*fits ' + outdir
    os.system(cmd)

#-----------------------------------------------------------------------------
#-- run_get_hrc_5_eng: extract hrc eng5 fits file and creates combined eng5
#-----------------------------------------------------------------------------

def run_get_hrc_5_eng(tstart, year, mon, day):
    """
     extract hrc eng5 fits file and creates combined eng5
    input:  tstart  --- time in seconds from 1998.1.1
            year    --- year
            mon     --- month
            day     --- day of the month
    output: hrc_5_eng0_<yyy><mm><dd>.fits
            hrc_rates_<yyy><mm><dd>.fits
            mcptot_a_stats.rdb  --- updated
            shield_a_stats.rdb  --- updated
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
    run_arc5gl(tstart, tstop, 'hrc5eng')
#
#--- merge all hrc4eng data into one fits file
#
    cmd  = 'ls hrcf*_5_eng0.fits.gz >dat.lis'
    os.system(cmd)

    cmd1 = "/usr/bin/env PERL5LIB="
    cmd2 = ' dmmerge infile=@dat.lis"[time='+ str(tstart) + ':' + str(tstop) + ',quality=0000000000000000000,mnf=0]" '
    cmd2 = cmd2 + 'outfile=hrc_5_eng0_'+ cyear + cmon + cday +'.fits mode=h'

    cmd  = cmd1 + cmd2
    bash(cmd,  env=ascdsenv)
#
#--- move the fits files to the saving direcotry
#
    cmd   = 'mv hrc_5_eng0*fits ' + outdir
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

def run_arc5gl(tstart, tstop, ftype):
    """
    run arc5gl
    input:  start   --- start time
            stop    --- stop time
            ftype   --- file type
    output: extracted fits files
    """

    fo    = open(zspace, 'w')

    line  = 'operation=retrieve\n'
    line  = line + 'tstart=' + str(tstart) + '\n'
    line  = line + 'tstop='  + str(tstop)  + '\n'
    line  = line + 'dataset=flight\n'
    line  = line + 'level=0\n'
    line  = line + 'detector=hrc\n'
    line  = line + 'subdetector=eng\n'
    line  = line + 'filetype=' + ftype + '\n'
    line  = line + 'go\n'

    fo.write(line)
    fo.close()

    cmd1 = "/usr/bin/env PERL5LIB="
    cmd2 = ' /proj/axaf/simul/bin/arc5gl -user isobe -script ' + zspace
    cmd  = cmd1 + cmd2
    bash(cmd,  env=ascdsenv)

#-----------------------------------------------------------------------------
#-- append_to_rdb: extract stat results and append the results to a rdb file -
#-----------------------------------------------------------------------------

def append_to_rdb(col, out_rdb, cyear, cmon, cday):
    """
    extract stat results and append the results to a rdb file
    input:  col     --- column name to be stat
            out_rdb --- output file name
            cyear   --- year
            cmon    --- month
            cday    --- day of the month
    output: out_rdb --- updated rdb file
    """
    ifits = 'hrc_rates_' + cyear + cmon + cday + '.fits'

    [nsamples, mean, sig, dmin, dmax] = get_stat(ifits, col)

    oline = cyear + cmon + cday + '\t' + nsamples + '\t' + mean + '\t' + sig + '\t' + dmin + '\t' + dmax + '\n'

    fo    = open(out_rdb, 'a')
    fo.write(oline)
    fo.close()

#-----------------------------------------------------------------------------
#-- get_stat: extract dmstat statistics                                    ---
#-----------------------------------------------------------------------------

def get_stat(ifits, col):
    """
    extract dmstat statistics
    input:  ifits       --- table fits file
            col         --- column name to be stat
    output: nsamples    --- numbers of good data points
            mean        --- mean
            sig         --- sigma
            dmin        --- min
            dmax        --- max
    """

    [dmin, dmax, mean, sig, med, nsamples] = fitsTableStat(ifits, col)
#
#--- adjust value to the older format. all values are in string format
#
    nsamples = str(int(nsamples))
    mean     = adjust_val(mean)
    sig      = adjust_val(sig)
    dmin     = adjust_val(dmin)
    dmax     = adjust_val(dmax)

    return [nsamples, mean, sig, dmin, dmax] 

#-----------------------------------------------------------------------------
#-- adjust_val: adjust value to match the format to the older data         ---
#-----------------------------------------------------------------------------

def adjust_val(val):
    """
    adjust value to match the format to the older data
    input:  val --- value
    output: val --- string form of the adjusted value
    """
    
    if val == 0:
        val = '0'
    elif int(val) == float(val):
        val = str(int(val))
    else:
        val = '%.11f' % (round(val, 11))

    return val

#-----------------------------------------------------------------------------
#-- fitsTableStat: find min, max, avg, std, med, and sample size of the colum of a table fits file
#-----------------------------------------------------------------------------

def fitsTableStat(infits, column, extension=1):

    """
    find min, max, avg, std, med, and sample size of the colum of a table fits file. 
    Input       infits  --- table fits file name
                column  --- name of the column(s). if there are more than one, must be 
                            in the form of list or tuple
                extension-- data extension #. default = 1
    Output      a list  [min, max, avg, std, med, sample size]
    """
    t     = pyfits.open(infits)
    tdata = t[extension].data
    t.close()

    data = tdata.field(column)

    results = []
    dmin = min(data)
    results.append(dmin)

    dmax = max(data)
    results.append(dmax)

    avg  = numpy.mean(data)
    results.append(avg)

    std  = numpy.std(data)
    results.append(std)

    med  = numpy.median(data)
    results.append(med)

    results.append(len(data))

    return results


#-----------------------------------------------------------------------------

if __name__ == '__main__':

    run_hrc_scripts()

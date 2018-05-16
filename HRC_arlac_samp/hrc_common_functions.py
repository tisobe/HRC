#!/usr/bin/env /proj/sot/ska/bin/python

#############################################################################################################
#                                                                                                           #
#           author: t. isobe (tisobe@cfa.harvard.edu)                                                       #
#                                                                                                           #
#           Last Update: May 15, 2017                                                                       #
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

ascdsenv = getenv('source /home/ascds/.ascrc -r release ;  source /home/mta/bin/reset_param ', shell='tcsh')
ciaoenv  = getenv('source /soft/ciao/bin/ciao.sh')


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


#-----------------------------------------------------------------------------------------
#-- find_fits_mean: find mean value of the column name given in the fits file           --
#-----------------------------------------------------------------------------------------

def find_fits_stats(fits, column, extension=1):
    """
    find mean value of the column name given in the fits file
    input:  fits        --- fits file name
            column      --- the column name
            extension   --- extension of the data block
    output: avg         --- mean value

    """

    t = pyfits.open(fits)
    tdata = t[extension].data
    t.close()

    data = tdata.field(column)
    avg  = numpy.mean(data)
    std  = numpy.std(data)
    med  = numpy.mediam(data)
    cnt  = len(data)

    return [avg, std, med, cnt]

#-----------------------------------------------------------------------------------------
#-- run_arc5gl: run arc5gl                                                              --
#-----------------------------------------------------------------------------------------

def run_arc5gl(start, stop, obsid = '', operation='retrieve', level ='1', filetype='evt1'):
    """
    run arc5gl 
    input:  start   --- starting time in the format of <yyyy>-<mm>-<dd>T<hh>:<mm>:<ss>
            stop    --- stoping time in the format of <yyyy>-<mm>-<dd>T<hh>:<mm>:<ss>
            obsid   --- obsid, if this is blank start and stop time will be used
            operation   --- retrieve or brwose  default: retrieve
            level       --- level               default: 1
            filetype    --- filetype            default: evt1
    output: the name of fits file or a list of the name of the fits files and the extracted file(s)
    """

    line = 'operation=' + operation + '\n'
    line = line + 'dataset=flight\n'
    line = line + 'detector=hrc\n'
    line = line + 'level='    + level    + '\n'
    line = line + 'filetype=' + filetype + '\n'

    if obsid == '':
        line = line + 'tstart=' + start  + '\n'
        line = line + 'tstop='  + stop   + '\n'
    else:
        line = line + 'obsid='  + str(obsid)  + '\n'
    line = line + 'go\n'

    fo   = open(zspace, 'w')
    fo.write(line)
    fo.close()

    cmd  = ' /proj/sot/ska/bin/arc5gl -user isobe -script ' + zspace +'> ztemp'
    os.system(cmd)

    mcf.rm_file(zspace)

    data = read_data('ztemp', remove=1)
    if len(data) > 2:
        fits = []
        for ent in data:
            mc = re.search('fits.gz', ent)
            if mc is not None:
                fits.append[ent]
    else:
        fits = ''
        for ent in data:
            mc = re.search('fits.gz', ent)
            if mc is not None:
                fits = ent

    return fits

#-----------------------------------------------------------------------------------------
#-- run_ascds: run the command in ascds environment                                     --
#-----------------------------------------------------------------------------------------

def run_ascds(cmd, clean =0):
    """
    run the command in ascds environment
    input:  cmd --- command line
    clean   --- if 1, it also resets parameters default: 0
    output: command results
    """
    if clean == 1:
        acmd = '/usr/bin/env PERL5LIB=""  source /home/mta/bin/reset_param ;' + cmd
    else:
        acmd = '/usr/bin/env PERL5LIB=""  ' + cmd
    
    try:
        bash(acmd, env=ascdsenv)
    except:
        try:
            bash(acmd, env=ascdsenv)
        except:
            pass

#-----------------------------------------------------------------------------------------
#-- run_ciao: running ciao comannds                                                    ---
#-----------------------------------------------------------------------------------------

def run_ciao(cmd, clean =0):
    """
    run the command in ciao environment
    input:  cmd --- command line
    clean   --- if 1, it also resets parameters default: 0
    output: command results
    """
    if clean == 1:
        acmd = '/usr/bin/env PERL5LIB=""  source /home/mta/bin/reset_param ;' + cmd
    else:
        acmd = '/usr/bin/env PERL5LIB="" LD_LIBRARY_PATH=""   ' + cmd
    
    try:
        bash(acmd, env=ciaoenv)
    except:
        try:
            bash(acmd, env=ciaoenv)
        except:
            pass


#-----------------------------------------------------------------------------------------
#--read_data: read data file                                                            --
#-----------------------------------------------------------------------------------------

def read_data(infile, remove=0):

    f    = open(infile, 'r')
    data = [line.strip() for line in f.readlines()]
    f.close()

    if remove == 1:
        mcf.rm_file(infile)

    return data

#-----------------------------------------------------------------------------------------
#-- get_data_from_db: extract observation information from the database                 --
#-----------------------------------------------------------------------------------------

def get_data_from_db(obsid):
    """
    extract observation information from the database
    input:  obsid   --- obsid
    output: tsec    --- the data of the observation in seconds from 1998.1.1
            line    --- a string of the information extract
                        <obsid> <target name> <obs date> <obs date in sec>
                        <target id> < sequence number>
    """

    try:
        dbase  = OcatDB(obsid)
        tname  = dbase.origValue('targname')
        target = clean_name(tname)              #--- limit target name to 14 characters
        inst   = dbase.origValue('instrument')
        odate  = dbase.origValue('soe_st_sched_date')
        tsec   = convert_date_to_sectime(odate)
        targid = dbase.origValue('targid')
        seqno  = dbase.origValue('seq_nbr')

        line   = str(obsid) + '\t' + target      + '\t' + str(odate) + '\t' + str(tsec) 
        line   = line       + '\t' + str(targid) + '\t' + str(seqno) + '\n'

        return [tsec, line]
    except:
        return NULL

#-----------------------------------------------------------------------------------------
#-- convert_date_to_sectime: convert time in <yyyy>-<mm>-<dd>T<hh>:<mm>:<ss> to seconds from 1998.1.1
#-----------------------------------------------------------------------------------------

def convert_date_to_sectime(odate):
    """
    convert time in <yyyy>-<mm>-<dd>T<hh>:<mm>:<ss> to seconds from 1998.1.1
    input:  odate   --- date in the format of <yyyy>-<mm>-<dd>T<hh>:<mm>:<ss>
    output: tsce    --- date in seconds from 1998.1.1
    """

    atemp  = re.split('\s+', str(odate))
    mon    = tcnv.changeMonthFormat(atemp[0])
    day    = int(float(atemp[1]))
    year   = int(float(atemp[2]))

    mc     = re.search('AM', atemp[3])
    if mc is not None:
        time  = atemp[3].replace('AM','')    
        btemp = re.split(':', time)
        hrs   = int(float(btemp[0]))
        mins  = int(float(btemp[1]))
    else:
        time  = atemp[3].replace('PM','')
        btemp = re.split(':', time)
        hrs   = int(float(btemp[0])) + 12
        mins  = int(float(btemp[1]))

    tsec  = tcnv.convertDateToTime2(year, mon, day, hours=hrs, minutes=mins)

    return tsec

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

def get_data(fits, col_list):

    hd    = pyfits.open(fits)
    tdata = hd[1].data

    save = []
    for col in col_list:
        exec "save.append(tdata['%s'])" % (col)

    return save

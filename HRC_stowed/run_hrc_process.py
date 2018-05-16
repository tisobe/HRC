#!/usr/bin/env /proj/sot/ska/bin/python

#############################################################################################
#                                                                                           #
#       run_hrc_process.py: create hrc evt1 file from evt0 file                             #
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
import random
import time
import operator
import math
import numpy
import astropy.io.fits  as pyfits
from datetime import datetime
import unittest

#
#--- from ska
#
from Ska.Shell import getenv, bash
ascdsenv = getenv('source /home/ascds/.ascrc -r release; source /home/mta/bin/reset_param', shell='tcsh')
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
import hrc_stowed_common_function as hcf        #---- contains hrc stored related functions
#
#--- temp writing file name
#
rtail  = int(time.time())
zspace = '/tmp/zspace' + str(rtail) + 'hrc'
#
#--- other settings
#
NULL      = 'null'

#-----------------------------------------------------------------------------------------------
#-- run_hrc_process: create hrc evt1 file from evt0 file                                      --
#-----------------------------------------------------------------------------------------------

def run_hrc_process(infits, inst, rsrfalv, wdthast):
    """
    create hrc evt1 file from evt0 file
    input:  infits  --- hrc evt0 fits file
            inst    --- instrument either hrc-i or hrc-s
            rsrfalv --- a value of estimated rsrfalv
            wdthast --- a value of estimated wdthast
    output: hrc evt1 fits file
    """
#
#--- find start time from the file
#
    start = hcf.read_header_value(infits, 'TSTART')
#
#--- set output file name
#
    outfits = infits.replace('evt0', 'evt1')
    outfits = outfits.replace('.gz', '')
#
#--- update header parameter values if needed
#
    add_header(infits, 'DETNAM',   value=inst,    datatype='string', comment='Detector')
    add_header(infits, 'WIDTHRES', value=wdthast, datatype='string', comment='tap-ringing correction')
    add_header(infits, 'RANGELEV', value=rsrfalv, datatype='short',  comment=' select between two possible calibration sets')

#
#--- set parameters depending on hrc i or hrc s
#
    mc   = re.search('hrc-i', inst)
    if mc is not None:
        paramlist = hrc_i_param(start)
    else:
        paramlist = hrc_s_param()
#
#--- set command and run hrc_process_events; first set ascds environment
#
    cmd1 = '/usr/bin/env PERL5LIB=""'
#
#--- process related environment settings
#
    cmd2 = ' punlearn ardlib;'
    cmd2 = cmd2 + ' pset hrc_process_events badfile=lev1_bad_evts.fits;'
    cmd2 = cmd2 + ' pset ardlib AXAF_HRC-I_BADPIX_FILE="' + house_keeping + 'hrci_bpix1.fits";'
    cmd2 = cmd2 + ' pset ardlib AXAF_HRC-S_BADPIX_FILE="' + house_keeping + 'hrcs_bpix1.fits";'
#
#--- actual command part
#
    cmd3 = " hrc_process_events "
    cmd3 = cmd3 + " infile="       + infits 
    cmd3 = cmd3 + " ampsatfile="   + paramlist[0]
    cmd3 = cmd3 + " ampsfcorfile=" + paramlist[1]
    cmd3 = cmd3 + " badpixfile="   + paramlist[2]
    cmd3 = cmd3 + " degapfile="    + paramlist[3]
    cmd3 = cmd3 + " evtflatfile="  + paramlist[4]
    cmd3 = cmd3 + " hypfile="      + paramlist[5]
    cmd3 = cmd3 + " obsfile="      + paramlist[6]
    cmd3 = cmd3 + " tapfile="      + paramlist[7]
    cmd3 = cmd3 + " gainfile="     + paramlist[8]
    cmd3 = cmd3 + " outfile="      + outfits
    cmd3 = cmd3 + " acaofffile=NONE"
    cmd3 = cmd3 + " clobber=yes"

    cmd  = cmd1 + cmd2 + cmd3

    try:
        bash(cmd, env=ascdsenv)
#
#--- zip the file
#
        cmd = 'gzip -f ' + outfits
        os.system(cmd)
#
#-- modify the output file name with "gz'
#
        name = outfits + '.gz'
    except:
        name = NULL

    return name

#-----------------------------------------------------------------------------------------------
#-- add_header: if a header parameter does not exist, add it                                 ---
#-----------------------------------------------------------------------------------------------

def add_header(infits, cname, value, datatype, comment):
    """
    if a header parameter does not exist, add it
    input:  infile      --- fits file name
            cname       --- parameter name
            value       --- parameter value
            datatype    --- data type
            comment     --- comment
    output: updated fits file
    """
    if hcf.read_header_value(infits, cname) == NULL:
        hcf.update_header(infits, cname, value=value, datatype=datatype, comment=comment)

#-----------------------------------------------------------------------------------------------
#-- hrc_i_param: set parameters for hrc i case                                               ---
#-----------------------------------------------------------------------------------------------

def hrc_i_param(start):
    """
    set parameters for hrc i case
    input:  start   --- start time in second from 1998.1.1
    output: [ampsatfile, ampsfcorfile, badpixfile, degapfile, evtflatfile, hypfile, obsfile, tapfile, gainfile]
    """
    ampsatfile   = calib_dir     + 'sattest/hrciD1999-07-22sattestN0002.fits'
    ampsfcorfile = calib_dir     + 'amp_sf_cor/hrciD1999-07-22amp_sf_corN0001.fits'
    badpixfile   = house_keeping + 'hrcf10702_000N001_bpix1.fits'
    degapfile    = calib_dir     + 'gaplookup/hrciD1999-07-22gaplookupN0004.fits'
    evtflatfile  = calib_dir     + 'eftest/hrciD1999-07-22eftestN0001.fits'
    hypfile      = calib_dir     + 'fptest/hrciD1999-07-22fptestN0003.fits'
    tapfile      = calib_dir     + 'tapringtest/hrciD1999-07-22tapringN0002.fits'
    obsfile      = house_keeping + 'obs.par'
#
#--- read gain selection list
#
    infile = house_keeping + 'Gain_files/gain_selection'
    data   = hcf.read_file_data(infile)
#
#--- select gain file name 
#
    for ent in data:
        atemp = re.split('\s+', ent)
        begin = int(float(atemp[0]))
        end   = int(float(atemp[1]))

        if start >=  begin and start < end:
            gainfile = calib_dir + 'gmap/' + atemp[2]
            break

    return [ampsatfile, ampsfcorfile, badpixfile, degapfile, evtflatfile, hypfile, obsfile, tapfile, gainfile]

#-----------------------------------------------------------------------------------------------
#-- hrc_s_param: set parameters for hrc s case                                               ---
#-----------------------------------------------------------------------------------------------

def hrc_s_param():
    """
    set parameters for hrc s case
    input:  none
    output: [ampsatfile, ampsfcorfile, badpixfile, degapfile, evtflatfile, hypfile, obsfile, tapfile, gainfile]
    """
    ampsatfile   = calib_dir     + 'sattest/hrcsD1999-07-22sattestN0002.fits'
    ampsfcorfile = calib_dir     + 'amp_sf_cor/hrcsD1999-07-22amp_sf_corN0001.fits'
    badpixfile   = house_keeping + 'hrcf10600_000N001_bpix1.fits'
    degapfile    = calib_dir     + 'gaplookup/hrcsD1999-07-22gaplookupN0003.fits'
    evtflatfile  = calib_dir     + 'eftest/hrcsD1999-07-22eftestN0001.fits'
    hypfile      = calib_dir     + 'fptest/hrcsD1999-07-22fptestN0004.fits'
    tapfile      = calib_dir     + 'tapringtest/hrcsD1999-07-22tapringN0002.fits'
    gainfile     = calib_dir     + 'gmap/hrcsD1999-07-22gainN0001.fits'
    obsfile      = house_keeping + 'obs_hrc_s.par'
    

    return [ampsatfile, ampsfcorfile, badpixfile, degapfile, evtflatfile, hypfile, obsfile, tapfile, gainfile]

#-----------------------------------------------------------------------------------------------
#-- update_gain_selection_file: check whether a new gain file is added. if so, update gain_selection file
#-----------------------------------------------------------------------------------------------

def update_gain_selection_file(): 
    """
    check whether a new gain file is added. if so, update gain_selection file
    input:  none, but read from <calib_dir>
    output: updated <house_keeping>/Gain_files/gain_selection
            also a copy of the gain file
    """
#
#--- read gain selection list
#
    infile = house_keeping + 'Gain_files/gain_selection'
    data   = hcf.read_file_data(infile)
#
#--- read gain file list from <calib_dir>
#
    gdata = get_gain_file_list()
#
#--- check whether a new file is added. if so go farther
#
    if len(gdata) > len(data):
        [dfile, dtime] = find_latest_file(data)
        [gfile, gtime] = find_latest_file(gdata)

        if gtime > dtime:
            fo  = open(infile, 'w')
            for ent in data[:-1]:
                fo.write(ent)
                fo.write('\n')
#
#--- update the last line
#
            atemp = re.split('\s+', data[-1])
            line  = atemp[0] + '\t' +  str(gtime)+ '\t' + atemp[2] + '\n'
            fo.write(line)
#
#--- and add the new line
#
            line  = str(gtime) + "\t1.0e12\t\t" + gfile + "\n"
            fo.write(line)

            fo.close()
#
#--- copy the file into Gain_files directory
#
            cmd = 'cp ' + calib_dir + 'gmap/' + gfile + ' ' + house_keeping + 'Gain_files/. '
            os.system(cmd)
            cmd = 'gzip -df ' + house_keeping + 'Gain_files/*.gz'
            os.system(cmd)

#-----------------------------------------------------------------------------------------------
#-- get_gain_file_list: read gain file list from <calib_dir>                                  --
#-----------------------------------------------------------------------------------------------

def get_gain_file_list():
    """
    read gain file list from <calib_dir>
    input:  none
    output: gdata   --- a list of gain files
    """

    cmd   = 'ls ' + calib_dir + 'gmap/hrci*sampgain*fits > ' + zspace
    os.system(cmd)
    adata = hcf.read_file_data(zspace, remove=1)
#
#--- get only file name
#
    gdata = []
    for ent in adata:
        atemp = re.split('hrciD', ent)
        name  = 'hrciD' + atemp[1]
        gdata.append(name)

    return gdata

#-----------------------------------------------------------------------------------------------
#-- find_latest_file: find the file name with the latest time stamp                           --
#-----------------------------------------------------------------------------------------------

def find_latest_file(dlist, precut='hrciD', postcut='sampgain'):
    """
    find the file name with the latest time stamp
        note: this works only on a gain file name as a default
    input:  dlist   --- a list of file name. it can be either a list of lists of
                        [<start time>, <stop time>, <file name>] of simply a list
                        of file names
            precut  --- a word preceeding the time stamp. default 'hrciD' for a gain file name
            postcut --- a word following the time stamp.  default 'sampgain' for a gain file name
    output: ofile   --- the name of file name with the latest time stamp
            stime   --- the time in seconds from 1998.1.1
    """
#
#-- there are two possible input; make sure to have a list of file names
#
    test = re.split('\s+', dlist[0])

    if len(test) > 1:
        tlist = []
        for ent in dlist:
            atemp = re.split('\s+', ent)
            tlist.append(atemp[2])
    else:
        tlist = dlist
#
#--- extract time part from the file name
#
    clist = []
    for ent in tlist:
        atemp = re.split(precut,  ent)
        btemp = re.split(postcut, atemp[1])
        clist.append(btemp[0])
#
#--- find the latest time
#
    clist.sort()
    cval = clist[-1]
#
#--- get the file name corresponding to the time
#
    for ent in tlist:
        mc = re.search(cval, ent)
        if mc is not None:
            ofile = ent
            break
#
#--- convert time into second from 1998.1.1
#
    cval  = cval + ',00:00:00'
    stime = hcf.convertto1998sec(cval)

    return [ofile, stime]

#-----------------------------------------------------------------------------------------
#-- TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST    ---
#-----------------------------------------------------------------------------------------

class TestFunctions(unittest.TestCase):
    """
    testing functions
    """

#------------------------------------------------------------

    def test_find_latest_file(self):

        infile = house_keeping + '/Gain_files/gain_selection'
        dlist  = hcf.read_file_data(infile)

        [ofile, stime] = find_latest_file(dlist)

        print "FILE NAME: " + ofile + '<--->' + str(stime)


        update_gain_selection_file()

        glist = get_gain_file_list()

        self.assertEquals(glist[-1], ofile)

#------------------------------------------------------------
    
    def test_hrc_i_param(self):

        start = 65750399                #--- Feb 1, 2000
        plist = hrc_i_param(start)

        for ent in plist:
            print str(ent)

        self.assertEquals(plist[0],  '/data/CALDB/sdp/data/chandra/hrc/sattest/hrciD1999-07-22sattestN0002.fits')

        self.assertEquals(plist[-1], '/data/CALDB/sdp/data/chandra/hrc/gmap/hrciD1999-10-04sampgainN0001.fits')

#------------------------------------------------------------

    def test_run_hrc_process(self):

        infits  = house_keeping + '/Test_data/hrcf567231386N001_evt0.fits.gz'
        outfits = house_keeping + '/Test_data/hrcf567231386N001_evt1.fits.gz'
        inst    = 'hrc-i'
        rsrfalv = 116
        wdthast = 2

        run_hrc_process(infits, inst, rsrfalv, wdthast)

        [cols, tbdata] = hcf.read_fits_file(outfits)

        crsv = list(tbdata['crsv'])
        self.assertEquals(crsv[:10], [13, 18, 25, 22, 58, 4, 21, 41, 30, 57])

        mcf.rm_file(outfits)

#-----------------------------------------------------------------------------------------------

#
if __name__ == "__main__":

    unittest.main()



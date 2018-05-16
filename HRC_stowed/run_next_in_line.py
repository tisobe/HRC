#!/usr/bin/env /proj/sot/ska/bin/python

#############################################################################################
#                                                                                           #
#       run_next_in_line.py:  extract hrc stowed event files and create evt 1 files         #
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
import run_hrc_process            as rhp        #---- contains run_hrc_process related function
#
#--- temp writing file name
#
rtail  = int(time.time())
zspace = '/tmp/zspace' + str(rtail) + 'hrc'
#
#--- hrc hk column list
#
hk_col_list = ['TIME','LLDIALV','RSRFALV','WDTHAST','ULDIALV','CBHUAST','CBHVAST','CBLUAST']
hk_col_list = hk_col_list + ['CBLVAST','IMTPAST','IMBPAST','IMHVLV','IMHBLV','SMTRATM','HVPSSTAT']
hk_col_list = hk_col_list + ['SPTPAST','SPBPAST','SPHVLV','SPHBLV','SCIDPREN', 'S2HVST', 'S2HVLV','QUALITY']
#
#--- set  a few things
#
NULL      = 'null'

#-----------------------------------------------------------------------------------------------
#-- extract_next_in_line_data: extract hrc stowed event files and create evt 1 files         ---
#-----------------------------------------------------------------------------------------------

def extract_next_in_line_data(begin, end, start, stop, outdir):
    """
    extract hrc stowed event files and create evt 1 files
    input:  begin   --- starting time in mm/dd/yy,hh:mm:ss
            end     --- ending time in mm/dd/yy,hh:mm:ss
            start   --- starting time in seconds from 1998.1.1
            stop    --- stopping time in seconds from 1998.1.1
            outdir  --- main output directory 
    output: evt0, evt1, ss0, and hk00 fits files corresponding to next_in_line condition
    """
    print "TIME: " + str(begin) +'<--->'+ str(end)
#
#--- select time spans set by SIM movement
#
    [start_list, stop_list]  = find_tscpos_positive_period(begin, end)
#
#--- keep the record
#
    fo = open('tscpos_positive', 'w')
    for m in range(0, len(start_list)):
        line = str(start_list[m]) +'\t' + str(stop_list[m])      + '\t'  
        line = line + str(hcf.covertfrom1998sec(start_list[m]))  + '\t' 
        line = line + str(hcf.covertfrom1998sec(stop_list[m]))   + '\n'
        fo.write(line)
    fo.close()
#
#--- check the next conditions for each SIM positive time span
#
    for k in range(0, len(start_list)):
#
#--- create hrc hk base fits files based on conditions. it creates hrc_i_115.fits etc
#--- which contains data appropriate for the conditions.
#
        comb_fits_list = create_hk_fits_with_extracted_data(start_list[k], stop_list[k])
#
#--- set a couple of parameters to be used. afits has a form of "hrc_i_115.fits" etc
#
        hk_chk = []
        en_chk = [[],[]]
        for afits in comb_fits_list:

            mc = re.search('hrc_i', afits)
            if mc is not None:
                hrc_type = 'IMAG'
                inst     = 'hrc-i'
            else:
                hrc_type = 'SPEC'
                inst     = 'hrc-s'
#
#--- set time span
#
            test_list  = create_data_period(start_list[k], stop_list[k], [afits])
            time_list  = []
            for ent in test_list:
                chk = 0
                for comp in hk_chk:
                    if ent == comp:
                        chk = 1
                        break
                if chk == 0:
                    hk_chk.append(ent)
                    time_list.append(ent)
            if len(time_list) == 0:
                continue
#
#--- keep the record
#
            tname = str(afits.replace('.fits',''))
            oname = tname + '_hk_time_period'
            fo = open(oname, 'a')
            for m in range(0, len(time_list[0])):
                line = str(time_list[0][m]) +'\t' + str(time_list[1][m])  + '\t'  
                line = line + str(hcf.covertfrom1998sec(time_list[0][m])) + '\t'  
                line = line + str(hcf.covertfrom1998sec(time_list[1][m])) + '\n'
                fo.write(line)
            fo.close()
#
#--- farther selet time span set by engineer data
#
            test_list  = get_time_period_engineer_data(time_list, hrc_type, start_list[k], stop_list[k])
            etime_list = [[],[]]
            for k in range(0, len(test_list[0])):
                chk = 0
                for comp in en_chk[0]:
                    if test_list[0][k] == comp:
                        chk = 1
                        break

                if chk == 0:
                    en_chk[0].append(test_list[0][k])
                    en_chk[1].append(test_list[1][k])
                    etime_list[0].append(test_list[0][k])
                    etime_list[1].append(test_list[1][k])
#
#--- keep the record
#
            oname = tname + '_eng_time_period'
            fo = open(oname, 'a')
            for m in range(0, len(etime_list[0])):
                line = str(etime_list[0][m]) +'\t' + str(etime_list[1][m]) + '\t'  
                line = line + str(hcf.covertfrom1998sec(etime_list[0][m])) + '\t'  
                line = line + str(hcf.covertfrom1998sec(etime_list[1][m])) + '\n'
                fo.write(line)
            fo.close()
#
#--- extract evt0, scinece and hrc hk data for given time spans
#--- science and hk are combined ones for the given time period by evt0
#
            [evt0_list, sci_list, hk_list] = extract_hrc_fits_files(etime_list)
#
#--- process evt0 file
#
            evt1_list = []
            for k in range(0, len(evt0_list)):
#
#--- estimate rsrfalv and wdthast from hk_file correspond to the evt0 fits file
#
                [rsrfalv, wdthast] = estimate_rsrfalv_wdthast(hk_list[k])
#
#--- evt1 could be NULL if the process failed
#
                evt1 = rhp.run_hrc_process(evt0_list[k], inst, rsrfalv, wdthast)
                if evt1 != NULL:
                    evt1_list.append(evt1)
#
#--- combine evt1 if their time difference is small (<180 sec)
#
            comb_list = combine_short_exposure_files(evt1_list)
#
#--- re-extract  hk and science supplemental fits files for the time period
#
            cmd = 'rm hrcf*hk0*fits* hrcf*ss0*fits* 2>/dev/null'
            os.system(cmd)

            dname = afits.replace('.fits', '')
            rdir  = outdir +  dname + '/'
            for efits in comb_list:
                atemp = re.split('_evt1', efits)
                head  = atemp[0]

                [fstart, fstop] = find_start_and_end(efits)
                if fstart == 0:
                    continue

                [ss_name, hk_name] = extract_supple_fits_files([[fstart], [fstop]], head)
#
#--- move the file to an appropriate directory
#
                ss_name2 = ss_name + '*'
                hk_name2 = hk_name + '*'
                #cmd   = 'mv ' + efits + '  ' + ss_name2 + ' ' + hk_name2 + ' '  + rdir 
                cmd   = 'mv ' + efits + '  ' + ss_name + ' ' + hk_name + ' '  + rdir + ' 2> /dev/null' 
                os.system(cmd)
#
#--- move all evt files into the directory, too
#
            tdir  = exc_dir + 'Temp_dir/'
            e0list= hcf.get_file_list(tdir, 'evt0.fits.gz')
            for e0 in e0list:
                cmd   = 'mv ' + e0 + ' ' + rdir
                os.system(cmd)

            cmd = 'gzip -f ' + tdir + '*evt1.fits 2> /dev/null'
            os.system(cmd)

            cmd = 'mv *evt1.fits ' + rdir + ' 2> /dev/null'
            os.system(cmd)

#-----------------------------------------------------------------------------------------------
#-- estimate_rsrfalv_wdthast: estimate RSRFALV and WDTHAST value from a hk fits file          --
#-----------------------------------------------------------------------------------------------

def estimate_rsrfalv_wdthast(fname):
    """
    estimate RSRFALV and WDTHAST value from a hk fits file
    input:  fname   --- hk fits file
    output: rval    --- estimated RSRFALV
            wval    --- esimtated WDTHAST
    """

    [cols, tbdata] = hcf.read_fits_file(fname)

    rdata = tbdata['rsrfalv']
    if len(rdata) > 0:
        ravg = rdata.mean()
        rval = int(2.0 * (ravg - 127) + 0.5)
    else:
        rval = 116

    wdata = tbdata['wdthast']
    if len(wdata) > 0:
        wval = int(wdata.mean())
    else:
        wval = 2

    return [rval, wval]

#-----------------------------------------------------------------------------------------------
#-- extract_hrc_fits_files: extract hrc evt0, secondary sci, hrchk files for given time period -
#-----------------------------------------------------------------------------------------------

def extract_hrc_fits_files(time_list):
    """
    extract hrc evt0, secondary sci, hrchk files for given time period
    input:  time_list   --- a list of lists of start and stop time
    output: hrc evt0 files
            hrc secondary science fits files
            hrck fits files
            evt0_list   --- a list of extracted hrc evt0 file names
            sci_list    --- a list of extracted secondary sci fits files
            hk_list     --- a list of extracted hrchk fits files
    """
#
#--- extreact hrc evt0 files
#
    sci_list = []
    hk_list  = []
    evt0_list = extract_fits_files(time_list, 'hrc', '0', 'evt0')

    if len(evt0_list) == 0:
        return [[], [], []]
#
#--- make sure that the files extracted are in the given time periods
#
    evt0_cleaned = []
    for evt0_file in evt0_list:
        #tstart    = hcf.read_header_value(evt0_file, 'tstart')
        #tstop     = hcf.read_header_value(evt0_file, 'tstop')
        [cols, fdata] = hcf.read_fits_file(evt0_file)
        tdata = fdata['time']
        tstart = tdata.min()
        tstop  = tdata.max()

        chk = 0
        for k in range(0, len(time_list[0])):

            cstart = time_list[0][k]
            cstop  = time_list[1][k]

            if ((tstart >= cstart) and (tstop <= cstop)):
                evt0_cleaned.append(evt0_file)
                chk = 1
                break 

        if chk == 0:
            continue
#
#--- extract corresponding ss and hk files
#
        atemp     = re.split('hrcf', evt0_file)
        btemp     = re.split('_', atemp[1])
        head      = 'hrcf' + btemp[0]

        [ss_name, hk_name] = extract_supple_fits_files([[cstart], [cstop]], head)
        sci_list.append(ss_name)
        hk_list.append(hk_name)

    return [evt0_cleaned, sci_list, hk_list]

#-----------------------------------------------------------------------------------------------
#-- extract_supple_fits_files: extract science and housekeeping fits files for given time period 
#-----------------------------------------------------------------------------------------------

def extract_supple_fits_files(time_list, head):
    """
    extract science and housekeeping fits files for given time period
    input:  time_list   --- a list of lists of start and stop time
            head        --- prefix of the fits files
    output: combined fits files of hrc ss0 and hk0
    """
#
#--- extreact secondary science files
#
    ss_name = head + '_comb_ss0.fits'
    extract_data_and_combine(time_list, 'hrc', '0', 'hrcss', ss_name)
#
#--- extreact hrc hk files
#
    hk_name = head + '_comb_hk0.fits'
#
#--- convert col list to string
#
    fcol_line = hk_col_list[0]
    for k in range(1, len(hk_col_list)):
        fcol_line = fcol_line + ', ' + str(hk_col_list[k])

    extract_data_and_combine(time_list, 'hrc', '0', 'hrchk', hk_name, fcol_line=fcol_line)

    return [ss_name, hk_name]

#-----------------------------------------------------------------------------------------------
#-- extract_data_and_combine: extract fits files in the time span(s) and combine all of them to one fits file 
#-----------------------------------------------------------------------------------------------

def extract_data_and_combine(time_list, detector, level, filetype, name, fcol_line=''):
    """
    extract fits files in the time span(s) and combine all of them to one fits file
    input:  time_list   --- a list of lists [[start_time>],[<stop_time>]]
            detector    --- detector name (e.g. hrc)
            level       --- level (e.g. 0, 1)
            filetype    --- file type (e.g. hrcss, hrchk)
            name        --- name of the resulted combined fits file
    output: name        --- resulted fits file
    """
#
#--- extreact files
#
    fits_list = extract_fits_files(time_list, detector, level, filetype)
#
#--- combine fits file to one
#
    hcf.combine_fits_files(fits_list, name, azip=0, fcol_line=fcol_line)
#
#--- remove indivisual fits files
#
    for ent in fits_list:
        mcf.rm_file(ent)

#-----------------------------------------------------------------------------------------------
#-- get_time_period_engineer_data: extract time period in which hrc engireer data gives  2preads=hrc type
#-----------------------------------------------------------------------------------------------

def get_time_period_engineer_data(time_list, hrc_type, start, stop):
    """
    extract time period in which hrc engireer data gives  2preads=hrc type
    input:  time_list   --- a list of lists of starting and stoppping time
                            from previous condition
            hrc_type    --- hrc type either IMG or SPEC
            start       --- start time of the entire period
            stop        --- stop time of the entire period
    output: rlist       --- a list of lists of starting and stoppping time
    """
#
#--- extract hrc engineer fits files
#
    fits_list = extract_fits_files(time_list, 'hrc', '0', 'hrc2eng', subdetector='eng')
#
#--- extract data needed
#
    dout = hcf.combine_and_select_data(fits_list, ['time','2PREADS'])
    if dout == NULL:
        return [[], []]
#
#--- select time with 2preads == hrc_type (either hrc i or hrc s coditino)
#
    odata = hcf.select_data_with_condition(dout ,'2PREADS', '==', hrc_type)
#
#--- get lists of start and stop time periods
#
    tlist = list(odata.field('TIME'))
    rlist = find_time_span(tlist, start, stop) 
#
#--- remove indivisual fits files
#
    for ent in fits_list:
        mcf.rm_file(ent)

    return rlist

#-----------------------------------------------------------------------------------------------
#-- extract_fits_files: extract fits file for given time period                              ---
#-----------------------------------------------------------------------------------------------

def extract_fits_files(time_list, detector, level, filetype, subdetector=''):
    """
    extract fits file for given time period
    input:  time_list   --- a list of lists of start and stop time
            detector    ---  detector name
            level       ---  level
            filetype    --- file type
            subdetector --- sub detector name
    output: extracted fits files in <exc_dir>/Temp_dir/
            fits_list   --- a list of extracted fits file names
    """

    if len(time_list) < 2:
        return []

    [start_list, stop_list] = time_list
#
#--- just in a case, start and stop time are not given in a list form 
#--- make them into a list form
#
    if not isinstance(start_list, list):
        start_list = [start_list]

    if not isinstance(stop_list, list):
        stop_list = [stop_list]

    fits_list = []                  #--- just in a case nothing extracted
    for m in range(0, len(start_list)):
        try:
            start = float(start_list[m])
            stop  = float(stop_list[m])
        except:
            break
#
#--- extract fits files
#
        dlist = hcf.run_arc5gl('retrieve', start, stop, detector=detector, level=level, filetype=filetype, subdetector=subdetector)
        if m == 0:
            fits_list = dlist
        else:
            fits_list = fits_list + dlist
#
#--- since there are occasional overlap in periods, make sure that all fits names are uniqe
#
    fits_list = hcf.remove_duplicate_from_list(fits_list)
    return fits_list

#-----------------------------------------------------------------------------------------------
#-- find_tscpos_positive_period: find time periods of  tscpos is located in positive position --
#-----------------------------------------------------------------------------------------------

def find_tscpos_positive_period(start, stop):
    """
    find time periods of  tscpos is located in positive position
    input:  start       --- starting time in mm/dd/yy,hh:mm:ss
            stop        --- stopping time in mm/dd/yy,hh:mm:ss
    output: positive_period a list of lists containing a list of time of period starting
            and a list of time of period ending. both in seconds from 1998.1.1
    """
#
#--- extract sim data
#
    dlist = hcf.run_arc5gl('retrieve', start, stop, detector='sim', level='0', filetype='sim')
#
#--- find time period when tscpos is in positive side
#
    positive_periods = create_data_period(start, stop, dlist, colname='TSCPOS', lgc='>=', val=0)
#
#-- clean the output directory
#
    for ent in dlist:
        mcf.rm_file(ent)

    return positive_periods

#-----------------------------------------------------------------------------------------------
#-- create_data_period: create a pair of time lists which indicate when tscpos is in positive side 
#-----------------------------------------------------------------------------------------------

def create_data_period(start, stop, dlist, colname='', lgc='', val=''):
    """
    create a pair of time lists which indicate when tscpos is in positive side
    input:  start       --- starting time in mm/dd/yy,hh:mm:ss or seconds from 1998.1.1
            stop        --- stopping time in mm/dd/yy,hh:mm:ss or seconds from 1998.1.1
            dlist       --- a list of fits file names
            colname     --- if given, the following two must be given. default: ''
            lgc         --- logical indicator such as "<=", "==". not effective unless colname is given
            val         --- value of the condition. not effective unless colname is given
    output: start_list  --- a list of starting time in format of seconds from 1998.1.1
            stop_list   --- a list of stopping tine in format of seconds from 1998.1.1
            
    """
#
#--- check time format
#
    if (not isinstance(start, float)) and (not isinstance(start, int)):
        start = hcf.convertto1998sec(start)
        stop  = hcf.convertto1998sec(stop)

    time_list = []
    chk = 0
    for fits in dlist:
#
#--- extract time and tscpos information from fits file
#
        [cols_in, hdata] = hcf.read_fits_file(fits)
#
#--- if an extra condition is given, run the next
#
        if colname != '':
            hdata = hcf.select_data_with_condition(hdata, colname, lgc, val)
        else:
            pass

        tlist = list(hdata.field('TIME'))

        if chk == 0:
            time_list =  tlist
            chk = 1
        else: 
            time_list = time_list + tlist
#
#--- create two lists: starting time and ending time
#
    time_list = sort_and_clean(time_list)

    t_lists   = find_time_span(time_list, start, stop)

    return t_lists

#-----------------------------------------------------------------------------------------------
#-- sort_and_clean: sort a list and remove duplicated entries numeric value only              --
#-----------------------------------------------------------------------------------------------

def sort_and_clean(tlist):
    """
    sort a list and remove duplicated entries, numeric value only
    input:  tlist   --- a list
    output: asave   --- the cleaned list
    """
    if len(tlist) == 0:
        return []

    atemp = []
    for ent in tlist:
        atemp.append(int(float(ent)))
    atemp.sort()
    prev  = atemp[0]
    asave = [prev]

    for k in range(1, len(atemp)):
        if atemp[k] == prev:
            continue
        else:
            asave.append(atemp[k])
            prev = atemp[k]

    return asave

#-----------------------------------------------------------------------------------------------
#-- find_time_span: find start and stop time of the period                                   ---
#-----------------------------------------------------------------------------------------------

def find_time_span(time_list, start, stop):
    """
    find start and stop time of the period. 
    input:  time_list       --- a list of time in seconds from 1998.1.1
            start           --- a start time of the entire period
            stop            --- a stop time of the entire period
    output: start_list      --- a list of starting time
            stop_list       --- a list of stopping time
    """
#
#--- if the list is empty or only one entry, just return empty lists
#
    tlen = len(time_list)

    if tlen < 2:
        return [[], []]

    elif tlen == 2:
        return [[time_list[0]], [time_list[1]]]

    else:
        time_list.sort()
        temp = []
        for ent in time_list:
            temp.append(int(float(ent)))
        time_list = temp

    prev       = time_list[0]
    start_list = [prev]
    stop_list  = []

    for i in range(1, tlen):
#
#--- if there is more than 5 min break, a new period start
#
        diff = time_list[i] - prev

        if diff > 300:
            stop_list.append(prev)
            start_list.append(time_list[i])

        prev = time_list[i]
#
#--- the last period ends at "stop" or time_list[-1] time
#
    if len(start_list) > len(stop_list):

        if time_list[-1] < stop:
            stop_list.append(time_list[-1])
        else:
            stop_list.append(stop)
#
#--- make sure that they are in the boundary
#
    [start_list, stop_list] = make_clean_cut(start_list, stop_list, start, stop)

    return [start_list, stop_list]

#-----------------------------------------------------------------------------------------------
#-- make_clean_cut: check the start and stop list are inside of the hard limits               --
#-----------------------------------------------------------------------------------------------

def make_clean_cut(start_list, stop_list, start, stop):
    """
    check the start and stop list are inside of the hard limits
    input:  start_list  --- a list of starting time
            stop_list   --- a list of stopping time
            start       --- hard starting limit
            stop        --- hard stopping limit
    output: start_list2 --- a list of start time
            stop_list2  --- a list of stopping time
    """

    start_list2 = []
    stop_list2  = []

    for k in range(0, len(start_list)):
#
#--- if start and stop time is same, drop from the list
#
        if start_list[k] == stop_list[k]:
            continue 
#
#--- starting point check
#
        if start_list[k] < start:
            if stop_list[k] > start:
                start_list2.append(start)
                stop_list2.append(stop_list[k])
            else:
                continue
#
#--- ending point check
#
        elif stop_list[k] >= stop:
            if start_list[k] < stop:
                start_list2.append(start_list[k])
                stop_list2.append(stop)
            break
#
#--- middle
#
        else:
            start_list2.append(start_list[k])
            stop_list2.append(stop_list[k])

    return [start_list2, stop_list2]

#-----------------------------------------------------------------------------------------------
#-- create_hk_fits_with_extracted_data: create fits files from extracted data for the given period 
#-----------------------------------------------------------------------------------------------

def create_hk_fits_with_extracted_data(start, stop):
    """
    create fits files from extracted data for the given period
    input:  start       --- starting time in the format in seconds from 1998.1.1
            stop        --- stoppping time in the format in  seconds from 1998.1.1
            conditions are read from <house_keeping><hrc_condition_file>
    output: fits_name   --- fits file (e.g.hrc_i_115.fits)
    """
    fits_list = []
    fdata  = extrct_hrchk_data(start, stop)
    if fdata == NULL:
        return fits_list
#
#--- the condition differs before and after 2012 Apr.
#
    if stop <= 449452797:
        hlist = ('hrc_i_90','hrc_i_115', 'hrc_s_125_1', 'hrc_s_125_hi_1', 'hrc_s_90_1', 'hrc_s_90_hi_1')
    else:
        hlist = ('hrc_i_90','hrc_i_115', 'hrc_s_125_2', 'hrc_s_125_hi_2', 'hrc_s_90_2', 'hrc_s_90_hi_2')

    for cfile in hlist:
#
#--- read the condition
#
        [select_list, fits_name] = read_condition(cfile)
#
#--- extract data meet the condition
#
        data  = extract_data_under_condition(fdata, select_list)
#
#--- create fits file
#
        if data != NULL:
            hcf.create_fits_file(data, hk_col_list, fits_name)
            fits_list.append(fits_name)

    return fits_list

#-----------------------------------------------------------------------------------------------
#-- read_condition: read data extraction condition file and also creates output fits file name -
#-----------------------------------------------------------------------------------------------

def read_condition(cfile):
    """
    read data extraction condition file and also creates output fits file name
    input:  cfile       --- condition file name
    output: condition   --- a list of lists containing column name and the value range
            fits        --- output fits file name
    """
#
#--- read a condition file
#
    ifile = house_keeping + 'Selection_coditions/' +  cfile
    data  = hcf.read_file_data(ifile)

    condition = []
    for ent in data:
        if ent[0] == '#':
            continue

        atemp = re.split('=', ent)
        condition.append(atemp)
#
#--- create output fits file name
#
    test  = str(cfile)
    test2 = test[-2:]           #--- checking the last two character

    if test2 == '_1' or test2 == '_2':
        test = test[:-2]

    fits = test + '.fits'

    return [condition, fits]

#-----------------------------------------------------------------------------------------------
#-- extract_data_under_condition: extract specified data from the hrc hk data set            ---
#-----------------------------------------------------------------------------------------------

def extract_data_under_condition(hdata, select_list):
    """
    extract specified data from the hrc hk data set
    input:  hdata       --- hrc hk data processed by extrct_hrchk_data
            select_list --- a list contains lists of [<column name>, <condition>]
                            used to select data
    output: hdata       --- selected data
    """
#
#--- the list select_list contains lists of [<name>, <condition>]
#
    for ent in select_list:
        name = ent[0]
        mc   = re.search(':', ent[1])
#
#-- numerical range case; need start and stop range
#
        if mc is not None:
            atemp = re.split(':', ent[1])
            cstart = int(float(atemp[0]))
            cstop  = int(float(atemp[1]))

            hdata  = hcf.select_data_with_condition(hdata, name, '>=', cstart)
            hdata  = hcf.select_data_with_condition(hdata, name, '<=', cstop)
#
#--- single condition case
#
        else:
            cond   = ent[1]
            hdata  = hcf.select_data_with_condition(hdata, name, '==', cond)

#
#--- check whether any data left in the datatable
#
        try:
            test = hdata.field('time')
            if len(test) == 0:
                return NULL
        except:
            return NULL
#
#--- check one more last time whether any data left in the datatable
#
    try:
        test = hdata.field(0)
    except:
        test = []

    if len(list(test)) > 0:
        return hdata
    else:
        return NULL

#-----------------------------------------------------------------------------------------------
#-- extrct_hrchk_data: extract hrchk fits file and get column data we need                    --
#-----------------------------------------------------------------------------------------------

def extrct_hrchk_data(start, stop):
    """
    extract hrchk fits file and get column data we need
    input:  start   --- starting time
            stop    --- stopping time
    output: out     --- pyfits data table format of extracted data from data with datamore: NEXT_IN_LINE 
    """
#
#--- extract hrchk fits files 
#
    dlist = hcf.run_arc5gl('retrieve', start, stop, dataset='flight', detector='hrc', level='0', filetype='hrchk')

    nil_list = []
    for fits in dlist:
#
#--- find datamode: we need only NEXT_IN_LINE case
#
        dmode = hcf.read_header_value(fits, 'DATAMODE')

        if dmode == 'NEXT_IN_LINE':
            nil_list.append(fits)

    if len(nil_list) > 0:
        dout = hcf.combine_and_select_data(nil_list, hk_col_list)
    else:
        dout = NULL
#
#--- empty out the temporary direcotry
#
    for ent in dlist:
        mcf.rm_file(ent)

    return dout

#-----------------------------------------------------------------------------------------------
#-- set_data_period: create a list of dates to be examined                                   ---
#-----------------------------------------------------------------------------------------------

def set_data_period(syear, smonth, eyear, emonth):
    """
    create a list of dates to be examined
    input:  syear   --- year of starting time
            smonth  --- month of starting time
            eyear   --- year of ending time
            emonth  --- month of ending time
    output: start_time  --- a list of starting time in form of mm/01/yy,00:00:00
            end_time    --- a list of ending time in format of mm/01/yy,00:00:00
                            the end time is set to the month after start_time
                            e.g if start_time is 05/01/14,00:00:00, endtime
                                is 06/01/14,00:00:00
    """

    start_time = []
    end_time   = []
#
#--- if the time period is not given, find starting and the ending time
#
    if syear == '':
        try:
            (syear, smonth, eyear, emonth) = find_time_period()
        except:
            exit(1)             #--- it seems that data are already up to date
#
#-- for the case the period is complete in the same year
#
    if syear == eyear:
        for month in range(smonth, emonth+1):
            start_time.append(hcf.conv_time_format(syear, month))
            end_time.append(hcf.conv_time_format(syear, month, next=1))
#
#--- for the case the period goes over more than one year
#
    elif syear < eyear:
        for lyear in range(syear, eyear+1):
            for month in range(1, 13):
                if (lyear == syear) and (month < smonth):
                    continue
                elif (lyear == eyear) and (month > emonth):
                    break

                start_time.append(hcf.conv_time_format(lyear, month))
                end_time.append(hcf.conv_time_format(lyear, month, next=1))

    else:
        print "Error in the time order"
        exit(1)

    return [start_time, end_time]


#-----------------------------------------------------------------------------------------------
#-- find_time_period: find time period from existing data directory name                     ---
#-----------------------------------------------------------------------------------------------

def find_time_period():
    """
    find time period from existing data directory name
    input:  none but read from <data_dir>/2*
    output: syear   --- starting year
            smon    --- starting month
            eyear   --- ending year
            emon    --- ending month
            if also return "False" if something wrong (e.g. wrong time order)
    """
#
#--- find the last data set extracted
#
    cmd  = 'ls -d ' + data_dir + '2* > ' + zspace
    os.system(cmd)
    data = hcf.read_file_data(zspace, remove=1)

    syear = 1999
    smon  = 1

    if len(data) > 0:
        for ent in data:
            atemp = re.split('\/', ent)
            val   = atemp[-1]
#
#--- extract year and month in digit from the directory name
#
            year = int(float(str(val[0:4])))
            lmon = val[4:]
    
            mon  = hcf.find_month(lmon)
#
#--- find the latest time
#
            if syear == year:
                if smon < mon:
                    syear = year
                    smon  = mon
            elif syear < year:
                syear = year
                smon  = mon
#
#--- find today's date
#
    today = time.localtime()

    eyear = int(float(today.tm_year))
    emon  = int(float(today.tm_mon))
#
#--- set the last month as we know this month's data are not filled yet
#
    emon -= 1

    if emon < 1:
        emon   = 12
        eyear -= 1
#
#--- quick error check
#
    if eyear < syear:
        return False

    elif eyear == syear:
        if smon >= emon:
            return False

    return [syear, smon, eyear, emon]

#-----------------------------------------------------------------------------------------------
#-- combine_short_exposure_files:  combine evt1 files if the time difference is less than < 180 sec
#-----------------------------------------------------------------------------------------------

def combine_short_exposure_files(evt1_list):
    """
    combine evt1 files if the time difference between two fits files is less than < 180 sec
    output: combined fits files if they meet the condition
    """
#
#--- check whether the file exists
#
    if len(evt1_list) < 2:
        return evt1_list

    else:
#
#-- find the first evt1 file start and stop time
#
        out_list = []
        prev  = evt1_list[0]
        [pstart, pstop] = find_start_and_end(prev)

        for i in range(1, len(evt1_list)):
#
#--- check the next evt1 file and check whether the time span is less than 180 sec
#--- if so, combine them
#
            [nstart, nstop] = find_start_and_end(evt1_list[i])

            if  (nstart - pstop) < 180:
                hcf.combine_fits_files([prev, evt1_list[i]], prev)
#
#--- remove the added evt1 file
#
                cmd = 'rm ' + evt1_list[i] 
                os.system(cmd)
            else:
                out_list.append(prev)

                prev     = evt1_list[i]
                pstart   = nstart
                pstop    = nstop 

        out_list.append(prev)

        return out_list

#-----------------------------------------------------------------------------------------------
#-- find_start_and_end: find start and stop time of the data table in the fits file           --
#-----------------------------------------------------------------------------------------------

def find_start_and_end(fits):
    """
    find start and stop time of the data table in the fits file
    input:  fits    --- the name of the fits file
    output: start   --- starting time in seconds from 1998.1.1
            stop    --- stopping time in seconds from 1998.1.1
    """

    [cols, tbdata] = hcf.read_fits_file(fits)
    time  = tbdata['time']
    try:
        start = numpy.min(time)
        stop  = numpy.max(time)
    except:
        start = 0
        stop  = 0

    return [start, stop]


#-----------------------------------------------------------------------------------------
#-- TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST    ---
#-----------------------------------------------------------------------------------------

class TestFunctions(unittest.TestCase):
    """
    testing functions
    """

#------------------------------------------------------------

    def test_create_hk_fits_with_extracted_data(self):

        start = '12/01/15,00:00:00'
        stop  = '01/01/16,00:00:00'
        fits_list = create_hk_fits_with_extracted_data(start, stop)

        self.assertEquals(fits_list[0], 'hrc_i_115.fits') 

        for ent in fits_list:
            mcf.rm_file(ent)

#------------------------------------------------------------

    def test_find_tscpos_positive_period(self):

        start_comp = [536504737, 536731844, 537190552]
        stop_comp  = [536681398, 536911687, 537235200]

        start = '01/01/15,00:00:00'
        stop  = '01/10/15,00:00:00'

        [b_list, e_list]  = find_tscpos_positive_period(start, stop)

        self.assertEquals(b_list, start_comp)
        self.assertEquals(e_list, stop_comp)

#------------------------------------------------------------

    def test_read_condition(self):

        cfile     = 'hrc_i_115'
        cond_comp = [['LLDIALV', '130:132'], ['RSRFALV', '184:186'], ['WDTHAST', '2:3'], \
                     ['ULDIALV', '254:255'], ['CBHUAST', '255:255'], ['CBHVAST', '255:255'],\
                     ['CBLUAST', '0:0'], ['CBLVAST', '0:0'], ['IMTPAST', '77:77'],\
                     ['IMBPAST', '89:89'], ['IMHVLV', '78:82'], ['IMHBLV', '85:88'],\
                     ['SMTRATM', '0:54'], ['HVPSSTAT', 'xxx1xxxx'],
        #             ['QUALITY', '0000000000000000000000000000000000000000000000000000000'],\
                     ['SCIDPREN', '0000xxxx00000001']]
        name      = 'hrc_i_115.fits'

        [condition, fits] = read_condition(cfile)

        self.assertEquals(condition, cond_comp)
        self.assertEquals(fits, name)

#------------------------------------------------------------
    
    def test_set_data_period(self):

        syear  = 2012
        smonth = 5
        eyear  = 2014
        emonth = 2

        [start_time, stop_time] = set_data_period(syear, smonth, eyear, emonth)

        self.assertEquals(start_time[0],  '05/01/12,00:00:00')
        self.assertEquals(start_time[7],  '12/01/12,00:00:00')
        self.assertEquals(start_time[-1], '02/01/14,00:00:00')

        self.assertEquals(stop_time[0],  '06/01/12,00:00:00')
        self.assertEquals(stop_time[7],  '01/01/13,00:00:00')
        self.assertEquals(stop_time[-1], '03/01/14,00:00:00')

#------------------------------------------------------------
    
    def test_find_time_period(self):

        #[syear, smon, eyear, emon] = find_time_period()
        alist = find_time_period()
        
        print "/ntest_find_time_period RESULTS:  " + str(alist)

#------------------------------------------------------------

    def test_extrct_hrchk_data(self):

        start = '12/22/15,00:00:00'
        stop  = '12/24/15,00:00:00'

        dout  = extrct_hrchk_data(start, stop)

        tcomp = [131.0, 185.0, 2.0, 255.0, 255.0]
        tout  = dout[-1]
        tlist = list(tout[1:6])

        self.assertEquals(tlist, tcomp)

#------------------------------------------------------------

    def test_estimate_rsrfalv_wdthast(self):

        hkfits = house_keeping + 'Test_data/hrcf567118652N001_hk0.fits.gz'

        [rsrfalv, wdthast] = estimate_rsrfalv_wdthast(hkfits)

        self.assertEquals(rsrfalv, 116)
        self.assertEquals(wdthast, 2)

#------------------------------------------------------------

    def test_extract_data_and_combine(self):

        time_list = [[553651195], [553694395]]
        detector  = 'hrc'
        level     = 0
        filetype  = 'hrchk'
        name      = 'test_hk.fits'
        comp      = [43, 43, 43, 43, 43]

        extract_data_and_combine(time_list, detector, level, filetype, name)

        [cols, tbdata] = hcf.read_fits_file(name)
        tout = tbdata['SPTPAST']
        self.assertEquals(list(tout[:5]), comp)

        mcf.rm_file(name)

#------------------------------------------------------------

    def test_extract_fits_files(self):

        time_list = [553651195, 553694395]
        detector  = 'hrc'
        level     = 0
        filetype  = 'hrchk'
        comp      = ['/data/aschrc6/wilton/isobe/Project1_new/Exc/Temp_dir/hrcf553641557N001_hk0.fits.gz', \
                     '/data/aschrc6/wilton/isobe/Project1_new/Exc/Temp_dir/hrcf553659400N001_hk0.fits.gz', \
                     '/data/aschrc6/wilton/isobe/Project1_new/Exc/Temp_dir/hrcf553692135N001_hk0.fits.gz']

        dlist = extract_fits_files(time_list, detector, level, filetype)
        
        self.assertEquals(dlist, comp)

        for ent in dlist:
            mcf.rm_file(ent)

#------------------------------------------------------------

    def test_get_time_period_engineer_data(self):

        time_list = [553651195, 553694395]
        hrc_type = "IMAG"
        start    = 553651195
        stop     = 553694395
        get_time_period_engineer_data(time_list, hrc_type, start, stop)

#------------------------------------------------------------

    def test_find_start_and_end(self):

        fits = house_keeping + 'Test_data/hrcf383576042N002_evt1.fits.gz'

        [start, stop]  = find_start_and_end(fits)

        self.assertEqual(int(start), 383576040)
        self.assertEqual(int(stop),  383588913)

#-----------------------------------------------------------------------------------------------
 
if __name__ == "__main__":
#
#--- if you like to specify the date period, give
#---  a year and starting yday and ending yday
#
    test = 0
    if len(sys.argv) == 2:
        if sys.argv[1] == 'test':
            test = 1
            del sys.argv[1:]
        else:
            exit(1)

    elif len(sys.argv) == 3:
        syear  = int(float(sys.argv[1]))
        smonth = int(float(sys.argv[2]))
        eyear  = syear
        emonth = eyear
    elif len(sys.argv) == 5:
        syear  = int(float(sys.argv[1]))
        smonth = int(float(sys.argv[2]))
        eyear  = int(float(sys.argv[3]))
        emonth = int(float(sys.argv[4]))
#
#--- if the date period is not specified,
#--- the period is set from the last entry date to
#--- the last day of the last month
#
    else:
        syear  = ''
        smonth = ''
        eyear  = ''
        emonth = ''

    if test == 0:
        lmon = tcnv.changeMonthFormat(smonth)
        lmon = lmon.upper()
        outdir = str(syear) + lmon

        ####extract_next_in_line_data(syear, smonth, eyear, emonth, outdir)
#following is not working --- you need to change the time formats!!
        extract_next_in_line_data(begin, end, start, stop, outdir)
    else:
        unittest.main()


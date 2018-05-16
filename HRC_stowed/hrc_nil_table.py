#!/usr/bin/env /proj/sot/ska/bin/python

#############################################################################################
#                                                                                           #
#       hrc_nil_table.py: create/update stat data tables                                    #
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
#--- month list
#
m_list = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
#
#--- status bit description
#
descriptor = ['Ringing corrected, V-axis     ',
              'Ringing corrected, U-axis     ',
              'Piled-up (tentative)          ',
              '          - spare -           ',
              'Shifted event time            ',
              'Next-in-line mode telemetry   ',
              'No-trigger, V-axis            ',
              'No-trigger, U-axis            ',
              'Center-blank, V-axis          ',
              'Center-blank, U-axis          ',
              'Width exceeded, V-axis        ',
              'Width exceeded, U-axis        ',
              'Antico-shield active          ',
              '          - spare -           ',
              'ULD exceeded                  ',
              'LLD not exceeded              ',
              'Event in Bad Region           ',
              'Tap total on U- or V-axis <= 0',
              'Bad CRSV (AV1 or AV3 > AV2)   ',
              'Bad CRSU (AU1 or AU2 > AU2)   ',
              'PHA-ratio test failed         ',
              'Sum of 6 taps = 0             ',
              'Grid-ratio test failed        ',
              'ADC sum on U- or V-axis = 0   ',
              'PI > 255                      ',
              'Out of time-sequence          ',
              'Flatness test failed, V-axis  ',
              'Flatness test failed, U-axis  ',
              'Saturation test failed, V-axis',
              'Saturation test failed, U-axis',
              'H-test failed, V-axis         ',
              'H-test failed, U-axis         ']
#
#--- set  a few things
#
NULL      = 'null'

#-----------------------------------------------------------------------------------------------
#-- hrc_nil_table: create data tables for hrc-i and hrc-s high precision data                 --
#-----------------------------------------------------------------------------------------------

def hrc_nil_table(year, month):
    """
    create data tables for hrc-i and hrc-s high precision data
    input:  year    --- year of the data 
            month   --- month of the data
    output: <hrc name>__stat_results    --- a table of statistics obtained
    """

    lmon = m_list[month-1]

    tdir  = data_dir + str(year) + lmon + '/'

    for sdir in ("hrc_i_115", "hrc_i_90", "hrc_s_125_hi", "hrc_s_90_hi"):
         print_stat_results(tdir, sdir)
#
#--- hrc_s_125 and hrc_s_90 have 3 sections to check
#
    for sdir in ("hrc_s_125", "hrc_s_90"):
         print_stat_results(tdir, sdir, chip_ind = 1)
#
#--- clean up the stat files
#
    remove_duplicate_stat_entries()

#-----------------------------------------------------------------------------------------------
#-- remove_duplicate_stat_entries: removing duplicated stat results from stat files          ---
#-----------------------------------------------------------------------------------------------

def remove_duplicate_stat_entries():
    """
    removing duplicated stat results from stat files
    input:  none, but read from <inst>_stat_results
    output: cleaned up <inst>_stat_results
    """
#
#--- clean stat data table
#
    cmd = 'ls ' + stat_dir + '*stat_results > ' + zspace
    os.system(cmd)

    data = hcf.read_file_data(zspace, remove=1)
    
    for sfile in data:
        hcf.remove_duplicate_from_file(sfile)
#
#--- clean deat time table
#    
    cmd = 'ls ' + stat_dir + '*dead_time  > ' + zspace
    os.system(cmd)

    data = hcf.read_file_data(zspace, remove=1)
    
    for sfile in data:
        hcf.remove_duplicate_by_column(sfile, 0)

#-----------------------------------------------------------------------------------------------
#-- print_stat_results: compute general statistiscs and dead_time correction and append to tables
#-----------------------------------------------------------------------------------------------

def print_stat_results(tdir, sdir, chip_ind = ''):
    """
    compute general statistiscs and dead_time correction and append to tables
    input:  tdir    --- the main directory e.g, 2010MAR
            sdir    --- the sub directory e.g., hrc_i_115
    output: <inst>_stat_results --- general stat result table. for hrc_s, _sec#_ is added
            <inst>_dead_time    --- dead time stat result table
    """
#
#--- hrc_s_125 and hrc_s_90 have three seciontions
#
    if chip_ind == '':
        chip_list = ['']
    else:
        chip_list = [1, 2, 3]
#
#--- find evt1 file names
#
    ldir  = tdir + sdir + '/'
    flist = hcf.get_file_list(ldir, 'evt1.fits.gz')

    for fits in flist:
        atemp   = re.split('_evt1', fits)
        ss_file = atemp[0] + '_comb_ss0.fits.gz'
        hk_file = atemp[0] + '_comb_hk0.fits.gz'
#
#--- hrc_s_125 and hrc_s_90 needs to go sec 1 to 3, but others just spit out one result
#
        tstart = hcf.read_header_value(fits, 'tstart')
        tstop  = hcf.read_header_value(fits, 'tstop')

        [tstat, vstat, sstat, astat]    = find_antico_rate(ss_file, tstart, tstop)

        [s2hvst_avg, s2hvlv_avg]        = find_hvlv(hk_file)
        [scint, sint_sig]               = find_scint(tstart, tstop)

        for chip_id in chip_list:

            [date_obs, tstart, tstop, avg, med, std, duration, tcnt, cnt_p_sec]   \
                                        = find_evt_stat(fits, chip_id)
#
#--- general statistics table output
#
            line = date_obs + '\t'
            line = line + str(tstart)     + '\t'
            line = line + str(tstop)      + '\t'
            line = line + str(duration)   + '\t'
            line = line + str(tcnt)       + '\t'
            line = line + str(cnt_p_sec)  + '\t'
            line = line + str(avg)        + '\t'
            line = line + str(med)        + '\t'
            line = line + str(std)        + '\t'
            line = line + str(tstat[0])   + '\t'
            line = line + str(tstat[1])   + '\t'
            line = line + str(tstat[2])   + '\t'
            line = line + str(vstat[0])   + '\t'
            line = line + str(vstat[1])   + '\t'
            line = line + str(vstat[2])   + '\t'
            line = line + str(sstat[0])   + '\t'
            line = line + str(sstat[1])   + '\t'
            line = line + str(sstat[2])   + '\t'
            line = line + str(astat[0])   + '\t'
            line = line + str(astat[1])   + '\t'
            line = line + str(astat[2])   + '\t'
            line = line + str(s2hvst_avg) + '\t'
            line = line + str(s2hvlv_avg) + '\t'
            line = line + str(scint)      + '\t'
            line = line + str(sint_sig)   + '\n'
#
#--- dead time correction table output
#
            if vstat[1] == 0:
                continue

            arate = cnt_p_sec / vstat[1]
            corr  = duration * arate
            tcorr = 1.0 / corr

            arate = '%4.4f' % arate
            corr  = '%4.1f' % corr
            tcorr = '%4.6f' % tcorr
            dura  = '%6.1f' % duration

            line2 = str(int(tstart))   + '\t'
            line2 = line2 + date_obs   + '\t'
            line2 = line2 + str(arate) + '\t'
            line2 = line2 + str(dura)  + '\t'
            line2 = line2 + str(corr)  + '\t'
            line2 = line2 + str(tcorr) + '\n'

            if chip_id == '':
                out  = stat_dir + sdir + '_stat_results'
                out2 = stat_dir + sdir + '_dead_time'
            else:
                out  = stat_dir + sdir + '_sec' + str(chip_id) + '_stat_results'
                out2 = stat_dir + sdir + '_sec' + str(chip_id) + '_dead_time'
    
            fo  = open(out,  'a')
            fo.write(line)
            fo.close()
    
            fo  = open(out2, 'a')
            fo.write(line2)
            fo.close()

#-----------------------------------------------------------------------------------------------
#-- find_evt_stat: find count rate, pha postion from hrc evt 1 file                          ---
#-----------------------------------------------------------------------------------------------

def find_evt_stat(efits, chip_id = ''):
    """
    find count rate, pha postion from hrc evt 1 file
    input   efits   --- hrc evt1 fits file name
            chip_id     --- chip_id for hrc_s_125 and hrc_s_90; default: 0
    output: date_obs    --- date of observation
            tstart      --- starting time in seconds from 1998.1.1
            tstop       --- stopping time in seconds from 1998.1.1
            avg         --- average of pha
            med         --- mediam of pha
            std         --- standard deviation of pha
            duration    --- duration of the evt1 data
            tcnt        --- the number of entries
            cnt_p_sec   --- counts per second
    """
#
#--- read values from header
#
    date_obs = hcf.read_header_value(efits, 'date-obs')
    tstart   = hcf.read_header_value(efits, 'tstart')
    tstop    = hcf.read_header_value(efits, 'tstop')
#
#--- read data 
#
    [cols, tbdata] = hcf.read_fits_file(efits)
#
#--- if chip_id is given, select out the data for the chip_id
#
    if chip_id != '':
        tbdata = hcf.select_data_with_condition(tbdata, 'chip_id', '==', chip_id)
#
#--- find stat of pha
#
    [avg, med, std, vmin, vmax, vcnt]    = get_basic_stat(tbdata['pha'])
#
#--- find time related quantities
#
    time     = tbdata['time']
    [tavg, tmed, tstd, tmin, tmax, tcnt] = get_basic_stat(tbdata['time'])
#
#--- compute duration etc
#
    duration = tmax - tmin
    if duration > 0:
        cnt_p_sec =  tcnt /duration
    else:
        cnt_p_sec = 0

    return [date_obs, tstart, tstop, avg, med, std, duration, tcnt,  cnt_p_sec]

#-----------------------------------------------------------------------------------------------
#-- get_basic_stat: compute basic statistics for a given numpy array                         ---
#-----------------------------------------------------------------------------------------------

def get_basic_stat(tarray):
    """
    compute basic statistics for a given numpy array
    input:  tarray  --- one dim numpy array
    output: avg     --- average 
            med     --- mediam
            std     --- standard deviation
            vmin    --- minimum value
            vmax    --- maximum value
            vcnt    --- numbers of entries
    """
    vcnt = tarray.shape[0]

    if vcnt > 0:
        avg  = numpy.mean(tarray)
        med  = numpy.median(tarray)
        std  = numpy.std(tarray)
        vmin = numpy.min(tarray)
        vmax = numpy.max(tarray)
    else:
        avg  = 0.0
        med  = 0.0
        std  = 0.0
        vmin = 0.0
        vmax = 0.0

    return [avg, med, std, vmin, vmax, vcnt]

#-----------------------------------------------------------------------------------------------
#-- find_antico_rate: find a ratio of valid/total of anti-co data                            ---
#-----------------------------------------------------------------------------------------------

def find_antico_rate(ss_file, tstart, tstop):
    """
    find a ratio of valid/total of anti-co data
    input:  ss_file --- science fits file name
            tstart  --- starting time in seconds from 1998.1.1
            tstop   --- stopping tine in seconds from 1998.1.1
    output: tstat   --- statistics of tlevart column
            vstat   --- statistics of vlevart column
            sstat   --- statistics of shevart column
            astat   --- statistics of antico values
            see get_basic_stat for output (a list of statistics)
    """
    [cols, tbdata] = hcf.read_fits_file(ss_file)
    
    tbdata  = hcf.select_data_with_condition(tbdata, 'time', ">=", tstart)
    tbdata  = hcf.select_data_with_condition(tbdata, 'time', "<=", tstop)
#
#--- remove outlyers; the first several entries are not reliable 
#
    tbdata  = hcf.select_data_with_condition(tbdata, 'tlevart', ">", 0)
    tbdata  = hcf.select_data_with_condition(tbdata, 'tlevart', "<", 500)
    tbdata  = hcf.select_data_with_condition(tbdata, 'vlevart', ">", 0)
    tbdata  = hcf.select_data_with_condition(tbdata, 'vlevart', "<", 300)

    tlevart = tbdata['tlevart']
    vlevart = tbdata['vlevart']
    shevart = tbdata['shevart']
#
#--- compute antico values
#
    antico = vlevart / tlevart      #--- array / array gives array 
#
#--- compute basic statistics for each column
#
    tstat = get_basic_stat(tlevart)
    vstat = get_basic_stat(vlevart)
    sstat = get_basic_stat(shevart)
    astat = get_basic_stat(antico)
#
#--- output are in a list: [avg, med, std, min, max, cnt]
#
    return [tstat, vstat, sstat, astat]

#-----------------------------------------------------------------------------------------------
#-- find_hvlv: find values of s2hvst and shvlv from hrc hk files                              --
#-----------------------------------------------------------------------------------------------

def find_hvlv(hk_file):
    """
    find values of s2hvst and shvlv from hrc hk files
    input:  hk_file --- a hrc hk fits file name
    output: s2hvst_avg  --- average of s2hvst
            s2hvlv_avg  --- average of shvlv
    """

    [cols, tbdata] = hcf.read_fits_file(hk_file)

    [s2hvst_avg, med, std, vmin, vmax, vcnt] = get_basic_stat(tbdata['s2hvst'])
    [s2hvlv_avg, med, std, vmin, vmax, vcnt] = get_basic_stat(tbdata['s2hvlv'])
#
#--- make sure that the values are integer
#
    s2hvst_avg = int(s2hvst_avg + 0.5)
    s2hvlv_avg = int(s2hvlv_avg + 0.5)

    return [s2hvst_avg, s2hvlv_avg]

#-----------------------------------------------------------------------------------------------
#-- find_scint: find integrated EPHIN flux                                                    --
#-----------------------------------------------------------------------------------------------

def find_scint(tstart, tstop):
    """
    find integrated EPHIN flux
    input: tstart   --- starting time in seconds from 1998.1.1
           tstop    --- stopping time in seconds from 1998.1.1
    output; avg     --- average of ephin flux
            std     --- standard deviation of flux
    """
#
#--- extract ephin data
#
    hcf.run_arc5gl('retrieve', tstart, tstop, dataset='flight',detector='ephin', level=1,filetype='ephrates')
#
#--- create extracted fits file list
#
    ldir  = exc_dir + 'Temp_dir/'
    flist = hcf.get_file_list(ldir, 'lc1.fits.gz')
#
#--- combine all fits files into one
#
    hcf.combine_fits_files(flist, 'ztemp.fits', azip=0)
    [cols, tbdata] = hcf.read_fits_file('ztemp.fits')
#
#--- limit data to between tstart and tstop
#
    tbdata  = hcf.select_data_with_condition(tbdata, 'time', ">=", tstart)
    tbdata  = hcf.select_data_with_condition(tbdata, 'time', "<=", tstop)
#
#--- compute the basic stat
#
    [avg, med, std, vmin, vmax, vcnt] = get_basic_stat(tbdata['scint'])

    cmd = 'rm ztemp.fits' + '*'
    os.system(cmd)

    if hcf.check_file_in_dir(ldir):
        cmd = 'rm  ' + ldir + '*'
        os.system(cmd)

    return [avg, std]

#-----------------------------------------------------------------------------------------------
#-- status_bit_report: create status bit statitics tables                                    ---
#-----------------------------------------------------------------------------------------------

def status_bit_report(year, month):
    """
    create status bit statitics tables
    input:  year    --- year  in yyyy form
            month   --- month in degit
    output: hrc_status_stat_<time samp> --- bit statistics of the corresponding evt1 file
            <inst name>_status_stat     --- bit statistics of all evt1 file of the <inst name>
    """

    #cmd  = 'ls ' +  data_dir + str(year) + str(m_list[month-1]) + '/*/hrcf*evt1.fits.gz > ' + zspace
    cmd  = 'ls ' +  data_dir + str(year) + str(m_list[month-1]) + '/*/* > ' + zspace
    os.system(cmd)

    data = hcf.read_file_data(zspace, remove=1)

    for efile in data:
        mc  = re.search('evt1.fits.gz', efile)
        if mc is None:
            continue
#
#--- find an output file names and the date of the data
#
        atemp = re.split('hrcf', efile)
        btemp = re.split('N', atemp[1])
        iname = atemp[0] + 'hrc_status_stat_' + str(btemp[0])  #--- name of the indivisula bits stat output
    
        atemp = re.split('\/', efile)
        hrc_n = stat_dir +  atemp[-2] + '_status_stat'          #--- name of the cumulative data file
        btemp = re.split('hrcf', atemp[-1])
        ctemp = re.split('N', btemp[1])
        date  = hcf.covertfrom1998sec2(int(float(ctemp[0])))

        run_status_bits(efile, iname, hrc_n, date)

#-----------------------------------------------------------------------------------------------
#-- run_status_bits: create two types of bit statistics tables                               ---
#-----------------------------------------------------------------------------------------------

def run_status_bits(efile, iname, hrc_n, date):
    """
    create two types of bit statistics tables
    input:  efile   --- a full path to the evt1 file
            iname   --- a bit stat data file (cumulative)
            hrc_n   --- a bit stat data file of the evt1 file
            date    --- the date of the evt1 file
    output: iname and hrc_n
    """
#
#--- read fits file
#
    [cols, tbdata] = hcf.read_fits_file(efile)

    status = tbdata['status']
#
#--- initialize status bit saver
#
    stat = []
    for i in range(0, 32):
        stat.append(0)
#
#--- read status bits and counts them
#
    k = 0
    for ent in status:
        for i in range(0, 32):
            if ent[i] == True:
                stat[i] += 1

    dat_len = len(status)
#
#--- print out two different files. one for indivisual output and other for cumulative output from all data
#
    line  = '\nChecking STATUS bits in file: {}'.format(efile)
    line  = line + '\n\nChecked {} events'.format(len(status))
    line  = line + '\n\nBit | Description                    | Number of Times Set (Fraction)'
    line  = line + '\n---------------------------------------------------------------------\n'

    line2 = date + '\t' + str(dat_len) 

    for i in range(0, 32):
        line  = line + ' {0:2d} | {1} | {2:12d} ({3:.3f})'.format(i+1,descriptor[i],stat[i],1.0*stat[i]/dat_len)
        line  = line + '\n'

        line2 = line2 + '\t' + str(stat[i])

    fo  = open(iname, 'w')
    fo.write(line)
    fo.close()
#
#--- append the result to the data table
#
    fo  = open(hrc_n, 'a')
    fo.write(line2)
    fo.write('\n')
    fo.close()
#
#--- clean up the file
#
    hcf.remove_duplicate_from_file(hrc_n)


#-----------------------------------------------------------------------------------------
#-- TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST    ---
#-----------------------------------------------------------------------------------------

class TestFunctions(unittest.TestCase):
    """
    testing functions
    """

#------------------------------------------------------------

    def test_find_evt_stat(self):
        comp = ['2010-02-26T12:54:03', 383576042.97157, 383591130.71608, 116.27673387990251, \
                 96.0, 78.538520763718381, 12873.232206702232, 45130, 3.5057240695544838]
#
#--- [date_obs, tstart, tstop, avg, med, std, duration, tcnt,  cnt_p_sec]
#
        efits = house_keeping + 'Test_data/hrcf383576042N002_evt1.fits.gz'
        out_list = find_evt_stat(efits)

        self.assertEquals(out_list, comp)

#------------------------------------------------------------

    def test_get_basic_stat(self):

        comp = [5.5, 5.5, 2.8722813232690143, 1, 10, 10]

        tarray   = numpy.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        out_list = get_basic_stat(tarray)
        
        self.assertEquals(out_list, comp)

#------------------------------------------------------------

    def test_find_antico_rate(self):
        comp1 = [346.35714285714283, 349.5, 22.5249294999916, 308.0, 401.0, 14]
        comp2 = [122.14285714285714, 119.5, 9.4706851679812694, 105.0, 137.0, 14]
        comp3 = [6362.3571428571431, 6374.0, 171.163083112008, 6085.0, 6626.0, 14]
        comp4 = [0.35312865980484759, 0.35618336886993607, 0.023742017572634166, \
                 0.30973451327433627, 0.38557993730407525, 14]


        ss_file = house_keeping + 'Test_data/hrcf383572828N002_comb_ss0.fits.gz'
        tstart = 383572600
        tstop  = 383576600

        [tstat, vstat, sstat, astat] = find_antico_rate(ss_file, tstart, tstop)

        self.assertEquals(tstat, comp1)
        self.assertEquals(vstat, comp2)
        self.assertEquals(sstat, comp3)
        self.assertEquals(astat, comp4)

#------------------------------------------------------------

    def test_find_hvlv(self):

        hk_file = house_keeping + 'Test_data/hrcf567118652N001_hk0.fits.gz'

        [s2hvst_avg, s2hvlv_avg] = find_hvlv(hk_file)

        self.assertEquals(s2hvst_avg, 7)
        self.assertEquals(s2hvlv_avg, 90)

#------------------------------------------------------------

    def test_find_scint(self):

        tstart = 383572600
        tstop  = 383576600

        [avg, std] = find_scint(tstart, tstop)

        avg = round(avg, 11)
        std = round(std, 13)
        self.assertEquals(avg, 1.34785830624)
        self.assertEquals(std, 0.0835570207673)

#------------------------------------------------------------

    def test_run_status_bits(self):

        fits = house_keeping + 'Test_data/hrcf383576042N002_evt1.fits.gz'
        name1 = './test1'
        name2 = './test2'
        date  = '9999-12-31T00:00:00'
        run_status_bits(fits, name1, name2, date)

        cmd   = 'cat ' + name2 + '> ' + zspace
        os.system(cmd)
        f     = open(zspace, 'r')
        out   = f.read()
        f.close()

        mcf.rm_file(zspace)
        mcf.rm_file(name1)
        mcf.rm_file(name2)

        test = re.split('\s+', out)
        self.assertEquals(test[1:4], ["45130", "6373", "5797"])


#-----------------------------------------------------------------------------------------------
 
if __name__ == "__main__":
#
#---  give a year and month to create data tables
#---  give "test" to run test
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

    if test == 0:
        hrc_nil_table(syear, smonth)
        status_bit_report(syear, smonth)
    else:
        unittest.main()

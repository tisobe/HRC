#!/usr/bin/env /proj/sot/ska/bin/python

#############################################################################################
#                                                                                           #
#   run_hrc_stowed_background_monthly_update.py: run monthly update for the hrc stowed      #
#                                                background data                            #
#                                                                                           #
#               author: t. isobe (tisobe@cfa.harvard.edu)                                   #
#                                                                                           #
#               last update: Jan 18, 2017                                                   #
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
import create_hrc_image_files     as chif       #---- contains functions to create image fits files
import hrc_plotting_maps          as hpm        #---- contains functions to create png maps 
import update_hrc_html_page       as uhhp       #---- contains functions to update html pages
import hrc_stowed_background      as hsb        #---- contains functions to extract and create hrc evt files
import hrc_plotting_routines      as hpr        #---- contains functions to plot trends
#
#--- temp writing file name
#
rtail  = int(time.time())
zspace = '/tmp/zspace' + str(rtail) + 'hrc'
#
#--- set  a few things
#
NULL      = 'null'

#-----------------------------------------------------------------------------------------------
#-- run_hrc_stowed_background_monthly_update: run monthly update for the hrc stowed background data 
#-----------------------------------------------------------------------------------------------

def run_hrc_stowed_background_monthly_update(lyear='', lmonth=''):
    """
    run monthly update for the hrc stowed background data
    input:  lyear, lmonth   --- the year and month of the data to be processed
                                if they are not given, the last month of year and year will be used
    output: newly proccessed data: <data_dir>/<yyyy><mmm>/*
            if it is jan or jul, cummulative maps for the year (and the entire period) are created
    """
#
#--- clean up exc dir
#
    cmd = 'rm -rf ' + exc_dir + 'Temp_dir/* 2> /dev/null'
    os.system(cmd)

    cmd = 'rm -rf ' + exc_dir + 'param 2> /dev/null'
    os.system(cmd)

    cmd = 'mkdir '  + exc_dir + 'param 2> /dev/null'
    os.system(cmd)
#
#--- if date (year/month) to be processed is not given, find the last month's year, month
#
    if lyear == '':
        ctime  = time.gmtime()      #--- today's date information
        lyear  = ctime.tm_year
        lmonth = ctime.tm_mon -1

        if lmonth < 1:
            lyear -= 1
            lmonth = 12
#
#--- check whether the data was already processed; if not, run the process
#
    smonth = tcnv.changeMonthFormat(lmonth)
    cdir   = data_dir + str(lyear) + smonth.upper()

    smail = 0
    if not os.path.isdir(cdir):
        try:
            hsb.hrc_stowed_background(lyear, lmonth, lyear, lmonth)
#
#--- if the processed month is june or december, update maps and html pages
#
            if (lmonth == 12) or (lmonth == 6):
#
#--- update total event number tables
#
                get_yearly_evt_count()

                chif.yearly_cummulative_image(lyear)
                hpm.plot_hrc_map(lyear)
                hpm.plot_hrc_map('total')       #---- creating cumulative maps
                hpm.plot_hrc_map('total', 1)    #---- creating cumulative histogram for the front page
                hpr.plot_hrc_trend_data()       #---- creating trending plots
                uhhp.create_html_pages(lyear)   #---- updating all html pages
                create_cumulative_image()

            if test_data_creation(lyear, lmonth):
                smail = 1
            else:
                smail = 3
        except:
            smail = 2
    else:
        smail = 0
#
#--- remove the temp file if it is still around
#
    cmd = 'rm /tmp/zspace*hrc* 2> /dev/null'
    os.system(cmd)

    cmd = 'rm -rf ' + exc_dir + 'Temp_dir/* 2> /dev/null'
    os.system(cmd)
#
#--- sending mail
#
    if smail == 1:
        message = 'hrc stowed proccess ran normally for the period: ' + str(lyear) + ' : ' + str(lmonth) + '\n'
        header  = '"Subject: HRC Stowed Data Process Complated"'

    elif smail == 2:
        message = 'hrc stowed proccess ran into problems with the data for the period: ' + str(lyear) + ' : ' + str(lmonth) + '\n'
        header  = '"Subject: HRC Stowed Data Process Failed"'

    elif smail == 3:
        message = 'hrc stowed proccess did not find any data for the period: ' + str(lyear) + ' : ' + str(lmonth) + '\n'
        header  = '"Subject: HRC Stowed Data Process: No Data Found"'

    if smail > 0:
        send_email(message, header)


#---------------------------------------------------------------------
#-- get_yearly_evt_count: counts the numbers of events in the yearly evt1 files
#---------------------------------------------------------------------

def get_yearly_evt_count():
    """
    counts the numbers of events in the yearly evt1 files
    input: none, but read from yearly evt1.fits files
    output: <data_dir>/<inst>/<inst>_yearly_evt_counts
    """
#
#--- find today's date
#
    out   = time.localtime()
    lyear = out.tm_year
    mon   = out.tm_mon
#
#--- if the date is the first half of the year, compute the last year
#
    if mon > 6:
        lyear += 1
#
#--- check each instrument
#
    for inst in inst_list:
        #print inst
        ofile = data_dir + inst + '/' + inst.lower() + '_yearly_evt_counts'
        fo    = open(ofile, 'w')

        for year in range(2000, lyear):
            #print '\t' + str(year)
            infile = data_dir + inst + '/'+ inst.lower() + '_' + str(year) + '_evt1.fits.gz'
#
#--- read the numbers of events
#
            try:
                out  = hscf.fitsTableStat(infile, 'time')
                cnt  = out[-1]
            except:
                cnt  = 0


            line = str(year) + '\t' + str(cnt) + '\n'
            fo.write(line)


#-----------------------------------------------------------------------------------------------
#-- send_email: sending out mail                                                              --
#-----------------------------------------------------------------------------------------------

def send_email(message, header, email = 'tisobe@cfa.harvard.edu'):
    """
    sending out mail
    input:  message --- content of the mail
            header  --- header of the mail
            email   --- a receiver's email address
    outout:  none but email is sent out        
    """

    fo = open(zspace, 'w')
    fo.write(message)
    fo.close()

    cmd = 'cat ' + zspace + ' | mailx -s ' + header + ' ' + email
    os.system(cmd)

    mcf.rm_file(zspace)

#-----------------------------------------------------------------------------------------------
#-- test_data_creation: check whether any data files are created                              --
#-----------------------------------------------------------------------------------------------

def test_data_creation(lyear, lmonth):
    """
    check whether any data files are created
    input:  lyear   --- year of the data
            lmonth  --- month of the data
    output: True or False
    """

    lmonth = tcnv.changeMonthFormat(lmonth)
    lmonth = lmonth.upper()
    cmd    = 'ls ' + data_dir + str(lyear) + str(lmonth) + '/*/*evt1* > ' + zspace
    os.system(cmd)
    data = hcf.read_file_data(zspace, remove=1)
    if len(data) > 0:
        return True
    else:
        return False

#-----------------------------------------------------------------------------------------------
#-- create_cumulative_image: create normalized combined image fits files                      --
#-----------------------------------------------------------------------------------------------

def create_cumulative_image():
    """
    create normalized combined image fits files
    input:  none but read from each directory
    output: combine image fits files, e.g. hrc_i_115_total_norm.fits.gz
            it also create png files
    """

    for hdir in ['Hrc_i_115', 'Hrc_s_125_hi', 'Hrc_s_125']:

        head  = hdir.lower()
#
#--- hrc_s_125 has three parts
#
        if hdir == 'Hrc_s_125':
            p_list = ['1', '2', '3']
        else:
            p_list = ['']

        for part in p_list:
#
#--- lev 1 image
#
            cmd   = 'ls ' + data_dir + hdir + '/' + head + '*_norm' + part + '.fits.gz > ' + zspace
            os.system(cmd)
#
#--- exclude the previous combined image file (with "_total_") and instmant (with "_instmap_")
#
            tlist = hcf.read_file_data(zspace, remove=1)
            slist = [x for x in tlist if 'total'   not in x]
            flist = [x for x in slist if 'instmap' not in x]
            add_image_fits_data(flist, part)
#
#--- lev 2 image
#
            cmd   = 'ls ' + data_dir + hdir + '/' + head + '_lev2*_norm' + part + '.fits.gz > ' + zspace
            os.system(cmd)
    
            tlist = hcf.read_file_data(zspace, remove=1)
            flist = [x for x in tlist if 'total' not in x]
            add_image_fits_data(flist, part)
#
#--- create map png files and histogram png files
#
    hpm.plot_hrc_map('total')

#-----------------------------------------------------------------------------------------------
#-- add_image_fits_data: combined image fits files to one                                    ---
#-----------------------------------------------------------------------------------------------

def add_image_fits_data(fits_list, part=''):
    """
    combined image fits files to one
    input:  fits_list   --- a list of fits file names
            part        --- a section number (for hrc_s_125)
    output: combined image fits file, e.g. hrc_i_115_total_norm.fits.gz
    """

    fnum = len(fits_list)
    if fnum == 0:
        return NULL
#
#--- if there is only one file, just copy the original to the new one
#
    elif fnum == 1:
        out = set_out_name(fits_list[0], part)
        cmd = 'cp ' + fits_list[0] + ' ' + out + '.gz'
        os.system(cmd)

    else:
#
#--- there are more than one fits file combine...
#
        add_two_image_fits(fits_list[0], fits_list[1])

        for k in range(2, fnum):
            cmd = 'mv zimg_out.fits zimg_in.fits'
            os.system(cmd)

            add_two_image_fits('zimg_in.fits',  fits_list[k])
            mcf.rm_file('zimg_in.fits')
#
#--- zip the output file and move to the output fits file location/name
#
        cmd = 'gzip -f zimg_out.fits'
        os.system(cmd)

        out = set_out_name(fits_list[0], part)
        cmd = 'mv zimg_out.fits.gz ' + out + '.gz'
        os.system(cmd)

#-----------------------------------------------------------------------------------------------
#-- set_out_name: create output fits file name for combined image file                       ---
#-----------------------------------------------------------------------------------------------

def set_out_name(infits, part=''):
    """
    create output fits file name for combined image file
    input:  infits  --- input file name
            part    --- a section number (for hrc_s_125)
    output: out     --- output fits file name e.g. hrc_i_115_total_norm.fits
                        with a full path
    """
#
#--- if the file is lev2 case
#
    mc = re.search('lev2', infits)
    if mc is not None:
        atemp = re.split('lev2', infits)
        out   = atemp[0] + 'lev2_total_norm' + part + '.fits'
#
#--- for the case fits is lev1
#
    else:
        for head in ['hrc_i_115', 'hrc_s_125_hi', 'hrc_s_125']:
            mc = re.search(head, infits)
            if mc is not None:
                atemp = re.split(head, infits)
                out   = atemp[0] + head +  '_total_norm' + part + '.fits'
                break

    return out

#-----------------------------------------------------------------------------------------------
#-- add_two_image_fits: combine two image fits file and take average                         ---
#-----------------------------------------------------------------------------------------------

def add_two_image_fits(fits1, fits2, outfile='zimg_out.fits'):
    """
    combine two image fits file and take average
    input:  fits1   --- an image fits file name
            fits2   --- an image fits file name
            outfile --- an output fits file name, default: zimg_out.fits
    output: outfile
    """
#
#--- open fits files
#
    f1 = pyfits.open(fits1)
    f2 = pyfits.open(fits2)
    h1 = f1[0]
    h2 = f2[0]
#
#--- combine image data and devide by 2 to take an average
#
    d1 = h1.data
    d2 = h2.data

    d3 = (d1 + d2) / 2.0
#
#--- write the result to the output file
#
    h1.data = d3
    mcf.rm_file(outfile)
    h1.writeto(outfile)
#
#--- update the fits file header dates
#
    db1 = h1.header['date-obs']
    db2 = h2.header['date-obs']
    de1 = h1.header['date-end']
    de2 = h2.header['date-end']
    ts1 = int(float(h1.header['tstart']))
    ts2 = int(float(h1.header['tstart']))
    te1 = int(float(h1.header['tstop']))
    te2 = int(float(h1.header['tstop']))

    if ts1 < ts2:
        begin  = db1
        tstart = ts1
    else:
        begin  = db2
        tstart = ts2

    if te1 < te2:
        end    = de2
        tstop  = te2
    else:
        end    = de1
        tstop  = te1

    f1.close()
    f2.close()

    ctime = hcf.covertfrom1998sec2(tcnv.currentTime(format='SEC1998'))

    hcf.update_header(outfile, 'DATE',     value=ctime,  comment='Date and time of file creation', exten=0)
    hcf.update_header(outfile, 'DATE-OBS', value=begin,  comment='TT, with clock correction if clockapp', exten=0)
    hcf.update_header(outfile, 'DATE-END', value=end,    comment='TT, with clock correction if clockapp', exten=0)
    hcf.update_header(outfile, 'TSTART',   value=tstart, comment='Data file start time', datatype='Real8', exten=0)
    hcf.update_header(outfile, 'TSTOP',    value=tstop,  comment='Data file stop time',  datatype='Real8', exten=0)


#-----------------------------------------------------------------------------------------------
 
if __name__ == "__main__":

#
#--- argument "img" is given, it just update combined images, but not process new month data
#
    if len(sys.argv) == 2:
        if sys.argv[1] == 'img':
            create_cumulative_image()
#
#--- for the case, the year and month of the data to be processed are given
#
    elif len(sys.argv) == 3:
        lyear  = int(float(sys.argv[1]))
        lmonth = int(float(sys.argv[2]))

        run_hrc_stowed_background_monthly_update(lyear, lmonth)
#
#--- if no date is given, process the last month's data
#
    else:
        run_hrc_stowed_background_monthly_update()
        fo.close()

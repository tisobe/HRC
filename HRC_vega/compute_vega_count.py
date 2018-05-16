#!/usr/bin/env /proj/sot/ska/bin/python

#############################################################################################################
#                                                                                                           #
#           compute_vega_count.py: update Vega Monitoring the UV/Ion Shield Health page                     #
#                                                                                                           #
#           author: t. isobe (tisobe@cfa.harvard.edu)                                                       #
#                                                                                                           #
#           Last Update: Mar 01, 2017                                                                       #
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
ciaoenv  = getenv('source /soft/ciao/bin/ciao.sh')

#
#--- reading directory list
#
path = '/data/aschrc6/wilton/isobe/Project8/Vega/Scripts/house_keeping/dir_list'

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

import mta_common_functions     as mcf
import convertTimeFormat        as tcnv
import find_coordinate          as fc
import extract_vega_data        as evd
import create_html_page         as chp

#
#--- temp writing file name
#
rtail  = int(time.time())
zspace = '/tmp/zspace' + str(rtail)

asec_p_pix = 0.1318         #---- arcsec per pixcel
pi         = 3.141592653

#-----------------------------------------------------------------------------------------
#-- update_vega_trend_page: update Vega Monitoring the UV/Ion Shield Health page        --
#-----------------------------------------------------------------------------------------

def update_vega_trend_page():
    """
    update Vega Monitoring the UV/Ion Shield Health page
    input:  none, but read from /data/mta4/obs_ss/sot_ocat.out
    output: <html_dir>/vega_vis_montior.html
    """
#
#--- create calibration vega observation list and check whether any new observations are added
#
    chk = evd.extract_calb_vega()
    if chk == 0:
        exit(1)
    else:
#
#--- save the older data
#
        cmd = 'ls ' + data_dir + 'hrc_*_results > ' + zspace
        os.system(cmd)
        data = read_data(zspace, remove=1)
        for ent in data:
            cmd = 'mv ' + ent + ' ' + ent + '~'
            os.system(cmd)
#
#--- if there are new observations, compute the statistics for all 
#--- this is becasue the older data are often re calibrated according to the new data
#
        table = house_keeping + 'hrc_i_list'
        create_data_tables(table)

        table = house_keeping + 'hrc_s_list'
        create_data_tables(table)
#
#--- update plots and html pages
#
        chp.create_html_and_plot()


    
#-----------------------------------------------------------------------------------------
#-- create_data_tables: read input table of observations and create data tables         --
#-----------------------------------------------------------------------------------------

def create_data_tables(table):
    """
    read input table of observations and create data tables
    input:  table   --- a file which ontains a list of observations
                        either obsid list or a list in which the first column is obsid
    output: <data_dir>/hrc_<det/pos>_results 
    """
    data = read_data(table)

    for ent in data:
        if mcf.chkNumeric(ent):
            obsid = ent.strip()
        else:
            atemp = re.split('\s+', ent)
            obsid = atemp[0]

        if mcf.chkNumeric(obsid) == False:
            continue

        print obsid

        fits  = run_arc5gl('retrieve', detector='hrc', level=2, filetype='evt2', obsid=obsid)

        out   = extract_count_stats(fits)
        if out[-1] <0:
            continue

        line  = str(obsid) + '\t'

        if float(obsid) < 1000:
            line = line + '\t'

        line  = line + str(fits)                 + '\t'
        line  = line + out[7]                    + '\t'
        line  = line + '%2.1f' % round(out[6],1) + '\t'
        line  = line + '%2.2f' % round(out[5],2) + '\t'
        line  = line + '%2.2f' % round(out[8],2) + '\t'
        line  = line + '%2.4f' % round(out[9],4) + '\n'

        if out[-1] == 10:
            outfile = data_dir + 'hrc_i_results'
        if out[-1] == 0:
            outfile = data_dir + 'hrc_s_0_results'
        if out[-1] == 1:
            outfile = data_dir + 'hrc_s_10_results'
        if out[-1] == 2:
            outfile = data_dir + 'hrc_s_25_results'
        if out[-1] == 3:
            outfile = data_dir + 'hrc_s_m10_results'
        if out[-1] == 4:
            outfile = data_dir + 'hrc_s_m25_results'

        fo = open(outfile, 'a')
        fo.write(line)
        fo.close()

        cmd = 'rm  -f ' + fits
        os.system(cmd)


#-----------------------------------------------------------------------------------------
#-- extract_count_stats: extract target area event counts                               --
#-----------------------------------------------------------------------------------------

def extract_count_stats(fits):
    """
    extract target area event counts
    input:  fits    --- event fits file
    output: c_cnt   --- event counts in the targeted area
            b_cnt   --- event counts in the background area
            b_avg   --- background event counts per pixel
            val     --- background event counts adjusted for the target area
            sval    --- sqrt of val
            u_cnt   --- background subtracted target area event counts
            expo    --- exposure time
            odate   --- observation date
            err     --- count error
            dtf     --- dead time correction
            pos     --- indicator of which type of observations
                            0   --- hrc s center
                            1   --- hrc s +10 arcmin
                            2   --- hrc s >+20 arcmin
                            3   --- hrc s -10 arcmin
                            4   --- hrc s <-20 arcmin
                            10  --- hrc i
                            -1  --- error case
    """
#
#--- get information from the fits file header
#
    [odate, obsid, expo, fyear, det, ra_pnt, dec_pnt] = get_info_from_header(fits)
#
#--- call simbad to get vega location at the observation date
#
    [ra, dec] = fc.find_coordinate('vega', fyear)
#
#--- get area setting parameters
#
    [radius, annula, annulb, pos] = get_area_param(fits, det, fyear, ra, dec, ra_pnt, dec_pnt)
#
#--- if something wrong, just return null data
#
    if pos == -1:
        return [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1]
#
#--- computer the size of the extracted areas
#
    area1  = pi * radius * radius
    area2  = pi * (annulb * annulb - annula * annula)
#
#--- create the target area and background area event fits files
#
    get_area(fits, ra, dec, radius, outfile='center_area.fits')
    get_area(fits, ra, dec, annula, annulb, outfile='bkg_area.fits')
#
#--- count events
#
    c_cnt = get_evt_cnt('center_area.fits')
    b_cnt = get_evt_cnt('bkg_area.fits')
    b_avg = b_cnt / area2

    val   = b_avg * area1
    sval  = math.sqrt(val)
    u_cnt = c_cnt - b_avg * area1
#
#--- get dtf for the observations
#
    #dsave = get_dtf(obsid, 'bkg_area.fits')
    dsave = get_dtf(obsid, 'center_area.fits')
    dtf   = numpy.mean(dsave)

    mcf.rm_file('center_area.fits')
    mcf.rm_file('bkg_area.fits')
#
#--- compute the error for the counts
#
    err = math.sqrt(c_cnt + b_avg**2 / b_cnt)

    return [c_cnt, b_cnt, b_avg, val, sval, u_cnt, expo, odate, err, dtf, pos]

#-----------------------------------------------------------------------------------------
#-- get_info_from_header: extract information from the fits header                      --
#-----------------------------------------------------------------------------------------

def get_info_from_header(fits):
    """
    extract information from the fits header
    input:  fits    --- fits file
    output: odate   --- observation date
            obsid   --- obsid
            expo    --- exposure time in seconds
            fyear   --- observation date in fractional year format
            det     --- detector name; either hrc-i or hrc-s
            ra_pnt  --- ra of pointing direction
            dec_pnt --- dec of pointing direction
    """

    hdr    = pyfits.getheader(fits, 1)
    odate  = hdr['DATE-OBS']
    obsid  = hdr['OBS_ID']
    expo   = hdr['EXPOSURE']
    expo   = float(expo)
    atemp  = re.split('-', odate)
    year   = float(atemp[0])
    mon    = float(atemp[1])
    fyear  = year + mon/12
    det    = hdr['DETNAM'].lower()
    ra_pnt = hdr['RA_PNT']
    ra_pnt = float(ra_pnt)
    dec_pnt= hdr['DEC_PNT']
    dec_pnt= float(dec_pnt)

    return [odate, obsid, expo, fyear, det, ra_pnt, dec_pnt]

#-----------------------------------------------------------------------------------------
#-- get_area_param: set area parameters                                                      --
#-----------------------------------------------------------------------------------------

def get_area_param(fits, det, fyear, ra, dec, ra_pnt, dec_pnt):
    """
    set area parameters
    input:  fits    --- fits file name
            det     --- detector name; hrc-i or hrc-s
            fyrea   --- observation date in fractional year format
            ra      --- target ra
            dec     --- target dec
            ra_pnt  --- pointing ra
            dec_pnt --- pointing dec
    output: radius  --- the radius of the center target area
            annula  --- the radius of the innter circle of annulus
            annulb  --- the radius of the outer circle of annulus
            pos     --- indicator of the case (see extract_count_stats help for more details)
    """
#
#--- convert  ra/dec to det coordinates
#
    [detx_p, dety_p] = cell_to_det(fits, ra_pnt, dec_pnt)
    [detx, dety]     = cell_to_det(fits, ra, dec)
#
#--- find the difference between Vega location and the pointing position
#
    detd = dety - dety_p
#
#--- for hrc-s case
#
    if det == 'hrc-s':
        diff = math.sqrt((ra - ra_pnt)**2 + (dec - dec_pnt)**2) * 60.0
        if diff < 5:
            pos = 0
            radius =  4.0 / asec_p_pix
            annula =  8.0 / asec_p_pix
            annulb = 12.0 / asec_p_pix
        elif diff > 7 and diff < 15:
            if detd > 0:
                pos = 1
            else:
                pos = 3
            radius = 16.0 / asec_p_pix
            annula = 32.0 / asec_p_pix
            annulb = 48.0 / asec_p_pix
        elif diff > 18:
            if detd > 0:
                pos = 2
            else:
                pos = 4
            radius =  64.0 / asec_p_pix
            annula = 100.0 / asec_p_pix
            annulb = 130.0 / asec_p_pix
        else:
            pos = -1
#
#--- for hrc-i case
#
    else:
        radius =  4.0 / asec_p_pix
        annula =  8.0 / asec_p_pix
        annulb = 12.0 / asec_p_pix
        pos    = 10

    return [radius, annula, annulb, pos]

#-----------------------------------------------------------------------------------------
#-- get_dtf: reate a list of dtf for the given time periods                            ---
#-----------------------------------------------------------------------------------------

def get_dtf(obsid, fits):
    """
    create a list of dtf for the given time periods
    input:  obsid   --- obsid of the observations
            fits    --- event fits file name
    output: dsave   --- a list of dtf values correspoinding to the event fits data
    """
#
#--- get fits time list
#
    t     = pyfits.open(fits)
    tdata = t[1].data
    t.close()

    ftime = tdata.field('time')
#
#--- get dtf time list and dtf rate list
#
    dtf   = run_arc5gl('retrieve', detector='hrc', level=1, filetype='dtf', obsid=obsid)

    t     = pyfits.open(dtf)
    tdata = t[1].data
    t.close()

    dtime = tdata.field('time')
    drate =  tdata.field('dtf')
#
#--- compare them and find dtf correspond to the time of fits data
#
    dlen  = len(dtime)
    dsave = []
    for ent in ftime:
        for k in range(0, dlen-1):
            if ent >= dtime[k] and ent < dtime[k+1]:
                dsave.append(drate[k])
                break

    mcf.rm_file(dtf)

    return dsave

#-----------------------------------------------------------------------------------------
#-- get_evt_cnt: count # of event in the fits file                                     ---
#-----------------------------------------------------------------------------------------

def get_evt_cnt(fits):
    """
    count # of event in the fits file
    input:  fits    --- event fits file name
    ouptut: cnt     --- # of event in the fits file
    """

    t     = pyfits.open(fits)
    tdata = t[1].data
    t.close()

    data  = tdata.field('time')
    cnt   = len(data)

    return cnt

#-----------------------------------------------------------------------------------------
#-- get_area: extract specified area from the event fits file and create a new fits file -
#-----------------------------------------------------------------------------------------

def get_area(fits, ra, dec, radius, annul=0, outfile='circle_area.fits'):
    """
    extract specified area from the event fits file and create a new fits file
    input:  fits    --- event fits file
            ra      --- ra of the center of the area
            dec     --- dec of the center of the area
            radius  --- a radius of the circular area or the inter radius of annulus
            annul   --- a outer radius of annuls. if it is 0, ignored.
            outfile --- a name of the created fits file
    output: outfile
    """

    [skyx, skyy] = get_sky_coords(fits, ra, dec)

    if annul == 0:
        cmd = 'dmcopy "' + fits + '[sky=circle(' + str(skyx) + ',' + str(skyy) + ',' + str(radius)
        cmd = cmd  + ')]" outfile=' + outfile + ' clobber=yes'
    else:
        cmd = 'dmcopy "' + fits + '[sky=annulus(' + str(skyx) + ',' + str(skyy) + ',' + str(radius)
        cmd = cmd  + ',' + str(annul) + ')]" outfile=' + outfile + ' clobber=yes'

    run_ascds(cmd)


#-----------------------------------------------------------------------------------------
#-- get_sky_coords: convert ra/dec to skyx/skyy                                         --
#-----------------------------------------------------------------------------------------

def get_sky_coords(fits, ra, dec):
    """
    convert ra/dec to skyx/skyy
    input:  fits    --- fits file name
            ra      --- ra
            dec     --- dec
    output: [skyx, skyy]    --- ra/dec in sky coordinates
    """

    cmd = 'dmcoords ' + fits  + ' opt=cel ' + ' ra=' + str(ra)  + ' dec=' + str(dec)
    cmd = cmd + ' verbose=1 > ' + zspace
    run_ascds(cmd)
    data = read_data(zspace, remove=1)

    for ent in data:
        mc  = re.search('SKY', ent)
        if mc is not None:
            atemp = re.split('\s+', ent)
            skyx  = atemp[-2]
            skyy  = atemp[-1]
            break

    return [skyx, skyy]



#-----------------------------------------------------------------------------------------
#-- run_arc5gl: run arc5gl                                                             ---
#-----------------------------------------------------------------------------------------

def run_arc5gl(operation, detector='hrc', subdetector='', level=1, filetype='evt1', obsid='', tstart='', tstop='', filename='', outfile='fits_list'):
    """
    run arc5gl
    input:  operation   --- retrieve/browse
            detector    --- detector
            subdetector --- subdetector
            level       --- level;      default: 1
            filetype    --- filetype;   default: evt1
            obsid       --- obsid;      default: ''
            tstart      --- start time; default: ''
            tstop       --- stop time;  default: ''
            filename    --- file name   default: ''
            outfile     --- output file name; default:  fits_list
            if the param is "", it is ignored
    output: extracted fits file
            outfile     --- a list of fits file names
            fits_list   --- returning the list of fits file name
    """

    line = 'operation=' + operation + '\n'
    line = line + 'dataset=flight\n'
    line = line + 'detector=' + detector + '\n'
    line = line + 'level='    + str(level) + '\n'
    line = line + 'filetype=' + filetype + '\n'
    if subdetector != '':
        line = line + 'subdetector=' + subdetector + '\n'
    if obsid != '':
        line = line + 'obsid='  + str(obsid)  + '\n'

    if tstart != '':
        line = line + 'tstart=' + str(tstart)  + '\n'
        line = line + 'tstop='  + str(tstop)   + '\n'
    if filename != '':
        line = line + 'filename='  + filename   + '\n'

    line = line + 'go\n'

    fo   = open(zspace, 'w')
    fo.write(line)
    fo.close()

    cmd  = '  /proj/sot/ska/bin/arc5gl  -user isobe -script ' +  zspace + ' > ' + outfile
    run_ascds(cmd, clean=0)
    mcf.rm_file(zspace)

    out = read_data(outfile)
    fits_list = []
    for ent in out:
        mc = re.search('.fits.gz', ent)
        if mc is not None:
            fits_list.append(ent)

    if len(fits_list) == 0:
        return 'na'
    elif len(fits_list) == 1:
        return fits_list[0]
    else:
        return fits_list

#-----------------------------------------------------------------------------------------
#-- read_data: read data file                                                           --
#-----------------------------------------------------------------------------------------

def read_data(infile, remove=0):

    f    = open(infile, 'r')
    data = [line.strip() for line in f.readlines()]
    f.close()

    if remove == 1:
        mcf.rm_file(infile)

    return data

#-----------------------------------------------------------------------------------------
#-- read_header: read header info                                                       --
#-----------------------------------------------------------------------------------------

def read_header(fits):
    """
    read header info
    input:  fits    --- fits file name
    output:     [target, ra_pnt, dec_pnt, roll_pnt, ra_targ, dec_targ, 
                ra_norm, dec_norm, roll_norm, detnam, obsid, date_obs, 
                date_end, tstart, tstop, grating]
    """

    hdr = pyfits.getheader(fits, 1)
    try:
        target = hdr['TARGET']
    except:
        target = hdr['OBJECT']

    ra_pnt    = hdr['RA_PNT']
    dec_pnt   = hdr['DEC_PNT']
    roll_pnt  = hdr['ROLL_PNT']
    ra_targ   = hdr['RA_TARG']
    dec_targ  = hdr['DEC_TARG']
    ra_norm   = hdr['RA_NOM']
    dec_norm  = hdr['DEC_NOM']
    roll_norn = hdr['ROLL_NOM']
    detnam    = hdr['DETNAM']
    obsid     = hdr['OBS_ID']
    date_obs  = hdr['DATE-OBS']
    date_end  = hdr['DATE-END']
    tstart    = hdr['TSTART']
    tstop     = hdr['TSTOP']
    grating   = hdr['GRATING']

    return [target, ra_pnt, dec_pnt, roll_pnt, ra_targ, dec_targ, ra_norm, dec_norm, roll_norm, detnam, obsid, date_obs, date_end, tstart, tstop, grating]

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
            bash(acmd, env=ascdsnv)
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
#-- cell_to_det: convert ra/dec to det coordinates                                      --
#-----------------------------------------------------------------------------------------

def cell_to_det(fits, ra, dec):
    """
    convert ra/dec to det coordinates
    input:  fits    --- event fits file
            ra      --- ra
            dec     --- dec
    output: detx    --- detx
            dety    --- dety
    """

    cmd = 'dmcoords ' + fits  + ' opt=cel '
    cmd = cmd + ' ra=' + str(ra)  + ' dec=' + str(dec)
    cmd = cmd + ' verbose=1 > ' + zspace
    run_ascds(cmd)
    data = read_data(zspace, remove=1)
    detx = 0
    dety = 0
    for ent in data:
        mc  = re.search('DET', ent)
        if mc is not None:
            atemp = re.split('\s+', ent)
            detx  = float(atemp[-2])
            dety  = float(atemp[-1])


    return [detx, dety]

#-----------------------------------------------------------------------------------------

if __name__ == '__main__':

    if len(sys.argv) == 1:
        update_vega_trend_page()

    elif len(sys.argv) == 2:
        table = sys.argv[1]
        table.strip()
        create_data_tables(table)
    else:
        print "please give an input table name"

#!/usr/bin/env /proj/sot/ska/bin/python

#############################################################################################
#                                                                                           #
#           create_hrc_image_files.py: create evt1 and evt2 image files                     #
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
import hrc_plotting_maps          as hpm        #---- contains plot_hrc_map tp creates exposure maps
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
#-- yearly_cummulative_image: create cumulative image fits files for the given year           --
#-----------------------------------------------------------------------------------------------

def yearly_cummulative_image(year):
    """
    create cumulative image fits files for the given year
    input:  year    --- the year you want to create the image fits files
    output: lev1, lev2, and inst map for the given year
    """
#
#--- hrc i 115 case
#
    hdir = 'hrc_i_115'
    run_map_fucntion(hdir, year) 
#
#--- hrc s 125 hi res case
#
    hdir = 'hrc_s_125_hi'
    run_map_fucntion(hdir, year) 
#
#--- hrc s 125 case: 3 sections
#
    hdir = 'hrc_s_125'
    for sec in range(1, 4):
        run_map_fucntion(hdir, year, sec)
#
#--- create histograms and maps of image files
#
    hpm.plot_hrc_map(year)

#-----------------------------------------------------------------------------------------------
#-- run_map_fucntion: run all functions to create lev1, lev2, and inst map for the given year --
#-----------------------------------------------------------------------------------------------

def run_map_fucntion(hdir, year,  sec=''):
    """
    run all functions to create lev1, lev2, and inst map for the given year
    input:  hdir    --- the name of the instrument, e.g. hrc_i_115
            year    --- year of the data
            sec     --- section name of the chip, if needed (for hrc_s_125)
    output: lev1, lev2, and instrument map fits files
    """
#
#--- set parameters for a given hrc
#
    out = set_parmeters(year, hdir, sec)

    if out == NULL:
        return NULL
    else:
        [odir, oname, lev1, lev2, instmap, sbin, factor, lev1_opt, lev2_opt] = out
#
#--- create image fits files
#
    try:
        chk = create_image_file(hdir, year, oname, lev1, lev2, instmap, sbin, factor, lev1_opt, lev2_opt, sec)

        if chk == NULL:
            print "No data for Year: " + + str(year) + ' of the instrument of: ' + str(hdir)
        else:
#
#--- for hrc_i case
#
            if hdir == 'hrc_i_115':
                cmd =  'gzip -f ' + oname + ' ' + lev1 + ' ' + lev2 + ' ' + instmap
#
#--- for hrc_s_125 sec1 and sec2
#
            elif (sec == 1) or (sec == 2):
                cmd =  'gzip -f ' +  lev1 + ' ' + lev2 
#
#--- for hrcs_125_hi and hrc_s_125 sec 3
#
            else:
                cmd =  'gzip -f ' + oname + ' ' + lev1 + ' ' + lev2 

            os.system(cmd)

    except:
        print "Something wrong or no data for Year: " + str(year) + ' of the instrument of: ' + str(hdir)

#-----------------------------------------------------------------------------------------------
#-- set_parmeters: set appropriate parameter values for a given instrument and the data year   -
#-----------------------------------------------------------------------------------------------

def set_parmeters(year, hdir, sec=''):
    """
    set appropriate parameter values for a given instrument and the data year
    input:  year    --- year of the data
            hdir    --- the name of the instrument, e.g. , hrc_i_115
            sec     --- section value (1, 2, or 3) for hrc_s_125 case   default: ""
    output: odir    --- output directory path
            oname   --- evt1 fits file name
            lev1    --- lev1 image fits file name
            lev2    --- lev2 image fits file name
            instmap --- the name of the instrument map fits file
            sbin    --- binning size (usually, 128, 265)
            factor  --- normalization factor based on exposure time and binning size
            lev1_opt--- optional command input for lev1 computation
            lev2_opt--- optional command input for lev2 computation

    """
#
#--- set image fits file names (with a full path)
#
    odir    = data_dir + hdir.replace('hrc', 'Hrc') + '/'
    oname   = odir + hdir + '_'      + str(year) + '_evt1.fits'                 #--- evt 1 file
    lev1    = odir + hdir + '_'      + str(year) + '_norm' + str(sec) +'.fits'  #--- lev 1 image file
    lev2    = odir + hdir + '_lev2_' + str(year) + '_norm' + str(sec) +'.fits'  #--- lev 2 image file
    instmap = odir + hdir + '_'      + str(year) + '_instmap_norm.fits'         #--- instrument map
#
#---  estimate total exposure times from dead time table for the year
#
    atime  = get_dead_time(hdir, year, sec=sec)

    if atime > 0.0:
#
#--- factor:    normalization factor
#--- sbin:      # of bins to use
#--- lev2_opt:  status setting for bit correction 
#
        factor = 1.0 /atime
#
#--- hrc i
#
        if hdir == 'hrc_i_115':
            sbin = 256
            factor /= (sbin * sbin)
            lev1_opt= '[events]'
            lev2_opt= '[events][status=xxxxxx00xxxx0xxx00000000xx000000]'
#
#--- hrc s; it has three secionts, but hi res case uses only middle one
#
        else:
            sbin = 128
            factor /= (sbin * sbin)

            if hdir == 'hrc_s_125_hi':
                lev1_opt= '[events][chip_id=2]'
                lev2_opt= '[events][chip_id=2][pha=0:254,status=xxxxxx00xxxx0xxx0000x000xx0000xx]'

            else:
                lev1_opt= '[events][chip_id=' + str(sec) + ']'
                lev2_opt= '[events][chip_id=' + str(sec) + '][pha=0:254,status=xxxxxx00xxxx0xxx0000x000xx0000xx]'
    else:
        return NULL
    
    return [odir, oname, lev1, lev2, instmap, sbin, factor, lev1_opt, lev2_opt]

#-----------------------------------------------------------------------------------------------
#-- create_image_file: create lev1 and lev2 image files from evt1 file                       ---
#-----------------------------------------------------------------------------------------------

def create_image_file(hdir, year,  oname, lev1, lev2, instmap, sbin, factor, lev1_opt, lev2_opt, sec=''):
    """
    create lev1 and lev2 image files from evt1 file
    input:  hdir    --- instrument name, e,g hrc_i_125
            year    --- the year that the data is processed
            oname   --- evt1 fits file name
            lev1    --- lev1 output fits file name
            lev2    --- lev2 output fits file name
            instmap --- the name of instrument map
            sbin    --- binning size (usually 128 or 256)
            factor  --- normalization factor
            lev2_opt--- option to be used to compute lev2 image file
    output: lev1, lev2, and instmap fits files
    """
#
#--- make a list of all evt1 file of the year
#
    e_list = get_evt1_list(year, hdir)
    if len(e_list) == 0:
        return NULL
#
#--- combined all evt1 files of the year
#
    if (sec == '') or (sec ==1):
        hcf.combine_fits_files(e_list, oname, azip=0)
#
#--- lev 1 image
#
    create_binned_image(oname, 'temp_img.fits', xsize=sbin, ysize=sbin, opt1=lev1_opt, opt2='image')
    os.system('cp temp_img.fits comb_temp_img.fits')
#
#--- normalized by factor
#
    img_calc('temp_img.fits', outfile=lev1, factor=factor)
#
#--- lev 2 image: bit corrected image
#
    create_binned_image(oname, 'temp_img.fits', xsize=sbin, ysize=sbin, opt1=lev2_opt, opt2='image')
#
#---  for the case of hrc i
#
    if hdir == 'hrc_i_115':
#
#--- compute a instrument map
#
        create_inst_map(oname, instmap)
#
#--- normalized with the instrument map
#
        img_calc('temp_img.fits', instmap,  outfile='temp_norm_img.fits', method='div')
#
#--- normalized by the factor
#
        img_calc('temp_norm_img.fits', outfile=lev2, factor=factor)
#
#--- for the case of hrc s
#
    else:
        img_calc('temp_img.fits', outfile=lev2, factor=factor)
#
#--- clean up!
#
    cmd = 'rm ' + exc_dir + '*.fits'
    os.system(cmd)

#-----------------------------------------------------------------------------------------------
#-- create_inst_map: create normalized instrument map for a given evt1 file                  ---
#-----------------------------------------------------------------------------------------------

def create_inst_map(oname, instmap):
    """
    create normalized instrument map for a given evt1 file
    input:  oname   --- evt1 fits file
            instmap --- the name of instrument map fits file
    output: instmap --- the normalized instrument map fits file
    """
#
#--- create instrument map
#
    run_mkinstmap(oname, 'instmap.fits')
#
#--- bin the instrument map to 256x256
#
    create_binned_image('instmap.fits', 'instmap_binned.fits', xs=1, xe=16384, ys=1, ye=16384)
#
#--- normalize the instrument map
#
    smax  = find_img_max('instmap_binned.fits')
    div   = 1.0 / smax

    img_calc('instmap_binned.fits', outfile=instmap, factor=div)

    cmd = 'rm instmap.fits  instmap_binned.fits'
    os.system(cmd)

#-----------------------------------------------------------------------------------------------
#-- run_mkinstmap: function to run mkinstmap                                                  --
#-----------------------------------------------------------------------------------------------

def run_mkinstmap(in_name, out_name, bsize = 16384):
    """
    function to run mkinstmap
    input:  in_name     --- input evt1 fits file name
            out_name    --- instrument map fits file name
    output: out_name    --- instrument map fits file
    """

    cmd1 = '/usr/bin/env PERL5LIB=""'

    cmd2 = ' punlearn ardlib;'
    cmd2 = cmd2 + ' pset hrc_process_events badfile=lev1_bad_evts.fits;'
    cmd2 = cmd2 + ' pset ardlib AXAF_HRC-I_BADPIX_FILE="' + house_keeping + 'hrci_bpix1.fits";' 
    cmd2 = cmd2 + ' pset ardlib AXAF_HRC-S_BADPIX_FILE="' + house_keeping + 'hrcs_bpix1.fits";' 
    cmd2 = cmd2 + ' mkinstmap obsfile=' + in_name 
    cmd2 = cmd2 + ' pixelgrid=1:16384:#' + str(bsize) +',1:16384:#' + str(bsize) 
    cmd2 = cmd2 + ' detsubsys="HRC-I;IDEAL"  outfile='    + out_name 
    cmd2 = cmd2 + ' mirror="HRMA;AREA=1" spectrumfile=NONE monoenergy=1.0  '
    cmd2 = cmd2 + ' grating=NONE maskfile=NONE clobber=yes'

    cmd  = cmd1 + cmd2

    bash(cmd, env=ascdsenv)

#-----------------------------------------------------------------------------------------------
#-- get_dead_time: find a total exposure time for the year of the instument from dead_time list 
#-----------------------------------------------------------------------------------------------

def get_dead_time(hdir, year, sec=''):
    """
    find a total exposure time for the year of the instument from dead_time list
    input:  hdir    --- the name of the instrument, e.g., hrc_i_115
            year    --- year of the data
            sec     --- section of the chip; used only in hrc_s_125
    output: asu     --- total exposure time in seconds
    """

    if sec == '':
        efile = stat_dir + hdir + '_dead_time'
    else:
        efile = stat_dir + hdir + '_sec' + str(sec)  + '_dead_time'

    data  = hcf.read_file_data(efile)
#
#--- set the time span to be Jan 1 to Dec 31
#
    start = str(year)   + ':001:00:00:00'
    start = tcnv.axTimeMTA(start)

    stop  = str(year+1) + ':001:00:00:00'
    stop  = tcnv.axTimeMTA(stop)
#
#--- add exposure time of the year
#
    asum  = 0.0
    for ent in data:
        atemp = re.split('\s+', ent)
        if atemp[0] == 'inf':
            continue

        stime = float(atemp[0])
        if stime < start:
            continue
        elif stime > stop:
            break
        try:
            val = float(atemp[4])
        except:
            continue

        asum += val

    return asum

#-----------------------------------------------------------------------------------------------
#-- get_evt1_list: make a list of evt1 for year and a specified instrument                    --
#-----------------------------------------------------------------------------------------------

def get_evt1_list(year, hdir):
    """
    make a list of evt1 for year and a specified instrument
    input:  year    --- year of the data
            hdir    --- instrument name, e.g., hrc_i_115
    output: out     --- a list of the names of evt1 files with a full path
    """

    cmd = 'ls ' + data_dir + str(year) + '*/' + hdir + '/hrcf*evt1.fits* > ' + zspace
    os.system(cmd)

    out = hcf.read_file_data(zspace, remove=1)

    return out

#-----------------------------------------------------------------------------------------------
#-- create_binned_image: function to run dmtool: dmcopy tailored to make an image file      ----
#-----------------------------------------------------------------------------------------------

def create_binned_image(image_in, image_out, xs='', xe='', xsize=256, ys='', ye='', ysize=256, opt1='', opt2=''):
    """
    function to run dmtool: dmcopy tailored to make an image file
    input:  image_in        --- input fits file
            image_out       --- output fits file
            xs              --- starting x coordinate
            xe              --- ending x coordinate
            xsize           --- binning size            defailt: 256
            ys              --- starting y coordinate
            ye              --- ending y coordinate
            ysize           --- binning size            defailt: 256
            opt1            --- option 1
            opt2            --- option 2
    output: image_out
    """

    cmd1 = '/usr/bin/env PERL5LIB=""'

    cmd2 = ' dmcopy "' + image_in 

    if opt1 != '':
        cmd2 = cmd2 + opt1

    cmd2 = cmd2 + '[bin chipx=' + str(xs) + ':' + str(xe) + ':' + str(xsize) + ', ' 
    cmd2 = cmd2 +     ' chipy=' + str(ys) + ':' + str(ye) + ':' + str(ysize) + ']" ' 
    cmd2 = cmd2 + ' outfile= '  + image_out 

    if opt2 != '':
        cmd2 = cmd2 + ' opt='  + opt2

    cmd2 = cmd2 + ' clobber=yes'

    cmd  = cmd1 + cmd2
    bash(cmd, env=ascdsenv)

#-----------------------------------------------------------------------------------------------
#-- img_calc: fuction to run dm tool dmimgcalc                                                --
#-----------------------------------------------------------------------------------------------

def img_calc(image1, image2='None', outfile='temp.fits', factor=1, method='add'):
    """
    fuction to run dm tool dmimgcalc
    input:  image1      --- input image fits file 1
            image2      --- input image fits file 2,        default: "None"
            outfile     --- output image fits file name,    default: "temp.fits"
            factor      --- a scaling factor                default:  1
            method      --- operation: add, div, etc        default: add
    output: <outfile>
    """

    cmd1 = '/usr/bin/env PERL5LIB=""'

    cmd2 = ' dmimgcalc ' + image1 + ' '  + image2 + ' '  + outfile + ' ' + method

    if factor != 1:
        cmd2 = cmd2 + ' weight=' + str(factor)

    cmd2 = cmd2 + ' clobber=yes'

    cmd  = cmd1 + cmd2

    bash(cmd, env=ascdsenv)

#-----------------------------------------------------------------------------------------------
#-- img_normalized: normalize fits image file to 1.                                           --
#-----------------------------------------------------------------------------------------------

def img_normalized(image_in, image_out):
    """
    normalize fits image file to 1.
    input:  image_in    --- input fits file name
            image_out   --- output fits file name
    output: image_out   --- normalized fits file
    """
    smax  = find_img_max(image_in)
    scale = 1.0 / smax

    img_scale(scale, image_in, image_out)

#-----------------------------------------------------------------------------------------------
#-- find_img_max: find the max count of the image fits file                                   --
#-----------------------------------------------------------------------------------------------

def find_img_max(image):
    """
    find the max count of the image fits file
    input:  image   --- the image fits file name
    output: smax    --- the max count of the file
    """

    hdu      = pyfits.open(image)
    scidata  = hdu[0].data
    smax     = numpy.max(scidata)
    hdu.close()

    return smax

#-----------------------------------------------------------------------------------------------
#-- img_scale: scale the counts of the image file by "scale"                                  --
#-----------------------------------------------------------------------------------------------

def img_scale(scale, image_in, image_out):
    """
    scale the counts of the image file by "scale"
    input:  scale   --- scaling factor (multiple)
            image_in    --- input fits file name
            image_out   --- output fits file name
    output: image_out   --- the scaled fits file 
    """
    hdu    = pyfits.open(image_in)
    hdata  = hdu[0].data
    hdata *= scale

    hdu[0].data = hdata
    hdu[0].scale('float')

    hdu.writeto(image_out)
    hdu.close()


#-----------------------------------------------------------------------------------------
#-- TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST    ---
#-----------------------------------------------------------------------------------------

class TestFunctions(unittest.TestCase):
    """
    testing functions
    """

#------------------------------------------------------------

    def test_create_binned_image(self):
        image_in = house_keeping + 'Test_data/hrcf383576042N002_evt1.fits.gz'
        image_out = 'test_img.fits'
        create_binned_image(image_in, image_out, opt1='[events]', opt2='image')

#------------------------------------------------------------

####    def test_find_img_max(self):

        smax = find_img_max('./test_img.fits')

        self.assertEquals(smax, 30)

#------------------------------------------------------------

####    def test_img_scale(self):

        scale = 2
        image_in  = './test_img.fits'
        image_out = './test_scaled.fits'
        img_scale(scale, image_in, image_out)

        smax = find_img_max('./test_scaled.fits')
        mcf.rm_file('./test_scaled.fits')

        self.assertEquals(smax, 60)

#------------------------------------------------------------

####    def test_img_calc(self):

        image1 = './test_img.fits'
        image2 = 'None'
        out_image = 'factored_image.fits'
        factor = 1.2e-5
        method = 'div'
        img_calc(image1, image2, out_image, factor=factor, method=method)

        smax = find_img_max('./factored_image.fits')
        mcf.rm_file('./test_scaled.fits')

        smax = round(smax, 5)
        self.assertEquals(smax, 0.00036)


        mcf.rm_file('./test_img.fits')

#------------------------------------------------------------

    def test_get_evt1_list(self):

        eout = '/data/aschrc6/wilton/isobe/Project1_new/Exc/Data/2000APR/hrc_i_115/hrcf072990712N004_evt1.fits.gz'
        year = 2000
        hdir ='hrc_i_115'
        tout = get_evt1_list(year, hdir)

        self.assertEquals(tout[0], eout)

#------------------------------------------------------------

    def test_set_parmeters(self):
       
        comp = '[events][status=xxxxxx00xxxx0xxx00000000xx000000]'

        hdir = 'hrc_i_115'
        year = 2000
        out = set_parmeters(year, hdir, sec='')
        self.assertEquals(out[-1], comp)

#------------------------------------------------------------

    def test_get_dead_time(self):

        hdir = 'hrc_i_115'
        year = 2000
        atime = get_dead_time(hdir, year)
        #print 'DEAD TIME: ' + str(atime)
        self.assertEquals(atime, 108088.6)

#------------------------------------------------------------

    def test_run_mkinstmap(self):

        in_name  = house_keeping + 'Test_data/hrcf383576042N002_evt1.fits.gz'
        out_name = 'instmap.fits'
        run_mkinstmap(in_name, out_name, bsize=256)

        #in_name = '/data/aschrc6/wilton/isobe/Project1_new/Exc/Data/Hrc_i_115/hrc_i_115_2000_evt1.fits'
        in_name = '/data/aschrc6/wilton/isobe/Project1/Hrc_i_115/hrc_i_115_2000_evt1.fits.gz'
        out_name = 'test_inst.fits'
        run_mkinstmap(in_name, out_name)

        smax = find_img_max(out_name)
        self.assertEquals(smax, 1)

        mcf.rm_file(out_name)

#------------------------------------------------------------

    def test_create_binned_image(self):

        image_in = house_keeping + 'Test_data/hrcf383576042N002_evt1.fits.gz'
        image_out = './test_binned.fits'
        sbin = 256
        factor = 1
        factor /= (sbin * sbin)
        lev2_opt= '[events][status=xxxxxx00xxxx0xxx00000000xx000000]'

        create_binned_image(image_in, image_out, xsize=sbin, ysize=sbin, opt1=lev2_opt)

        smax  = find_img_max(image_out)
        self.assertEquals(smax, 25)

        mcf.rm_file(image_out)

#------------------------------------------------------------

    def test_create_inst_map(self):

            oname = hrc_i_dir + 'hrc_i_115_2000_evt1.fits.gz'
            instmap = 'inst_map.fits'
            create_inst_map(oname, instmap)

#------------------------------------------------------------

    def test_combine_evt1_files(self):

#
#--- make a list of all evt1 file of the year
#
        year = 2000
        hdir = 'hrc_i_115'
        e_list = get_evt1_list(year, hdir)
        if len(e_list) == 0:
            print "Something wrong in getting files"
#
#--- combined all evt1 files of the year
#
        oname = 'test_combined_evt1.fits'
        hcf.combine_fits_files(e_list, oname)

        tstart = hcf.read_header_value(oname, 'tstart')
        tstop  = hcf.read_header_value(oname, 'tstop')

        tstart = int(float(tstart))
        tstop  = int(float(tstop))

        self.assertEquals(tstart, 65961070)
        self.assertEquals(tstop,  93190294)

        mcf.rm_file(oname)


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
            year = int(float(sys.argv[1]))

    else:
        exit(1)

    if test == 0:
        yearly_cummulative_image(year)
    else:
        unittest.main()


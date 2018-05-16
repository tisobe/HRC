#!/usr/bin/env /proj/sot/ska/bin/python

#############################################################################################################
#                                                                                                           #
#           extract_vega_data.py: extract vega calibration data                                             #
#                                                                                                           #
#           author: t. isobe (tisobe@cfa.harvard.edu)                                                       #
#                                                                                                           #
#           Last Update: Feb 28, 2017                                                                       #
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

#
#--- temp writing file name
#
rtail  = int(time.time())
zspace = '/tmp/zspace' + str(rtail)

asec_p_pix = 0.1318         #---- arcsec per pixcel
pi         = 3.141592653

#-----------------------------------------------------------------------------------------
#-- extract_calb_vega: create a list of calibration observations of Vega from the database 
#-----------------------------------------------------------------------------------------

def extract_calb_vega():
    """
    create a list of calibration observations of Vega from the database
    input:  none
    output: <house_keeping>/hrc_i   --- a list of vega observation on hrc i
            <house_keeping>/hrc_s   --- a list of vega observation on hrc s
            chk                     --- if htere are new entries, return 1, otherwise 0
    """
#
#--- read the past data
#
    ifile_i = house_keeping + 'hrc_i_list'
    hrc_i_p = read_data(ifile_i)
    ifile_s = house_keeping + 'hrc_s_list'
    hrc_s_p = read_data(ifile_s)
#
#--- read database
#
    data = read_data('/data/mta4/obs_ss/sot_ocat.out')
#
#--- extract Vega calibration data
#
    hrc_i = []
    hrc_s = []
    for ent in data:
        mc1  = re.search('VEGA',        ent)
        mc1a = re.search('Vega',        ent)
        mc2  = re.search('CAL',        ent)
        mc3  = re.search('archived',    ent)
        mc4  = re.search('HRC-I',       ent)
        mc5  = re.search('HRC-S',       ent)
        if (mc1 is not None) or (mc1a is not None):
            if mc2 is not None:
                if mc3 is not None:
                    atemp = re.split('\^', ent)
                    obsid = atemp[1].strip()
                    obsid = int(float(obsid))
                    if mc4 is not None:
                        hrc_i.append(obsid)
                    elif mc5 is not None:
                        hrc_s.append(obsid)
#
#--- check whether there are any new vega calibration data
#
    hrc_i = set(hrc_i)
    hrc_s = set(hrc_s)
    chk = 0
    if hrc_i != set(hrc_i_p):
        chk = 1
    if hrc_s != set(hrc_s_p):
        chk = 1
#
#--- if so, update the list
#
    if chk > 0:
        fo = open(ifile_i, 'w')
        for ent in list(hrc_i):
            fo.write(str(ent))
            fo.write('\n')
        fo.close()
    
        fo = open(ifile_s, 'w')
        for ent in list(hrc_s):
            fo.write(str(ent))
            fo.write('\n')
        fo.close()

    return chk

#-----------------------------------------------------------------------------------------
#-- read_data: read data file                                                           --
#-----------------------------------------------------------------------------------------

def read_data(infile, remove=0):

    try:
        f    = open(infile, 'r')
        data = [line.strip() for line in f.readlines()]
        f.close()
    
        if remove == 1:
            mcf.rm_file(infile)
    
        return data
    except:
        return []


#-----------------------------------------------------------------------------------------

if __name__ == "__main__":

    extract_calb_vega()

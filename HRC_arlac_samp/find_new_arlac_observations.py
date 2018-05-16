#!/usr/bin/env /proj/sot/ska/bin/python

#############################################################################################################
#                                                                                                           #
#           find_new_arlac_observations.py: find new ar lac observations                                    #
#                                                                                                           #
#           author: t. isobe (tisobe@cfa.harvard.edu)                                                       #
#                                                                                                           #
#           Last Update: May 25, 2017                                                                       #
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
path = '/data/aschrc6/wilton/isobe/Project8/ArLac/Scripts3/house_keeping/dir_list'

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
#
#--- temp writing file name
#
rtail  = int(time.time())
zspace = '/tmp/zspace' + str(rtail)

#-----------------------------------------------------------------------------------------
#-- find_new_arlac_observations: create a list of calibration observations of ArLac from the database 
#-----------------------------------------------------------------------------------------

def find_new_arlac_observations(ret='0'):
    """
    create a list of calibration observations of ArLac from the database
    input:  none
    output: <house_keeping>/hrc_i   --- a list of arlac observation on hrc i
            <house_keeping>/hrc_s   --- a list of arlac observation on hrc s
            [<hrc_i list>, <hrc_s list>] --- a list of lists of  new hrc i and hrc s obsids
    """
#
#--- read the past data
#
    ifile_i = house_keeping + 'hrc_i_list'
    htemp   = read_data(ifile_i)
    hrc_i_p = []
#
#--- convert them into int
#
    for ent in htemp:
        hrc_i_p.append(int(float(ent)))

    ifile_s = house_keeping + 'hrc_s_list'
    htemp   = read_data(ifile_s)
    hrc_s_p = []
    for ent in htemp:
        hrc_s_p.append(int(float(ent)))

    hrc_i_p = set(hrc_i_p)
    hrc_s_p = set(hrc_s_p)
#
#--- read database
#
    data = read_data('/data/mta4/obs_ss/sot_ocat.out')
#
#--- extract ArLac calibration data
#
    hrc_i = []
    hrc_s = []
    for ent in data:
        mc1  = re.search('ARLAC',     ent)
        mc1a = re.search('ArLac',     ent)
        mc1b = re.search('AR LAC',    ent)
        mc2  = re.search('CAL',       ent)
        mc3  = re.search('archived',  ent)
        mc4  = re.search('HRC-I',     ent)
        mc5  = re.search('HRC-S',     ent)
        if (mc1 is not None) or (mc1a is not None) or (mc1b is not None):
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
#--- check whether there are any new arlac calibration data
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

    if chk == 0:
        return [[], []]
    else:
        diff_i = hrc_i - hrc_i_p
        diff_s = hrc_s - hrc_s_p

        if len(diff_i) > 0:
            ofile = house_keeping + 'hrci_new'
            fo    = open(ofile, 'w')
            for ent in diff_i:
                line = str(ent) + '\n'
                fo.write(line)
            fo.close()

        if len(diff_s) > 0:
            ofile = house_keeping + 'hrcs_new'
            fo    = open(ofile, 'w')
            for ent in diff_s:
                line = str(ent) + '\n'
                fo.write(line)
            fo.close()

        if ret == 1:
            return [list(diff_i), list(diff_s)]
        

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

    if len(sys.argv) == 2:
        ret = sys.argv[1]
        [hrci_list, hrcs_list] = find_new_arlac_observations(ret)
    
        if len(hrci_list) > 0:
            print 'HRC I NEW ArLac Observations:'
            for ent in hrci_list:
                print ent
     
        if len(hrcs_list) > 0:
            print 'HRC S NEW ArLac Observations:'
            for ent in hrcs_list:
                print ent
    else:
        find_new_arlac_observations()
    

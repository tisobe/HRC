#!/usr/bin/env /proj/sot/ska/bin/python

#############################################################################################################
#                                                                                                           #
#           create_hz43_center_location_list.py: extract coordinates from manually process data list        #
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

ascdsenv = getenv('source /home/ascds/.ascrc -r release; source /home/mta/bin/reset_param ', shell='tcsh')
ciaoenv  = getenv('source /soft/ciao/bin/ciao.sh')


#
#--- reading directory list
#
path = '/data/aschrc6/wilton/isobe/Project8/HZ43/Scripts/house_keeping/dir_list'

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
import filter_evt_file      as fef
import find_hz43_data       as fhd
import hrc_common_functions as hcf
#
#--- temp writing file name
#
rtail  = int(time.time())
zspace = '/tmp/zspace' + str(rtail)

NULL   = 'NULL'

tpath  = '/data/hrc/HAT/20170510.hz43_monitor_qe/'
dlist  = ['hrci_hz43_zero.tar', 'hrcs_hz43_zero.tar']
olist  = ['hrc_i_coords', 'hrc_s_coords']

#-----------------------------------------------------------------------------------------
#-- create_hz43_center_location_list: extract coordinates from manually process data list 
#-----------------------------------------------------------------------------------------

def create_hz43_center_location_list():
    """
    extract coordinates from manually process data list
    input:  none but read from 'hrci_hz43_zero.tar, hrcs_hz43_zero.tar
    output: <house_keeping>/hrc_i_coords, <hosue_keeping>/hrc_s_coords
            both has <obsid> <ra> <dec> <skyx> <skyy> 
    """

    for k in range(0, 2):
#
#--- find out which obsids are already processed
#
        cfile = house_keeping + olist[k]
        f     = open(cfile, 'r')
        data  = [line.strip() for line in f.readlines()]
        f.close()
        obs_list = []
        for ent in data:
            atemp = re.split('\s+', ent)
            obs_list.append(atemp[0])

        fo  = open('./zdata', 'w')
#
#--- check which data are available
#
        cmd  = 'cp ' + tpath + dlist[k] + ' . '
        os.system(cmd)
        cmd  = 'tar xf ' + dlist[k]
        os.system(cmd)
        cmd  = 'ls src*.reg > ' + zspace
        os.system(cmd)

        data = read_data(zspace, remove=1)

        chk  = 0
        for ent in data:
            mc = re.search('src_zero_', ent)
            if mc is not None:
                cut = 'src_zero_'
            else:
                cut = 'src'
            atemp = re.split(cut, ent)
            btemp = re.split('.reg', atemp[1])
            obsid = btemp[0]
#
#--- check this obsid is in the saved data list
#
            if obsid in obs_list:
                continue

            f     = open(ent, 'r')
            out   = f.read()
            f.close()
            atemp = re.split('\(', out)
            btemp = re.split(',', atemp[1])
            ra    = btemp[0]
            dec   = btemp[1]

#
#--- extract evt1 data
#
            evt1 = hcf.run_arc5gl(0, 0,  obsid = obsid, operation='retrieve', level ='1',   filetype='evt1')
            if evt1 == "":
               continue
#
#--- convert coordinates from cel to sky
#
            cmd = 'dmcoords ' + evt1 + ' opt=cel ra=' + ra + ' dec=' + dec + ' verbose=1 > ' + zspace
            run_ciao(cmd)
#
#--- extract sky coordindats
#
            info = read_data(zspace, remove=1)
            atemp = re.split('\s+', info[-2])
            skyx  = atemp[2]
            skyy  = atemp[3]
            
            line  = obsid + '\t' + ra + '\t' + dec + '\t' + skyx + '\t' + skyy + '\n'
            fo.write(line)
            chk   = 1
            mcf.rm_file(evt1)

        fo.close()
        if chk > 0:
            cmd = 'cat ./zdata >> ' + cfile
            os.system(cmd)


        cmd = 'rm -rf *.reg *.tar ./zdata'
        os.system(cmd)

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

    f= open(infile, 'r')
    data = [line.strip() for line in f.readlines()]
    f.close()
    
    if remove == 1:
        mcf.rm_file(infile)

    return data



#-----------------------------------------------------------------------------------------

if __name__ == '__main__':

    create_hz43_center_location_list()


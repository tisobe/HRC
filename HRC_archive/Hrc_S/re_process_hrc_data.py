#!/usr/bin/env /proj/sot/ska/bin/python

#################################################################################################
#                                                                                               #
#           re_process_hrc_data.py: a control script to run reprocess csh scripts               #
#                                                                                               #
#           author: t. isobe (tisobe@cfa.harvard.edu)                                           #
#                                                                                               #
#           Last Update: Apr 09, 2018                                                           #
#                                                                                               #
#################################################################################################

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
#
#--- set ciao environment 
#
ciaoenv  = getenv('source /soft/ciao/bin/ciao.csh; source /home/mta/bin/reset_param; setenv PFILES "${PDIRS}"; set path=(/soft/ciao/bin/ $path);', shell='tcsh')
#
#--- reading directory list
#
path = '/data/aschrc6/wilton/isobe/Project9/Scripts/Hrc_S/house_keeping/dir_list'
#path = '/data/aschrc6/wilton/isobe/Project9/Scripts/Hrc_I/house_keeping/dir_list'

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

import OcatSQL                  as sql
from   OcatSQL                  import OcatDB
#
#--- temp writing file name
#
rtail  = int(time.time())
zspace = '/tmp/zspace' + str(rtail)

#-----------------------------------------------------------------------------------------
#-- run_process: a control script to run reprocess csh scripts                         ---
#-----------------------------------------------------------------------------------------

def run_process(hrc):
    """
    a control script to run reprocess csh scripts 
    input:  hrc: either "hrc_i" or "hrc_s"
    output: hrc_i_list  --- a list of hrc i obsids which need to be re-processed
            hrc_s_list  --- a list of hrc s obsids which need to be re-processed
            <data_dir>/<obsid>    --- re-processed data direcotry
    """
#
#--- set conditions for either hrc-i or hrc s
#
    if hrc == 'hrc_i':
        out_list = 'hrc_i_list'
        data_dir = '/data/hrc/i/'
        inst     = 'i'
    else:
        out_list = 'hrc_s_list'
        data_dir = '/data/hrc/s/'
        inst     = 's'
#
#--- find un process data 
#
    hlist = find_un_processed_data(inst)

    if hrc == 'hrc_i':
        print "HRC I : " + str(hlist)
    else:
        print "HRC S : " + str(hlist)
    
    for obsid in hlist:
        fo  = open(out_list, 'w')
        fo.write(str(obsid))
        fo.write('\n')
        fo.close()
#
#--- extract fits data needed for analysis
#
        extract_hrc_data(obsid, data_dir)

        if hrc == 'hrc_i':
            cmd = 'csh -f ' + bin_dir + 'repro_all_new.csh hrc_i_list'
        else:
            cmd = 'csh -f ' + bin_dir + 'repro_all_S_new.csh hrc_s_list'

        try:
            run_ciao(cmd)
            cdir = data_dir + '/' + str(obsid)
            if os.path.isdir(cdir):
                cmd = 'chgrp -R hat ' + cdir 
                os.system(cmd)
                cmd = 'chmod -R 774 ' + cdir 
                os.system(cmd)
        except:
            pass

        rm_file(out_list)

#-----------------------------------------------------------------------------------------
#-- find_un_processed_data: find hrc obsids which need to be reprocessed                --
#-----------------------------------------------------------------------------------------

def find_un_processed_data(inst):
    """
    find hrc obsids which need to be reprocessed
    input: inst     --- insturment indicator "i" or "s"
    output: uhrc    --- a list of obsids of either hrc i or hrc s
    """
#
#--- extract all hrc obsid listed in database
#
    infile = '/data/mta4/obs_ss/sot_ocat.out'
    data   = read_data(infile)

    h_list = []
    h_dict = {}
    for ent in data:
        atemp = re.split('\^', ent)
        if inst == 'i':
            mc = re.search('HRC-I', atemp[12])
        else:
            mc = re.search('HRC-S', atemp[12])

        if mc is not None:
            atemp = re.split('\^\s+', ent)
            atemp[0].strip()
            try:
                val   = int(float(atemp[1]))
            except:
                continue
            h_list.append(val)
            h_dict[val] = check_status(ent)

        else:
            continue

    uhrc = clean_the_list(h_list, h_dict, inst)

    return uhrc

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

def check_status(line):

    mc1 = re.search('archived',   line)
    mc2 = re.search('observed',   line)
    mc3 = re.search('scheduled',  line)
    mc4 = re.search('unobserved', line)
    mc5 = re.search('canceled',   line)
    mc6 = re.search('discarded',  line)
    if mc1 is not None:
        return 'archived'
    elif (mc2 is not None) and (mc4 is None):
        return 'observed'
    elif mc3 is not None:
        return 'scheduled'
    elif mc4 is not None:
        return 'unobserved'
    elif mc5 is not None:
        return 'canceled'
    elif mc6 is not None:
        return 'discarded'
    else:
        return 'tbd'

#-----------------------------------------------------------------------------------------
#-- find_processed_data: find the hrc obsids which are already re-processed             --
#-----------------------------------------------------------------------------------------

def find_processed_data(inst):
    """
    find the hrc obsids which are already re-processed
    input:  inst    --- instrument designation: "i" or "s"
    output: out     --- a list of obsids
    """
    if inst == 'i':
        data_dir = '/data/hrc/i/'
    else:
        data_dir = '/data/hrc/s/'

    cmd  = 'ls -d ' + data_dir + '/* > ' + zspace
    os.system(cmd)
    data = read_data(zspace, remove=1)

    out = []
    for ent in data:
        atemp = re.split('\/', ent)
        try:
            val   = int(float(atemp[-1]))
        except:
            continue
        if chkNumeric(val):
            out.append(val)
#
#--- remove duplicate
#
    oset = set(out)
    out  = list(oset)

    return out

#-----------------------------------------------------------------------------------------
#-- clean_the_list: select out obsids which need re-process                            ---
#-----------------------------------------------------------------------------------------

def clean_the_list(current, cdict,  inst):
    """
    select out obsids which need re-process
    input:  current --- a list of all hrc obsid found in database
            inst    --- instrument designation; "i" or "s"
    output: good    --- a list of obsids to be re-processed
            <house_keeping>/cancelled_list  --- a list of observations canceleld or discarded
    """
#
#--- read the past cancelled obsids
#
    rfile  = house_keeping + 'cancelled_list'
    remove = set(read_data(rfile, num =1))
#
#--- find obsids already re-processed
#
    phrc = find_processed_data(inst)
    uhrc = set(current) - set(phrc)
    uhrc = list(uhrc - remove)
#
#--- select out obsids which need to be reprocessed
#
    good = []
    bad  = []
    for obsid in uhrc:
        try:
            status = find_status(obsid)
        except:
            try:
                status = cdict[obsid]
            except:
                status = 'tbd'

        if status in ['archived', 'observed']:
        #if status == 'archived':
            good.append(obsid)

        elif status in ['unobserved', 'scheduled']:
        #elif status in ['unobserved', 'scheduled', 'observed']:
            continue

        elif status in ['canceled', 'discarded']:
            bad.append(obsid)

        else:
            continue
#
#--- update cancelled_list if there are new cancelled observations
#
    sbad = set(bad)
    nbad = (set(bad) - remove)

    if len(bad) > 0:
        ncancel = list(remove) + list(nbad)
        fo = open(rfile, 'a')
        for ent in  ncancel:
            fo.write(str(ent))
            fo.write('\n')

        fo.close()

    return good

#-----------------------------------------------------------------------------------------
#-- find_status: find the status of the observations from the database                  --
#-----------------------------------------------------------------------------------------

def find_status(obsid):
    """
    find the status of the observations from the database
    input:  obsid   --- obsid
    output: status  --- status
    """
    try:
        dbase  = OcatDB(obsid)
        status = dbase.origValue('status')
        return status
    except:
        return ""

#-----------------------------------------------------------------------------------------
#-- chkNumeric: check the entry is numeric. If so return True, else False              ---
#-----------------------------------------------------------------------------------------

def chkNumeric(elm):

    """
    check the entry is numeric. If so return True, else False.
    """
    
    try:
        test = float(elm)
    except:
        return False
    else:
        return True

#-----------------------------------------------------------------------------------------
#-- read_data: read data file                                                           --
#-----------------------------------------------------------------------------------------

def read_data(infile, remove=0, num =0):

    f    = open(infile, 'r')
    data = [line.strip() for line in f.readlines()]
    f.close()

    if num == 1:
        temp = []
        for ent in data:
            temp.append(int(float(ent)))
        data = temp

    if remove == 1:
        rm_file(infile)

    return data

#-----------------------------------------------------------------------------------------
#-- rm_file: check whether the file exists and if it does, remove it                    --
#-----------------------------------------------------------------------------------------

def rm_file(file):
    """
    remove file
    Input:  file --- a name of file to be removed
    Output: none
    """
    chk = chkFile(file)
    if chk > 0:
        cmd = 'rm -rf ' + file
        os.system(cmd)

#------------------------------------------------------------------------------------------
#-- chkFile: check whether a file/directory exits in the directory given                ---
#------------------------------------------------------------------------------------------

def chkFile(inline, name = 'NA'):

    """
    check whether a file/directory exits in the directory given, 
    Input: a file/directory name with a full path   or a directory path and a file/directory name
    """
#
#--- if the second element is not given, assume that the first element contains 
#--- a full path and file/directory name
#
    if name == 'NA':
        cmd =  inline
    else:
        cmd = inline + '/' + name
    
    chk  = os.path.isfile(cmd)
    chk2 = os.path.isdir(cmd)
    
    if (chk == True) or (chk2 == True):
        return 1
    else:
        return 0

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

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

def extract_hrc_data(obsid, data_dir):
    """
    extract fits data needed for analysis
    input:  obsid       --- obsid
            data_dir    --- data directory
    output: <obsid>/primary/*fits etc
    """
#
#--- extract fits data
#
    line = 'operation=retrieve\n'
    line = line + 'dataset=flight\n'
    line = line + 'level=1\n'
    line = line + 'detector=hrc\n'
    line = line + 'obsid=' + str(obsid) + '\n'
    line = line + 'go\n'

    fo   = open('zline', 'w')
    fo.write(line)
    fo.close()

    cmd  = ' /proj/sot/ska/bin/arc5gl  -user isobe -script zline > zout'
    os.system(cmd)
#
#--- create directories and move the data into them
#
    cmd  = 'mkdir primary secondary'
    os.system(cmd)

    cmd  = 'mv *dtf1*fits* *fov*fits* ./primary/.'
    os.system(cmd)

    cmd  = 'mv *bpix1*fits* *evt1*fits* *msk1*fits* *mtl1*fits* *std_dtfstat1.fits* *std_flt1.fits* ./secondary/.'
    os.system(cmd)


    line = 'operation=retrieve\n'
    line = line + 'dataset=flight\n'
    line = line + 'level=1\n'
    line = line + 'detector=pcad\n'
    line = line + 'subdetector=aca\n'
    line = line + 'obsid=' + str(obsid) + '\n'
    line = line + 'go\n'

    fo   = open('zline', 'w')
    fo.write(line)
    fo.close()

    cmd  = ' /proj/sot/ska/bin/arc5gl  -user isobe -script zline > zout'
    os.system(cmd)
 
    cmd  = 'mv *asol*fits* ./primary/.'
    os.system(cmd)

    cmd  = 'rm -rf *fits* zline zout'
    os.system(cmd)

    hdir = data_dir + '/' + str(obsid)
    if os.path.isdir(hdir):
        cmd = 'rm -rf ' + hdir + '/*'
        os.system(cmd)
    else:
        cmd = 'mkdir ' + hdir 
        os.system(cmd)

    cmd = 'chmod 774 primary/* secondary/*'
    os.system(cmd)

    cmd = 'mv primary secondary ' + hdir + '/.'
    os.system(cmd)


#-----------------------------------------------------------------------------------------

if __name__ == '__main__':

    if len(sys.argv) > 1:
        hrc = sys.argv[1].strip()
    else:
        hrc = 'hrc_s'

    run_process(hrc)


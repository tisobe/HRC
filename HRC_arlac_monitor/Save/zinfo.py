#!/usr/bin/env /proj/sot/ska/bin/python

#############################################################################################################
#                                                                                                           #
#               hrc_gain_find_ar_lac.py:  find  new AL Lac observations and put in a list                   #
#                                                                                                           #
#               author: t. isobe  (tisobe@cfa.harvard.edu)                                                  #
#                                                                                                           #
#               Last Update: Sep 24, 2014                                                                   #
#                                                                                                           #
#############################################################################################################

import os
import sys
import re
import string
import random
import operator
import numpy
import unittest
from astropy.table import Table
from Ska.DBI import DBI

#
#--- reading directory list
#
path = '/data/mta/Script/HRC/Gain/house_keeping/dir_list_py'

f    = open(path, 'r')
data = [line.strip() for line in f.readlines()]
f.close()

for ent in data:
    atemp = re.split(':', ent)
    var  = atemp[1].strip()
    line = atemp[0].strip()
    exec "%s = %s" %(var, line)
#
#--- append a path to a private folder to python directory
#
sys.path.append(bin_dir)
sys.path.append(mta_dir)
#
#--- converTimeFormat contains MTA time conversion routines
#
import convertTimeFormat    as tcnv
import mta_common_functions as mcf

#
#--- sql realated settings
#
db_user   = 'browser'
db_server = 'ocatsqlsrv'
file      = bdata_dir + '/.targpass'
db_passwd = mcf.get_val(file)

#---------------------------------------------------------------------------------------------------
#-- hrc_gain_find_ar_lac: find  new AL Lac observations and put in a list                         --
#---------------------------------------------------------------------------------------------------

def hrc_gain_find_ar_lac():

    """
    find  new AL Lac observations and put in a list
    Input: none, but the data will be read from mp_reports and also hrc_obsid_list in <house_keeping>
    Output: "./candidate_list"  which lists obsids of new AR Lac observations
            candidate_list      it also returns the same list
    """

    hrc_list       = hrc_gain_find_hrc_obs()
    candidate_list = hrc_gain_test_obs(hrc_list)

    return candidate_list

#---------------------------------------------------------------------------------------------------
#-- hrc_gain_find_ar_lac: find  new AL Lac observations and put in a list                         --
#---------------------------------------------------------------------------------------------------

def hrc_gain_find_hrc_obs():

    """
    select out the current hrc observations and create test candidate list
    Input: none, but the data will be read from mp_reports and also hrc_obsid_list in <house_keeping>
    Output: new_obs   --- recently observed HRC obsid list
    """

#
#--- read obsid list of AR Lac we already checked
#
    file = house_keeping + '/hrc_obsid_list'
    f    = open(file, 'r')
    obsid_list = [line.strip() for line in f.readlines()]
    f.close()
#
#--- find  HRC events from a recent mp_reports
#
    page = '/data/mta_www/mp_reports/events/mta_events.html'
    f    = open(page, 'r')
    data = [line.strip() for line in f.readlines()]
    f.close()
    data_list = []
    for ent in data:
        m1 = re.search('HRC',   ent)
        m2 = re.search('Obsid', ent)
        if (m1 is not None) and (m2 is not None):
            atemp = re.split('\/', ent)
            data_list.append(atemp[1])
#
#--- select out obsids which we have not checked before
#
    new_obs = []
    for ent in data_list:
        chk = 0
        for comp in  obsid_list:
            if ent == comp:
                chk = 1
                continue
        if chk > 0:
            continue
        new_obs.append(ent)

    return new_obs

#---------------------------------------------------------------------------------------------------
#-- hrc_gain_test_obs: find  new AL Lac observations from a hrc obsid list                      ----
#---------------------------------------------------------------------------------------------------

def hrc_gain_test_obs(new_obs, test=''):

    """
    find  new AL Lac observations from a hrc obsid list
    Input: new_obs  --- a list of hrc obsids
           test     --- a test indicator. if it is other than "", test will run
    Output: "./candidate_list"  which lists obsids of new AR Lac observations
            candidate_list      it also returns the same list
    """

    if test == "":
        f1    = open('./candidate_list', 'w')

        file  = house_keeping + 'hrc_obsid_list'
        file2 = house_keeping + 'hrc_obsid_list~'
        cmd   = 'cp -f ' + file + ' ' + file2
        os.system(cmd)
        f2    = open(file, 'a')

    candidate_list = []
    for obsid in new_obs:
#
#--- open sql database and extract data we need
#
        db  = DBI(dbi='sybase', server=db_server, user=db_user, passwd=db_passwd,  database='axafocat')

        cmd = 'select obsid,targid,seq_nbr,targname,grating,instrument from target where obsid='  + obsid 
        query_results = db.fetchall(cmd)
        if len(query_results):
            query_results = Table(query_results)

        line       = query_results['targname'].data
        targname   = line[0]
#
#--- if the observation is AR Lac, write it down in candidate_list
#
        m1 = re.search('arlac', targname.lower())
        if m1 is not None:
            line = obsid + '\n'
            candidate_list.append(obsid)

            if test == '':
                f1.write(line)
                f2.write(line)

    if test == '':
        f1.close()
        f2.close()

    return candidate_list

#-----------------------------------------------------------------------------------------
#-- TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST    ---
#-----------------------------------------------------------------------------------------

class TestFunctions(unittest.TestCase):
    """
    testing functions
    """

    def test_hrc_gain_test_obs(self):

        page = house_keeping + '/Test_prep/candidate'
        f    = open(page, 'r')
        data_list = [line.strip() for line in f.readlines()]
        f.close()

        test_candidates = ['14313', '14314', '14315', '14316']

        candidates = hrc_gain_test_obs(data_list, test='test')

        self.assertEquals(candidates, test_candidates)


#--------------------------------------------------------------------

if __name__ == '__main__':

    unittest.main()

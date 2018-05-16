#!/usr/bin/env /proj/sot/ska/bin/python

#################################################################################################
#                                                                                               #
#           find_coordinate.py: find coordinates of a given object                              #
#                                                                                               #
#           author: t. isobe (tisobe@cfa.harvard.edu)                                           #
#                                                                                               #
#           Last Update: Feb 10, 2017                                                           #
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
#--- reading directory list
#
path = '/data/aschrc6/wilton/isobe/Project7/Scripts/house_keeping/dir_list'

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

import mta_common_functions as mcf
import convertTimeFormat    as tcnv

#
#--- temp writing file name
#
rtail  = int(time.time())
zspace = '/tmp/zspace' + str(rtail)

#-----------------------------------------------------------------------------------------
#-- find_coordinate: find coordinates of the target                                    ---
#-----------------------------------------------------------------------------------------

def find_coordinate(target, cyear=''):
    """
    find coordinates of the target 
    input:  target  --- target name
            cyear   --- the year of the observed. it can be in year date. 
                        if it is not given, no corrections for the propoer motion
    output: coordinate  [ra, dec] in decimal format.
    """

#
#--- set several initial values
#
    target = target.strip()
    target = target.replace(' ', '\ ')
    tra    = 'na'
    tdec   = 'na'
    pra    = 0
    pdec   = 0
    schk   = 0
#
#--- call simbad to get the coordinate data from simbad site
#
    cmd    = 'lynx -source http://simbad.u-strasbg.fr/simbad/sim-id\?output.format=ASCII\&Ident=' + target + '>' + zspace
    os.system(cmd)

    f      = open(zspace, 'r')
    data   = [line.strip() for line in f.readlines()]
    f.close()
    mcf.rm_file(zspace)
#
#--- read the coordinates and the proper motion
#
    tchk = 0
    pchk = 0
    for ent in data:
        mc1 = re.search('Coordinates', ent)
        mc2 = re.search('Proper', ent)
        if tchk == 0 and mc1 is not None:
            try:
                atemp = re.split(':', ent)
                btemp = re.split('\s+', atemp[1])
                ahr   = float(btemp[1])
                amin  = float(btemp[2])
                asec  = float(btemp[3])
                tra   = 15.0 * (ahr + amin / 60.0 + asec / 3600.0)
     
                deg   = float(btemp[4])
                dmin  = float(btemp[5])
                dsec  = float(btemp[6])
                sign  = 1
                if deg < 0:
                    sign = -1
     
                tdec  = abs(deg) + dmin / 60.0 + dsec / 3600.0
                tdec *= sign
                tchk += 1
            except:
                schk = 1
                break


        if pchk == 0 and mc2 is not None:
            try:
                atemp = re.split(':', ent)
                btemp = re.split('\s+', atemp[1])
                pra   = btemp[1]
                pdec  = btemp[2]
                if mcf.chkNumeric(pra):
                    pra = float(pra) / 3600.0 / 1000.0
                else:
                    pra  = 0.0
                if mcf.chkNumeric(pdec):
                    pdec = float(pdec) / 3600.0 / 1000.0
                else:
                    pdec = 0.0
    
                pchk += 1
            except:
                pass
        if tchk == 1 and pchk == 1:
            break

#
#--- if the year is given, correct for the proper motion
#
    if schk == 0 and mcf.chkNumeric(cyear):
        dyear = float(cyear) - 2000.0
        try:
            tra  += dyear * pra
            tdec += dyear * pdec
        except:
            pass
        

    return [tra, tdec]

#-----------------------------------------------------------------------------------------

if __name__ == "__main__":

    if len(sys.argv) == 2:
        target = sys.argv[1].strip()
        year   = ''
    elif len(sys.argv) == 3:
        target = sys.argv[1].strip()
        year   = float(sys.argv[2].strip())

    [ra, dec]  = find_coordinate(target, year)
    print 'RA: ' + str(ra) + '   DEC: ' + str(dec)


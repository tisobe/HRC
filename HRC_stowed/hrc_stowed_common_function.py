#!/usr/bin/env /proj/sot/ska/bin/python

#############################################################################################
#                                                                                           #
#   hrc_stowed_common_function.py: collection of function used in hrc stowed computation    #
#                                                                                           #
#               author: t. isobe (tisobe@cfa.harvard.edu)                                   #
#                                                                                           #
#               last update: Jan 17, 2017                                                   #
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
from os import listdir
from os.path import isfile, join
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
#
#--- temp writing file name
#
rtail   = int(time.time())
zspace  = '/tmp/zspace' + str(rtail) + 'hrc'
zspace2 = zspace + 'x' 
#
#--- a couple of things needed
#
dare   = mcf.get_val('.dare',   dir = bindata_dir, lst=1)
hakama = mcf.get_val('.hakama', dir = bindata_dir, lst=1)
#
#--- other settings
#
NULL   = 'null'
#
#--- month list
#
m_list = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
#
#--- set base time at 1998.1.1
#
BTFMT    = '%m/%d/%y,%H:%M:%S'
FMT2     = '%Y-%m-%d,%H:%M:%S'
FMT3     = '%Y-%m-%dT%H:%M:%S'
basetime = datetime.strptime('01/01/98,00:00:00', BTFMT)
#
#--- set epoch
#
Epoch    = time.localtime(0)
Ebase_t  = time.mktime((1998, 1, 1, 0, 0, 0, 5, 1, 0))

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

def run_test():

    cmd = 'ls *fits.gz > ./ztemp2'
    os.system(cmd)

    f     = open('./ztemp2', 'r')
    flist = [line.strip() for line in f.readlines()]
    f.close()

    cond_list = [['QUALITY', '==', '0000000000000000000'], ['MNF', '==', 0]]

    outfits = './test.fits'

    fits_merge(flist, outfits, cond_list)

#-----------------------------------------------------------------------------------------------
#-- fits_merge: combine fits files in a list and extract data matched for the condition       --
#-----------------------------------------------------------------------------------------------

def fits_merge(flist, outfits,  codition_lists = ''):
    """
    combine fits files in a list and extract data matched for the condition
    this is a substitue for ciao dmmerge
    input:  flist           --- a list of fits file names
            outfits         --- output fits file name
            condition_lists --- a list of lists of condition in the form of:
                                [<col name>, <condition>, <value>] 
                                    e.g., ['QUALITY', '==', '0000000000000000000']
    output: outfits         --- combined fits file
    """
#
#--- find the first fits file; must be unzipped
#
    mc  = re.search('gz', flist[0])
    mcf.rm_file('temp.fits')
    mcf.rm_file('temp.fits.gz')
    if mc is not None:
        cmd = 'cp ' + flist[0] + ' temp.fits.gz'
        os.system(cmd)
        cmd = 'gzip -d temp.fits.gz'
        os.system(cmd)
    else:
        cmd = 'cp ' + flist[0] + ' temp.fits'
        os.system(cmd)
#
#--- combine all fits files
#
    for k in range(1, len(flist)):    
            mcf.rm_file('zout.fits')
            appendFitsTable('temp.fits', flist[k], 'zout.fits')
            cmd = 'mv zout.fits temp.fits'
            os.system(cmd)
#
#--- if there is no selection condition, just move the temp file to output fits file
#
    if codition_lists == '':
        cmd = 'mv temp.fits ' + outfits
        os.system(cmd)
#
#--- if the selection conditions are given, select data accordingly
#
    else:
        [cols, dout] = read_fits_file('temp.fits')
        for condition in codition_lists:
            col = condition[0]
            cnd = condition[1]
            val = condition[2]

            dout = select_data_with_condition(dout, col, cnd, val)

        hdu    = pyfits.open('temp.fits')
        hdu[1].data = dout

        hdu.writeto(outfits)

    mcf.rm_file('temp.fits')

#-------------------------------------------------------------------------------------------------------
#-- appendFitsTable: Appending one table fits file to the another                                    ---
#-------------------------------------------------------------------------------------------------------

def appendFitsTable(file1, file2, outname, extension = 1):

    """
    Appending one table fits file to the another
    the output table will inherit column attributes of the first fits table
    Input:  file1   --- fits table
            file2   --- fits table (will be appended to file1)
            outname --- the name of the new fits file
    Output: a new fits file "outname"
    """

    t1 = pyfits.open(file1)
    t2 = pyfits.open(file2)
#
#-- find numbers of rows (two different ways as examples here)
#
    nrow1 = t1[extension].data.shape[0]
    nrow2 = t2[extension].header['naxis2']
#
#--- total numbers of rows to be created
#
    nrows = nrow1 + nrow2
    hdu   = pyfits.BinTableHDU.from_columns(t1[extension].columns, nrows=nrows)
#
#--- append by the field names
#
    for name in t1[extension].columns.names:
        hdu.data.field(name)[nrow1:] = t2[extension].data.field(name)
#
#--- write new fits data file
#
    hdu.writeto(outname)

    t1.close()
    t2.close()
#
#--- attach the header
#
    header = t1[extension].header
    data, temp = pyfits.getdata(outname, extension, header=True)
    mcf.rm_file('ztemp.fits')
    pyfits.writeto('ztemp.fits', data, header)

    update_header_time('ztemp.fits', 'ztemp2.fits')
    mcf.rm_file('ztemp.fits')

    cmd = 'mv ztemp2.fits ' + outname
    os.system(cmd)

#-----------------------------------------------------------------------------------------------
#-- read_fits_file: read table fits data and return col names and data                       ---
#-----------------------------------------------------------------------------------------------

def read_fits_file(fits):
    """
    read table fits data and return col names and data
    input:  fits    --- fits file name
    output: cols    --- column name
            tbdata  --- table data
                        to get a data for a <col>, use:
                        data = list(tbdata.field(<col>))
    """
    hdulist = pyfits.open(fits)
#
#--- get column names
#
    cols_in = hdulist[1].columns
    cols    = cols_in.names
#
#--- get data
#
    tbdata  = hdulist[1].data

    hdulist.close()

    return [cols, tbdata]

#-----------------------------------------------------------------------------------------------
#-- read_header_value: read fits header value for a given parameter name                      --
#-----------------------------------------------------------------------------------------------

def read_header_value(fits, name):
    """
    read fits header value for a given parameter name
    input:  fits    --- fits file name
            name    --- parameter name
    output: val     --- parameter value
                        if the parameter does not exist, reuturn "NULL"
    """

    hfits = pyfits.open(fits)
    hdr   = hfits[1].header
    try:
        val   = hdr[name.lower()]
    except:
        val   = NULL

    hfits.close()

    return val

#-----------------------------------------------------------------------------------------------
#-- update_header: update header vlue of a fits file                                         ---
#-----------------------------------------------------------------------------------------------

def update_header(infile, cname, value=NULL, datatype='string', comment='', exten=1):
    """
    update header vlue of a fits file
    input:  infile      --- fits file name
            cname       --- parameter name
            value       --- parameter value
            datatype    --- data type
            comments    --- comment
            dpos        --- extension position  default: 1
    """

    hdulist = pyfits.open(infile, mode='update')

    header  = hdulist[exten].header
    if datatype == 'string':
        val = str(value)
        header[cname] = (val, comment)
    else:
        try:
            val = int(value)
        except:
            val = str(value)

        header[cname] = (val, comment)

    hdulist[exten].header = header

    hdulist.flush()
    hdulist.close()

#-----------------------------------------------------------------------------------------------
#-- create_fits_file: create a new fits file from a given data set                           ---
#-----------------------------------------------------------------------------------------------

def create_fits_file(data, col_list, fits_name):
    """
    create a new fits file from a given data set
    input:  data        --- data to be inserted to fits file
                            data can be a list of lists, or numpy data 
            col_list    --- a list of column names
            fits_name   --- output fits file name
    output: fits_name   --- fits file written
    Note: data must be numeric values only
    """
    col_l = []
    j     = 0
    for name in col_list:

        try:
            out = data.field(name)          #--- numpy array case
        except:
            out = data[j]                   #--- list case
            j += 1

        cols = set_format_for_col(name, out)
        col_l.append(cols)

    dcol  = pyfits.ColDefs(col_l)
    tbhdu = pyfits.BinTableHDU.from_columns(dcol)

    mcf.rm_file(fits_name)
    tbhdu.writeto(fits_name)

#-----------------------------------------------------------------------------------------------
#-- extract_col_data: extract pyfits table data from needed column and return a table data     --
#-----------------------------------------------------------------------------------------------

def extract_col_data(tbdata, col_list):
    """
    extract pyfits table data from needed column and return a table data
    input:  tbdata      --- pyfits table data
            col_list    --- a list of column names to be extracted
    output: dout        --- a table data with the extracted column data
    """
#
#--- extract column requested
#
    out = []
    for col_name in col_list:
        try:
            cdata = tbdata[col_name]
            col   = set_format_for_col(col_name, cdata)
            out.append(col)
        except:
            pass
#
#--- recreate table data format
#
    dcol  = pyfits.ColDefs(out)
    dbhdu = pyfits.BinTableHDU.from_columns(dcol)
    dout  = dbhdu.data

    return dout

#-----------------------------------------------------------------------------------------------
#-- set_format_for_col: find a format of the input data and set column object                 --
#-----------------------------------------------------------------------------------------------

def set_format_for_col(name, cdata):
    """
    find a format of the input data and set column object
    input:  name    --- column name
            cdata   --- column data in numpy array form
    output: column object
    """
    test = str(list(cdata)[0])
    mc1  = re.search('True',  test)
    mc2  = re.search('False', test)
#
#--- check whether the value is numeric
#
    if mcf.chkNumeric(test):
        ft = 'E'
#
#--- check whether the value is logical
#
    elif  (mc1 is not None) or (mc2 is not None):
        tcnt = len(list(cdata)[0])
        ft   = str(tcnt) + 'L'
#
#--- all others are set as character set
#
    else:
        tcnt = len(test)
        ft   = str(tcnt) + 'A'

    return pyfits.Column(name=name, format=ft, array=cdata)

#-----------------------------------------------------------------------------------------------
#-- combine_and_select_data: for a given list of fits files, combine them and extract columns needed 
#-----------------------------------------------------------------------------------------------

def combine_and_select_data(fits_list, selected_cols=''):
    """
    for a given list of fits files, combine them and extract columns needed
    input:  fits_list       --- a list of fits files
            selected_cols   --- a list of columns to extracted
    output: dout            --- fits data objects with the selected data
    """
#
#--- extract data 
#
    chk   = 0
    for fits in fits_list:
        [cols, tbdata] = read_fits_file(fits)
        if chk == 0:
            ctbdata = tbdata
            chk = 1
        else:
            ctbdata = combine_tabledata(ctbdata, tbdata, cols)

    if selected_cols == '':
        return ctbdata
#
#--- extract only columns we need
#
    else:
        try:
            dout = extract_col_data(ctbdata, selected_cols)
        except:
            dout = NULL
     
        return dout

#-----------------------------------------------------------------------------------------------
#-- combine_tabledata: combine two pyfits table data into one                                 --
#-----------------------------------------------------------------------------------------------

def combine_tabledata(tbdata1, tbdata2, col):
    """
    combine two pyfits table data into one
    input:  tbdata1 --- table data
            tbdata2 --- table data
            col     --- table data column name
    output: dout    --- combined table data
    """

    nrow1 = tbdata1.shape[0]
    nrow2 = tbdata2.shape[0]
    nrows = nrow1 + nrow2
    hdu = pyfits.BinTableHDU.from_columns(tbdata1.columns, nrows=nrows)

    for colname in col:
        try:
            hdu.data[colname][nrow1:] = tbdata2[colname]
        except:
            pass

    data = hdu.data

    return data

#-----------------------------------------------------------------------------------------------
#-- combine_fits_files: combine all fits file in a list to one                                --
#-----------------------------------------------------------------------------------------------

def combine_fits_files(fits_list, out, azip=1, fcol_line=''):
    """
    combine all fits file in a list to one
    input:  fits_list   --- a list of fits files
            out         --- output fits file name
            azip        --- if it is 1, the final file is zipped
            fcol_line   --- a string contain a list of column names (e.g. "time,lldialv,rsrfalv")
                            default: ""
    output: temp.fits   --- a combined fits file
    """
#
#--- check whether there is any file
#
    temp_list = []
    for ent in fits_list:
        if isfile(ent):
            temp_list.append(ent)

    fits_list = temp_list

    flen = len(fits_list)

    if flen == 0:
        return NULL
#
#--- if there in only one file, check gipping is asked or not. 
#
    elif flen == 1:
        mc = re.search('gz', fits_list[0])
        mtemp = fits_list[0]
        if mc is not None:
            if azip != 1:
                cmd = 'gzip -d ' + fits_list[0]
                os.system(cmd)
                mtemp = mtemp.replace('.gz','')
        else:
            if azip == 1:
                cmd = 'gzip -f ' + fits_list[0]
                os.system(cmd)
                mtemp = mtemp + '.gz'

        cmd = 'mv  ' + mtemp + ' ' + out
        os.system(cmd)
    else:
#
#--- combine the first two
#
        combine_two_fits_files_ascds(fits_list[0], fits_list[1], fcol_line=fcol_line)
#
#--- combine the rest
#
        for k in range(2,len(fits_list)):
            cmd = 'mv ztemp99.fits zxc99.fits'
            os.system(cmd)
            combine_two_fits_files_ascds('zxc99.fits', fits_list[k], fcol_line=fcol_line)
#
#--- change the fits name to the specified one
#
        out = out.replace('.gz', '')
#
#--- update header time
#
        update_header_time('ztemp99.fits', out)
#
#--- if requested, gip the file
#
        if azip == 1:
            out = out.replace('.gz', '')
            cmd = 'gzip -f ' + out + ' 2>/dev/null'
            os.system(cmd)

    mcf.rm_file('zxc99.fits')
    mcf.rm_file('ztemp99.fits')

#-----------------------------------------------------------------------------------------------
#-- append_two_fits_files: append one fits file to the other                                  --
#-----------------------------------------------------------------------------------------------

def append_two_fits_files(fits1, fits2, out, fcol_line=''):
    """
    append one fits file to the other
    input:  fits1   --- a fits file
            fits2   --- a fits file
            out     --- name of the combined fits file
    output: out     --- a combined fits file
    """

    chk1 = 0
    chk2 = 0
    try:
        t1 = pyfits.open(fits1)
    except:
        chk1 = 1
    try:
        t2 = pyfits.open(fits2)
    except:
        chk2 = 1

    if chk1 == 0 and chk2 == 0:

        combine_two_fits_files_ascds(fits1, fits2, fcol_line)
#
#--- update header
#
        update_header_time('ztemp99.fits', out)
        mcf.rm_file('ztemp99.fits')

        return 1

    elif chk1 == 1 and chk2 == 1:
        return 0
    elif chk2 == 1:
        cmd = 'mv ' + fits1 + ' ' + out
        os.system(cmd)
        return 1
    elif chk1 == 1:
        cmd = 'mv ' + fits2 + ' ' + out
        os.system(cmd)
        return 1
    else:
        return 0

    try:
        t1.close()
    except:
        pass
    try:
        t2.close()
    except:
        pass

#-----------------------------------------------------------------------------------------------
#-- update_header_time: update time part of the header                                        --
#-----------------------------------------------------------------------------------------------

def update_header_time(fits, out):
    """
    update time part of the header
    input:  fits    --- input fits file name
            out     --- output fits file name
    output: out     --- updated fits file
    """
#
#--- read the fits file
#
    hdu = pyfits.open(fits)
#
#--- update header
#
    header  = hdu[1].header
    tlist   = hdu[1].data['time']
    tstart  = numpy.min(tlist)
    tstop   = numpy.max(tlist)
    begin   = covertfrom1998sec2(tstart)
    end     = covertfrom1998sec2(tstop)
    current = covertfrom1998sec2(tcnv.currentTime(format='SEC1998'))

    header['tstart']   = (tstart, 'Data file start time')
    header['tstop']    = (tstop,  'Data file stop time')
    header['date']     = (current,'Date and time of file creation')
    header['date-obs'] = (begin,  'TT, with clock correction if clockapp')
    header['date-end'] = (end,    'TT, with clock correction if clockapp')

    hdu[1].header = header

    mcf.rm_file(out)
    hdu.writeto(out)
    hdu.close()

#-----------------------------------------------------------------------------------------------
#-- combine_two_fits_files_ascds: combine two fits files with dmtool dmmerge                 ---
#-----------------------------------------------------------------------------------------------

def combine_two_fits_files_ascds(fits1, fits2, outfile='ztemp99.fits',  fcol_line=''):
    """
    combine two fits files with dmtool dmmerge
    input:  fits1       --- first fits file name
            fits2       --- second fits file name
            outfile     --- the name of the combined fits file. default: ztemp99.fits
            fcol_line   --- a string of a list of column name to be combined
    output: ztemp99.fits  --- combined fits file
    """
    cmd1 = '/usr/bin/env PERL5LIB=""'
    cmd2 = ' dmmerge "' + str(fits1) +', ' +str(fits2) 
    cmd2 = cmd2 + '" outfile=' + outfile + '  outBlock="" columnList="' + str(fcol_line) + '" clobber="yes"'

    cmd  = cmd1 + cmd2

    bash(cmd, env=ascdsenv)
    
#-----------------------------------------------------------------------------------------------
#-- select_data_with_condition: select data accorinding to a condition given to the column "col" 
#-----------------------------------------------------------------------------------------------

def select_data_with_condition(tbdata, col, cond, val):
    """
    select data accorinding to a condition given to the column "col"
    input:  tbdata  --- table data
            col     --- column name
            cond    --- boolen  e.g. "<" ">" "<=" ">=" and "=="
    output: odata   --- updated table data
    """
    chk = 0
    if cond == '<':
        mask = tbdata[col] <  val

    elif cond == '>':
        mask = tbdata[col] >  val

    elif cond == '<=':
        mask = tbdata[col] <= val

    elif cond == '>=':
        mask = tbdata[col] >= val

    else:
#
#--- often the data needs to be tested agaist an array of condition (e.g. [True, True, False])
#--- so we should check whether it is the case
#
        if (isinstance(val, float)) or (isinstance(val, int)):
            mask = tbdata[col] == val
        else:
            if val == 'SPEC' or val == 'IMAG':
                mask = tbdata.field(col) == val
            else:
                odata = select_data_with_logical_mask(tbdata, col, val)
                chk = 1

    if chk == 0:
        odata = tbdata[mask]

    return odata

#-----------------------------------------------------------------------------------------------
#-- create_array_mask: create a mask for the array  ---- CURRENTLY NOT USED!!                ---
#-----------------------------------------------------------------------------------------------

def create_logic_mask(cbdata, mask):
    """
    create a mask for the array
    input:  cbdata  --- table data for the column
            mask    --- original mask with possibly array of coditions
    output: test    --- array of True and False for a mask
    """
#
#--- convert mask format from xxx1xxx (1 is True 0 is False and x is no checking)
#
    poslist   = []
    condition = []
    for i in range(0, len(mask)):
        ent = mask[i]
        if ent == 'x':
            continue
        else:
            poslist.append(i)
            if ent == '1':
                condition.append(True)
            else:
                condition.append(False)

#
#--- check each row of the column data to see whether the condition matches
#
    lcondition = []
    for ent in cbdata:
        lkey = True
        for m in range(0, len(poslist)):
            k = poslist[m]
            if ent[k] != condition[m]:
                lkey = False
                break
        lcondition.append(lkey)


    return lcondition

#-----------------------------------------------------------------------------------------------
#-- select_data_with_logical_mask: for given column name and logical mask, select out data and return data table
#-----------------------------------------------------------------------------------------------

def select_data_with_logical_mask(tbdata, col, mask):
    """
    for given column name and logical mask, select out data and return data table
    input:  tbdata  --- table data
            col     --- column name to be examined
            mask    --- original mask with possibly array of coditions
    output: ntbdata --- table data with selected data rows
    """
#
#--- find conlumn names
#
    col_names = tbdata.columns.names
#
#--- convert mask format from xxx1xxx (1 is True 0 is False and x is no checking)
#
    poslist   = []
    condition = []
    for i in range(0, len(mask)):
        ent = mask[i]
        if ent == 'x':
            continue
        else:
            poslist.append(i)
            if ent == '1':
                condition.append(True)
            else:
                condition.append(False)

#
#--- check each row of the column data to see whether the condition matches
#
    save = []
    for rdata in tbdata:
        ent = rdata[col]
        lkey = True
        for m in range(0, len(poslist)):
            k = poslist[m]
            if ent[k] != condition[m]:
                lkey = False
                break
        if lkey == True:
            save.append(rdata)
#
#--- recreate data table
#
    out = []
    for col in col_names:
        vdata = []
        for ent in save:
            val = ent.field(col)
            vdata.append(val)

        vcol = set_format_for_col(col, vdata)
        out.append(vcol)

    dcol    = pyfits.ColDefs(out)
    dbhdu   = pyfits.BinTableHDU.from_columns(dcol)
    dout    = dbhdu.data

    return dout

#-------------------------------------------------------------------------------------------------------
#-- fitsTableStat: find min, max, avg, std, and mediam of the column                                  --
#-------------------------------------------------------------------------------------------------------

def fitsTableStat(file, column, extension=1):

    """
    find min, max, avg, std, and mediam of the column. 
    Input   file--- table fits file name
    column  --- name of the column(s). if there are more than one, must be 
    in the form of list or tuple
    extension-- data extension #. default = 1
    Output  a list or a list of lists of [min, max, avg, std, med, sample size]
    """
    
    t = pyfits.open(file)
    tdata = t[extension].data
    t.close()
    
    if isinstance(column, list) or isinstance(column, tuple):
    
        results = []
        for ent in column:
            data = tdata.field(ent)
            line = []
            dmin = min(data)
            line.append(dmin)
            dmax = max(data)
            line.append(dmax)
            avg  = numpy.mean(data)
            line.append(avg)
            std  = numpy.std(data)
            line.append(std)
            med  = numpy.median(data)
            line.append(med)
            line.append(len(data))
     
            if len(column) > 1:
                results.append(line)
            else:
                results = line
                break
    
    else:
        data = tdata.field(column)
        results = []
        dmin = min(data)
        results.append(dmin)
        dmax = max(data)
        results.append(dmax)
        avg  = numpy.mean(data)
        results.append(avg)
        std  = numpy.std(data)
        results.append(std)
        med  = numpy.median(data)
        results.append(med)
        results.append(len(data))
    
    return results


#-----------------------------------------------------------------------------------------------
#-- run_arc5gl: extract data from archive using arc5gl                                        --
#-----------------------------------------------------------------------------------------------

def run_arc5gl(operation, start='', stop='', dataset='flight', detector='hrc', level='0', filetype='hrchk',filename='', subdetector=''):
    """
    extract data from archive using arc5gl (fomerlly arc4gl)
    input:  operation   --- operation command either, browse or retrieve
            start       --- starting time 
            stop        --- stoping time
            dataset     --- dataset name.       default = flight
            detector    --- detector name       default = hrc  
            level       --- level               default = o
            filetype    --- filetype            default = hrchk
            filename    --- file name           default = ''
            subdetector --- sub detector name   default = ''
    output: extracted data (if the operation is retrieve)
            data        --- a list of the data 
    """
    temp_dir = exc_dir + 'Temp_dir/'
#
#--- modify the time format acceptable by arc5gl
#
    start = convert_time_format(start)
    stop  = convert_time_format(stop)
#
#---  arc5gl command part which used more than once
#
    mline = 'dataset = '  + dataset    + '\n'
    mline = mline + 'detector = ' + detector   + '\n'

    if subdetector != '':
        mline = mline + 'subdetector = ' + subdetector   + '\n'

    mline = mline + 'level = '    + str(level) + '\n'
    mline = mline + 'filetype = ' + filetype   + '\n'

    if start != '':
        mline = mline + 'tstart=' + str(start) + '\n'

    if stop  != '':
        mline = mline + 'tstop='  + str(stop)  + '\n'

    if filename  != '':
        mline = mline + 'filename'+ str(filename) + '\n'

    mline = mline + 'go\n'

    cmd1 = '/usr/bin/env PERL5LIB=""'
    cmd2 = ' /proj/axaf/simul/bin/arc5gl -user ' + dare + ' -script ' + zspace
#
#--- set operation : browse/retrieve
#
    if operation == 'browse':
        line = 'operation = browse\n'
        cmd = cmd1 + cmd2 + ' > ' + zspace2

    elif operation == 'retrieve':
        line = 'operation = retrieve\n'
        cmd = 'cd ' + temp_dir + '; ' + cmd1 + cmd2 + '> ' + zspace2

    else:
        return []

    line = line + mline
#
#--- write arc5gl command in a file
#
    fo = open(zspace, 'w')
    fo.write(line)
    fo.close()
#
#--- run arc5gl command
#
    chk = 0
    try:
        bash(cmd, env=ascdsenv)
    except:
        try:
            bash(cmd, env=ascdsenv)
        except:
            chk == 1
     
    mcf.rm_file(zspace)
#
#--- read and/or clean up the list of the data obtained
#
    if chk == 0:
        data = clean_up_fits_list(operation, zspace2)
    else:
        data = []

    return data

#-----------------------------------------------------------------------------------------------
#-- convert_time_format: convert time format to be acceptable by arc5gl                      ---
#-----------------------------------------------------------------------------------------------

def convert_time_format(stime):
    """
    convert time format to be acceptable by arc5gl
    input:  stime   --- time
    output  simte   --- modified (if needed) time
    note: only time format does not accept is mm/dd/yy,hh:mm:ss which is aceptable in ar4cgl
          so convert that into an acceptable format: yyyy-mm-ddThh:mm:ss
    """
#
#--- if the time is seconds from 1998.1.1, just pass
#
    if isinstance(stime, float) or isinstance(stime, int):
        return stime
#
#--- time word cases; only mm/dd/yy,hh:mm:ss is modified
#
    mc = re.search('\,', stime)
    if mc is not None:
        atemp = re.split('\,', stime)
        btemp = re.split('\/', atemp[0])
        mon   = btemp[0]
        day   = btemp[1]
        yr    = int(float(btemp[2]))
        if yr > 90:
            yr += 1900
        else:
            yr += 2000
        stime = str(yr) + '-' + mon + '-' + day + 'T' + atemp[1]

    return stime

#-----------------------------------------------------------------------------------------------
#-- clean_up_fits_list: clean up the fits list                                                --
#-----------------------------------------------------------------------------------------------
    
def clean_up_fits_list(operation, fname):
    """
    clean up the fits list 
    input:  operation   --- browse or retrieve
            fname       --- the name of the file to be read
    output: data        --- a list of fits files
    """

    tdir = exc_dir + 'Temp_dir/'

    try:
        out = read_file_data(fname, remove=1)
    except:
        return []
#
#--- extract fits file names
#
    data = []
    for ent in out:
        mc = re.search('fits', ent)
        if mc is not None:
            atemp = re.split('\s+', ent)
            if operation == 'retrieve':
                name = tdir + atemp[0]
            else:
                name = atemp[0]
            data.append(name)

    data.sort()

    return data

#-----------------------------------------------------------------------------------------------
#-- read_file_data: read the content of the file and return it                                --
#-----------------------------------------------------------------------------------------------

def read_file_data(infile, remove=0):
    """
    read the content of the file and return it
    input:  infile  --- file name
            remove  --- if 1, remove the input file after read it, default: 0
    output: out     --- output
    """

    f   = open(infile, 'r')
    out = [line.strip() for line in f.readlines()]
    f.close()

    if remove == 1:
        cmd = 'rm ' + infile
        os.system(cmd)

    return out

#-----------------------------------------------------------------------------------------------
#-- clean_dir: empty out the directory content                                                --
#-----------------------------------------------------------------------------------------------

def clean_dir(tdir):
    """
    empty out the directory content
    input:  tdir    --- a directory path
    output: none
    """

    chk = 0
    if os.listdir(tdir):
        chk = 1

    if chk == 1:
        cmd   = 'rm -rf ' +  tdir + '/*'
        os.system(cmd)

#-----------------------------------------------------------------------------------------------
#-- get_file_list: create a file name list in a given directory                               --
#-----------------------------------------------------------------------------------------------

def get_file_list(tdir, tail= ''):
    """
    create a file name list in a given directory
    input:  tdir    --- a directory path
            tial    --- a sufix of the file name. default is blank
    output: data    --- a list of file names
    """
#
#--- check whether the given directory exists
#
    if os.path.isdir(tdir) and check_file_in_dir(tdir):
#
#--- check whether any files in that directory
#
        cmd = 'ls -a  ' + tdir + '* > '+  zspace
        os.system(cmd)
        f    = open(zspace, 'r')
        test = f.read()
        f.close()
        mcf.rm_file(zspace)

        if len(str(test)) > 0:
            chk = 0
            if tail == '':
                chk = 1
            else:
                mc   = re.search(tail, test)
                if mc is not None:
                    chk = 1 

            if chk == 1:
                cmd  = 'ls ' + tdir + '*' +  tail + '> ' + zspace
                os.system(cmd)
            
                data = read_file_data(zspace, remove=1)
            else:
                data = []
        else:
            data =  []
    else:
        data = []

    return data

#-----------------------------------------------------------------------------------------------
#-- check_file_in_dir: check whether there are any files in the directory                     --
#-----------------------------------------------------------------------------------------------

def check_file_in_dir(path):
    """
    check whether there are any files in the directory
    input:  path    --- a path to the directory
    output: True/False
    """

    return any(isfile(join(path, i)) for i in listdir(path))

#-----------------------------------------------------------------------------------------------
#-- conv_time_format: change time forat to arc5gl usable form                                ---
#-----------------------------------------------------------------------------------------------

def conv_time_format(year, month, next = 0):
    """
    change time forat to arc5gl usable form
    input:  year    --- year
            month   --- month
            next    --- if > 0, set to the 1st of the next month
    output: date    --- date in the format of "mm/01/yy,00:00:00"
    """
    if next > 0:
        month += 1
        if month > 12:
            month = 1
            year += 1

    lyear = str(year)
    lmon  = str(month)
    if month < 10:
        lmon = '0' + lmon

    date = lmon + '/01/' + lyear[2] + lyear[3] + ',00:00:00'

    return date

#-----------------------------------------------------------------------------------------------
#-- convertto1998sec: convert time format from mm/dd/yy,hh:mm:ss to seconds from 1998.1.1    ---
#-----------------------------------------------------------------------------------------------

def convertto1998sec(ftime):
    """
    convert time format from mm/dd/yy,hh:mm:ss to seconds from 1998.1.1
    input   ftime   --- time in mm/dd/yy,hh:mm:ss or yyyy-mm-dd,hh:mm:ss
    output  stime   --- time in seconds from 1998.1.1
    """
    mc = re.search('\/', str(ftime))
    mc2= re.search('\:', str(ftime))
    if mc is None and mc2 is None:
        if isinstance(ftime, int):
            return ftime
        elif isinstance(ftime, float):
            return ftime
        else:
            try:
                sec1998 = int(float(ftime))
                return sec1998
            except:
                return NULL
    else:
#
#--- base time 1998 Jan 1, 00:00:00 (see top setting section)
#
#--- remove any empty space
#
        ftime = ftime.replace('\s+', '')

        mc    = re.search('-', ftime)
#
#--- for yyyy-mm-dd,hh:mm:ss
#
        if mc is not None:
            ftime = datetime.strptime(ftime, FMT2)
#
#--- for mm/dd/yy,hh:mm:ss
#
        else:
            ftime = datetime.strptime(ftime, BTFMT)

        tdel  = ftime - basetime

        sec1998 = 86400 * tdel.days + tdel.seconds

        return sec1998


#-----------------------------------------------------------------------------------------------
#-- covertfrom1998sec: convert second from 1998.1.1 to mm/dd/yy,hh:mm:ss format               --
#-----------------------------------------------------------------------------------------------

def covertfrom1998sec(stime):
    """
    convert second from 1998.1.1 to mm/dd/yy,hh:mm:ss format
    input:  stime   --- second from 1998.1.1
            etime   --- time in mm/dd/yy,hh:hh:ss
    """

    etime = Ebase_t + stime
    etime = time.localtime(etime)
    etime = time.strftime(BTFMT, etime)

    return etime
        

#-----------------------------------------------------------------------------------------------
#-- covertfrom1998sec2: convert second from 1998.1.1 to yyyy-mm-ddThh:mm:ss format            --
#-----------------------------------------------------------------------------------------------

def covertfrom1998sec2(stime):
    """
    convert second from 1998.1.1 to yyyy-mm-ddThh:mm:ss format 
    input:  stime   --- second from 1998.1.1
            etime   --- time in yyyy-mm-ddThh:mm:ss
    """

    etime = Ebase_t + stime
    etime = time.localtime(etime)
    etime = time.strftime(FMT3, etime)

    return etime
        
#-----------------------------------------------------------------------------------------------
#-- find_month: convert month format from digit to letter or letter to digit                  --
#-----------------------------------------------------------------------------------------------

def find_month(mon):
    """
    convert month format from digit to letter or letter to digit
    input:  mon     --- either digit month or letter month
    output: lmon    --- either digit month or letter month
    """

    if mcf.chkNumeric(mon):
        mon = int(float(mon))
        lmon = m_list[mon-1]
        return lmon
    else:
        mon = mon.lower()
        for i in range(0, 12):
            if mon == m_list[i].lower():
                lmon = i + 1
                break
        return lmon
            
#-----------------------------------------------------------------------------------------------
#-- remove_duplicate_from_file: remove duplcated lines from a given file                      --
#-----------------------------------------------------------------------------------------------

def remove_duplicate_from_file(efile):
    """
    remove duplcated lines from a given file
    input:  efile   --- an input file name
    output: efile   --- the cleaned file
    """
    try:
        rdata = read_file_data(efile)
    except:
        return False
#
#--- make sure that there are no duplicated line
#
    cleaned = remove_duplicate_from_list(rdata)
#
#--- print out the cleaned data
#
    fo  = open(efile, 'w')
    for ent in cleaned:
        fo.write(ent)
        fo.write('\n')
    fo.close()

    return True

#-----------------------------------------------------------------------------------------------
#-- remove_duplicate_from_list: remove duplicated input from a list                          ---
#-----------------------------------------------------------------------------------------------

def remove_duplicate_from_list(alist):
    """
    remove duplicated input from a list
    alist   --- an input list
    cleaned --- a cleaned list
    """

    test = len(alist)
    if test > 0:
        alist.sort()
        prev    = alist[0]
        cleaned = [prev]
        for i in range(1,len(alist)):
            if alist[i] != prev:
                cleaned.append(alist[i])
                prev = alist[i]
            else:
                continue
    else:
        cleaned = []

    return cleaned

#-----------------------------------------------------------------------------------------------
#-- remove_duplicate_by_column: sort a list by a specific coloum position and remove duplicated lines   
#-----------------------------------------------------------------------------------------------

def remove_duplicate_by_column(alist, colpos, sep = ''):
    """
    sort a list by a specific coloum position and remove duplicated lines
    input:  alist   --- a list of lists or strings
            colpos  --- a position of column or for the case of string, the column position after separated
            sep     --- a separator, if it is not given, the function will guess
    output: cleaned --- a cleaned list
    """

    test = len(alist)
    if test < 2:
        return alist

    olist = sorted_by_column(alist, colpos, sep)

    prev    = olist[0]
    cleaned = [prev]
    for i in range(1,len(olist)):
        if olist[i] != prev:
            cleaned.append(olist[i])
            prev = olist[i]
        else:
            continue

    return cleaned


#-----------------------------------------------------------------------------------------------
#-- sorted_by_column: sort a list or a string separated by a separtor by a specified column   --
#-----------------------------------------------------------------------------------------------

def sorted_by_column(alist, colpos, sep = '', vald=1):
    """
    sort a list or a string separated by a separtor by a specified column
    input:  alist   --- a list of lists or a string
            colpos  --- a position of column; if the case of a string, the psition after 
                        splited by the separator
            sep     --- a column sepatator. if it is not given, it will geuss from the line
            vald    --- indicator of whether the column value is digit. defalut: 1 (digit)
    output: olist   --- sorted list
    """

    test   = len(alist)
    if test < 2:
        return alist

    adict    = {}
    ind_list = []
    if isinstance(alist[0], list):
        chk = 0
    else:
        chk = 1
        if sep == '':
            sep = find_separator(alist[0])

    for ent in alist:
        if chk == 0:
            val = float(ent[colpos])
            adict[val] = ent
            ind_list.append(val)
        else:
            atemp = re.split(sep, ent)
            if vald == 1:
                try:
                    val = float(atemp[colpos])
                except:
                    continue
            else:
                val = atemp[colpos]

            adict[val] = ent
            ind_list.append(val)

    ind_list.sort()

    olist = []
    for oind in ind_list:
        olist.append(adict[oind])

    return olist

#-----------------------------------------------------------------------------------------------
#-- find_separator: find a separator                                                          --
#-----------------------------------------------------------------------------------------------

def find_separator(line):
    """
    find a separator 
    input:  line    --- an example line to be examined
    output: sep     --- a separator
    """

    mc1 = re.search('\,', line)
    mc2 = re.search(':',  line)
    mc3 = re.search(';',  line)
    mc4 = re.search('\t', line)

    if mc1 is not None:
        sep = '\,'
    elif mc2 is not None:
        sep = ':'
    elif mc3 is not None:
        sep = ';'
    elif mc4 is not None:
        sep = '\t+'
    else:
        sep = '\s+|\t+'

    return sep


#-----------------------------------------------------------------------------------------
#-- TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST    ---
#-----------------------------------------------------------------------------------------

class TestFunctions(unittest.TestCase):
    """
    testing functions
    """

#------------------------------------------------------------

    def test_fitsfile_manupulation(self):

        a1       = numpy.array([1.0, 2.0, 3.0, 4.0])
        a2       = numpy.array([5.0, 6.0, 7.0, 8.0])
        a3       = numpy.array([8.0, 9.0, 10.0, 11.0])
        a4       = [1.0, 2.0, 3.0, 4.0, 1.0, 2.0, 3.0, 4.0]
        a5       = [5.0, 6.0]
        data     = [a1, a2, a3]
        col_list = ['a1', 'a2', 'a3']
        fname    = 'test.fits'
#
#--- creating fits file
#
        create_fits_file(data, col_list, fname)
#
#--- reading the fits file
#
        [cols, odata] = read_fits_file(fname)

        o1 = odata.field(cols[0])
        o2 = odata.field(cols[1])

        self.assertEquals(list(a1), list(o1))
        self.assertEquals(list(a2), list(o2))
#
#--- extracting a part of the data
#
        p1  = extract_col_data(odata, ['a3'])
        out = p1.field(0)
        self.assertEquals(list(a3), list(out))

        tfits     = house_keeping + 'Test_data/hrcf567118652N001_hk0.fits.gz'
        tcol_list = ['TIME','LLDIALV','RSRFALV','WDTHAST']
        [dummy, pdata] = read_fits_file(tfits)
        p1  = extract_col_data(pdata, tcol_list)
        tout = list(p1['RSRFALV'])
        ttest = tout[0:10]
        tcomp = [185., 185., 185., 185., 185., 185., 185., 185., 185., 185.]
        self.assertEquals(ttest, tcomp)

#
#--- combine two table data
#
        out  = combine_tabledata(odata, odata, cols)
        test = out.field('a1')
        self.assertEquals(list(test), a4)


        tfile = house_keeping + 'Test_data/hrcf567118652N001_hk0.fits.gz'
        [cols, tbdata] = read_fits_file(tfile)
        tfile = house_keeping + 'Test_data/hrcf567189107N001_hk0.fits.gz'
        [cols2, tbdata2] = read_fits_file(tfile)

        ctbdata = combine_tabledata(tbdata, tbdata2, cols)
        p1  = extract_col_data(ctbdata, tcol_list)
        tout = list(p1['RSRFALV'])
        ttest = tout[-10:]
        tcomp = [185.0, 185.0, 185.0, 185.0, 185.0, 185.0, 185.0, 185.0, 185.0, 185.0]
        self.assertEquals(ttest, tcomp)
#
#--- selecting data with condition on one column
#
        out  = select_data_with_condition(odata, 'a1', '<=', 2.0)
        test = out.field('a2')
        self.assertEquals(list(test), a5)

        out  = select_data_with_condition(odata, 'a1', '==', 2.0)
        test = out.field('a2')
        self.assertEquals(list(test), [6.0])

        mcf.rm_file(fname)

        tfile = house_keeping + 'Test_data/hrcf567118652N001_hk0.fits.gz'
        [cols, odata] = read_fits_file(tfile)

        comp    =  [False, True,False, True, True,False,False,False]
        lpart   = '1xxxxxxx'
        lpart   = 'xxx1xxxx'
        out  = select_data_with_condition(odata, 'HVPSSTAT', '==', lpart)
        test = out.field('HVPSSTAT')
        self.assertEquals(list(test[0]), comp)
#
#--- header update
#
        cmd = 'cp ' +  house_keeping + 'Test_data/hrcf567118652N001_hk0.fits.gz ./ztest.fits.gz'
        os.system(cmd)

        infile   = './ztest.fits.gz'
        cname    = 'TESTNAME'
        datatype = 'string'
        value    = 'aaaaa'
        comment  = 'this is a test param value'
        update_header(infile, cname, value=value, datatype=datatype, comment=comment)

        cname2   = 'TESTINT'
        datatype = 'short'
        value2   =  1256
        comment  = 'this is a test param int value'
        update_header(infile, cname2, value=value2, datatype=datatype, comment=comment)

        val      = read_header_value(infile, cname)
        self.assertEquals(val, value)

        val      = read_header_value(infile, cname2)
        self.assertEquals(val, value2)

        mcf.rm_file(infile)

        #fits1 = house_keeping + 'Test_data/hrcf567118652N001_hk0.fits.gz'
        #fits2 = house_keeping + 'Test_data/hrcf567189107N001_hk0.fits.gz'
        #out   = './test_comb.fits'
        #append_two_fits_files(fits1, fits2, out)
        #hfits = pyfits.open(out)
        #hdr   = hfits[1].header
        #print str(hdr)

#------------------------------------------------------------

    def test_combine_fits_files(self):

        fits1 = house_keeping + 'Test_data/hrcf567118652N001_hk0.fits.gz'
        fits2 = house_keeping + 'Test_data/hrcf567189107N001_hk0.fits.gz'
        fits_list = [fits1, fits2]

        out   = exc_dir + 'test.fits'

        combine_fits_files(fits_list, out)


#------------------------------------------------------------

    def test_run_arc5gl(self):
        comp_list = ['hrcf536454389N001_hk0.fits', 'hrcf536471051N001_hk0.fits', 'hrcf536502671N001_hk0.fits', \
                     'hrcf536502933N001_hk0.fits', 'hrcf536514315N001_hk0.fits', 'hrcf536541375N001_hk0.fits', \
                     'hrcf536559349N001_hk0.fits', 'hrcf536589099N001_hk0.fits']
        comp_list2= ['hrcf536454389N001_hk0.fits.gz', 'hrcf536471051N001_hk0.fits.gz', 'hrcf536502671N001_hk0.fits.gz', \
                     'hrcf536502933N001_hk0.fits.gz', 'hrcf536514315N001_hk0.fits.gz', 'hrcf536541375N001_hk0.fits.gz', \
                     'hrcf536559349N001_hk0.fits.gz', 'hrcf536589099N001_hk0.fits.gz']

        start = '01/01/15,00:00:00'
        stop  = '01/03/15,00:00:00'
        #data  = run_arc5gl('browse', start, stop,  dataset='flight', detector='hrc', level='0', filetype='hrchk')
        
        #self.assertEquals(data, comp_list)

        #data  = run_arc5gl('retrieve',start, stop,  dataset='flight', detector='hrc', level='0', filetype='hrchk')

        #self.assertEquals(data, comp_list2)
        #tdir = exc_dir + 'Temp_dir'
        #clean_dir(tdir)

#------------------------------------------------------------

    def test_conv_time_format(self):

        year  = 2014
        month = 5
        out   = conv_time_format(year, month, next = 0)

        self.assertEquals(out, '05/01/14,00:00:00')
        
        year  = 2014
        month = 5
        out   = conv_time_format(year, month, next = 1)
        self.assertEquals(out, '06/01/14,00:00:00')
        
#------------------------------------------------------------

    def test_convertto1998sec(self):

        ftime = '01/01/15,00:00:00'
        out  = convertto1998sec(ftime)

        self.assertEquals(out, 536457600)

        ftime = '2015-01-01,00:00:00'
        out  = convertto1998sec(ftime)

        self.assertEquals(out, 536457600)

#------------------------------------------------------------

    def test_covertfrom1998sec(self):

        stime = 544194930
        ftime = covertfrom1998sec(stime)

        self.assertEquals(ftime, '03/31/15,14:15:30')

#------------------------------------------------------------

    def test_dir_file(self):

        test_list = ['hrc_i_115', 'hrc_i_90', 'hrc_s_125_1', 'hrc_s_125_2', 'hrc_s_125_hi_1', \
               'hrc_s_125_hi_2', 'hrc_s_90_1', 'hrc_s_90_2', 'hrc_s_90_hi_1', 'hrc_s_90_hi_2'] 

        cmd = 'mkdir ./Test_in'
        os.system(cmd)
        cmd = 'cp ' + house_keeping + 'Selection_coditions/hrc_* ./Test_in/.'
        os.system(cmd)
#
#--- testing get_file_list
#
        out = get_file_list('./Test_in')
        self.assertEquals(out, test_list)
#
#--- testing clean_dir
#
        clean_dir('./Test_in')

        chk = 0
        if os.listdir('./Test_in'):
            chk = 1

        self.assertEquals(chk, 0)

        cmd =  'rm -rf ./Test_in'
        os.system(cmd)

#------------------------------------------------------------

    def test_find_month(self):

        lmon = 5
        out  = find_month(lmon)
        self.assertEquals(out, 'MAY')

        lmon = 'May'
        out  = find_month(lmon)
        self.assertEquals(out, 5)

#------------------------------------------------------------

    def test_read_header_value(self):

        fits = house_keeping + 'Test_data/hrcf567118652N001_hk0.fits.gz'
        name = 'DATAMODE'
        val = read_header_value(fits, name)
        
        self.assertEquals(val, 'NEXT_IN_LINE')

#------------------------------------------------------------

    def test_remove_duplicate_from_list(self):

        test   = [1, 3, 2, 6, 4, 5]
        result = remove_duplicate_from_list(test)

        self.assertEquals(result, [1, 2, 3, 4, 5, 6])

#------------------------------------------------------------

    def test_remove_duplicate_by_column(self):

        alist  = ('3 1 3234 5', '2 23342 3 1', '1 33 22 33')
        out   = remove_duplicate_by_column(alist, 0)
        self.assertEquals(out,  ['1 33 22 33', '2 23342 3 1', '3 1 3234 5'])

        alist  = ([3,1,3234,5], [2,23342,3,1], [1,33,22,33])
        out   = remove_duplicate_by_column(alist, 0)
        self.assertEquals(out,  [[1, 33, 22, 33], [2, 23342, 3, 1], [3, 1, 3234, 5]])

#------------------------------------------------------------

    def test_combine_and_select_data(self):

        hk_col_list = ['TIME','LLDIALV','RSRFALV','WDTHAST','ULDIALV','CBHUAST','CBHVAST','CBLUAST', 'HVPSSTAT']

        lpart   = numpy.array( [False, True,False, True, True,False,False,False])
        compare =  [131.0, 185.0, 2.0, 255.0, 97.0, 97.0, 90.0, lpart]

        tdir  = house_keeping + 'Test_data/' 
        flist = []
        for ent in ('hrcf567118652N001_hk0.fits.gz', 'hrcf567189107N001_hk0.fits.gz'):
            efile = tdir + ent
            flist.append(efile)

        dout = combine_and_select_data(flist, hk_col_list)

        alist = list(dout[0])
        tlist = alist[1:8]
        self.assertEquals(tlist, compare[:-1])
        llist = alist[-1]
        llist = list(llist)
        clist = list(compare[-1])
        self.assertEquals(llist, clist)

        #test = dout.field('HVPSSTAT')
        #print "I AM HERE TEST: "  + str(test)

#------------------------------------------------------------

    def test_convert_time_format(self):

        stime = '01/01/15,00:00:00'
        otime = convert_time_format(stime)

        self.assertEquals(otime, '2015-01-01T00:00:00')

        stime = 553651195.0
        otime = convert_time_format(stime)
        self.assertEquals(otime, stime)

#------------------------------------------------------------
#
#    def test_find_word_with_wild_card(self):
#
#        ex_word = '*a.*b'
#        tlist   =['a.b', 'cad.b', 'cca.dadb', 'sscab.dd']
#
#        out = find_word_with_wild_card(ex_word, tlist)
#
#        print "I AM HERE: " + str(out)
#
#-----------------------------------------------------------------------------------------------
#
if __name__ == "__main__":

    run_test()


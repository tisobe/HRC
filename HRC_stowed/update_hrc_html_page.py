#!/usr/bin/env /proj/sot/ska/bin/python

#############################################################################################
#                                                                                           #
#       update_hrc_html_page.py: update HRC Stowed background html pages                    #
#                                                                                           #
#               author: t. isobe (tisobe@cfa.harvard.edu)                                   #
#                                                                                           #
#               last update: Jan 17, 2017                                                   #
#                                                                                           #
#############################################################################################

import os
import os.path
import sys
import re
import string
import random
import time
from datetime import datetime
import numpy
import locale

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
#
#--- temp writing file name
#
rtail  = int(time.time())
zspace = '/tmp/zspace' + str(rtail) + 'hrc'
#

template_list = ['hrc_i_115_year_template', 'hrc_s_125_hi_year_template', 'hrc_s_125_year_template']

locale.setlocale(locale.LC_ALL, 'en_US')
#--- set base time at 1998.1.1
#
BTFMT    = '%m/%d/%y,%H:%M:%S'
basetime = datetime.strptime('01/01/98,00:00:00', BTFMT)
Ebase_t  = time.mktime((1998, 1, 1, 0, 0, 0, 5, 1, 0))
#
#-----------------------------------------------------------------------------------------------
#-- create_html_pages: update HRC Stowed background html pages                               ---
#-----------------------------------------------------------------------------------------------

def create_html_pages(year):
    """
    update HRC Stowed background html pages
    input:  year    --- indicator of which year to be updated. if it is "all", update all years
    output: updated/newly created html pages
    """
#
#--- create page for indivisual years
#
    stime = tcnv.currentTime()
    tyear = stime[0]
    if year == 'all':
        eyear = stime[0] + 1
        for syear in range(2000, eyear):
            create_page(syear, tyear)
    
    elif year == "":
        stime = tcnv.currentTime()
        eyear = tyear
        create_page(eyear, tyear)
    
    else:
#
#--- if the current year is the same as the given year, it is mid year update; so tyear = year+1
#
        if tyear == year:
            tyear += 1
        create_page(year, tyear)
#
#--- crate main pages: hrc_i_115_main.html, hrc_s_125_hi_main.html, and hrc_s125_main.html
#
    print_main_page(tyear)
#
#--- update hrc i image correction page
#
    update_image_correction_page(tyear)
#
#--- update inserting html pages
#
    create_slide_map_pages()
#
#--- copy data files to html dir data holder
#
    copy_data()
#
#--- if a new gain file is added, update html gain list and data_selection.html page
#
    updata_data_selection_page()

#-----------------------------------------------------------------------------------------------
#-- create_page: create the html page for indivisual year for all three setting              ---
#-----------------------------------------------------------------------------------------------

def create_page(year, tyear):
    """
    create the html page for indivisual year for all three setting
    input:  year    --- the year that the page is created
            tyear   --- the current year
    output: the html pages of the given year
    """
    tdir = house_keeping + 'Templates/'
#
#--- first check whether the data actuall exists
#
    file_chk = check_data_exisitance(year)
#
#--- create pages for hrc_i_115, hrc_s_125_hi, and hrc_s_125
#
    for k in range(0, 3):
        if file_chk[k] > 0:
            template = template_list[k]
            outfile  = html_dir + 'Yearly/' +  template
#
#--- if there is no data, use a no_data_template 
#
        else:
            if k == 0:
                template = 'no_data_template1'
            elif k == 1:
                template = 'no_data_template2'
            else:
                template = 'no_data_template3'

            outfile = template_list[k]
            outfile = html_dir + 'Yearly/' + outfile
#
#--- create previous/next year button
#
        pdirection  = create_direct_button(year, tyear, k)
#
#--- set outputfile name
#
        infile  = tdir + template
        subl    = str(year) + '.html'
        outfile = outfile.replace('_template', subl)
#
#--- read the template file
#
        f    = open(infile, 'r')
        data = f.read()
        f.close()
#
#--- add (substitute) the data
#
        data = data.replace("#YEAR#", str(year))
        data = data.replace("#DLINE#", pdirection)

        fo   = open(outfile, 'w')
        fo.write(data)
        fo.close()

#-----------------------------------------------------------------------------------------------
#-- create_direct_button: create previous year and next year link buttons                     --
#-----------------------------------------------------------------------------------------------

def create_direct_button(year, tyear, k):
    """
    create previous year and next year link buttons
    input:  year    --- the year of the page
            tyear   --- the current year
            k       --- the indicator of hrc_i_115, hrc_s_125_hi, or hrc_s_125 (in order: 0, 1, 2)
    output: line    --- html line with previous and next year buttons
    """

    if k == 0:
        inst = "hrc_i_115"
    elif k == 1:
        inst = "hrc_s_125_hi"
    else:
        inst = "hrc_s_125"

    if year > 2000:
        lyear = year -1

        pline = '<a href="' + inst + '_year' + str(lyear) + '.html">&lt;&lt;Prev</a>'
    else:
        pline = "&#160;"

    if year == tyear:
        nline = "&#160;"
    else:
        nyear = year + 1
        nline = '<a href="' + inst + '_year' + str(nyear) + '.html">Next&gt;&gt;</a>'

    line = '<th>' +  pline + '</th>\n<th>' + nline + '</th>'

    return line


#-----------------------------------------------------------------------------------------------
#-- check_data_exisitance: check whether the input data exist for the given year              --
#-----------------------------------------------------------------------------------------------

def check_data_exisitance(year):
    """
    check whether the input data exist for the given year
    input:  year    --- the year you want to check
    output: [ans1, ans2, ans3] --- a list of indicator 0 (no data) or 1 (yes) for
                                   hrc_i_115, hrc_s_125_hi , and hrc_s_125
    """

    file1 = data_dir + 'Hrc_i_115/hrc_i_115_' + str(year) + '_evt1.fits.gz'
    file2 = data_dir + 'Hrc_s_125_hi/hrc_s_125_hi_' + str(year) + '_evt1.fits.gz'
    file3 = data_dir + 'Hrc_s_125/hrc_s_125_' + str(year) + '_evt1.fits.gz'

    ans1 = 0
    if os.path.isfile(file1):
        ans1 = 1

    ans2 = 0
    if os.path.isfile(file2):
        ans2 = 1
    ans3 = 0
    if os.path.isfile(file3):
        ans3 = 1

    return [ans1, ans2, ans3]

#-----------------------------------------------------------------------------------------------
#-- print_main_page: create three main map pages                                              --
#-----------------------------------------------------------------------------------------------

def print_main_page(tyear):
    """
    create three main map pages
    input:  tyear   --- the current year
    output: hrc_stowed_position_study.html/hrc_i_115_main.html / hrc_s_125_hi_main.html / hrc_s_125_main.html
    """

    tdir = house_keeping + 'Templates/'
#
#--- update the top page
#
    file = tdir + 'hrc_stowed_position_study_template'
    f    = open(file, 'r')
    data = f.read()
    current_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())
    data = data.replace("#UPDATE#", current_time)

    file = html_dir + 'hrc_stowed_position_study.html'
    fo   = open(file, 'w')
    fo.write(data)
    fo.close()
#
#--- print three stowed map pages
#

    for m in range(0, 3):
        tname = template_list[m]
#
#--- for hrc s 125 has three sections, but the exposure time are all same
#
        if m == 2:
            dname = tname.replace('year_template', 'sec1_dead_time')
        else:
            dname = tname.replace('year_template', 'dead_time')
#
#--- read cumulative dead time corrected exposure time for each year
#
        [y_list, e_list] = read_dead_time_list(dname, tyear)
#
#--- read total events counts
#
        ecnt_list        = read_tot_evnt_cnt(y_list, tname)
#
#--- create atable row entry for each year
#
        ycnt = tyear - 2000
        line = ''
        for k in range(0, ycnt):
            hpart = tname.replace('year_template', '')
            try:
                ypart = str(y_list[k])
            except:
                break

            html  = './Yearly/'        + hpart + 'year'  + ypart + '.html'
            hdir  = tname.replace('_year_template', '')
            hdir  = hdir.replace('hrc', 'Hrc')

            if m == 2:
                slist = ['1', '2', '3']
            else:
                slist = ['']
#
#--- every 10 yrs, add scale of the map
#
            if k % 10  == 0:
                scale = './Maps/' + hdir + '/scale.png'
                line  = line + '<tr><th colspan=6><img src="' + scale + '"></th></tr>'

            line = line + '<tr><th><a href="' + html + '">' + str(y_list[k]) + '</a></th>\n'
            line = line + '<th>' + str(e_list[k]) +  '<br />' + str(ecnt_list[k])  + '</th>\n'
#
#--- evt1 with peak region enhanced / evt1 with the uniform scale, / evt2 with peak region / evt2 with the uniform scale
#
            line = line + make_table_line(html, hdir, hpart, ypart, slist, '',      '_b')
            line = line + make_table_line(html, hdir, hpart, ypart, slist, '',      '')
            line = line + make_table_line(html, hdir, hpart, ypart, slist, 'lev2_', '_b')
            line = line + make_table_line(html, hdir, hpart, ypart, slist, 'lev2_', '')
#
#--- add the cumulative maps
#
        esum = 0
        for ent in e_list:
            ent = ent.replace(',', '')
            esum += int(float(ent))

        esum = locale.format("%d", esum, grouping=True)

        ecnt = 0
        for ent in ecnt_list:
            ent = ent.replace(',', '')
            ecnt += int(float(ent))

        ecnt = locale.format("%d", ecnt, grouping=True)

        ypart = 'total'
        html  = './Yearly/'        + hpart + 'total.html'
        hdir  = tname.replace('_year_template', '')
        hdir  = hdir.replace('hrc', 'Hrc')

        line  = line + '\n<tr><th colspan=5>Cumulative Maps</th></tr>\n'

        scale = './Maps/' + hdir + '/scale.png'
        line  = line + '<tr><th colspan=6><img src="' + scale + '"></th></tr>'

        line = line + '<tr><th><a href="' + html + '">Total</a></th>\n'
        line = line + '<th>' + str(esum)+ '<br />' + str(ecnt) + '</th>\n'

        line = line + make_table_line(html, hdir, hpart, ypart, slist, '',      '_b')
        line = line + make_table_line(html, hdir, hpart, ypart, slist, '',      '')
        line = line + make_table_line(html, hdir, hpart, ypart, slist, 'lev2_', '_b')
        line = line + make_table_line(html, hdir, hpart, ypart, slist, 'lev2_', '')
#
#--- read the template
#
        mfile = tdir + hpart + 'main_template'
        f     = open(mfile, 'r')
        fdata = f.read()
        f.close()
#
#--- add (substitute) table rows to the template 
#
        fdata = fdata.replace('#TABLE#', line)
#
#--- save the result
#
        outname = html_dir + hpart + 'main.html'
        fo      = open(outname, 'w')
        fo.write(fdata)
        fo.close()


#-----------------------------------------------------------------------------------------------
#-- make_table_line: make hrc main page table entries  to display exposure maps               --
#-----------------------------------------------------------------------------------------------

def make_table_line(html, hdir, hpart, ypart, slist,  lev, part):
    """
    make hrc main page table entries  to display exposure maps
    input:  html    --- a html page linked from this line
            hdir    --- a which directory the maps are kept (e.g. Hrc_i_115)
            hpart   --- a header part of the map (e.g., hrc_i_115)
            ypart   --- a year of the map
            lev     --- a lev of the map, either "" or "lev2_"
            part    --- an indicator of whether peak area or standard cropping. _b for the peak area cropping
    """

    line = '<th><a href="' + html  + '">'
    for sec in slist:
        ifile = './Maps/' + hdir + '/Simage/' + hpart +  ypart + '_thumb' + sec + part + '.png'
        line  = line + '<img src="' + ifile  + '">'
    line = line + '</a></th>\n'

    return line
#-----------------------------------------------------------------------------------------------
#-- update_image_correction_page: update hrci_image_correction.html page                      --
#-----------------------------------------------------------------------------------------------

def update_image_correction_page(year):
    """
    update hrci_image_correction.html page
    input;  year    --- the year
    output: <html_dir>/hrci_image_correction.html
    """
    [y_list, e_list] = read_dead_time_list('hrc_i_115_dead_time', year)
    evt_cnt          = read_tot_evnt_cnt(y_list, 'hrc_i_115_year_template')

    line = ''
    for k in range(0, len(y_list)):
        line = line + '<tr>\n'
        line = line + '<th>' + str(y_list[k]) + '</th>\n'
        line = line + '<td style="text-align:center"><a href="./Data/Hrc_i_115/hrc_i_115_' + str(y_list[k]) 
        line = line + '_evt1.fits.gz">hrc_i_115_' + str(y_list[k]) + ' evt1 file </a></td>\n'
        line = line + '<td style="text-align:center">' + str(e_list[k]) + '</td>\n'
        line = line + '<td style="text-align:center">' + str(evt_cnt[k]) + '</td>\n'
        line = line + '</tr>\n\n'

    file = house_keeping + 'Templates/hrci_image_correction_template'
    f    = open(file, 'r')
    data = f.read()
    f.close()

    out  = data.replace("#TABLE#", line)
    
    file = html_dir + 'hrci_image_correction.html'
    fo   = open(file, 'w')
    fo.write(out)
    fo.close()

#-----------------------------------------------------------------------------------------------
#-- read_dead_time_list: a name of dead time list and return dead time corrected exposure time -
#-----------------------------------------------------------------------------------------------

def read_dead_time_list(dead_list_name, tyear):
    """
    read dead time list 
    input:  dead_list_name      --- a name of dead time list and return dead time corrected exposure time
            tyear               --- current year
    output: [y_list, e_list]    --- a list of years and a list of dead time corrected exposure time
    """

    y_list = []
    e_list = []

    ifile = data_dir + 'Stats/' + dead_list_name
    data  = hcf.read_file_data(ifile)

#
#--- accumulate dead time corrected exposure time for each year
#
    prev = 2000
    esum = 0
    for ent in data:
        atemp = re.split('\s+', ent)
        btemp = re.split('-', atemp[1])
        year  = int(float(btemp[0]))

        try:
            val = int(float(atemp[4]))
        except:
            continue

        if year == prev:
            esum += val
        else:
            y_list.append(prev)
            dval = locale.format("%d", esum, grouping=True)
            e_list.append(dval)
            prev = year
            esum = val
#
#--- save the last entry year data (it can be a partial year)
#
    if val != 0:
        y_list.append(prev)
        dval = locale.format("%d", esum, grouping=True)
        e_list.append(dval)

    return [y_list, e_list]

#-----------------------------------------------------------------------------------------------
#-- read_tot_evnt_cnt: create a list of the total numbers of events for the year              --
#-----------------------------------------------------------------------------------------------

def read_tot_evnt_cnt(y_list, tname):
    """
    create a list of the total numbers of events for the year
    (start from 2000)
    input:  y_list  ---- a list of year which we want to extract data
            tname   ---- instrument template name
    output: ecnt    ---- a list of total event counts starting from year 2000
    """

    atemp  = re.split('_year', tname)
    infile = data_dir + atemp[0].capitalize() + '/' + atemp[0] + '_yearly_evt_counts'
    f      = open(infile, 'r')
    data   = [line.strip() for line in f.readlines()]
    f.close()

    ecnt   = []
    for ent in data:
        atemp = re.split('\s+', ent)
        year  = int(float(atemp[0]))
        if year in y_list:
            val   = int(float(atemp[1]))
            dval  = locale.format("%d", val, grouping=True)
            ecnt.append(dval)

    return ecnt

#-----------------------------------------------------------------------------------------------
#-- update_stat_data_table: create html pages to display stat results                        ---
#-----------------------------------------------------------------------------------------------

def update_stat_data_table():
    """
    create html pages to display stat results
    input:  none, but read from <data_dir>/Stats/<head>_stat_results
    output: <html_dir><head>_stats.html
    Table columns:
            0   --- Date
            1   --- Tstart
            2   --- Tstop
            3   --- Duration
            4   --- Total Counts
            5   --- Counts per Sec
            6   --- Mean of PHA Peak Position
            7   --- Median of PHA Peak Position
            8   --- Sigma of PHA Peak Position
            9   --- Mean of MCP Valid Count Rate
            10  --- Median of MCP Valid Count Rate
            11  --- Sigma of MCP Valid Count Rate
            12  --- Mean of MCP Total Count Rate
            13  --- Median of MCP Total Count Rate
            14  --- Sigma of MCP Total Count Rate
            15  --- Mean of MCP Shield Count Rate
            16  --- Median of MCP Shield Count Rate
            17  --- Sigma of MCP Shield Count Rate
            18  --- Mean of Valid/Total MCP Rate
            19  --- Median of Valid/Total MCP Rate
            20  --- Sigma of  Valid/Total MCP Rate
            21  --- S2HVST
            22  --- S2HVLV
            23  --- EPHIN Integrated Flux
    """
#
#--- the table will display values with either integer or float.
#--- the following indicies indicate how to round the value. if it is 0, integer is used
#---- pos : [0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]
#
    r_ind = [0, 0, 0, 0, 5, 2, 2, 0, 3, 3, 0, 3, 3, 0, 3, 3, 0, 3, 3, 3, 3, 0, 0, 4]
#
#--- read a template file
#
    t_file   = house_keeping + 'Templates/stat_table_template'
    f        = open(t_file, 'r')
    template = f.read()
    f.close()
#
#--- go around three hrc device. hrc_s_125 has 3 sections
#
    for hname in ['hrc_i_115', 'hrc_s_125_hi', 'hrc_s_125']:
        if hname == 'hrc_s_125':
            s_list = ['_sec1_', '_sec2_', '_sec3_']
        else:
            s_list = ['_']

        for sec in s_list:
#
#--- make data table
#
            ename = 'Stats/' + hname + sec + 'stat_results'
            fname = data_dir + ename
            data  = hcf.read_file_data(fname)
            line = ''
            for ent in data:
                val = re.split('\s+', ent)
                line = line + '\t<tr>\n'
                line = line + '\t\t<td>' + val[0] + '</td>\n'
                for k in range(1, 24):

                    try:
                        fval = float(val[k])
                        if r_ind[k] == 0:
                            fval = str(int(fval))
                        else:
                            fval = str(round(fval, r_ind[k]))
                    except:
                        fval = 'NaN'

                    line = line + '\t\t<td>' + fval + '</td>\n'
#
#--- set output file name
#
            oname  = html_dir + hname + sec +  'stats.html'
            title  = hname + sec
            title  = title.replace('_', ' ')
            title2 = title.upper()
#
#--- a link to the ascii data
#
            adata  = '<a href="./Data/' + ename + '">Open Plain Text Version</a>'
#
#--- today's date
#
            tdate = time.strftime("%a, %d %b %Y", time.gmtime())
#
#--- insert the information to the template and print it out
#
            out = template
            out = out.replace('#HTITLE#', title)
            out = out.replace('#MTITLE#', title2)
            out = out.replace('#ADATA#',  adata)
            out = out.replace('#TABLE#',  line)
            out = out.replace('#UPDATE#', tdate)

            fo  = open(oname, 'w')
            fo.write(out)
            fo.close()

#-----------------------------------------------------------------------------------------------
#-- create_slide_map_pages: create three sliding insert html pages                            --
#-----------------------------------------------------------------------------------------------

def create_slide_map_pages():
    """
    create three sliding insert html pages
    input:  none
    output: hrc_i_115_slide.html, hrc_s_125_slide.html, and hrc_s_125_hi_slide.html in <html_dir>
    """
    
    for name in ['Hrc_i_115', 'Hrc_s_125_hi', 'Hrc_s_125']:
        lname = name.lower()
#
#--- first check whether there are any years which do not have any data
#
        cmd   = 'ls ' + data_dir + name + '/*evt1* > ' +zspace
        os.system(cmd)

        data = hcf.read_file_data(zspace, remove=1)
        y_list = []
        for ent in data:
            atemp = re.split('_', ent)
            y_list.append(int(float(atemp[-2])))

        tarray = set(numpy.array(y_list))
        y_list = list(tarray)
        y_list.sort()

        stime = tcnv.currentTime()
        tyear = stime[0]
        lyear = tyear -1

        ymax  = lyear
        if ymax < y_list[-1]:
            ymax = y_list[-1]
#
#--- pcant is used to set the width of the table in pixel unit
#--- each year, 70 pixels are asigned for the width
#
        pcnt  = 80 * (ymax - 2000)

        outname = html_dir + lname + '_slide.html'

        line = '<!DOCTYPE html> \n'
        line = line + '<html> \n'
        line = line + '<head>\n\t <title>HRC Stowed Background Study ' + name 
        line = line + ' Map Display</title>\n</head>\n'
        line = line + '<body>\n'
        line = line + '<table border=0 style="width:' + str(pcnt) + 'px">\n'
        line = line + '<tr> \n'

        for kyear in range(2000, ymax+1):
            line = line + '<th>' +  str(kyear) + '</th>\n'

        line = line + '</tr>\n'
        line = line + '<tr>\n'

        kmax = ymax - 2000
        for k in range(0, kmax+1):
            kyear = 2000 + k
            if kyear in y_list:
#
#--- hrc s 125 has three sections
#
                if name == 'Hrc_s_125':
                    line = line + '<th style="width:70px">\n'
                    line = line + '\t<a href="./Yearly/' + lname + '_year' + str(kyear) + '.html"'
                    line = line + ' target="blank">\n'
                    for sec in range(1,4):
                        line = line + '\t\t<img src="./Maps/' + name + '/Simage/' + lname+ '_' 
                        line = line +  str(kyear) + '_thumb' + str(sec) + '.png">\n'
                    line = line + '\n\t</a>\n</th>\n'
#
#--- hrc s 125 hi and hrc i 115 cases
#
                else:
                    if name == 'Hrc_s_125_hi':
                        line = line + '<th style="width:70px">\n'
                    else:
                        line = line + '<th>\n'
                    line = line + '\t<a href="./Yearly/' + lname + '_year' + str(kyear) + '.html"'
                    line = line + ' target="blank">\n'
                    line = line + '\t\t<img src="./Maps/' + name + '/Simage/' + lname+ '_' 
                    line = line +  str(kyear) + '_thumb.png">\n\t</a>\n</th>\n'
#
#--- if there is no data for the year, say so
#
            else:
                line = line + '<th>No<br /> Data</th>\n'

        line = line + '</tr>\n'
        line = line + '</table>\n'
        line = line + '</body>\n'
        line = line + '</html>\n'

        fo = open(outname, 'w')
        fo.write(line)
        fo.close()

#-----------------------------------------------------------------------------------------------
#-- copy_data: copy stat results files and data files from the original location to html data directory 
#-----------------------------------------------------------------------------------------------

def copy_data():
    """
    copy stat results files and data files from the original location to html data directory
    input: origial data in <data_dir>/Stats/ and <data_dir>/<hrc inst>/*fits.gz
    output: stat result files and data files in <html_dir>/Data/...
    """
#
#--- copy stat files
#
    cmd = 'cp ' + data_dir + 'Stats/* ' + html_dir + 'Data/Stats/.'
    os.system(cmd)
#
#--- copy data fits files
#
    for hdir in ['Hrc_i_115', 'Hrc_s_125', 'Hrc_s_125_hi']:

        odir = data_dir + hdir + '/'
        hdir = html_dir + 'Data/' + hdir + '/'

        cmd   = 'ls ' + odir + '*.fits.gz > ' + zspace
        os.system(cmd)
        data  = hcf.read_file_data(zspace, remove=1)
#
#--- remove a directory path to the file listed
#
        data1 = get_the_last_part(data)

        cmd   = 'ls ' + hdir + '*.fits.gz > ' + zspace
        os.system(cmd)
        data  = hcf.read_file_data(zspace, remove=1)
        data2 = get_the_last_part(data)
#
#--- check which files are not in the hrml data directory
#
        missing = list(set(data1) - set(data2))
#
#--- copy only the missing files
#
        for ent in missing:
            cmd = 'cp ' +  odir + '/' + ent + ' ' + hdir + '/.'
            os.system(cmd)

#-----------------------------------------------------------------------------------------------
#-- get_the_last_part: remove a directory path and return a list of file names                --
#-----------------------------------------------------------------------------------------------

def get_the_last_part(data):
    """
    remove a directory path and return a list of file names
    input   data    --- a list of data with a full path
    outpu   clpped  --- a list of data withtout a directory path
    """

    if len(data) == 0:
        return []

    mc = re.search('\/', data[0])
    if mc is None:
        return data
    else:
        clipped = []
        for ent in data:
            atemp = re.split('\/', ent)
            clipped.append(atemp[-1])

        return clipped

#-----------------------------------------------------------------------------------------------
#-- update_data_selection_page: update html gain file information and data_selection.html page --
#-----------------------------------------------------------------------------------------------

def updata_data_selection_page():
    """
    update html gain file information and data_selection.html page
    input:  none but read two gain_selection files (from <house_keeping> and <html_dir>
    output: <html_file>/Gain_files/gain_selection and gain fits file 
            <html_file>/data_selection.html
    """
#
#--- check two gain list files are different. if so, new gain file was found
#
    gain_file1 = house_keeping + 'Gain_files/gain_selection'
    gain_file2 = html_dir      + 'Data_save/Gain_files/gain_selection'

    cmd  =  'diff ' + gain_file1 + ' ' + gain_file2 + ' >' + zspace
    os.system(cmd)

    if mcf.isFileEmpty(zspace) == 0:
        mcf.rm_file(zspace)

    else:
        mcf.rm_file(zspace)
#
#--- find which gain files are added 
#
        data   = hcf.read_file_data(gain_file1)
        flist1 = []
        for ent in data:
            atemp = re.split('\s+', ent)
            flist1.append(atemp[-1])

        data2 = hcf.read_file_data(gain_file2)
        flist2 = []
        for ent in data2:
            atemp = re.split('\s+', ent)
            flist2.append(atemp[-1])
#
#--- add new gain file to <html_dir>/Data_save/Gain_files
#
        missing = list(set(flist1) - set(flist2))
        
        for ent in missing:
            cmd = 'cp -f  ' + house_keeping + 'Gain_files/' + ent + ' ' + html_dir + 'Data_save/Gain_files/.'
            os.system(cmd)

        cmd = 'cp  -f ' + gain_file1 + ' ' + gain_file2
        os.system(cmd)
#
#--- update data_selection.html page
#
        line = ''
#
#--- create the table of the gain file list
#
        for ent in data:
            atemp = re.split('\s+', ent)
            line = line + '<tr>\n\t<th>'
            if atemp[0] != '0':
                dtime1 = gformat_time(atemp[0])
                line   = line + 'After ' + str(dtime1)
    
            if atemp[1] != '1.0e12':
                dtime2 = gformat_time(atemp[1])
                line   = line + ' Before ' + str(dtime2)

            line = line + '</th>\n\t<td><a href="./Data_sve/Gain_file/' + atemp[2] + '">' + atemp[2] + '</a></td>\n</tr>\n'
#
#--- find the current date
#
        current_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())
#
#--- read the template and substitute the table and update date
#
        template   = house_keeping + 'Templates/data_selection_template'
        f    = open(template, 'r')
        data = f.read()
        f.close()

        data = data.replace("#TABLE#", line)
        data = data.replace("#UPDATE#", current_time)
#
#--- print out the file
#
        out  = html_dir + 'data_selection.html'
        fo   = open(out, 'w')
        fo.write(data)
        fo.close()


#-----------------------------------------------------------------------------------------------
#-- gformat_time: convert time in seconds from 1998.1.1 to <mmm> <dd>, <yyyy> format         ---
#-----------------------------------------------------------------------------------------------

def gformat_time(stime):
    """
    convert time in seconds from 1998.1.1 to <mmm> <dd>, <yyyy> format
    input:  stime   --- time in seconds from 1998.1.1
    output: etime   --- time in the format of <mmm> <dd>, <yyyy>
    """
    etime = Ebase_t + float(stime)
    etime = time.localtime(etime)
    etime = time.strftime('%b %d, %Y', etime)

    return etime

#-----------------------------------------------------------------------------------------------
 
if __name__ == "__main__":

    test = 0
    if len(sys.argv) == 1:
        create_html_pages('')
#
#--- for the case, you want to update all pages
#
    elif len(sys.argv) == 2:
        if sys.argv[1] == 'all':
            create_html_pages('all')
            update_stat_data_table()
#
#--- for the case, you want to update only a specified year
#
        else:
            year = int(float(sys.argv[1]))
            create_html_pages(year)
            update_stat_data_table()
#
#--- for the case,  you want to update only main map paged
#
    elif len(sys.argv) == 3:
        if sys.argv[1] == 'bmp':
            tyear = int(float(sys.argv[2]))
            print_main_page(tyear)

    else:
        exit(1)

#!/usr/bin/env /proj/sot/ska/bin/python

#################################################################################
#                                                                               #
#   create_html_page.py: create hz43.html page                                  #
#                                                                               #
#       author: t. isobe (tisobe@cfa.harvard.edu)                               #
#                                                                               #
#       last update: May 10, 2017                                               #
#                                                                               #
#################################################################################

import os
import sys
import re
import string
import time

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
#--- append a path to a private folder to python directory
#
sys.path.append(bin_dir)
sys.path.append(mta_dir)
#
#--- converTimeFormat contains MTA time conversion routines
#
import convertTimeFormat    as tcnv
import mta_common_functions as mcf

f_list  = ['samp_n_s', 'samp_p_s', 'pi_n_s', 'pi_p_s']
f_list2 = ['samp_n_i', 'samp_p_i', 'pi_n_i', 'pi_p_i']

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

def create_html_page():
#
#--- read fitting results and extract slope data; put in a dictionary form
#
    [s_dict, s_dict2]  = read_slope()
#
#--- create the table part
#
    text = "<table border=1 cellpadding=2>\n"
    text = text + "<tr>\n"
    text = text + "<th>&#160;</th><th colspan=4>Scaled Sum Amp</th><th colspan=4>PI</th>\n"
    text = text + "</tr>\n"
    text = text + "<tr>\n"
    text = text + "<th>&#160;</th><th colspan=2>Negative</th><th colspan=2>Positive</th>\n"
    text = text + "<th colspan=2>Negative</th><th colspan=2>Positive</th>\n"
    text = text + "<tr>\n"
#
#--- HRC I  ---- only position 3 (70-90 A) will be displayed
#
    for k in range(3, 4):
        start = k * 20 + 10
        stop  = start + 20

        text2 = text  + "<tr>\n"
        ###text2 = text2 + "<th>" +  str(start) + " - " + str(stop) + " &#8491</th>\n"
        text2 = text2 + "<th>60 - 80&#8491</th>\n"
        
        for m in range(0, 4):
#
#--- HRC I have a bit strange slope value saving; so we need to use a trick
#
            if len(s_dict[f_list2[m]]) > 1:
                kpos = 1
            else:
                kpos = 0
            text2 = text2 + create_table_cell(k, s_dict, s_dict2,  f_list2, m, kpos)

        text2 = text2 + '</tr>\n'
    text2 = text2 + "</table>\n"
#
#--- HRC S
#
    for k in range(5, 18):
        start = k * 10
        stop  = start + 10

        text = text + "<tr>\n"
        text = text + "<th>" +  str(start) + " - " + str(stop) + " &#8491</th>\n"
        
        for m in range(0, 4):
            text = text + create_table_cell(k, s_dict, s_dict2,  f_list, m)
    
        text = text + '</tr>\n'
    text = text + "</table>\n"

#
#--- read the template
#
    f = open("/data/aschrc6/wilton/isobe/Project8/HZ43/Scripts/house_keeping/hz43_template", 'r')
    
    html = f.read()
    f.close()
#
#--- modified time
#
    update = time.strftime("%b %d, %Y", time.localtime())

#
#--- insert the table and update the html page
#
    html = html.replace('#TABLE2#', text2) 
    html = html.replace('#TABLE#',  text) 
    html = html.replace('#UPDATE#', update)

    text = s_dict['samp_s_center'][0] + '<br />' + s_dict2['samp_s_center'][0]
    html = html.replace('#HRCSSAMP#', text)

    text = s_dict['samp_i_center'][0] + '<br />' + s_dict2['samp_i_center'][0]
    html = html.replace('#HRCISAMP#', text)

    text = s_dict['pi_s_center'][0] + '<br />' + s_dict2['pi_s_center'][0]
    html = html.replace('#HRCSPI#', text)

    text = s_dict['pi_i_center'][0] + '<br />' + s_dict2['pi_i_center'][0]
    html = html.replace('#HRCIPI#', text)

    fo = open("/proj/web-cxc-dmz/htdocs/contrib/cxchrc/HRC_trendings/HZ43/hz43.html", 'w')
    fo.write(html)
    fo.close()

#-----------------------------------------------------------------------------------------------
#-- create_table_cell: create a row of the table for the main page                            --
#-----------------------------------------------------------------------------------------------

def create_table_cell(kpos, s_dict, s_dict2,  flist, m, kalt=''):
    """
    create a row of the table for the main page
    input:  kpos    --- cell position
            s_dict  --- a dictionary which contains the slope value
            s_dict2 --- a dictionary which contains the second slope value
            flist   --- a list of the names of the dict index
            m       --- a position of the index name in flist
            kalt    --- there is an occasion that the slope position and cell position are different;
                        if that is the case, use this
    output: text    --- a row fo the table
    """

    if kalt == '':
        kadj = kpos -5                          #--- data table start the index of 0 which correpsond HRC S 5th
    else:
        kadj = kalt                             #--- hrc i case
    try:
        out  =  s_dict[flist[m]][kadj]          #---- this tests whether the value exists before doing anything else

        text = "<td><a href=\"" + html_top + 'Plots/'+ flist[m] + "_" + str(kpos) + ".html\">" 
        text = text + '<img src="' + html_top + 'Plots/' + 'Thumb_plots/' + flist[m] + "_" 
        text = text + str(kpos) + '_thumb_plot.png"></a></td>\n'
        text = text + "<td><a href=\"" + html_top + 'Plots/' + flist[m] + "_" + str(kpos) + ".html\">" 
        text = text +  s_dict[flist[m]][kadj]  + '<br />'

        if s_dict2[flist[m]][kadj] != 0:
            text = text +  s_dict2[flist[m]][kadj] 

        text = text + "</a></td>"
    except:
        text = "<td>No Plot</td><td>NA</td>"

    return text

#-----------------------------------------------------------------------------------------------
#-- read_slope: read the slope for each category                                              --
#-----------------------------------------------------------------------------------------------

def read_slope():
    """
    read the slope for each category
    input:  none, but read from <data_dir><samp/pi>_<i/s>_<p/n>_fitting_resluts
    output: save    --- a dictionary with index of <samp/pi>_<i/s>_<p/n> with a list of slope values
                        the list can contain only 1 value or as many as 20 slope values.
            save2   --- a dictionary with index of <samp/pi>_<i/s>_<p/n> with a list of the secnd slope values
    """
    
    save  = {}
    save2 = {}
    for head in ['samp', 'pi']:
        for inst in ['s', 'i']:
        #for inst in ['s']:
            for sign in ['n', 'p']:
                infile = data_dir + head + '_' + inst + '_' + sign + '_fitting_results'
                name   = head + '_' + sign + '_' + inst
                [s_list, s_list2] = read_fitting_results(infile)
                save[name]  = s_list
                save2[name] = s_list2

    for head in ['samp', 'pi']:
        for inst in ['s','i']:
            infile = data_dir + head + '_' + inst + '_center_fitting_results'
            [s_list, s_list2]  = read_fitting_results(infile)
            name = head + '_' + inst + '_center'
            save[name]  = s_list
            save2[name] = s_list2

    return [save, save2]

#-----------------------------------------------------------------------------------------------
#-- read_fitting_results: read the file given and make a list of slope with the error         --
#-----------------------------------------------------------------------------------------------

def read_fitting_results(infile):
    """
    read the file given and make a list of slope with the error 
    input:  infile  --- the data file name
    output: s_list  --- a list of <slope>+/-<error>
    """

    try:
        f    = open(infile, 'r')
        data = [line.strip() for line in f.readlines()]
        f.close()
    except:
        data = []

    s_list1 = []
    s_list2 = []
    for ent in data:
        atemp = re.split('\s+', ent)
        if len(atemp) == 6:
            slope  = atemp[1] + '+/-' + atemp[2]
            slope2 = atemp[4] + '+/-' + atemp[5]
        else:
            slope  = atemp[2] + '+/-' + atemp[3]
            slope2 = atemp[5] + '+/-' + atemp[6]

        s_list1.append(slope)
        s_list2.append(slope2)

    return [s_list1, s_list2]

#-----------------------------------------------------------------------------------------------

if __name__ == "__main__":

    create_html_page()

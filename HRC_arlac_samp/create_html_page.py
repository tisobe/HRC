#!/usr/bin/env /proj/sot/ska/bin/python

#################################################################################
#                                                                               #
#   create_html_page.py: create arlac.html page                                 #
#                                                                               #
#       author: t. isobe (tisobe@cfa.harvard.edu)                               #
#                                                                               #
#       last update: Jun 01, 2017                                               #
#                                                                               #
#################################################################################

import os
import sys
import re
import string
import time

import matplotlib as mpl

if __name__ == '__main__':
    mpl.use('Agg')

from pylab import *
import matplotlib.pyplot       as plt
import matplotlib.font_manager as font_manager
import matplotlib.lines        as lines

import mpld3
from mpld3 import plugins, utils


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

f_list  = ['samp_s', 'pi_s']
f_list2 = ['samp_i', 'pi_i']

hrci_offset = [[0,0], [2,0], [0,2], [-2,0], [0,-2], [2,2], [-2,2], [-2,-2], [2,-2], [4,0], [0,4],[-4,0],\
               [0,-4], [6,0], [0,6], [-6,0], [0,-6], [10,10], [-10,10], [-10,-10], [10,-10]]

hrcs_offset = [[0,0], [2,0], [0,2], [-2,0], [0,-2], [2,2], [-2,2], [-2,-2], [2,-2], [4,0], [-4,0], [4,2],\
               [-4,2], [-4, -2], [4,-2], [6,0], [-6,0], [10,2], [-10,2], [-10,-2],[10,-2]]

point_list  = [['blue', 'o'],    ['blue', '*'],    ['blue', '^'],    ['blue', 's'],\
               ['red', 'o'],     ['red', '*'],     ['red', '^'],     ['red', 's'],\
               ['green', 'o'],   ['green', '*'],   ['green', '^'],   ['green', 's'],\
               ['fuchsia', 'o'], ['fuchsia', '*'], ['fuchsia', '^'], ['fuchsia', 's'],\
               ['lime', 'o'],    ['lime', '*'],    ['lime', '^'],    ['lime', 's'],\
               ['yellow', 'o'],  ['yellow', '*'],  ['yellow', '^'],  ['yellow', 's']]

#-----------------------------------------------------------------------------------------------
#-- create_html_page: create ArLac trend html top age                                         --
#-----------------------------------------------------------------------------------------------

def create_html_page(chk=''):
    """
    create ArLac trend html top age
    input:  none
    output: <html_dir>/arlac_trends.html
    """
#
#--- check whether new data came in
#
    ifile = house_keeping + 'hrci_new'
    sfile = house_keeping + 'hrcs_new'

    if chk == '':
        if (not os.path.isfile(ifile)) and (not os.path.isfile(sfile)):
            exit(1)

#
#--- read fitting results and extract slope data; put in a dictionary form
#
    [s_dict, s_dict2]  = read_slope()
#
#--- create hrci/hrcs plots inserts (interactive plots etc)
#
    [hrci_text, hrcs_text] = create_inserting_texts()
#
#--- create the table part
#
    text = "<table border=1 cellpadding=2>\n"
    text = text + "<tr>\n"
    text = text + "<th>Y/Z Offsets<br />(arcmin)</th><th colspan=2>Scaled Sum Amp</th><th colspan=2>PI</th>\n"
    text = text + "</tr>\n"
    text = text + "<tr>\n"
    text = text + "<tr>\n"
#
#--- hrc i table
#
    text2 = text
    for k in range(0, 21):

        text2 = text2 + "<tr>\n"
        text2 = text2 + "<th>" + str(hrci_offset[k]) + "</th>\n"
        
        for m in range(0, 2):
            text2 = text2 + create_table_cell(k, s_dict, s_dict2,  f_list2, m)

        text2 = text2 + '</tr>\n'
    text2 = text2 + "</table>\n"
#
#--- hrc s table
#
    for k in range(0, 21):
        text = text + "<tr>\n"
        text = text + "<th>[" + str(hrcs_offset[k]) + "</th>\n"
        
        for m in range(0, 2):
            text = text + create_table_cell(k, s_dict, s_dict2,  f_list, m)
    
        text = text + '</tr>\n'
    text = text + "</table>\n"

#
#--- read the template
#
    hfile = house_keeping + 'arlac_template'
    f = open(hfile, 'r')
    
    html = f.read()
    f.close()
#
#--- modified time
#
    update = time.strftime("%b %d, %Y", time.localtime())

#
#--- insert the table and update the html page
#
    html = html.replace('#HRCITEXT#', hrci_text) 
    html = html.replace('#HRCSTEXT#', hrcs_text) 
    html = html.replace('#TABLE2#', text2) 
    html = html.replace('#TABLE#',  text) 
    html = html.replace('#UPDATE#', update)

    outname = web_dir + 'arlac_trends.html'
    fo = open(outname, 'w')
    fo.write(html)
    fo.close()
#
#--- remove the new data files
#
    mcf.rm_file(ifile)
    mcf.rm_file(sfile)

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

    try:
        out  =  s_dict[flist[m]][kpos]          #---- this tests whether the value exists before doing anything else

        text = "<td><a href=\"" + html_top + 'Plots3/'+ flist[m] + "_" + str(kpos) + ".html\">" 
        text = text + '<img src="' + html_top + 'Plots3/' + 'Thumb_plots/' + flist[m] + "_" 
        text = text + str(kpos) + '_thumb_plot.png"></a></td>\n'
        text = text + "<td><a href=\"" + html_top + 'Plots3/' + flist[m] + "_" + str(kpos) + ".html\">" 
        text = text +  s_dict[flist[m]][kpos]  + '<br />'

        if s_dict2[flist[m]][kpos] != 0:
            text = text +  s_dict2[flist[m]][kpos] 

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
                infile = data_dir + head + '_' + inst + '_fitting_results'
                name   = head +  '_' + inst
                [s_list, s_list2] = read_fitting_results(infile)
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


#---------------------------------------------------------------------------------------------------
#-- create_inserting_texts: create interactive sky coordinates plots and y/z offset plots link    --
#---------------------------------------------------------------------------------------------------

def create_inserting_texts():
    """
    create interactive sky coordinates plots and y/z offset plots link
    input:  none
    output: [hrci_text, hrcs_text]: texts which contains the codes for these plot sections
    """
#
#--- hrci_offset_surface.png coordinate info
#
    hrci_cx = 255       #--- center point x pixel coordinate
    hrci_cy = 210       #--- center point y pixel coordinate
    ixstep  = 20        #--- step size in x
    iystep  = 20        #--- step size in y
#
#--- hrcs_offset_surface.png coordinate info
#
    hrcs_cx = 253
    hrcs_cy = 96 
    sxstep  = 20
    systep  = 20
#
#--- create interactive plots for hrc i
#
    infile = data_dir + 'samp_list_i'
    [xdata, ydata, obsids, inst] = read_data_for_plots(infile)
    xname  = 'Sky X'
    yname  = 'Sky Y'
    fig1  = create_interactive_plot(xdata, ydata, obsids, inst, xname, yname)

    [xdata, ydata, obsids, inst] = read_data_for_plots(infile, posx=24, posy=25)
    xname  = 'Det X'
    yname  = 'Det Y'
    fig1b = create_interactive_plot(xdata, ydata, obsids, inst, xname, yname)
#
#--- create interactive plots for hrc s
#
    infile = data_dir + 'samp_list_s'
    [xdata, ydata, obsids, inst] = read_data_for_plots(infile)
    xname  = 'Sky X'
    yname  = 'Sky y'
    fig2  = create_interactive_plot(xdata, ydata, obsids, inst, xname, yname)

    [xdata, ydata, obsids, inst] = read_data_for_plots(infile, posx=24, posy=25)
    xname  = 'Det X'
    yname  = 'Det Y'
    fig2b = create_interactive_plot(xdata, ydata, obsids, inst, xname, yname)
#
#--- hrc i
#
    out = '<p style="text-align:left;;">'
    out = out + 'The following plot shows the locations of all observations in <b>sky coordinates</b>. '
    out = out + ' The markers and the colors indicate the location of the observation in Y/Z offset sarface. '
    out = out + ' If you hover the mouse over the marker, it will display the distribution ofsamp data.</p>\n'

    out = out + '<p style="text-align:left;;">'
    out = out + 'If you like to see a detail, click the magnifier icon at the bottom left, and select ' 
    out = out + 'the area you want to see. Then click the cross icon before hover the mouse. If you '
    out = out + 'like to go back to the entier view, click the house icon.</p> \n'
#
#---- hrc i sky coordinate map
#
    out = out + '<a name="hrci_sky"></a>\n'
    out = out + '<div style="padding-left:80px;">\n'
#
#--- convert the plot into html language
#
    text = mpld3.fig_to_html(fig1)

    text = text.replace('None', '')
    out = out + text
    out = out + '</div>\n\n'

    out = out + '<p style="text-align:left;">'
    out = out + '<a href="#top">Back to Top</a></p>\n'
#
#--- det coordinates plot
#
    out = out + '<p style="text-align:left;;">'
    out = out + 'The following plot is similar to the above, but shows the locations of '
    out = out + 'all observations in <b>det coordinates</b>.</p> '
#
#---- hrc i det coordinate map
#
    out = out + '<a name="hrci_det"></a>\n'
    out = out + '<div style="padding-left:80px;">\n'
#
#--- convert the plot into html language
#
    text = mpld3.fig_to_html(fig1b)

    text = text.replace('None', '')
    out = out + text
    out = out + '</div>\n\n'

    out = out + '<p style="text-align:left;">'
    out = out + '<a href="#top">Back to Top</a></p>\n'


    out = out + '<p style="text-align:left;padding-bottom:30px;">'
    out = out + 'The following plot shows  21 <b>Y/Z Offset coordinates</b> which all observations '
    out = out + 'fall in. If you click the marker on the location, it will open the trend page of the location. '
    out = out + 'The left plot will open the scaled samp trend and the righ plot will open the PI trend page.</p>'
#
#--- hrc i y/z offset map: samp
#
    out = out + '<a name="hrci_off"></a>\n'
    out = out + '<table border=0 cellpadding=2>\n'
    out = out + '<tr><th>Scaled Samp</th><th>PI</th></tr>\n'
    out = out + '<tr><td>\n'
#
#--- create a map to the link position
#
    out = out + '<map name="hrci_samp_map">\n'
    for k in range(0, 21):
        out = out + '<area alt="circle" shape="circle" coords="' 

        out = out + str(hrci_cx + hrci_offset[k][0] * ixstep) + ','
        out = out + str(hrci_cy - hrci_offset[k][1] * iystep) + ',20 '

        hlink = html_top + 'Plots3/' + 'samp_i_' + str(k) + '.html'
        out = out + '" href="' + hlink + '">\n'
    out = out + '</map>\n'
    out = out + '<img src="./Plots3/hrci_offset_surface.png" alt="hrci yzoffset" usemap="#hrci_samp_map">'
#
#--- hrc i y/z offset map: pi 
#
    out = out + '</td><td>\n'

    out = out + '<map name="hrci_pi_map">\n'
    for k in range(0, 21):
        out = out + '<area alt="circle" shape="circle" coords="' 

        out = out + str(hrci_cx + hrci_offset[k][0] * ixstep) + ','
        out = out + str(hrci_cy - hrci_offset[k][1] * iystep) + ',20 '

        hlink = html_top + 'Plots3/' + 'pi_i_' + str(k) + '.html'
        out = out + '" href="' + hlink + '">\n'
    out = out + '</map>\n'
    out = out + '<img src="./Plots3/hrci_offset_surface.png" alt="hrci yzoffset" usemap="#hrci_pi_map">'

    out = out + '</td></tr>\n'
    out = out + '</table>\n'

    hrci_text = out
#
#--- hrc s
#
    out = '<p style="text-align:left;">'
    out = out + 'The following plot shows the locations of all observations in <b>sky coordinates.</b> '
    out = out + ' The markers and the colors indicate the location of the observation in Y/Z offset sarface. '
    out = out + ' If you hover the mouse over the marker, it will display the distribution ofsamp data.</p>'

    out = out + '<p style="text-align:left;;">'
    out = out + 'If you like to see a detail, click the magnifier icon at the bottom left, and select ' 
    out = out + 'the area you want to see. Then click the cross icon before hover the mouse. If you '
    out = out + 'like to go back to the entier view, click the house icon.</p> \n'

#
#--- hrc s sky coordindate map
#
    out = out + '<a name="hrcs_sky"></a>\n'
    out = out + '<div style="padding-left:80px;">\n'
    text = mpld3.fig_to_html(fig2)
    text = text.replace('None', '')
    out = out + text
    out = out + '</div>\n\n'

    out = out + '<div style="padding-bottom:50px;"></div>\n'

    out = out + '<p style="text-align:left;padding-bottom:30px;">'
    out = out + '<a href="#top">Back to Top</a></p>\n'
#
#--- det coordinate
#
    out = out + '<p style="text-align:left;">'
    out = out + 'The following plot is similar to the above, but shows the locations of '
    out = out + 'all observations in <b>det coordinates</b>.</p> '

#
#--- hrc s det coordindate map
#
    out = out + '<a name="hrcs_det"></a>\n'
    out = out + '<div style="padding-left:80px;">\n'
    text = mpld3.fig_to_html(fig2b)
    text = text.replace('None', '')
    out = out + text
    out = out + '</div>\n\n'

    out = out + '<div style="padding-bottom:50px;"></div>\n'

    out = out + '<p style="text-align:left;padding-bottom:30px;">'
    out = out + '<a href="#top">Back to Top</a></p>\n'



    out = out + '<p style="text-align:left;;">'
    out = out + 'The following plot shows  21 <b>Y/Z Offset coordinates</b> which all observations '
    out = out + 'fall in. If you click the marker on the location, it will open the trend page of the location. '
    out = out + 'The left plot will open the scaled samp trend and the righ plot will open the PI trend page.</p>'
#
#--- hrc s y z offset map: smap
#
    out = out + '<a name="hrcs_off"></a>\n'
    out = out + '<table border=0 cellpadding=2>\n'
    out = out + '<tr><th>Scaled Samp</th><th>PI</th></tr>\n'
    out = out + '<tr><td>\n'

    out = out + '<map name="hrcs_samp_map">\n'
    for k in range(0, 21):
        out = out + '<area alt="circle" shape="circle" coords="' 

        out = out + str(hrcs_cx + hrcs_offset[k][0] * sxstep) + ','
        out = out + str(hrcs_cy - hrcs_offset[k][1] * systep) + ',20 '

        hlink = html_top + 'Plots3/' + 'samp_s_' + str(k) + '.html'
        out = out + '" href="' + hlink + '">\n'
    out = out + '</map>\n'
    out = out + '<img src="./Plots3/hrcs_offset_surface.png" alt="hrcs yzoffset" usemap="#hrcs_samp_map">'
#
#--- hrcs y z offset map: pi
#
    out = out + '</td><td>\n'

    out = out + '<map name="hrcs_pi_map">\n'
    for k in range(0, 21):
        out = out + '<area alt="circle" shape="circle" coords="' 

        out = out + str(hrcs_cx + hrcs_offset[k][0] * sxstep) + ','
        out = out + str(hrcs_cy - hrcs_offset[k][1] * systep) + ',20 '

        hlink = html_top + 'Plots3/' + 'pi_s_' + str(k) + '.html'
        out = out + '" href="' + hlink + '">\n'
    out = out + '</map>\n'
    out = out + '<img src="./Plots3/hrcs_offset_surface.png" alt="hrcs yzoffset" usemap="#hrcs_pi_map">'

    out = out + '</td></tr>\n'
    out = out + '</table>\n'

    hrcs_text = out


    return [hrci_text, hrcs_text]


#---------------------------------------------------------------------------------------------------
#-- read_data_for_plots: read data and classify which y/z offset group that observation belongs to -
#---------------------------------------------------------------------------------------------------

def read_data_for_plots(infile, posx=6, posy=7):
    """
    read data and classify which y/z offset group that observation belongs to
    input:  infile      --- data file name
            posx        --- position of x data; defalut 6 for sky x
            posy        --- position of y data; defalut 7 for sky y
    output: xdata       --- a list of lists of x data; it has 21 sub lists
            ydata       --- a list of lists of y data; it has 21 sub lists
            obisds      --- a list of lists of obsids; it has 21 sub lists
            inst        --- instrument; either i or s
    """
    mc = re.search('_i', infile)
    if mc is not None:
        offset = hrci_offset
        inst = 'i'
    else:
        offset = hrcs_offset
        inst = 's'

    data  = read_data(infile)
    if data == []:
        return False

    obsids = []
    xdata  = []
    ydata  = []
    for k in range(0, 21):
        obsids.append([])
        xdata.append([])
        ydata.append([])

    for ent in data:
        atemp = re.split('\s+', ent)
        xt = float(atemp[posx])
        yt = float(atemp[posy])
#
#--- find y/z offsets
#
        yo = int(round(float(atemp[10])))
        zo = int(round(float(atemp[11])))
#
#--- find which group this data belongs to
#
        pos = 0
        for m in range(0, 21):
            if (yo == offset[m][0]) and (zo == offset[m][1]):
                pos = m
                break

        obsids[pos].append(atemp[1])
        xdata[pos].append(xt)
        ydata[pos].append(yt)

    return [xdata, ydata, obsids, inst]

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

def read_data(infile, remove=0):

    try:
        f    = open(infile, 'r')
        data = [line.strip() for line in f.readlines()]
        f.close()

        if remove == 1:
            mcf.rm_file(infile)
    except:
        data = []

    return data

#---------------------------------------------------------------------------------------------------
#-- create_interactive_plot: create interactive plots                                             --
#---------------------------------------------------------------------------------------------------

def create_interactive_plot(xdata, ydata, obsids,  inst, xname, yname):

    """
    create interactive plots
    input:  xdata   --- a list of x data
            ydata   --- a list of y data
            obsids  --- a list of obsid
            inst    --- instrument; either i or s
    Output: fig     --- plotting code; this will be modified into html data
    """
#
#--- set css for the plot
#
    css = """
            body{
                width:600px;
                height:300px;
            }
            p{
                text-align:center;
            }
    """
    fsize = 12 
    lsize = 2
    resolution = 100.0
    connect = 0
#
#--- close everything opened before
#
    plt.close('all')
#
#--- set font size
#
    mpl.rcParams['font.size'] = fsize
    props = font_manager.FontProperties(size=9)
#
#---- the first panel
#
#--- set plotting range
#
    test  = xdata[0]
    test2 = ydata[0]
    for m in range(1, 21):
        test  = test  + xdata[m]
        test2 = test2 + ydata[m]
    [xmin, xmax] = set_range(test)
    [ymin, ymax] = set_range(test2)

    fig, ax  = plt.subplots()
    ax.set_autoscale_on(False)
    ax.set_xbound(xmin,xmax)
    ax.set_xlim(xmin=xmin, xmax=xmax, auto=False)
    ax.set_ylim(ymin=ymin, ymax=ymax, auto=False)
    ax.set_axis_bgcolor('#F5F5DC')
    ax.yaxis.set_label_coords(-0.11, 0.5)
    ax.xaxis.set_label_coords(0.5, -0.07)

    if inst == 'i':
        offset = hrci_offset
    else:
        offset = hrcs_offset
#
#--- plot data
#
    xset = xdata[0]
    yset = ydata[0]
    oset = obsids[0]
    for k in range(1,21):
        xset = xset + xdata[k]
        yset = yset + ydata[k]
        oset = oset + obsids[k]

    for k in range(0,21):
        pv = ax.plot(xdata[k], ydata[k], color=point_list[k][0], marker=point_list[k][1], markersize=7, lw = connect, label=offset[k])
        [info_list, html_list] = create_info_link(obsids[k], 'samp')
        plugins.connect(fig, mpld3.plugins.PointHTMLTooltip(pv[0], info_list, css=css, hoffset=-150))

    mpl.rcParams['font.size'] = 7 
    props = font_manager.FontProperties(size=9)
#
#--- add legend
#
    legend(loc=1, bbox_to_anchor=(1.05, 1))
#
#--- label axes
#
    plt.xlabel(xname, size=fsize)
    plt.ylabel(yname, size=fsize)
#
#--- set the size of the plotting area in inch
#
    fig = matplotlib.pyplot.gcf()
    if inst == 'i':
        fig.set_size_inches(10.0, 10.0)
    else:
        fig.set_size_inches(10.0, 5.0)

    plt.close('all')

    return fig


#---------------------------------------------------------------------------------------------------
#-- set_y_range: set min and max range for the plotting                                           --
#---------------------------------------------------------------------------------------------------

def set_range(data):
    """
    set min and max range for the plotting
    input:  data    --- data list
    output: dmin    --- min of the plotting range
            dmax    --- max of the plotting range
    """
    dmin = min(data)
    dmax = max(data)
    diff = dmax - dmin

    dmax = dmax + 0.3 * diff
    dmin = dmin - 0.1 * diff 
    if diff > 1:
        dmax = int(10 * dmax) + 5
        dmin = int(10 * dmin) - 5
        dmax = 0.1 * dmax
        dmin = 0.1 * dmin

    return [dmin, dmax]

#---------------------------------------------------------------------------------------------------
#-- create_info_link: create a list of links used by an interactive plot                          --
#---------------------------------------------------------------------------------------------------

def create_info_link(obsid_list, col):
    """
    create a list of links used by an interactive plot
    input:  obsid_list  --- a list of obsids
    col --- pi or samp
    output: outlist --- a list of html address to a distribution plot
    """
    
    html_plot = html_top + '/Plots3/Indivisual_Plots/'
    outlist  = []
    hlist= []
    for obsid in obsid_list:
        ofile = html_plot + str(obsid)   + '_' + col +  '_list_vfit.png'
        olink = '<p><img src="' + ofile + '" width= 600px></p>'
    
        hlist.append(ofile)
        outlist.append(olink)
    
    return [outlist, hlist]


#-----------------------------------------------------------------------------------------------

if __name__ == "__main__":

    if len(sys.argv) == 2:
        chk = sys.argv[1]
    else:
        chk = ''
    create_html_page(chk)

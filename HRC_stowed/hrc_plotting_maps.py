#!/usr/bin/env /proj/sot/ska/bin/python

#############################################################################################
#                                                                                           #
#       hrc_plotting_maps.py: create hrc maps and histgram                                  #
#                                                                                           #
#               author: t. isobe (tisobe@cfa.harvard.edu)                                   #
#                                                                                           #
#               last update: Jan 08, 2018                                                   #
#                                                                                           #
#############################################################################################

import os
import sys
import re
import string
import time
import math
import numpy
import astropy.io.fits  as pyfits

import matplotlib as mpl

if __name__ == '__main__':
    mpl.use('Agg')

from pylab import *
import matplotlib.pyplot       as plt
import matplotlib.font_manager as font_manager
import matplotlib.lines        as lines

#
#--- from ska
#
from Ska.Shell import getenv, bash
ascdsenv = getenv('source /home/ascds/.ascrc -r release; source /home/mta/bin/reset_param', shell='tcsh')


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
import hrc_stowed_common_function as hcf
#
#--- temp writing file name
#
import random
rtail   = int(time.time())
zspace  = '/tmp/zspace' + str(rtail)

#-----------------------------------------------------------------------------------------------
#-- plot_hrc_map: create maps and histograms of the hrc_i and hrc_s image data                --
#-----------------------------------------------------------------------------------------------

def plot_hrc_map(year, chk=0):
    """
    create maps and histograms of the hrc_i and hrc_s image data
    input:  year    --- the year of the data to be used
            chk     --- indicator of different plotting range defalut 0 (don't use the sepcial range)
    output: <head_name>_<year>_norm.png
            <head_name>_lev2_<year>_norm.png
            hist_<part>_<year>_norm.png
            hist_<part>_lev2_<year>_norm.png
    """
#
#--- Hrc_i_115  
#
    hdir = 'Hrc_i_115'
    xmin    = 1.0e-7
    #xmax    = 8.0e-7
    xmax    = 9.0e-7
    ymin    = 0
    ymax    = 300
    if chk == 1:
        #xmin    = 2.5e-7
        #xmax    = 4.5e-7
        xmin    = 3.0e-7
        xmax    = 7.0e-7
        ymax    = 200
    binsize = 0.2e-8
    hrc     = 'i'
    try:
        run_map_scripts(hdir, year,  xmin, xmax, ymin, ymax, binsize,  hrc=hrc, chk=chk)
    except:
        pass
#
#--- Hrc_s_125_hi
#
    hdir = 'Hrc_s_125_hi'
    xmin    = 0.2e-6
    xmax    = 3.0e-6 
    ymin    = 0
    ymax    = 300

    binsize = 0.1e-7
    hrc     = 's'
    xp      = [0.5e-6, 1.0e-6, 1.5e-6, 2.0e-6, 2.5e-6, 3.0e-6]
    xl      = ['0.5e-6', '1.0e-6', '1.5e-6', '2.0e-6', '2.5e-6', '3.0e-6']
    if chk == 1:
        xmin    = 0.3e-6
        xmax    = 1.5e-6 
        ymax    = 250
        xp      = [0.5e-6, 1.0e-6, 1.5e-6]
        xl      = ['0.5e-6', '1.0e-6', '1.5e-6']
    try:
        run_map_scripts(hdir, year,  xmin, xmax, ymin, ymax, binsize, hrc=hrc, xp=xp, xl=xl, chk=chk)
    except:
        pass
#
#--- Hrc_s_125
#
    hdir = 'Hrc_s_125'
    xmin    = 1.0e-6
    xmax    = 8.0e-6
    ymin    = 0
    ymax    = 300
    binsize = 0.5e-7
    hrc     = 'p'
    xp      = [2.0e-6, 4.0e-6, 6.0e-6]
    xl      = ['2.0', '4.0', '6.0']
    if chk == 1:
        #xmin    = 2.0e-6
        #xmax    = 4.0e-6
        #ymax    =270
        #xp      = [2.0e-6, 3.0e-6, 4.0e-6]
        #xl      = ['2.0', '3.0', '4.0']
        xmin    = 3.0e-6
        xmax    = 5.0e-6
        ymax    =270
        xp      = [3.0e-6, 4.0e-6, 5.0e-6]
        xl      = ['3.0', '4.0', '5.0']
#
#--- hrc_s_125 has three sections
#
    try:
        for part in range(1, 4):
            run_map_scripts(hdir, year,  xmin, xmax, ymin, ymax, binsize, hrc=hrc,xp=xp, xl=xl, part=part,chk=chk)
    except:
        pass

#-----------------------------------------------------------------------------------------------
#-- run_map_scripts: setup the file names and plots histograms and image maps                ---
#-----------------------------------------------------------------------------------------------

def run_map_scripts(hdir, year,  xmin, xmax, ymin, ymax, binsize, hrc='i',xp='', xl='',  part='', chk=0):
    """
    setup the file names and plots histograms and image maps
    input:  hdir    --- a directory name in which the data is kept. e.g. Hrc_i_115
            year    --- a year of the data
            xmin    --- min of x plotting range
            xmax    --- max of x plotting range
            ymin    --- min of y plotting range
            ymax    --- max of y plotting range
            binsize --- size of histogram bin 
            hrc     --- indicator of which instrument is assigned. 'i' or 's'
            xp      --- x axis tick positions default = ''
            xl      --- x axis tick values    default = ''
            part    --- section of hrc_s_125  default = '' for hrc_i and hrc_s_125_hi
            chk     --- indicator of making a special hist plot for the front page default=0 (not making)
    output: <head>_<year>_normal.png
            <head>_<year>_normal_b.png
            <head>_<year>_lev2_normal.png
            <head>_<year>_lev2_normal_b.png
            hist_<year>_normal.png
            hist_<year>_normal_b.png
            hist_<year>_lev2_normal.png
            hist_<year>_lev2_normal_b.png
            for hrc_s_125, section 1, 2, 3 are also created
    """
    head = hdir.lower()
#
#--- input file names
#
    file1   = data_dir + hdir + '/' + head + '_' + str(year) + '_norm' + str(part) + '.fits.gz'
    file2   = data_dir + hdir + '/' + head + '_lev2_' + str(year) + '_norm' + str(part) + '.fits.gz'
#
#--- output file names
#
    img1    = map_dir  + hdir + '/' + head + '_' + str(year) + '_norm' + str(part) + '.png'
    img1b   = map_dir  + hdir + '/' + head + '_' + str(year) + '_norm' + str(part) + '_b.png'
    img2    = map_dir  + hdir + '/' + head + '_lev2_' + str(year) + '_norm' + str(part) + '.png'
    img2b   = map_dir  + hdir + '/' + head + '_lev2_' + str(year) + '_norm' + str(part) + '_b.png'

    thumb1  = map_dir  + hdir + '/Simage/' + head + '_' + str(year) + '_thumb' + str(part) + '.png'
    thumb1b = map_dir  + hdir + '/Simage/' + head + '_' + str(year) + '_thumb' + str(part) + '_b.png'
    thumb2  = map_dir  + hdir + '/Simage/' + head + '_lev2_' + str(year) + '_thumb' + str(part) + '.png'
    thumb2b = map_dir  + hdir + '/Simage/' + head + '_lev2_' + str(year) + '_thumb' + str(part) + '_b.png'

    hist1   = img1.replace(head, 'hist')
    hist2   = img2.replace(head, 'hist')
    if chk == 1:
        hist1   = hist1.replace(str(year),'total_f') 
        hist2   = hist2.replace(str(year),'total_f') 
#
#--- lev1 images
#
    try:
        plot_image_ascds(file1, img1, xmin, xmax, hrc)
        [lpos, tpos] = plot_hist(file1, xmin, xmax, ymin, ymax, binsize, hist1, hrc, xp, xl, part)
        plot_image_ascds(file1, img1b, lpos, tpos, hrc)
#
#--- if there are no data, just copy "no_data" png file to the position
#
    except:
        copy_no_data_image(img1, img1b, img2, img2b, thumb1, thumb1b, thumb2, thumb2b, hist1, hist2)
        return  False
#
#--- lev2 images
#
    plot_image_ascds(file2, img2, xmin, xmax, hrc)
    [lpos, tpos] = plot_hist(file2, xmin, xmax, ymin, ymax, binsize, hist2, hrc, xp, xl, part)
    plot_image_ascds(file2, img2b, lpos, tpos, hrc)
#
#--- thumb nail images
#
    for cset in [[img1, thumb1], [img1b, thumb1b],[img2, thumb2],[img2b, thumb2b]]:

        if hrc == 'i':                                                      #--- hrc_i_115
            cmd = 'convert ' + cset[0] + ' -resize 100x100 '  + cset[1]
        else:
            if part != '':                                                  #--- hrc_s_125
                cmd = 'convert ' + cset[0] + ' -resize 100x150 '  + cset[1]
            else:                                                           #--- hrc_s_125_hi
                cmd = 'convert ' + cset[0] + ' -resize 150x150 '  + cset[1]

        os.system(cmd)

#-----------------------------------------------------------------------------------------------
#-- copy_no_data_image: copy no_data png file to the images                                   --
#-----------------------------------------------------------------------------------------------

def copy_no_data_image(img1, img1b, img2, img2b, thumb1, thumb1b, thumb2, thumb2b, hist1, hist2):
    """
    copy no_data png file to the images
    input: img1, img1b, img2, img2b, thumb1, thumb1b, thumb2, thumb2b, hist1, hist2 --- image file names
    output: all above files have no_data png image
    """

    no_data = house_keeping + 'Templates/no_data.png'

    for ent in [img1, img1b, img2, img2b,  hist1, hist2]:
        cmd = 'cp ' + no_data + ' ' + ent
        os.system(cmd)

    no_data = house_keeping + 'Templates/no_data_small.png'

    for ent in  [thumb1, thumb1b, thumb2, thumb2b]:
        cmd = 'cp ' + no_data + ' ' + ent
        os.system(cmd)


#-----------------------------------------------------------------------------------------------
#-- plot_image_ascds: create an hrc exposure map                                              --
#-----------------------------------------------------------------------------------------------

def plot_image_ascds(img_fits, outfile, xmin, xmax, hrc='i'):
    """
    create an hrc exposure map
    input:  img_fits    --- image fits file name
            outfile     --- output png file name
            xmin        --- a lower limit of the data
            xmax        --- an upper limit of the data
    output: outfile png file
    """
#
#--- set a few parameters
#
    scale = 'linear'
    color = 'sls'
    size  = '125x125'
#
#--- adding a max value reference point so that the image scale correctly
#
    if hrc == 'i':
        mpx = 63
        mpy = 63
    else:
        mpx = 128
        mpy = 31

    cmd1 = '/usr/bin/env PERL5LIB=""'
#
#--- trim data data around the main peak
#
    cmd2 = ' dmimgthresh infile=' + str(img_fits) + ' outfile=trimed_out.fits '
    cmd2 = cmd2 + ' cut='+ str(xmin) + ':' + str(xmax) 
    cmd2 = cmd2 + ' value=' + str(xmin) + ' clobber=yes'

    cmd  = cmd1 + cmd2 
    bash(cmd, env=ascdsenv)

    add_max_point_at_edge('trimed_out.fits', mpx, mpy, xmax)
#
#--- convert the fits image to an eps image
#
    #cmd2 = ' dmimg2jpg trimed_out.fits  ' 
    cmd2 = ' dmimg2jpg zmax_img.fits  ' 
    cmd2 = cmd2 + " greenfile='' bluefile='' regionfile='' outfile='foo.jpg'"
    cmd2 = cmd2 + " scalefunction='" + str(scale) + "' " 
    cmd2 = cmd2 + " psfile='foo.ps' lut=')lut." + str(color) + "' "
    cmd2 = cmd2 + "  showgrid=no   clobber='yes' showaimpoint='no' "

    cmd  = cmd1 + cmd2

    bash(cmd, env=ascdsenv)
#
#--- convert to png
#
    cmd2 = "echo ''|gs -sDEVICE=ppmraw  -r" + str(size) + "   -q -NOPAUSE "
    cmd2 = cmd2 + " -sOutputFile=-  ./foo.ps| pnmtopng > " + str(outfile)

    cmd  = cmd1 + cmd2
    bash(cmd, env=ascdsenv)
#
#--- trim the edges (see: http://www.bsu.edu/libraries/wiki/index.php?title=ImageMagick)
#
    cmd = 'mogrify -trim ' + str(outfile)           #---- removing the white space
    os.system(cmd)
#
#--- farther trim "black" part
#
    if hrc == 'i':
        cmd = ' mogrify -crop 780x780+100+325 '  + str(outfile)
    elif hrc == 'p':
        cmd = ' mogrify -crop 108x1220+470+130 ' + str(outfile)
    else:
        cmd = ' mogrify -crop 208x1000+420+192 ' + str(outfile)

    os.system(cmd)
#
#--- clean up
#
    os.system("rm trimed_out.fits zmax_img.fits  foo.jpg foo.ps")

#-----------------------------------------------------------------------------------------------
#-- add_max_point_at_edge: add a reference point to an image fits file                        --
#-----------------------------------------------------------------------------------------------

def add_max_point_at_edge(fits, x, y, val, out='zmax_img.fits'):
    """
    add a reference point to an image fits file
    input:  fits    --- image fits file name
            x       --- x position (image coordinates)
            y       --- y position (image coordinates)
            val     --- reference value
            out     --- output image fits file name default: zmax_iimg.fitx
    output: out     --- image fits file with the reference point
    """
    hd = pyfits.open(fits)
    data = hd[0].data
    data[x, y] = val
    mcf.rm_file(out)
    hd.writeto(out)

#-----------------------------------------------------------------------------------------------
#-- plot_hist: create histogram plot                                                          --
#-----------------------------------------------------------------------------------------------

def plot_hist(img_fits, xmin, xmax, ymin, ymax, binsize, outname, hrc='i',xp ='', xl = '',part=''):
    """
    create histogram plot
    input:  img_fits    --- image fits file
            xmin        --- min of plotting range
            xmax        --- max of plotting range
            ymin        --- min of y 
            ymax        --- max of y
            binsize     --- bin size
            outname     --- name of output file name
            hrc         --- indicator of instrument, "i" or "s"
            xp          --- a list of x tick positions 
            xl          --- a list of x tick values
            part        --- section of hrc (for hrc_s_125), 1, 2, or 3
    output: png plot in outname
    """
#
#--- read fits data
#
    image_data = pyfits.getdata(img_fits)
#
#--- create a bin list and convert to an array
#
    nbins      = int((xmax - xmin) / binsize) + 1
    bin_list   = []
    for k in range (0, nbins):
        bin_list.append(xmin + binsize * k)
    bin_array = numpy.array(bin_list)
#
#--- clean up the plotting device
#
    plt.close('all')
#
#---- set a few parameters
#
    mpl.rcParams['font.size']   = 12
    mpl.rcParams['font.weight'] = 'bold'
#
#--- set plotting range
#
    plt.axis([xmin, xmax, ymin, ymax])
#
#--- create histogram plot
#
    plt.hist(image_data.flat, bin_array, color='green' )
#
#--- for the case  x axis are munaully set.
#
    if xp != '':
        plt.xticks(xp, xl)
#
#--- hrc_s_125 has bit more setting to do
#
    if hrc =='p':
        if part == 1:
            plt.ylabel('Number of Occurance',  fontweight='bold')
            pass
        elif part == 2:
            plt.xlabel('Cnts/Sec/Pix (x10^-6)', fontweight='bold')

        create_plot_file(outname, xsize=3.2)
    else:
        plt.xlabel('Counts per Sec per Pixel', fontweight='bold')
        plt.ylabel('Number of Occurance',      fontweight='bold')

        create_plot_file(outname)
#
#--- trim white space
#
    if hrc == 'p':
        if part == 1:
            cmd = 'mogrify -crop 900x2400+0+180 '   + str(outname)
        else:
            cmd = 'mogrify -crop 790x2400+115+180 ' + str(outname)
    else:
        cmd = 'mogrify -crop 2600x2450+140+160 '    + str(outname)

    os.system(cmd)
#
#--- now find the range of the main peak
#
    return find_main_hist_range(image_data, bin_array)

#-----------------------------------------------------------------------------------------------
#-- find_main_hist_range: find the range of the main peak of the histogram                    --
#-----------------------------------------------------------------------------------------------

def find_main_hist_range(image_data, bin_array):
    """
    find the range of the main peak of the histogram
    input:  image_data  --- pyfits data of image data
            bin_array   --- array of bins
    output: bin_list[lpos]  --- the lower boundary of the main peak
            bin_list[tpos]  --- the upper boundary of the main peak
    """
    out      = numpy.histogram(image_data, bin_array)
    clist    = list(out[0])
    bin_list = list(bin_array)
#
#--- first find the max position
#
    pmax  = 0
    vmax  = 0
    top   = len(clist)
    for k in range(0, top):
        if clist[k] > vmax:
            vmax = clist[k]
            pmax = k
        else:
            continue
#
#--- find the 30 %  height point in left
#
    hlim = 0.30 * vmax
    for i in range(0, pmax):
        if clist[i] >= hlim:
            lstart = i
            break
        else:
            lstart = pmax
#
#--- if we find zero count bin between max position and 30% postion, the data are 
#--- too small to make a peak area selection; so skip it
#
    chk = 0
    for k in range (lstart, pmax):
        if clist[k] == 0:
            chk = 1
            break
#
#--- find the first zero position to left 
#
    if chk == 1:
        lpos = 0
        tpos = top
    else:   
        for i in range(lstart, -1, -1):
            try:
                if clist[i] == 0:
                    lpos = i
                    break
                else:
                    lpos = 0
            except:
                lpos = 0
#
#--- find the first zero position to right
#
        for i in range(pmax, top):
            try:
                if clist[i] == 0:
                    tpos = i
                    break
                else:
                    tpos = top
            except:
                tpos = top
#
#--- return the range 
#
    return [bin_list[lpos], bin_list[tpos]]

#-----------------------------------------------------------------------------------------------
#-- create_plot_file: save the plot in a file                                                 --
#-----------------------------------------------------------------------------------------------

def create_plot_file(outname, xsize=10.0, ysize=9.0, resolution=300, fmt='png'):
    """
    save the plot in a file
    input:  outname     --- output file name
    xsize       --- x axis size in inch     default: 10 in
    ysize       --- y axis size in inch     defalut:  8 in
    resolution  --- resolution of the plot  defalut: 300 ppi
    fmg         --- format of the plot      defalut: png
    output: outname     --- plot in <fmt>
    """
#
#--- set the size of the plotting area in inch
#
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(xsize, ysize)
#
#--- save the plot in png format
#
    plt.savefig(outname, format=fmt, dpi=resolution)


#-----------------------------------------------------------------------------------------------

if __name__ == "__main__":

    if len(sys.argv) == 2:
        if sys.argv[1] == 'total':
            plot_hrc_map('total')
        else:
            year = int(float(sys.argv[1]))
            plot_hrc_map(year)
    elif len(sys.argv) == 3:
            plot_hrc_map('total', chk=1)
    else:
        print "Please supply year you want to create the map"

#!/usr/bin/perl

##############################################################
#
# Checks HRC EVT Header for RANGELEV and WIDTHRES keywords, 
# and sets them to the proper values if the don't exist
#
# initialize CIAO before running
#
# Takes filename as argument 
#
# JPB, 8 Feb 2007
#
#   t. isobe Mar 09, 2017
#
###############################################################

$outdir = '/data/hrc/i/';

$obsfile=$ARGV[0];  #for other sources

print "$obsfile\n";

$obsfile = "$outdir"."$obsfile";

open(OBS, $obsfile);
while (<OBS>){
    @tmp=split;
    $obsid=trim($tmp[0]);
   
    print "$obsid\t";

    $blah=`ls $outdir$obsid/secondary/hrcf*N*evt1.fits*`;
    $file=trim($blah);

    print "$file\n";

    $rangelev=`dmlist ${file} head | grep RANGELEV`;
    print "rangelev = $rangelev\n";
    if (trim($rangelev) =~ /^$/) {
	$date=`dmlist ${file} head | grep DATE-OBS | cut -c27-36`;
	@pieces=split /-/, $date;
	$yr=trim($pieces[0]);
	$mo=trim($pieces[1]);
	$day=trim($pieces[2]);
	print "yr=${yr}\tmon=${mo}\tday=${day}\n";
	if ($yr > 1999 || ($yr == 1999 && $mo == 12 && $day > 6)) {
	    print `dmhedit infile=${file} filelist=none operation=add key=RANGELEV value=115\n`;
	    print "setting rangelev to 115\n";
	} else {
	    print `dmhedit infile=${file} filelist=none operation=add key=RANGELEV value=90\n`;
	    print "setting rangelev to 90\n";
	}
    }

    $widthres=`dmlist ${file} head | grep WIDTHRES`;
    print "widthres=$widthres\n";
    if (trim($widthres) =~ /^$/) {
	$date=`dmlist ${file} head | grep DATE-OBS | cut -c27-36`;
	@pieces=split /-/, $date;
	$yr=trim($pieces[0]);
	$mo=trim($pieces[1]);
	$day=trim($pieces[2]);
	print "yr=${yr}\tmon=${mo}\tday=${day}\n";
	if ($yr > 2000 || ($yr == 2000 && $mo > 10) || ($yr == 2000 && $mo == 10 && $day > 5)) {
	    print `dmhedit infile=${file} filelist=none operation=add key=WIDTHRES value=2\n`;
	    print "setting wedthres to 2\n";
	} else {
	    print `dmhedit infile=${file} filelist=none operation=add key=WIDTHRES value=3\n`;
	    print "setting widthres to 3\n";
	} 
    }

}

close(OBS);
    
sub trim {
    my @out = @_;
    for (@out){
        s/^\s+//;
        s/\s+$//;
    }
    return wantarray ? @out: $out[0];
}

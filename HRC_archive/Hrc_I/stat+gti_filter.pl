#!/usr/bin/perl

##############################################################
#
# Filter new evt1 on new GTI
#
# initialize CIAO before running
#
# Takes aoN as argument, e.g.
# gti_filter.pl ao7 
#
# JPB, 9 Feb 2007
#
#   t. isobe Sep 27, 2017
#
###############################################################

$outdir = '/data/hrc/i/';

$obsfile=$ARGV[0];

print "$obsfile\n";

$obsfile = "$outdir"."$obsfile";

open(OBS, $obsfile);
while (<OBS>){
    @tmp=split;
    $obsid=trim($tmp[0]);
   
    print "$obsid\n";
 
    $gti=`ls $outdir${obsid}/secondary/*std_flt1*`;
    $gti=trim($gti);
    
    $evt1="$outdir${obsid}/secondary/hrcf${obsid}_evt1_new.fits";
    $stflt="$outdir${obsid}/secondary/hrcf${obsid}_evt1_new_stflt.fits";
    $outfile="$outdir${obsid}/analysis/hrcf${obsid}_evt2.fits";
    
    print `rm -rf param; mkdir param`;
    print `source /home/mta/bin/reset_param`;

    print `dmcopy \"${evt1}[status=xxxxxx00xxxx0xxx00000000x0000000]\" ${stflt} verbose=3 clobber=yes`;
    print `dmcopy \"${stflt}[EVENTS][\@${gti}]\" ${outfile} verbose=3 clobber=yes`;
}

sub trim {
    my @out = @_;
    for (@out){
        s/^\s+//;
        s/\s+$//;
    }
    return wantarray ? @out: $out[0];
}


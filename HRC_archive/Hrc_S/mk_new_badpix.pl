#!/usr/bin/perl

##############################################################
#
# Make a new observation-specific bad pixel file to using the 
# latest CALDB degap.
#
# initialize CIAO before running
#
# JPB, 23 June 2007
#
#   t. isobe Sep 27, 2017
#
###############################################################

$outdir = '/data/hrc/s/';

$obsfile=$ARGV[0];

print "$obsfile\n";

open(OBS, $obsfile);
while (<OBS>){
    @tmp=split;
    $obsid=trim($tmp[0]);

    $blah=`ls $outdir$obsid/secondary/*evt1.fits*`;
    $file=trim($blah);

    $dir="$outdir${obsid}/secondary/";

    print "\n\n$obsid\t$file";

    print `rm -rf param; mkdir param`;
    print `source /home/mta/bin/reset_param`;

    print `dmmakepar ${file} ${dir}hrcf${obsid}_obs.par clobber=yes mode=h verbose=3`;
   
    print `punlearn hrc_build_badpix`;
    print `pset hrc_build_badpix infile='CALDB(DETECTOR=HRC-S)'`;
    print `pset hrc_build_badpix degapfile='CALDB(DETECTOR=HRC-S)'`;
    print `pset hrc_build_badpix outfile=${dir}hrcf${obsid}_new_bpix1.fits`;
    print `pset hrc_build_badpix obsfile=${dir}hrcf${obsid}_obs.par`;
    print `pset hrc_build_badpix logfile=${dir}hrc_build_badpix.log`;
    print `pset hrc_build_badpix verbose=3`;
    print `pset hrc_build_badpix clobber=yes`;
    print `hrc_build_badpix mode=h`;
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

#!/usr/bin/perl

##############################################################
#
# Uses hrc_dtfstats to recompute deadtime statistics, 
# updates relevant keywords in evt1 header
#
# initialize CIAO before running
#
# JPB, 12 Feb 2007
#
#   t.isobe Mar 09, 2017
#
###############################################################

$outdir = '/data/hrc/s/';

$obsfile=$ARGV[0];

print "$obsfile\n";

open(OBS, $obsfile);
while (<OBS>){

    @tmp=split;
    $obsid=trim($tmp[0]);
   
    print "$obsid\n";

    $blah = `ls $outdir${obsid}/primary/*dtf1*`;
    $dtf = trim($blah);

    print "$dtf\n\n";

    $olddtcor = `dmkeypar $outdir${obsid}/secondary/hrcf${obsid}_evt1_new.fits DTCOR echo+`;
    $olddtor = trim($olddtcor);

    $gti=`ls $outdir${obsid}/secondary/*std_flt1.fits*`;
    $outfile="$outdir${obsid}/analysis/new_dtfstat.fits";

    print `punlearn hrc_dtfstats`;
    print `pset hrc_dtfstats infile=${dtf}`;
    print `pset hrc_dtfstats outfile=${outfile}`;
    print `pset hrc_dtfstats gtifile=${gti}`;
    print `pset hrc_dtfstats verbose=3`;
    print `pset hrc_dtfstats clobber=yes`;
    print `pset hrc_dtfstats mode=h`;
    print `hrc_dtfstats`;


    $foo=`dmlist \"${outfile}[cols DTCOR]\" data,clean`;
    my ($dtcor) = ($foo =~ /\#\s+DTCOR\s+([0-9\.\-]+)/s);
    print "$olddtcor\t$dtcor\t";

    $evt1="$outdir${obsid}/analysis/hrcf${obsid}_evt2.fits";

    print `punlearn dmkeypar`;
    $tmp= `dmkeypar $evt1 ONTIME echo+`;
    $ontime=trim($tmp);
    print "$ontime\t";
    
    $livetime=$ontime * $dtcor;
    $exposure=$ontime * $dtcor;
    
    print "$livetime\n";
    
    print `punlearn dmhedit`;
    print `dmhedit ${evt1} filelist=\"\" op=add key=LIVETIME value=${livetime}`;
    print `dmhedit ${evt1} filelist=\"\" op=add key=EXPOSURE value=${exposure}`;
    print `dmhedit ${evt1} filelist=\"\" op=add key=DTCOR value=${dtcor}`;
    
}

sub trim {
    my @out = @_;
    for (@out){
        s/^\s+//;
        s/\s+$//;
    }
    return wantarray ? @out: $out[0];
}

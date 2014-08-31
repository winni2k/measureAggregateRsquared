#!/usr/bin/env perl

use strict;
use warnings;
use File::Path qw(make_path);
use FindBin qw($Bin);

use Test::More tests=>1;
use Test::Files;

my $sampDir = "$Bin/../samples/v2.0";
my $tempDir = "$Bin/results/10-chr20";
if(!-d $tempDir){
    make_path($tempDir);
}
my $truth = "$Bin/../samples/05/MARv0.1.EUR.snps";
system("$Bin/../bin/measureAggregateRsquared.dbg --validation $sampDir/ALL.20130523.snps_indels.CGI.genotypes.nonGL.inRef.inGL.head1000.chr20.gen.gz --imputed $sampDir/ALL.20130523.snps_indels.CGI.genotypes.nonGL.inRef.inGL.inChip.impute2.notInChip.head500.chr20.gen.gz --sample $sampDir/CGI2.ILLU1M.sample --freq $sampDir/ALL.chr20.unionAC10NM.founders.B7.R1.P8.M30.K100.W500kb.inInt.inTruth.GL.head1000.freqs --bin $sampDir/olivierBins.txt --output $tempDir/test");
compare_ok($truth, "$tempDir/test.EUR.snps", "chr20 results OK");

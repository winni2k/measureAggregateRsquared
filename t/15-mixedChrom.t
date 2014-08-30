#!/usr/bin/env perl

use strict;
use warnings;
use File::Path qw(make_path);
use FindBin qw($Bin);

use Test::More tests=>1;
use Test::Files;

my $sampDir = "$Bin/../samples/v2.0";
my $tempDir = "$Bin/results/15-mixedChrom";
if(!-d $tempDir){
    make_path($tempDir);
}
my $truth = "$Bin/../samples/v1.0/test.EUR.snps";
my $cmd = "$Bin/../bin/measureAggregateRsquared.dbg --validation $sampDir/mixedChrom/ALL.20130523.snps_indels.CGI.genotypes.nonGL.inRef.inGL.head1000.randChr.gen.gz --imputed $sampDir/mixedChrom/ALL.20130523.snps_indels.CGI.genotypes.nonGL.inRef.inGL.inChip.impute2.notInChip.head500.randChr.gen.gz --sample $sampDir/CGI2.ILLU1M.sample --freq $sampDir/ALL.chr20.unionAC10NM.founders.B7.R1.P8.M30.K100.W500kb.inInt.inTruth.GL.head1000.freqs --bin $sampDir/olivierBins.txt --output $tempDir/test.mixedChrom";
print STDOUT "Running: $cmd\n";
system($cmd);
compare_ok($truth, "$tempDir/test.mixedChrom.EUR.snps", "mixedChrom results OK");

#!/usr/bin/env perl

use strict;
use warnings;
use File::Path qw(make_path);
use FindBin qw($Bin);

use Test::More tests => 4;
use Test::Files;

my $sampDir = "$Bin/../samples/05";
my $tempDir = "$Bin/results/05-chr20Simple";
if ( !-d $tempDir ) {
    make_path($tempDir);
}
my $truth = "$Bin/../samples/05/MARv0.1.EUR.snps";
my $cmd =
"$Bin/../bin/measureAggregateRsquared.dbg --validation $sampDir/ALL.20130523.snps_indels.CGI.genotypes.nonGL.inRef.inGL.head1000.499.chr20.gen.gz --imputed $sampDir/ALL.20130523.snps_indels.CGI.genotypes.nonGL.inRef.inGL.inChip.impute2.notInChip.head500.499.chr20.gen.gz --sample $sampDir/CGI2.ILLU1M.sample --freq $sampDir/ALL.chr20.unionAC10NM.founders.B7.R1.P8.M30.K100.W500kb.inInt.inTruth.GL.head1000.499.freqs --bin $sampDir/olivierBins.txt --output $tempDir/test\n";
print STDOUT "Running: $cmd";
system($cmd);

my $truthALL = "$Bin/../samples/05/MARv0.7.EUR.all";
compare_ok( "$tempDir/test.EUR.snps", $truth, "simple chr20 results OK" );
compare_ok( "$tempDir/test.EUR.snps.map", "$truthALL.map",
    "simple chr20 map OK" );

compare_ok( "$tempDir/test.EUR.all", "$truth", "simple chr20 ALL OK" );
compare_ok( "$tempDir/test.EUR.all.map", "$truthALL.map",
    "simple chr20 ALL map OK" );

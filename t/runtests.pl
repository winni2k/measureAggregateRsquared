#!/usr/bin/env perl

use strict;
use warnings;

use TAP::Harness;
my $harness = TAP::Harness->new();
my $tests   = qx/find . -name "*.t"/;
my @tests   = split( /\s+/, $tests );
$harness->runtests(@tests);


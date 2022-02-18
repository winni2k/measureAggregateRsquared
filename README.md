measureAggregateRsquared
========================

measureAggregateRsquared measures aggregate r squared from imputed gen files

Synopsis
--------

```bash
measureAggregateRsquared --validation truth.gen.gz --imputed imputed.gen.gz \
--sample truth_and_imputed.samples --freq allele_frequencies_of_imputed_sites.freq \
--bin allele_frequency_bins.txt --output output_base
```

Make sure the truth and imputed gen files contain the same samples in the same order, which is defined in the .samples file. 

Authorship
----------
This code was written by Olivier Delaneau and Warren Kretzschmar.  The maintainer of this code is Warren Kretzschmar.

Bug reports and questions
------------------------
Please raise an issue on the [github page](https://github.com/winni2k/measureAggregateRsquared/issues).

Input file formatting
---------------------

The `--validation` and `--imputed` input files are [Impute2 .gen](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#gensample) files.  
Below are a set of examples for the other input files.

### .samples file

comparison is by population, multiple populations allowed...

    ID_1 ID_2 missing pop
    0 0 0 D
    NA07346 NA07346 0 EUR
    NA11832 NA11832 0 EUR

### frequencies file

Population is first line. After that, each line corresponds to the allele frequency in 
that population in the truth.gen file.  Multiple columns, one for each population
allowed.

    EUR
    0.2214
    0.02241
    0.3206

### bins file

Each line defines a boundary of bins.  
Whether or not the first boundary is included can be changed using the "--discard-monomorphic" flag (I think).

    0.000
    0.005
    0.010

### output files

Each output file will be written to:
* `<outfix>.<POP>.snps`
* `<outfix>.<POP>.complexs`
* `<outfix>.<POP>.all`

, where the suffix includes the type of variant being assessed.

Each output file has the following header attached:
"Bin_frequency r_square num_genotypes freq_validation freq_imputation"


Testing
-------

The test suite requires an installed version of CPAN. To install the perl dependencies for the test runners:

    make test-setup

To run the tests:

    make test

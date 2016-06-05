measureAggregateRsquared
========================

measureAggregateRsquared measures aggregate r squared from imputed gen files

Here is a sample command line:

```bash
measureAggregateRsquared --validation truth.gen.gz --imputed imputed.gen.gz \
--samples truth_and_imputed.samples --freq allele_frequencies_of_imputed_sites.freq \
--bin allele_frequency_bins.txt --output output_base
```

Make sure the truth and imputed gen files contain the same samples in the same order, which is defined in the .samples file. 

## sample .samples file
comparison is by population, multiple populations allowed...

    ID_1 ID_2 missing pop
    0 0 0 D
    NA07346 NA07346 0 EUR
    NA11832 NA11832 0 EUR

## sample frequencies file
Population is first line. After that, each line corresponds to the allele frequency in 
that population in the truth.gen file.  Multiple columns, one for each population
allowed.

    EUR
    0.2214
    0.02241
    0.3206

## sample bins txt
Each line defines a boundary of bins.  
Whether or not the first boundary is included can be changed using the "--discard-monomorphic" flag (I think).

    0.000
    0.005
    0.010

output files are written to the output base


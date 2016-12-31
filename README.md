# Abundance Model Accounting for Detection and Availability from Bird Point Counts

## General description

Point counts at ~79 sites with most visited on two occasions in 2015 and 2 additional occasions in 2016 (4 closed periods for an open population model). At each visit, three time bins of 3 min, 2 min, and 5 min were used with only new birds recorded in each bin. Birds were counted in three distance bins of approximately 50, 100, and 150 m. 

## Use

The `Abundance_Kery_Royle.R` is a stand alone script based on a script from Marc Kery and Andy Royle's Hierarchical Modeling Book (Volume 1, 9.3), but it does not work with this data currently. Error: `Node inconsistent with parents`. Potentially the starting values don't work with the data-model combination. GitHub issue: [https://github.com/djhocking/Bird_Communities/issues/2](https://github.com/djhocking/Bird_Communities/issues/2)

The updated way to run the scripts in this repo will be:

1. Data prep 

2. JAGS analysis

3. Output organization

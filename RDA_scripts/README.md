# RDA

These scripts are used to detect SNPs under selection by using the R package 'vegan' to perform an RDA. This workflow is based on the procedures detailed in this [vignette](https://popgen.nescent.org/2018-03-27_RDA_GEA.html#data-packages)

* [run\_RDA\_sims.sh](#run_rda_simssh)
* [proc\_sims\_RDA.R](#proc_sims_rdar)
* [combine\_RDA\_tables.R](#combine_rda_tablesr)
---


## run\_RDA\_sims.sh

This script is used to run the other scripts on a set of many simulations.

### Usage

`$ run_RDA_sims.sh`

### Inputs

This script takes no command line arguments, instead the variables at the top of the script should be adjusted depending on the situation. These variables are as follows:

* myDir  : by default, the directory the script is run in. This is where the program looks for the VCF files to analyze
* begin  : first sim to analyze
* end    : last sim to analyze
* simType: type of sim, the script uses this to find the names of the files it processes 

These scripts require two compressed VCFs for each simulation, one pruned to contain quasi-independent alleles only and one with data for the entire set of alleles. They must be named as follows:

* [number]\_[simType]\_VCFallFILT.vcf.gz
* [number]\_PRUNED\_[simType]\_VCFallFILT.vcf.gz

For example, to analyze simulation 10900 of type "Invers" you would need two files named "10900\_Invers\_VCFallFILT.vcf.gz" and "10900\_PRUNED\_Invers\_VCFallFILT.vcf.gz"

### Outputs

Outputs are contained in three directories that the script will create under 'myDir'

* log\_files : captures standard out for each run of proc\_sims\_RDA.R
* plots     : each run of proc\_sims\_RDA.R creates a plot that represents the effects of the first two RDAs
* results   : tables of RDA loading values. One table for each sim, contains both pruned and unpruned results
---


## proc\_sims\_RDA.R

This Rscript processes a single VCF, which it identifies based on the inputs it is given by run\_RDA\_sims.sh. It outputs a table of loadings into a directory called "raw\_data", which must then be processed by combine\_RDA\_tables.R. It can be run parallel with other instances of proc\_sims\_RDA.R, but all instances must finish before the raw data tables can be combined. run\_RDA\_sims.sh runs two instances of this script at a time by default

---


## combine\_RDA\_tables.R

This Rscript takes a starting sim and an ending sim from run\_RDA\_sims.sh and processes all sims from the start to the end. It finds the pruned and unpruned tables for each sim, combines the data, and renames the columns.    
 

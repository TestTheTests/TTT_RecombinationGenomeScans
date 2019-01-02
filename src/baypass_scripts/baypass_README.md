# Run Baypass

The baypass pipeline is run using run_baypass.sh. This script does not take any variables from the command line. If everything is set up correctly, you should be able to just change a few variables at the top of the script and run it.

The first few lines of the script:
```
mypath="/shared_lab/20181217_TTT_RecombinationGenomeScans/src/baypass_scripts
ncpus=60
nsims=210
start=10900
finish=$(($start + $nsims-1))
```
'ncpus' is a parameter passed directly to baypass, it can be changed depending on the resources available. 'mypath' is just the absoulute path to the folder containing the baypass scripts. 'nsims' and 'start' are used to generate the sequence of simulations to run the scripts on. If there are missing sims (ex: sims go 10900, 10901, 10903) then nsims should not correspond to the actual number of sims, instead to the difference between the first and last sim.

After checking these variables to make sure they look right, you can just run the script:

`$ run_baypass.sh`

It shouldn't matter where you run the script from as long as your 'mypath' variable is correct and the directory structure of the scripts and the 'results_final' is as expected.


## Directory Structure

    |------results_final
    |       	|-------xxxx_Invers_indFILT.txt
    |			|------xxxx_Invers_ScanResults.txt
    |			|------xxxx_Invers_VCFallFILT.vcf.gz
    |
    |-----src	
	|---------------run_baypass     
 			|---converted_files     (created by script)
   			|	|------xxxx_Invers_VCFallFILT.covar
 			|	|------xxxx_Invers_VCFallFILT.geno
 			|	|------xxxx_Invers_VCFallFILT_PRUNED.covar
			|	|------xxxx_Invers_VCFallFILT_PRUNED.geno
 			|-----baypasss_results     (created by script)
 			|	. (many files, see baypass 2.1 documnetation for more info)
 			|	.
 			|	. 
 			|-----log_files            (created by script)
 			|	|------xxxx_baypass_err.txt
 			|	|------xxxx_baypass_log.txt					
 	     		|-----run_baypass.sh
 			|----vcf2baypass.pl
 			|----pruneSNPs.pl
 			|----getTableOfBaypassResults.pl					


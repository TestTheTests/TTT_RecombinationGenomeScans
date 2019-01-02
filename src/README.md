# TestTheTests

### To run the script for simulation and all analysis in serial (one after another)
```
conda activate msprime-env
cd TTT_RecombinationGenomeScans
./src/repInversion_serial.sh > output.txt
```
### To run script for all the steps after the neutral mutations have been overlayed on the tree sequence:
```
conda activate msprime-env
cd TTT_RecombinationGenomeScans
./src/repInversion_serial_onlyR.sh > output.txt
```

katie@Katies-MacBook-Air:~$ conda activate msprime-env
(msprime-env) katie@Katies-MacBook-Air:~$ conda list
```
# packages in environment at /anaconda3/envs/msprime-env:
#
# Name                    Version                   Build  Channel
attrs                     18.2.0                    <pip>
attrs                     18.2.0                     py_0    conda-forge
blas                      1.0                         mkl
bzip2                     1.0.6                         1    conda-forge
ca-certificates           2018.11.29           ha4d7672_0    conda-forge
certifi                   2018.11.29            py36_1000    conda-forge
cycler                    0.10.0                    <pip>
gsl                       2.2.1                h002c638_3
h5py                      2.8.0            py36h097b052_4    conda-forge
hdf5                      1.10.3               hc401514_2    conda-forge
intel-openmp              2019.0                      118
jsonschema                3.0.0a3               py36_1000    conda-forge
kiwisolver                1.0.1                     <pip>
libffi                    3.2.1                hfc679d8_5    conda-forge
libgfortran               3.0.1                h93005f0_2
libopenblas               0.3.3                hdc02c5d_3
matplotlib                3.0.2                     <pip>
mkl                       2018.0.3                      1
mkl_fft                   1.0.6                    py36_0    conda-forge
mkl_random                1.0.2                    py36_0    conda-forge
msprime                   0.6.1            py36hcb787e7_0    conda-forge
ncurses                   6.1                  hfc679d8_1    conda-forge
numpy                     1.15.4           py36h6a91979_0
numpy-base                1.15.4           py36h8a80b8c_0
openblas                  0.3.3                ha44fe06_1    conda-forge
openssl                   1.0.2p               h470a237_1    conda-forge
pandas                    0.23.4           py36hf8a1672_0    conda-forge
pip                       18.1                  py36_1000    conda-forge
pyparsing                 2.3.0                      py_0    conda-forge
pyrsistent                0.14.6           py36h470a237_0    conda-forge
pyslim                    0.1                       <pip>
python                    3.6.6                h5001a0f_0    conda-forge
python-dateutil           2.7.5                      py_0    conda-forge
python.app                1.2                    py36_200    conda-forge
pytz                      2018.7                     py_0    conda-forge
readline                  7.0                  haf1bffa_1    conda-forge
scikit-allel              1.1.10                    <pip>
setuptools                40.6.2                   py36_0    conda-forge
six                       1.11.0                py36_1001    conda-forge
sqlite                    3.25.3               hb1c47c0_0    conda-forge
svgwrite                  1.2.1                     <pip>
svgwrite                  1.2.1                      py_0    conda-forge
tk                        8.6.9                ha92aebf_0    conda-forge
wheel                     0.32.2                   py36_0    conda-forge
xz                        5.2.4                h470a237_1    conda-forge
zlib                      1.2.11               h470a237_3    conda-forge
```

Session Info for R
```
sessionInfo()
R version 3.4.0 (2017-04-21)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS  10.13.4

Matrix products: default
BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] tools     grid      stats     graphics  grDevices utils     datasets  methods  
[9] base     

other attached packages:
 [1] rehh_2.0.2          rehh.data_1.0.0     related_1.0         bigsnpr_0.2.1      
 [5] bigstatsr_0.2.3     LEA_1.8.1           lfmm_0.0            pcadapt_3.1.0      
 [9] adegenet_2.1.1      ade4_1.7-10         OutFLANK_0.2        qvalue_2.8.0       
[13] data.table_1.10.4-3 bindrcpp_0.2        RColorBrewer_1.1-2  vcfR_1.6.0         
[17] ROCR_1.0-7          gplots_3.0.1        dbplyr_1.2.0        RSQLite_2.0        
[21] forcats_0.2.0       stringr_1.2.0       dplyr_0.7.4         purrr_0.2.4        
[25] readr_1.1.1         tidyr_0.7.2         tibble_1.3.4        ggplot2_2.2.1      
[29] tidyverse_1.2.1    

loaded via a namespace (and not attached):
 [1] colorspace_1.3-2    RcppEigen_0.3.3.3.1 seqinr_3.4-5        deldir_0.1-14      
 [5] rprojroot_1.3-2     rstudioapi_0.7      bit64_0.9-7         RSpectra_0.12-2    
 [9] lubridate_1.7.1     xml2_1.1.1          codetools_0.2-15    splines_3.4.0      
[13] mnormt_1.5-5        memuse_4.0-0        knitr_1.18          RcppRoll_0.2.2     
[17] jsonlite_1.5        broom_0.4.3         cluster_2.0.6       shiny_1.0.5        
[21] compiler_3.4.0      httr_1.3.1          backports_1.1.2     assertthat_0.2.0   
[25] Matrix_1.2-12       lazyeval_0.2.1      cli_1.0.0           htmltools_0.3.6    
[29] igraph_1.1.2        coda_0.19-1         gtable_0.2.0        glue_1.2.0         
[33] reshape2_1.4.3      gmodels_2.16.2      Rcpp_0.12.15        cellranger_1.1.0   
[37] spdep_0.7-4         gdata_2.18.0        ape_5.0             nlme_3.1-131       
[41] iterators_1.0.9     pinfsc50_1.1.0      psych_1.7.8         rvest_0.3.2        
[45] mime_0.5            gtools_3.5.0        LearnBayes_2.15     MASS_7.3-48        
[49] scales_0.5.0        hms_0.4.0           parallel_3.4.0      expm_0.999-2       
[53] yaml_2.1.16         memoise_1.1.0       stringi_1.1.6       foreach_1.4.4      
[57] permute_0.9-4       caTools_1.17.1      boot_1.3-20         spData_0.2.6.9     
[61] rlang_0.1.6         pkgconfig_2.0.1     bitops_1.0-6        evaluate_0.10.1    
[65] lattice_0.20-35     bindr_0.1           htmlwidgets_1.0     labeling_0.3       
[69] bit_1.1-12          plyr_1.8.4          magrittr_1.5        R6_2.2.2           
[73] DBI_0.7             haven_1.1.0         foreign_0.8-69      mgcv_1.8-22        
[77] sp_1.2-5            modelr_0.1.1        crayon_1.3.4        xgboost_0.6-4      
[81] KernSmooth_2.23-15  plotly_4.7.1        rmarkdown_1.8       readxl_1.0.0       
[85] blob_1.1.0          vegan_2.4-6         digest_0.6.15       xtable_1.8-2       
[89] httpuv_1.3.5        munsell_0.4.3       viridisLite_0.2.0  

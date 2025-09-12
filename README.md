This repository includes the following folders and files corresponding to the manuscripts \
[1] "Resampling-based confidence intervals and bands for the average treatment effect in observational studies with competing risks" \
and [2] "Using propensity score matching for inference about the average treatment effect in competing-risks data" \
by Rühl, J. and Friedrich, S.:

- ./Results/ \
  A folder containing (interim) results of the simulations (Section 4) in manuscript [1].
  
  - ATE_true.Rda \
    An Rda file containing results for the true average treatment effect considered in the simulations. \
    This file is produced by l. 22-88 in the R script 'simu_masterscript.R'. 

  - res_[effect]ATE_[scenario]_n[n].Rda \
    (effect: adv/no/-, \
     scenario: noCens/lowCens/highCens/lowTreatProb/highTreatProb/lowVarCov/highVarCov/typeII, \
     n: 50/75/100/200/300) \
     Rda files containing results that summarize the outcomes of the simulations for each scenario. \
     These files are produced by l. 100-150 in the R script 'simu_masterscript.R'.  

  - total_coverages.Rda \
    An Rda file containing results for the coverages of the simulated confidence intervals and bands. \
    This file is produced by l. 158-244 in the R script 'simu_masterscript.R'.

- ./Results_PSM/ \
  A folder containing (interim) results of the simulations (Section X) in manuscript [2].
  
  - res_[effect]ATE_[scenario]_n[n].Rda \
    (effect: adv/no/-, \
     scenario: noCens/lowCens/highCens/lowTreatProb/highTreatProb/lowVarCov/highVarCov/typeII, \
     n: 50/75/100/200/300) \
     Rda files containing results that summarize the outcomes of the simulations for each scenario. \
     These files are produced by l. 74-124 in the R script 'simu_masterscript_PSM.R'.  

  - total_coverages.Rda \
    An Rda file containing results for the coverages of the simulated confidence intervals and bands. \
    This file is produced by l. 129-219 in the R script 'simu_masterscript_PSM.R'.

- ./ATESurvival_1.0.tar.gz \
  An R package for the derivation of confidence intervals and bands for the average treatment effect for survival data using the classical bootstrap, an 
  influence function approach, the wild bootstrap, or a double-resampling approach for PS-matched data. \
  To install the package, use the following command: \
  `install.packages("ATESurvival_1.0.tar.gz", repos = NULL, source = TRUE)` \
  Note that on Windows systems, Rtools is required (see https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html).

- ./hd_analysis.R \
  An R script that performs the analysis of the hd data reported in manuscript [1], reproducing Table 3 and Figure 8.

- ./hd_analysis_PSM.R \
  An R script that performs the analysis of the hd data reported in manuscript [2], reproducing Table X and Figure X.

- ./simu_functions.R \
  An R script that defines the functions used for the simulations reported in manuscript [1]. \
  The R script 'simu_masterscript.R' is based on these functions.

- ./simu_functions_PSM.R \
  An R script that defines the functions used for the simulations reported in manuscript [2]. \
  The R script 'simu_masterscript_PSM.R' is based on these functions.

- ./simu_masterscript.R \
  An R script that performs the simulations reported in manuscript [1], reproducing Figures 1 - 7. \
  The simulations were run in parallel on a Linux server with 16 cores. \
  Replication on a Windows system occasionally yielded slightly different confidence intervals/bands, but the differences should be neglegible. \
  Interim results are saved in the folder ./Results as the execution of the complete simulation study takes several days. To check reproducibility, one might reduce the number of iterations by choosing a smaller number for the parameter 'iter' of the function 'run' (l. 146 in the script 'simu_masterscript.R').

- ./simu_masterscript_PSM.R \
  An R script that performs the simulations reported in manuscript [2], reproducing Figures X - X. \
  The simulations were run in parallel on a Linux server with 16 cores. \
  Replication on a Windows system occasionally yielded slightly different confidence intervals/bands, but the differences should be neglegible. \
  Interim results are saved in the folder ./Results_PSM as the execution of the complete simulation study takes several days. To check reproducibility, one might reduce the number of iterations by choosing a smaller number for the parameter 'iter' of the function 'run' (l. 120 in the script 'simu_masterscript_PSM.R').

---

To reproduce the simulation results presented in manuscript [1] (Section 4), install the package 'ATESurvival' and run the script 'simu_masterscript.R'. \
Interim results are stored in the 'Results' folder. \
To reproduce the analysis of the real data application (Section 5), install the package 'ATESurvival' and run the script 'hd_analysis.R'.

To reproduce the simulation results presented in manuscript [2] (Section X), install the package 'ATESurvival' and run the script 'simu_masterscript_PSM.R'. \
Interim results are stored in the 'Results_PSM' folder. \
To reproduce the analysis of the real data application (Section X), install the package 'ATESurvival' and run the script 'hd_analysis_PSM.R'.

The code was created and evaluated in R using the following software:
```
R version 4.4.1 (2024-06-14 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26100)

Matrix products: default

locale:
[1] LC_COLLATE=German_Germany.utf8  LC_CTYPE=German_Germany.utf8   
[3] LC_MONETARY=German_Germany.utf8 LC_NUMERIC=C                   
[5] LC_TIME=German_Germany.utf8    

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] patchwork_1.3.0           ggplot2_3.5.1            
 [3] ATESurvival_1.0           doRNG_1.8.6              
 [5] rngtools_1.5.2            doParallel_1.0.17        
 [7] iterators_1.0.14          foreach_1.5.2            
 [9] survival_3.7-0            prodlim_2024.06.25       
[11] riskRegression_2023.12.21

loaded via a namespace (and not attached):
 [1] gtable_0.3.6        xfun_0.49           htmlwidgets_1.6.4  
 [4] lattice_0.22-6      numDeriv_2016.8-1.1 vctrs_0.6.5        
 [7] tools_4.4.1         generics_0.1.3      sandwich_3.1-1     
[10] tibble_3.2.1        fansi_1.0.6         cluster_2.1.6      
[13] pkgconfig_2.0.3     Matrix_1.7-0        data.table_1.16.2  
[16] checkmate_2.3.2     lifecycle_1.0.4     farver_2.1.2       
[19] compiler_4.4.1      stringr_1.5.1       MatrixModels_0.5-3 
[22] munsell_0.5.1       codetools_0.2-20    SparseM_1.84-2     
[25] quantreg_5.99       htmltools_0.5.8.1   htmlTable_2.4.3    
[28] Formula_1.2-5       pillar_1.9.0        MASS_7.3-60.2      
[31] cmprsk_2.2-12       rms_6.8-2           Hmisc_5.2-0        
[34] multcomp_1.4-26     rpart_4.1.23        nlme_3.1-166       
[37] parallelly_1.39.0   lava_1.8.0          timereg_2.0.6      
[40] tidyselect_1.2.1    digest_0.6.37       polspline_1.1.25   
[43] mvtnorm_1.3-2       stringi_1.8.4       future_1.34.0      
[46] dplyr_1.1.4         listenv_0.9.1       splines_4.4.1      
[49] fastmap_1.2.0       grid_4.4.1          colorspace_2.1-1   
[52] cli_3.6.3           magrittr_2.0.3      base64enc_0.1-3    
[55] utf8_1.2.4          TH.data_1.1-2       future.apply_1.11.3
[58] withr_3.0.2         foreign_0.8-87      mets_1.3.4         
[61] scales_1.3.0        backports_1.5.0     rmarkdown_2.29     
[64] globals_0.16.3      nnet_7.3-19         gridExtra_2.3      
[67] zoo_1.8-12          evaluate_1.0.1      knitr_1.49         
[70] rlang_1.1.4         Rcpp_1.0.13         glue_1.7.0         
[73] rstudioapi_0.17.1   R6_2.5.1
```

The simulations were run in parallel on a Linux server with the subsequent software versions:
```
R version 4.4.3 (2025-02-28)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 22.04.5 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
[1] ATESurvival_1.0           doRNG_1.8.2              
[3] rngtools_1.5.2            doParallel_1.0.17        
[5] iterators_1.0.14          foreach_1.5.2            
[7] survival_3.5-5            riskRegression_2023.03.22

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.0     timeDate_3043.102    dplyr_1.1.4         
 [4] fastmap_1.1.0        TH.data_1.1-0        pROC_1.18.0         
 [7] caret_6.0-90         digest_0.6.34        rpart_4.1.24        
[10] lifecycle_1.0.4      conquer_1.2.1        cluster_2.1.8.1       
[13] magrittr_2.0.3       compiler_4.4.3       rlang_1.1.0         
[16] Hmisc_4.6-0          tools_4.4.3          utf8_1.2.4          
[19] data.table_1.16.4    knitr_1.36           timereg_2.0.1       
[22] htmlwidgets_1.5.4    plyr_1.8.6           RColorBrewer_1.1-3  
[25] multcomp_1.4-17      polspline_1.1.19     foreign_0.8-88      
[28] withr_2.5.0          purrr_0.3.4          numDeriv_2016.8-1.1 
[31] stats4_4.4.3         nnet_7.3-20          grid_4.4.3          
[34] fansi_1.0.6          latticeExtra_0.6-29  mets_1.2.9          
[37] colorspace_2.1-0     future_1.23.0        ggplot2_3.4.4       
[40] globals_0.14.0       scales_1.3.0         MASS_7.3-65         
[43] cli_3.6.2            mvtnorm_1.2-4        rms_6.2-0           
[46] generics_0.1.3       rstudioapi_0.15.0    future.apply_1.8.1  
[49] reshape2_1.4.4       stringr_1.4.0        splines_4.4.3       
[52] matrixStats_1.2.0    base64enc_0.1-3      vctrs_0.6.5         
[55] sandwich_3.0-1       Matrix_1.7-3         SparseM_1.81        
[58] Formula_1.2-4        htmlTable_2.3.0      listenv_0.8.0       
[61] jpeg_0.1-9           gower_0.2.2          cmprsk_2.2-10       
[64] recipes_0.1.17       glue_1.6.2           parallelly_1.28.1   
[67] codetools_0.2-19     lubridate_1.8.0      stringi_1.7.12      
[70] gtable_0.3.4         munsell_0.5.0        tibble_3.2.1        
[73] pillar_1.9.0         htmltools_0.5.2      quantreg_5.86       
[76] ipred_0.9-12         lava_1.6.10          R6_2.5.1            
[79] lattice_0.22-5       png_0.1-7            backports_1.4.1     
[82] class_7.3-23         MatrixModels_0.5-0   Rcpp_1.0.14          
[85] gridExtra_2.3        nlme_3.1-167         prodlim_2019.11.13  
[88] checkmate_2.3.1      xfun_0.32            zoo_1.8-9           
[91] pkgconfig_2.0.3      ModelMetrics_1.2.2.2
``` 

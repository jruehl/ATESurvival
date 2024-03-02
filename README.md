This repository includes the following folders and files corresponding to the manuscript \
"Resampling-based confidence intervals and bands for the average treatment effect in observational studies with competing risks" \
by Rühl, J. and Friedrich, S.:

- ./Results/ \
  A folder containing (interim) results of the simulations (Section 4 in the manuscript).
  
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

- ./ATESurvival_1.0.tar.gz \
  An R package for the derivation of confidence intervals and bands for the average treatment effect for survival data using the classical bootstrap, an 
  influence function approach and the wild bootstrap. \
  To install the package, use the following command: \
  `install.packages("ATESurvival_1.0.tar.gz", repos = NULL, source = TRUE)` \
  Note that on Windows systems, Rtools is required (see https://cran.r-project.org/bin/windows/Rtools/rtools43/rtools.html).

- ./hd_analysis.R \
  An R script that performs the analysis of the hd data reported in the manuscript, reproducing Table 3 and Figure 8.

- ./simu_functions.R \
  An R script that defines the functions used for the simulations reported in the manuscript. \
  The R script 'simu_masterscript.R' is based on these functions.

- ./simu_masterscript.R \
  An R script that performs the simulations reported in the manuscript, reproducing Figures 1 - 7. \
  The simulations were run in parallel on a Linux server with 16 cores. \
  Replication on a Windows system occasionally yielded slightly different confidence intervals/bands, but the differences should be neglegible. \
  Interim results are saved in the folder ./Results as the execution of the complete simulation study takes several days. To check reproducibility, one might reduce the number of iterations by choosing a smaller number for the parameter 'iter' of the function 'run' (l. 146 in the script 'simu_masterscript.R').

---

To reproduce the simulation results presented in the manuscript (Section 4), install the package 'ATESurvival' and run the script 'simu_masterscript.R'. \
Interim results are stored in the 'Results' folder. \
To reproduce the analysis of the real data application (Section 5), install the package 'ATESurvival' and run the script 'hd_analysis.R'.

The code was created and evaluated in R using the following software: \
R version 4.1.2 (2021-11-01) \
Platform: x86_64-w64-mingw32/x64 (64-bit) \
Running under: Windows 10 x64 (build 19045)

Matrix products: default

locale: \
[1] LC_COLLATE=German_Germany.1252  LC_CTYPE=German_Germany.1252 \
[3] LC_MONETARY=German_Germany.1252 LC_NUMERIC=C \
[5] LC_TIME=German_Germany.1252    

attached base packages: \
[1] parallel  stats     graphics  grDevices utils     datasets  methods \
[8] base

other attached packages: \
 [1] patchwork_1.1.2           ggplot2_3.4.1            
 [3] ATESurvival_1.0           doRNG_1.8.2              
 [5] rngtools_1.5.2            doParallel_1.0.17        
 [7] iterators_1.0.14          foreach_1.5.2            
 [9] survival_3.5-5            prodlim_2019.11.13       
[11] riskRegression_2023.03.22

loaded via a namespace (and not attached): \
 [1] splines_4.1.2       Formula_1.2-4       latticeExtra_0.6-30
 [4] globals_0.16.2      timereg_2.0.4       numDeriv_2016.8-1.1
 [7] pillar_1.9.0        backports_1.4.1     lattice_0.20-45    
[10] quantreg_5.94       glue_1.6.2          digest_0.6.31      
[13] RColorBrewer_1.1-3  checkmate_2.1.0     colorspace_2.0-3   
[16] sandwich_3.0-2      cmprsk_2.2-11       rms_6.3-0          
[19] htmltools_0.5.5     Matrix_1.5-3        pkgconfig_2.0.3    
[22] SparseM_1.81        listenv_0.8.0       mvtnorm_1.2-3      
[25] scales_1.2.1        jpeg_0.1-10         lava_1.7.0         
[28] MatrixModels_0.5-1  htmlTable_2.4.1     tibble_3.2.1       
[31] mets_1.3.1          generics_0.1.3      withr_2.5.0        
[34] TH.data_1.1-2       nnet_7.3-18         cli_3.6.1          
[37] magrittr_2.0.3      deldir_1.0-6        polspline_1.1.22   
[40] future_1.33.0       fansi_1.0.3         parallelly_1.36.0  
[43] nlme_3.1-160        MASS_7.3-58.1       foreign_0.8-83     
[46] tools_4.1.2         data.table_1.14.6   lifecycle_1.0.3    
[49] multcomp_1.4-20     stringr_1.5.0       interp_1.1-3       
[52] munsell_0.5.0       cluster_2.1.4       compiler_4.1.2     
[55] rlang_1.1.1         grid_4.1.2          rstudioapi_0.14    
[58] htmlwidgets_1.5.4   base64enc_0.1-3     gtable_0.3.1       
[61] codetools_0.2-18    R6_2.5.1            gridExtra_2.3      
[64] zoo_1.8-11          knitr_1.44          dplyr_1.1.3        
[67] fastmap_1.1.1       future.apply_1.10.0 utf8_1.2.2         
[70] Hmisc_4.7-2         stringi_1.7.12      Rcpp_1.0.11        
[73] vctrs_0.6.3         rpart_4.1.19        png_0.1-8          
[76] tidyselect_1.2.0    xfun_0.40


The simulations were run in parallel on a Linux server with the subsequent software versions: \
R version 4.3.2 (2023-10-31) \
Platform: x86_64-pc-linux-gnu (64-bit) \
Running under: Ubuntu 22.04.3 LTS \

Matrix products: default \
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 \
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0

locale: \
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages: \
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages: \
[1] ATESurvival_1.0           doRNG_1.8.2              
[3] rngtools_1.5.2            doParallel_1.0.17        
[5] iterators_1.0.14          foreach_1.5.2            
[7] survival_3.5-5            riskRegression_2023.03.22

loaded via a namespace (and not attached): \
 [1] tidyselect_1.1.1     timeDate_3043.102    dplyr_1.0.8         
 [4] fastmap_1.1.0        TH.data_1.1-0        pROC_1.18.0         
 [7] caret_6.0-90         digest_0.6.33        rpart_4.1.23        
[10] lifecycle_1.0.4      conquer_1.2.1        cluster_2.1.5       
[13] magrittr_2.0.3       compiler_4.3.2       rlang_1.1.0         
[16] Hmisc_4.6-0          tools_4.3.2          utf8_1.2.4          
[19] data.table_1.14.8    knitr_1.36           timereg_2.0.1       
[22] htmlwidgets_1.5.4    plyr_1.8.6           RColorBrewer_1.1-3  
[25] multcomp_1.4-17      polspline_1.1.19     foreign_0.8-86      
[28] withr_2.5.0          purrr_0.3.4          numDeriv_2016.8-1.1 
[31] stats4_4.3.2         nnet_7.3-19          grid_4.3.2          
[34] fansi_1.0.6          latticeExtra_0.6-29  mets_1.2.9          
[37] colorspace_2.1-0     future_1.23.0        ggplot2_3.4.4       
[40] globals_0.14.0       scales_1.3.0         MASS_7.3-60         
[43] cli_3.6.1            mvtnorm_1.2-4        rms_6.2-0           
[46] generics_0.1.1       rstudioapi_0.15.0    future.apply_1.8.1  
[49] reshape2_1.4.4       stringr_1.4.0        splines_4.3.2       
[52] matrixStats_1.1.0    base64enc_0.1-3      vctrs_0.5.2         
[55] sandwich_3.0-1       Matrix_1.6-4         SparseM_1.81        
[58] Formula_1.2-4        htmlTable_2.3.0      listenv_0.8.0       
[61] jpeg_0.1-9           gower_0.2.2          cmprsk_2.2-10       
[64] recipes_0.1.17       glue_1.6.2           parallelly_1.28.1   
[67] codetools_0.2-19     lubridate_1.8.0      stringi_1.7.12      
[70] gtable_0.3.4         munsell_0.5.0        tibble_3.2.1        
[73] pillar_1.9.0         htmltools_0.5.2      quantreg_5.86       
[76] ipred_0.9-12         lava_1.6.10          R6_2.5.1            
[79] lattice_0.22-5       png_0.1-7            backports_1.4.1     
[82] class_7.3-22         MatrixModels_0.5-0   Rcpp_1.0.8          
[85] gridExtra_2.3        nlme_3.1-163         prodlim_2019.11.13  
[88] checkmate_2.3.1      xfun_0.32            zoo_1.8-9           
[91] pkgconfig_2.0.3      ModelMetrics_1.2.2.2

The package 'ATESurvival' in this repository includes the function 'ATE' that computes the average treatment effect together with confidence intervals and bands for time-to-event data with competing risks based on the influence function and/or the wild bootstrap.
To install the package, use the following command: install.packages("ATESurvival_1.0.tar.gz", repos = NULL, source = TRUE)  
Note that on Windows systems, Rtools is required (see https://cran.r-project.org/bin/windows/Rtools/rtools43/rtools.html).

The R scripts 'simu_masterscript' (using 'simu_functions') and 'hd_analysis' reproduce the simulation results and the data analysis in the manuscript 'Resampling-based confidence intervals and bands for the average treatment effect in observational studies with competing risks' (Rühl, Friedrich), respectively.

---

R version 4.1.2 (2021-11-01)

Platform: x86_64-w64-mingw32/x64 (64-bit)

Running under: Windows 10 x64 (build 19045)

Packages:
- Rcpp (v1.0.11)
- RcppArmadillo (v0.12.4.1.0)
- RcppClock (v1.1)
- riskRegression (v2023.03.22)
- prodlim (v2019.11.13)
- survival (v3.5-5)
- doParallel (v1.0.17)
- rngtools (v1.5.2)
- doRNG (v1.8.2)

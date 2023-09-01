The package 'ATESurvival' in this repository includes the function 'ATE_IF_WBS' that computes the average treatment effect and corresponding confidence intervals/bands for time-to-event data with competing risks, using influence functions and/or the wild bootstrap. Besides, the function 'EBS' yields the ATE for data that have been drawn with replacement (i.e. in the context of bootstrapping).
The R scripts 'simulations' and 'hd_analysis' reproduce the simulation results and the data analysis in the manuscript 'Resampling-based confidence intervals and bands for the average treatment effect in observational studies with competing risks' (Rühl, Friedrich), respectively.

---

R version 4.1.2 (2021-11-01)

Platform: x86_64-w64-mingw32/x64 (64-bit)

Running under: Windows 10 x64 (build 19045)

Packages:
- Rcpp (v1.0.10)
- RcppArmadillo (v0.11.4.2.1)
- RcppClock (v1.1)
- riskRegression (v2023.03.22)
- prodlim (v2019.11.13)
- survival (v3.5-5)
- doParallel (v1.0.17)
- rngtools (v1.5.2)
- doRNG (v1.8.2)

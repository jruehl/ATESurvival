# This script includes the functions that produce the results of the simulation study.

# generate data ################################################################

# generate data
generate_data <- function(n, 
                          alpha_0 = 0, beta_0, gamma_0 = 0, delta_0 = 0, 
                          sd_cov = 1, 
                          true_ATE = FALSE, 
                          cens = !true_ATE, cens_par = 200, m = NULL){
  
  # input:
  ## n: sample size
  ## alpha_0: parameter for treatment probability
  ## beta_0: parameter for treatment effect on type I events
  ## gamma_0: parameter for treatment effect on type II events (if m is NULL)
  ## delta_0: parameter for treatment effect on censoring (if m is NULL)
  ## sd_cov: standard deviation of the covariates
  ## true_ATE: logical value indicating whether the true risk difference is to be
  ##   determined based on the generated data
  ##   (i.e. implement random treatment allocation & no censoring)
  ## cens: logical value indicating whether the data should be censored
  ##   (if FALSE, the true risk difference is determined for random treatment
  ##    allocation & no censoring)
  ## cens_par: parameter for censoring time distribution (if m is NULL)
  ## m: number of events to observe before censoring is imputed in case of
  ##   type II censoring
  
  
  # output: data set with the subsequent columns:
  ## time: (observed) event time
  ## event: event type (either 0,1,2, with 0 indicating censored observations)
  ## A: treatment group (0 or 1)
  ## z1-z12: covariates (z1-z6: normally distributed, z7-z12: binary)
  
  # covariates
  Z <- cbind(matrix(rnorm(6*n, 0, sd_cov), nrow=n), 
             matrix(rbinom(6*n, 1, 0.5), nrow=n))
  colnames(Z) <- c("z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12")
  # treatment
  if(true_ATE){
    p_treat <- 0.5
  }else{
    arg_treat <- alpha_0 + colSums(rep(c(1,-1,0,0,0,1),2)*log(2) * t(Z))
    p_treat <- exp(arg_treat) / (1 + exp(arg_treat))
  }
  A <- rbinom(n, 1, p_treat)
  # type I events
  arg_T1 <- beta_0*A + colSums(rep(c(1,0,1,0,0,1),2)*log(2) * t(Z))
  T1 <- rweibull(n, shape = 2, scale = 1/sqrt(1/100*exp(arg_T1)))
  
  if(is.null(m)){ # competing risks and censoring
    # type II events
    arg_T2 <- gamma_0*A + colSums(rep(c(-1,0,0,0,1,1),2)*log(2) * t(Z))
    T2 <- rweibull(n, shape = 2, scale = 1/sqrt(1/100*exp(arg_T2)))
    # censoring
    if(!cens){
      C <- Inf
    }else{
      arg_C <- delta_0*A + colSums(rep(c(-1,0,0,1,0,-1),2)*log(2) * t(Z))
      C <- rweibull(n, shape = 2, scale = 1/sqrt(1/cens_par*exp(arg_C)))
    }
    # save data
    data <- data.frame(time = pmin(T1, T2, C), 
                       event = ifelse(pmin(T1, T2, C) == T1, 1, 
                                      ifelse(pmin(T1, T2, C) == T2, 2, 0)),
                       A = as.factor(A),
                       Z)
  }else{ # standard survival data with type II censoring & staggered entry
    if(!cens){
      time <- T1
      event <- 1
    }else{
      # staggered entry
      entry <- runif(n, 0, min(qweibull(m/n, 2, 1/sqrt(1/100*exp(arg_T1)))))
      # type II censoring
      c <- sort(entry + T1)[m]
      event <- as.numeric(entry + T1 <= c)
      time <- ifelse(event == 1, T1, c - entry)
    }
    # save data
    data <- data.frame(time = time,
                       event = event,
                       A = as.factor(A),
                       Z)
  }
  
  return(data)
  
}


# run simulations ##############################################################

# run simulations and save summarized results
run <- function(n, 
                alpha_0 = 0, beta_0, gamma_0 = 0, delta_0 = 0, 
                sd_cov = 1, 
                cens = TRUE, cens_par = 200, m = NULL,
                t = NULL, 
                ATE_true, 
                EBS_iter = 1e3, BS_iter = 1e4, 
                iter = 5000, seed = 1234, cores = detectCores()-1){
  # input:
  ## n: sample size
  ## alpha_0: parameter for treatment probability
  ## beta_0: parameter for treatment effect on type I events
  ## gamma_0: parameter for treatment effect on type II events (if m is NULL)
  ## delta_0: parameter for treatment effect on censoring (if m is NULL)
  ## sd_cov: standard deviation of the covariates
  ## cens: logical value indicating whether the data should be censored
  ## cens_par: parameter for censoring time distribution (if m is NULL)
  ## m: number of events to observe before censoring is imputed in case of
  ##   type II censoring
  ## t: time point(s) to evaluate
  ## ATE_true: true ATE
  ## EBS_iter: number of samples to use for the EBS approach
  ## BS_iter: number of random multipliers to use for the IF/WBS approach
  ## iter: number of Monte Carlo repetitions
  ## seed: seed for reproducibility
  ## cores: number of CPU cores for parallel computations
  ##   (use 16 cores to reproduce the results from the manuscript)
  
  # output: list with (summarized) simulation results
  
  
  # set time points if t is not specified
  if(is.null(t) & is.null(m)){ # competing risks
    t <- c(1,3,5,7,9)
  }else if(is.null(t)){ # type II censoring
    if(beta_0 == -2){
      t <- c(2,4,6,8,10)
    }else if(beta_0 == 0){
      t <- c(1,2,3,4,5)
    }else if(beta_0 == 2){
      t <- c(0.5,1,1.5,2,2.5)
    }
  }
  
  # create matrices to store coverages
  cov_CI_EBS <- 
    cov_CI_IF <- 
    cov_CI_WBS_calc <- cov_CI_WBS_Lin <- cov_CI_WBS_Bey <- cov_CI_WBS_Weird <- 
    matrix(NA, nrow=iter, ncol=length(t))
  cov_CB_EBS <- 
    cov_CB_IF <- 
    cov_CB_WBS_Lin_calc <-cov_CB_WBS_Lin <- 
    cov_CB_WBS_Bey_calc <- cov_CB_WBS_Bey <- 
    cov_CB_WBS_Weird_calc <- cov_CB_WBS_Weird <- 
    rep(NA, iter)
  
  res <- list()
  set.seed(seed)
  
  for(i in 1:iter){
    
    res[[i]] <- NA
    
    try({
      
      # generate data
      if(is.null(m)){ # competing risks
        # ensure that at least 10 events of each type are observed
        while(TRUE){ 
          data <- generate_data(
            n = n, 
            alpha_0 = alpha_0, beta_0 = beta_0, 
            gamma_0 = gamma_0, delta_0 = delta_0, 
            sd_cov = sd_cov, 
            cens = cens, cens_par = cens_par
          )
          if(min(dim(data[data$event == 1,])[1], 
                 dim(data[data$event == 2,])[1]) >= 10) break
        }
        # determine transition frequencies at each time point
        event_prob <- sapply(t, 
                             function(t){
                               c(sum(data$event == 0 & data$time <= t),
                                 sum(data$event == 1 & data$time <= t),
                                 sum(data$event == 2 & data$time <= t))/n
                             })
        row.names(event_prob) <- c("0","1","2")
      }else{ # type II censoring
        data <- generate_data(n = n, beta_0 = beta_0, m = m)
        # determine transition frequencies at each time point
        event_prob <- sapply(t, 
                             function(t){
                               c(sum(data$event == 0 & data$time <= t),
                                 sum(data$event == 1 & data$time <= t))/n
                             })
        row.names(event_prob) <- c("0","1")
      }
      colnames(event_prob) <- paste0("t = ", t)
      # determine treatment frequency
      treat_prob <- sum(data$A == 1)/n
      
      
      # perform EBS
      EBS_time <- proc.time()
      cl <- makeCluster(cores)
      registerDoParallel(cl)
      if(is.null(m)){ # competing risks
        res_EBS <- 
          foreach(j = 1:EBS_iter, .combine = cbind, 
                  .packages = c('riskRegression', 'ATESurvival'), 
                  .options.RNG = seed+i) %dorng% { 
                    tryCatch({
                      data_bs <- data[sample(1:n, n, replace = TRUE),]
                      invisible(capture.output(
                        csc_bs <- CSC(
                          formula = 
                            list(Hist(time, event) ~ 
                                   A+z1+z2+z3+z4+z5+z6+z7+z8+z9+z10+z11+z12,
                                 Hist(time, event) ~ 
                                   A+z1+z2+z3+z4+z5+z6+z7+z8+z9+z10+z11+z12),
                          data_bs[order(data_bs$time),])
                      ))
                      EBS(
                        Z = list(csc_bs$models$`Cause 1`$x,
                                 csc_bs$models$`Cause 2`$x),
                        event = data_bs$event[order(data_bs$time)],
                        time = sort(data_bs$time),
                        t = t,
                        beta = coefficients(csc_bs),
                        index_A = c(
                          c(which(colnames(csc_bs$models$`Cause 1`$x) == 
                                    paste0("A", levels(data$A)[2])), NA)[1],
                          c(which(colnames(csc_bs$models$`Cause 2`$x) == 
                                    paste0("A", levels(data$A)[2])), NA)[1])
                      )
                    }, error = function(e){rep(NA, length(t))})
                  }
      }else{ # type II censoring
        res_EBS <- 
          foreach(j = 1:EBS_iter, .combine = cbind, 
                  .packages = c('survival', 'ATESurvival'), 
                  .options.RNG = seed+i) %dorng% { 
                    tryCatch({
                      data_bs <- data[sample(1:n, n, replace = TRUE),]
                      invisible(capture.output(
                        cox_bs <- coxph(
                          Surv(time, event) ~ 
                            A+z1+z2+z3+z4+z5+z6+z7+z8+z9+z10+z11+z12,
                          data_bs[order(data_bs$time),]
                        )
                      ))
                      EBS(
                        Z = list(model.matrix(cox_bs)),
                        event = data_bs$event[order(data_bs$time)],
                        time = sort(data_bs$time),
                        t = t,
                        beta = list(cox_bs$coefficients),
                        index_A = c(which(colnames(model.matrix(cox_bs)) == 
                                            paste0("A", levels(data$A)[2])), 
                                    NA)[1]
                      )
                    }, error = function(e){rep(NA, length(t))})
                  }
      }
      stopCluster(cl)
      CI_EBS <- apply(res_EBS, 1, function(res_EBS_t){
        res_EBS_t <- res_EBS_t[is.finite(res_EBS_t) & !is.na(res_EBS_t)]
        if(length(res_EBS_t) == 0) return(list(rep(NA, 2), 0))
        return(list(c(max(quantile(res_EBS_t, 0.025), -1), 
                      min(quantile(res_EBS_t, 0.975), 1)), 
                    length(res_EBS_t)))
      })
      q_EBS_sup <- quantile(apply(abs(res_EBS - rowMeans(res_EBS, na.rm = TRUE)) / 
                                    sqrt(apply(res_EBS, 1, var, na.rm = TRUE)), 
                                  2, max, na.rm = TRUE), 0.95)
      EBS_time <- ifelse(any(!is.na(unlist(lapply(CI_EBS, `[[`, 1)))), 
                         (proc.time() - EBS_time)[3]*1e3, 
                         NA)
      
      
      # prepare IF/WBS: fit (cause-specific) Cox models
      Cox_time <- proc.time()
      if(is.null(m)){ # competing risks
        csc <- CSC(formula = list(Hist(time, event) ~ 
                                    A+z1+z2+z3+z4+z5+z6+z7+z8+z9+z10+z11+z12,
                                  Hist(time, event) ~ 
                                    A+z1+z2+z3+z4+z5+z6+z7+z8+z9+z10+z11+z12),
                   data[order(data$time),])
      }else{ # type II censoring
        cox <- coxph(Surv(time, event) ~ 
                       A+z1+z2+z3+z4+z5+z6+z7+z8+z9+z10+z11+z12,
                     data[order(data$time),])
      }
      Cox_time <- (proc.time() - Cox_time)[3]*1e3
      
      res_IF_WBS <- list(rep(NA, length(t)),
                         matrix(NA, nrow=2, ncol=length(t)), 
                         matrix(NA, nrow=2, ncol=length(t)),
                         matrix(NA, nrow=8, ncol=length(t)), 
                         matrix(NA, nrow=12, ncol=length(t)))
      # perform IF/WBS
      IF_time <- NA
      WBS_time <- NA
      try({
        if(is.null(m)){ # competing risks
          res_IF_WBS <- ATE_IF_WBS(
            Z = list(csc$models$`Cause 1`$x,
                     csc$models$`Cause 2`$x),
            event = data$event[order(data$time)], 
            time = sort(data$time),
            t = t,
            beta = coefficients(csc),
            index_A = c(c(which(colnames(csc$models$`Cause 1`$x) == 
                                  paste0("A", levels(data$A)[2])), NA)[1],
                        c(which(colnames(csc$models$`Cause 2`$x) == 
                                  paste0("A", levels(data$A)[2])), NA)[1]),
            bs_iter = BS_iter
          )
        }else{ # type II censoring
          res_IF_WBS <- ATE_IF_WBS(
            Z = list(model.matrix(cox)),
            event = data$event[order(data$time)], 
            time = sort(data$time),
            t = t,
            beta = list(cox$coefficients),
            index_A = c(which(colnames(model.matrix(cox)) == 
                                paste0("A", levels(data$A)[2])), NA)[1],
            bs_iter = BS_iter
          )
        }
        IF_time <- unname(
          Cox_time + sum(clock$timer[clock$ticker %in% c("ATE","IF")]/1e6)
        )
        WBS_time <- unname(
          Cox_time + sum(clock$timer[clock$ticker %in% c("ATE","WBS")]/1e6)
        )
      }, silent = TRUE)
      
      
      # save results
      CI <- rbind(matrix(unlist(lapply(CI_EBS, `[[`, 1)), nrow=2),
                  res_IF_WBS[[2]],
                  res_IF_WBS[[4]])
      row.names(CI) <- c(
        "EBS (lower)", "EBS (upper)",
        "IF (lower)", "IF (upper)",
        "WBS - calculated (lower)", "WBS - calculated (upper)",
        "WBS - Lin (lower)", "WBS - Lin (upper)",
        "WBS - Beyersmann (lower)", "WBS - Beyersmann (upper)",
        "WBS - Weird (lower)", "WBS - Weird (upper)"
      )
      CB <- rbind(pmax(as.numeric(res_IF_WBS[[1]]) - 
                         q_EBS_sup * apply(res_EBS, 1, sd, na.rm = TRUE), -1), 
                  pmin(as.numeric(res_IF_WBS[[1]]) + 
                         q_EBS_sup * apply(res_EBS, 1, sd, na.rm = TRUE), 1),
                  res_IF_WBS[[3]],
                  res_IF_WBS[[5]])
      row.names(CB) <- c(
        "EBS (lower)", "EBS (upper)",
        "IF (lower)", "IF (upper)",
        "WBS - Lin calculated (lower)", "WBS - Lin calculated (upper)",
        "WBS - Lin (lower)", "WBS - Lin (upper)",
        "WBS - Beyersmann calculated (lower)", "WBS - Beyersmann calculated (upper)",
        "WBS - Beyersmann (lower)", "WBS - Beyersmann (upper)",
        "WBS - Weird calculated (lower)", "WBS - Weird calculated (upper)",
        "WBS - Weird (lower)", "WBS - Weird (upper)"
      )
      colnames(CI) <- colnames(CB) <- paste0("t = ", t)
      
      res[[i]] <- list(CI = CI,
                       CB = CB,
                       event_distribution = event_prob,
                       treatment_probability = treat_prob,
                       EBS_valid = unlist(setNames(lapply(CI_EBS, `[[`, 2), 
                                                   paste0("t = ", t))),
                       times = c("EBS"=EBS_time, "IF"=IF_time, "WBS"=WBS_time))
      
      cov_CI_EBS[i,] <- CI[1,] <= ATE_true & ATE_true <= CI[2,]
      cov_CI_IF[i,] <- CI[3,] <= ATE_true & ATE_true <= CI[4,]
      cov_CI_WBS_calc[i,] <- CI[5,] <= ATE_true & ATE_true <= CI[6,]
      cov_CI_WBS_Lin[i,] <- CI[7,] <= ATE_true & ATE_true <= CI[8,]
      cov_CI_WBS_Bey[i,] <- CI[9,] <= ATE_true & ATE_true <= CI[10,]
      cov_CI_WBS_Weird[i,] <- CI[11,] <= ATE_true & ATE_true <= CI[12,]
      cov_CB_EBS[i] <- all(CB[1,] <= ATE_true & ATE_true <= CB[2,])
      cov_CB_IF[i] <- all(CB[3,] <= ATE_true & ATE_true <= CB[4,])
      cov_CB_WBS_Lin_calc[i] <- all(CB[5,] <= ATE_true & ATE_true <= CB[6,])
      cov_CB_WBS_Lin[i] <- all(CB[7,] <= ATE_true & ATE_true <= CB[8,])
      cov_CB_WBS_Bey_calc[i] <- all(CB[9,] <= ATE_true & ATE_true <= CB[10,])
      cov_CB_WBS_Bey[i] <- all(CB[11,] <= ATE_true & ATE_true <= CB[12,])
      cov_CB_WBS_Weird_calc[i] <- all(CB[13,] <= ATE_true & ATE_true <= CB[14,])
      cov_CB_WBS_Weird[i] <- all(CB[15,] <= ATE_true & ATE_true <= CB[16,])
      
    })
    
    if(i%%(iter/20) == 0){print(paste0(100* i/iter, "%"))}
    
  }
  
  # summarize & save results
  return(list(
    n = n,
    effect = ifelse(beta_0 == 0, "no", 
                    ifelse(beta_0 == 2, "yes", 
                           ifelse(beta_0 == -2, "adverse", NA))),
    treatProb = ifelse(alpha_0 == 0, "normal", 
                       ifelse(alpha_0 == -2, "low", 
                              ifelse(alpha_0 == 2, "high", NA))),
    sdCovariates = ifelse(sd_cov == 1, "1", 
                          ifelse(sd_cov == 0.5, "0.5", 
                                 ifelse(sd_cov == 2, "2", NA))),
    ATE_true = unlist(setNames(ATE_true, paste0("t = ", t))),
    censoring = ifelse(!is.null(m), "typeII", 
                       ifelse(cens == FALSE, "no", 
                              ifelse(cens_par == 200, "low", 
                                     ifelse(cens_par == 50, "high", NA)))),
    results = res,
    errors = sum(is.na(res)),
    mean_event_distribution = 
      apply(sapply(res, 
                   function(i){
                     if(length(i) > 1){
                       i[[3]]
                     }else{
                       matrix(NA, nrow=ifelse(is.null(m),3,2),ncol=length(t))}
                   }, simplify = "array"), 
            c(1,2), mean, na.rm = TRUE),
    mean_treatment_probability = 
      mean(sapply(res, 
                  function(i){
                    if(length(i) > 1){i[[4]]}else{NA}
                  }), na.rm=TRUE),
    mean_EBS_valid = 
      rowMeans(sapply(res, 
                      function(i){
                        if(length(i) > 1){i[[5]]}else{rep(NA, length(t))}
                      }), na.rm = TRUE),
    mean_time = 
      rowMeans(sapply(res, 
                      function(i){
                        if(length(i) > 1){i[[6]]}else{rep(NA, 3)}
                      }), na.rm = TRUE),
    mean_width_CI = 
      matrix(rowMeans(sapply(res, 
                             function(i){
                               if(length(i) > 1){
                                 i[[1]][seq(2,12,2),] - i[[1]][seq(1,11,2),]
                               }else{
                                 rep(NA, 30)
                               }
                             }), na.rm = TRUE), 
             nrow=6, dimnames=list(c("EBS",
                                     "IF",
                                     "WBS - calculated",
                                     "WBS - Lin",
                                     "WBS - Beyersmann",
                                     "WBS - Weird"), 
                                   paste0("t = ", t))),
    mean_width_CB = 
      matrix(rowMeans(sapply(res, 
                             function(i){
                               if(length(i) > 1){
                                 i[[2]][seq(2,16,2),] - i[[2]][seq(1,15,2),]
                               }else{
                                 rep(NA, 40)
                               }
                             }), na.rm = TRUE), 
             nrow=8, dimnames=list(c("EBS",
                                     "IF",
                                     "WBS - Lin calculated",
                                     "WBS - Lin",
                                     "WBS - Beyersmann calculated",
                                     "WBS - Beyersmann",
                                     "WBS - Weird calculated",
                                     "WBS - Weird"), 
                                   paste0("t = ", t))),
    coverage_CI = matrix(100*c(colMeans(cov_CI_EBS, na.rm = TRUE), 
                               colMeans(cov_CI_IF, na.rm = TRUE), 
                               colMeans(cov_CI_WBS_calc, na.rm = TRUE), 
                               colMeans(cov_CI_WBS_Lin, na.rm = TRUE), 
                               colMeans(cov_CI_WBS_Bey, na.rm = TRUE), 
                               colMeans(cov_CI_WBS_Weird, na.rm = TRUE)), 
                         byrow=TRUE, nrow=6, 
                         dimnames=list(c("EBS",
                                         "IF",
                                         "WBS - calculated",
                                         "WBS - Lin",
                                         "WBS - Beyersmann",
                                         "WBS - Weird"), 
                                       paste0("t = ", t))),
    coverage_CB = 100*c("EBS" = mean(cov_CB_EBS, na.rm = TRUE),
                        "IF" = mean(cov_CB_IF, na.rm = TRUE), 
                        "WBS - Lin calculated" = 
                          mean(cov_CB_WBS_Lin_calc, na.rm = TRUE), 
                        "WBS - Lin" = 
                          mean(cov_CB_WBS_Lin, na.rm = TRUE), 
                        "WBS - Beyersmann calculated" = 
                          mean(cov_CB_WBS_Bey_calc, na.rm = TRUE), 
                        "WBS - Beyersmann" = 
                          mean(cov_CB_WBS_Bey, na.rm = TRUE), 
                        "WBS - Weird calculated" = 
                          mean(cov_CB_WBS_Weird_calc, na.rm = TRUE),
                        "WBS - Weird" = 
                          mean(cov_CB_WBS_Weird, na.rm = TRUE))
  ))
  
}


# find true ATE ################################################################

library(riskRegression)
library(prodlim)

# competing risks
ATE_true <- matrix(nrow = 1000, ncol = 5)
for(i in 1:1000){
  data <- generate_data(n=100000, beta_0=-2, true_ATE = TRUE) 
  # data <- generate_data(n=100000, beta_0=-2, sd_cov=0.5, true_ATE = TRUE)
  # data <- generate_data(n=100000, beta_0=-2, sd_cov=2, true_ATE = TRUE)
  # data <- generate_data(n=100000, beta_0=2, true_ATE = TRUE)
  # data <- generate_data(n=100000, beta_0=2, sd_cov=0.5, true_ATE = TRUE)
  # data <- generate_data(n=100000, beta_0=2, sd_cov=2, true_ATE = TRUE)
  ATE_true[i,] <- diff(predictRisk(prodlim(Hist(time, event) ~ A, data=data), 
                                  newdata=data.frame(A=0:1), cause=1, times=c(1,3,5,7,9)))
  if(i %% 100 == 0) print(paste0(i/10, "%"))
}
apply(ATE_true, 2, median, na.rm = TRUE)
# beta_0 = -2, sd_cov = 1:    c(-0.050, -0.221, -0.311, -0.339, -0.341)
# beta_0 = -2, sd_cov = 0.5:  c(-0.033, -0.208, -0.345, -0.398, -0.406)
# beta_0 = -2, sd_cov = 2:    c(-0.097, -0.198, -0.226, -0.233, -0.233)
# beta_0 = 2, sd_cov = 1:     c(0.216, 0.406, 0.365, 0.315, 0.285)
# beta_0 = 2, sd_cov = 0.5:   c(0.192, 0.494, 0.428, 0.345, 0.301) 
# beta_0 = 2, sd_cov = 2:     c(0.207, 0.265, 0.252, 0.236, 0.223)

library(survival)

# type II censoring
t <- c(2,4,6,8,10)
# t <- c(0.5, 1, 1.5, 2, 2.5)
ATE_true <- matrix(nrow = 1000, ncol = length(t))
for(i in 1:1000){
  data <- generate_data(n=100000, beta_0=-2, true_ATE = TRUE, m = 50000)
  # data <- generate_data(n=100000, beta_0=2)
  ATE_true[i,] <- diff(predictRisk(survfit(Surv(time, event) ~ A, data=data), 
                                   newdata=data.frame(A=0:1), times=t))
  if(i %% 100 == 0) print(paste0(i/10, "%"))
}
apply(ATE_true, 2, median, na.rm = TRUE)
# beta_0 = -2:  c(-0.149, -0.317, -0.400, -0.424, -0.415)
# beta_0 = 2:   c(0.084, 0.219, 0.321, 0.384, 0.415)


# execution ####################################################################

library(doParallel)
library(rngtools)
library(doRNG)
library(ATESurvival)

cat("Scenario: no censoring:\n")
for(n in c(50, 75, 100, 200, 300)){
  cat(paste0("n = ", n, "\n"))
  for(beta_0 in c(-2, 0, 2)){
    cat(paste0("n = ", n, ", beta_0 = ", beta_0, ":\n"))
    ATE_true = switch(as.character(beta_0),
                     "-2" = c(-0.050,-0.221,-0.311,-0.339,-0.341),
                     "0" = rep(0,5),
                     "2" = c(0.216,0.406,0.365,0.315,0.285))
    assign(paste0("res_", 
                  ifelse(beta_0 == 0, "no", 
                         ifelse(beta_0 == 2, "", 
                                ifelse(beta_0 == -2, "adv", NA))), 
                  "ATE_noCens_n", n, "_SSS"), 
           run(n = n, beta_0 = beta_0, cens = FALSE, ATE_true = ATE_true))
  }
  cat("\n")
}

cat("Scenario: low censoring:\n")
for(n in c(50, 75, 100, 200, 300)){
  cat(paste0("n = ", n, "\n"))
  for(beta_0 in c(-2, 0, 2)){
    cat(paste0("n = ", n, ", beta_0 = ", beta_0, ":\n"))
    ATE_true = switch(as.character(beta_0),
                     "-2" = c(-0.050,-0.221,-0.311,-0.339,-0.341),
                     "0" = rep(0,5),
                     "2" = c(0.216,0.406,0.365,0.315,0.285))
    assign(paste0("res_", 
                  ifelse(beta_0 == 0, "no", 
                         ifelse(beta_0 == 2, "", 
                                ifelse(beta_0 == -2, "adv", NA))), 
                  "ATE_lowCens_n", n, "_SSS"), 
           run(n = n, beta_0 = beta_0, ATE_true = ATE_true))
  }
  cat("\n")
}

cat("Scenario: high censoring:\n")
for(n in c(50, 75, 100, 200, 300)){
  cat(paste0("n = ", n, "\n"))
  for(beta_0 in c(-2, 0, 2)){
    cat(paste0("n = ", n, ", beta_0 = ", beta_0, ":\n"))
    ATE_true = switch(as.character(beta_0),
                     "-2" = c(-0.050,-0.221,-0.311,-0.339,-0.341),
                     "0" = rep(0,5),
                     "2" = c(0.216,0.406,0.365,0.315,0.285))
    assign(paste0("res_", 
                  ifelse(beta_0 == 0, "no", 
                         ifelse(beta_0 == 2, "", 
                                ifelse(beta_0 == -2, "adv", NA))), 
                  "ATE_highCens_n", n, "_SSS"), 
           run(n = n, beta_0 = beta_0, cens_par = 50, ATE_true = ATE_true))
  }
  cat("\n")
}

cat("Scenario: low treatment probability:\n")
for(n in c(50, 75, 100, 200, 300)){
  cat(paste0("n = ", n, "\n"))
  for(beta_0 in c(-2, 0, 2)){
    cat(paste0("n = ", n, ", beta_0 = ", beta_0, ":\n"))
    ATE_true = switch(as.character(beta_0),
                     "-2" = c(-0.050,-0.221,-0.311,-0.339,-0.341),
                     "0" = rep(0,5),
                     "2" = c(0.216,0.406,0.365,0.315,0.285))
    assign(paste0("res_", 
                  ifelse(beta_0 == 0, "no", 
                         ifelse(beta_0 == 2, "", 
                                ifelse(beta_0 == -2, "adv", NA))), 
                  "ATE_lowTreatProb_n", n, "_SSS"), 
           run(n = n, alpha_0 = -2, beta_0 = beta_0, ATE_true = ATE_true))
  }
  cat("\n")
}

cat("Scenario: high treatment probability:\n")
for(n in c(50, 75, 100, 200, 300)){
  cat(paste0("n = ", n, "\n"))
  for(beta_0 in c(-2, 0, 2)){
    cat(paste0("n = ", n, ", beta_0 = ", beta_0, ":\n"))
    ATE_true = switch(as.character(beta_0),
                     "-2" = c(-0.050,-0.221,-0.311,-0.339,-0.341),
                     "0" = rep(0,5),
                     "2" = c(0.216,0.406,0.365,0.315,0.285))
    assign(paste0("res_", 
                  ifelse(beta_0 == 0, "no", 
                         ifelse(beta_0 == 2, "", 
                                ifelse(beta_0 == -2, "adv", NA))), 
                  "ATE_highTreatProb_n", n, "_SSS"), 
           run(n = n, alpha_0 = 2, beta_0 = beta_0, ATE_true = ATE_true))
  }
  cat("\n")
}

cat("Scenario: low variance of the covariates:\n")
for(n in c(50, 75, 100, 200, 300)){
  cat(paste0("n = ", n, "\n"))
  for(beta_0 in c(-2, 0, 2)){
    cat(paste0("n = ", n, ", beta_0 = ", beta_0, ":\n"))
    ATE_true = switch(as.character(beta_0),
                     "-2" = c(-0.033,-0.208,-0.345,-0.398,-0.406),
                     "0" = rep(0,5),
                     "2" = c(0.192,0.494,0.428,0.345,0.301))
    assign(paste0("res_", 
                  ifelse(beta_0 == 0, "no", 
                         ifelse(beta_0 == 2, "", 
                                ifelse(beta_0 == -2, "adv", NA))), 
                  "ATE_lowVarCov_n", n, "_SSS"), 
           run(n = n, beta_0 = beta_0, sd_cov = 0.5, ATE_true = ATE_true))
  }
  cat("\n")
}

cat("Scenario: high variance of the covariates:\n")
for(n in c(50, 75, 100, 200, 300)){
  cat(paste0("n = ", n, "\n"))
  for(beta_0 in c(-2, 0, 2)){
    cat(paste0("n = ", n, ", beta_0 = ", beta_0, ":\n"))
    ATE_true = switch(as.character(beta_0),
                     "-2" = c(-0.097,-0.198,-0.226,-0.233,-0.233),
                     "0" = rep(0,5),
                     "2" = c(0.207,0.265,0.252,0.236,0.223))
    assign(paste0("res_", 
                  ifelse(beta_0 == 0, "no", 
                         ifelse(beta_0 == 2, "", 
                                ifelse(beta_0 == -2, "adv", NA))), 
                  "ATE_highVarCov_n", n, "_SSS"), 
           run(n = n, beta_0 = beta_0, sd_cov = 2, ATE_true = ATE_true))
  }
  cat("\n")
}

cat("Scenario: type II censoring:\n")
for(n in c(50, 75, 100, 200, 300)){
  cat(paste0("n = ", n, "\n"))
  for(beta_0 in c(-2, 0, 2)){
    cat(paste0("n = ", n, ", beta_0 = ", beta_0, ":\n"))
    ATE_true = switch(as.character(beta_0),
                      "-2" = c(-0.149,-0.317,-0.400,-0.424,-0.415),
                      "0" = rep(0,5),
                      "2" = c(0.084,0.219,0.321,0.384,0.415))
    assign(paste0("res_", 
                  ifelse(beta_0 == 0, "no", 
                         ifelse(beta_0 == 2, "", 
                                ifelse(beta_0 == -2, "adv", NA))), 
                  "ATE_typeII_n", n, "_SSS"), 
           run(n = n, beta_0 = beta_0, m = floor(n/2), ATE_true = ATE_true))
  }
  cat("\n")
}

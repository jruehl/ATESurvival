# This script includes the functions to produce the results of the simulation study.

# generate data ################################################################

# function to generate data
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
  ## sd_cov: standard deviation of the standard normal covariates
  ## true_ATE: logical value indicating whether the true risk difference is to be
  ##   determined based on the generated data
  ##   (i.e. implement random treatment allocation & no censoring)
  ## cens: logical value indicating whether the data should be censored
  ## cens_par: parameter for censoring time distribution (if m is NULL)
  ## m: in case of type II censoring: number of events to observe before censoring
  ##   is imputed (else: NULL)
  
  
  # output: data set with the subsequent columns:
  ## time: (observed) event time
  ## event: event type (either 0,1,2, with 0 indicating censored observations)
  ## A: treatment group (0 or 1)
  ## z1-z12: covariates (z1-z6: normally distributed, z7-z12: binary)
  
  # generate covariates
  Z <- cbind(matrix(rnorm(6*n, 0, sd_cov), nrow=n), 
             matrix(rbinom(6*n, 1, 0.5), nrow=n))
  colnames(Z) <- c("z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12")
  # determine treatment group
  if(true_ATE){
    p_treat <- 0.5
  }else{
    arg_treat <- alpha_0 + colSums(rep(c(1,-1,0,0,0,1),2)*log(2) * t(Z))
    p_treat <- exp(arg_treat) / (1 + exp(arg_treat))
  }
  A <- rbinom(n, 1, p_treat)
  # generate type I event times
  arg_T1 <- beta_0*A + colSums(rep(c(1,0,1,0,0,1),2)*log(2) * t(Z))
  T1 <- rweibull(n, shape = 2, scale = 1/sqrt(1/100*exp(arg_T1)))
  
  if(is.null(m)){ # competing risks scenarios
    # generate type II event times
    arg_T2 <- gamma_0*A + colSums(rep(c(-1,0,0,0,1,1),2)*log(2) * t(Z))
    T2 <- rweibull(n, shape = 2, scale = 1/sqrt(1/100*exp(arg_T2)))
    # generate censoring times
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
  }else{ # standard survival scenario
    if(!cens){ # setting without censoring
      time <- T1
      event <- 1
    }else{ # setting with staggered entry & type II censoring
      # implement staggered entry (calendar time scale)
      entry <- runif(n, 0, min(qweibull(m/n, 2, 1/sqrt(1/100*exp(arg_T1)))))
      # generate type II censoring times (calendar time scale)
      c <- sort(entry + T1)[m]
      event <- as.numeric(entry + T1 <= c)
      # transform observed event times to study time scale
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

# function to run simulations and save summarized results
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
  ## sd_cov: standard deviation of the standard normal covariates
  ## cens: logical value indicating whether the data should be censored
  ## cens_par: parameter for censoring time distribution (if m is NULL)
  ## m: in case of type II censoring: number of events to observe before censoring
  ##   is imputed (else: NULL)
  ## t: time point(s) to evaluate
  ## ATE_true: true average treatment effect
  ## EBS_iter: number of samples to use for the EBS approach
  ## BS_iter: number of random multipliers to use for the IF/WBS approach
  ## iter: number of Monte Carlo repetitions
  ## seed: seed for reproducibility
  ## cores: number of CPU cores for parallel computations
  
  # output: list with (summarized) simulation results
  
  
  # set time points if t is not specified
  if(is.null(t)){
    if(is.null(m)){ # competing risks scenarios
      t <- c(1,3,5,7,9)
    }else if(beta_0 == -2){
      t <- c(2,4,6,8,10)
    }else if(beta_0 == 0){
      t <- c(1,2,3,4,5)
    }else if (beta_0 == 2){
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
    cov_CB_WBS_Lin_calc <- cov_CB_WBS_Lin <- 
    cov_CB_WBS_Bey_calc <- cov_CB_WBS_Bey <- 
    cov_CB_WBS_Weird_calc <- cov_CB_WBS_Weird <- 
    rep(NA, iter)
  
  # prepare results list
  res <- list()
  set.seed(seed, kind = "Mersenne-Twister") # for reproducibility
  
  for(i in 1:iter){
    
    res[[i]] <- NA
    
    try({
      
      # generate data
      if(is.null(m)){ # competing risks scenarios
        while(TRUE){ 
          data <- generate_data(
            n = n, 
            alpha_0 = alpha_0, beta_0 = beta_0, 
            gamma_0 = gamma_0, delta_0 = delta_0, 
            sd_cov = sd_cov, 
            cens = cens, cens_par = cens_par
          )
          # ensure that at least 10 events of each type are observed
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
      }else{ # setting with staggered entry & type II censoring
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
      if(.Platform$OS.type == "windows"){
        cl <- makeCluster(cores)
        registerDoParallel(cl)
      }else{
        registerDoParallel(cores = cores)
      }
      if(is.null(m)){ # competing risks scenarios
        res_EBS <- 
          foreach(j = 1:EBS_iter, .combine = cbind, 
                  .packages = c('riskRegression', 'ATESurvival'), 
                  .options.RNG = seed+i) %dorng% { 
                    tryCatch({
                      # draw from data with replacement
                      data_bs <- data[sample(1:n, n, replace = TRUE),]
                      # fit cause-specific Cox models
                      invisible(capture.output(
                        csc_bs <- CSC(
                          formula = 
                            list(Hist(time, event) ~ 
                                   A+z1+z2+z3+z4+z5+z6+z7+z8+z9+z10+z11+z12,
                                 Hist(time, event) ~ 
                                   A+z1+z2+z3+z4+z5+z6+z7+z8+z9+z10+z11+z12),
                          data_bs[order(data_bs$time),])
                      ))
                      # calculate average treatment effect
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
      }else{ # setting with staggered entry & type II censoring
        res_EBS <- 
          foreach(j = 1:EBS_iter, .combine = cbind, 
                  .packages = c('survival', 'ATESurvival'), 
                  .options.RNG = seed+i) %dorng% { 
                    tryCatch({
                      # draw from data with replacement
                      data_bs <- data[sample(1:n, n, replace = TRUE),]
                      # fit Cox model
                      invisible(capture.output(
                        cox_bs <- coxph(
                          Surv(time, event) ~ 
                            A+z1+z2+z3+z4+z5+z6+z7+z8+z9+z10+z11+z12,
                          data_bs[order(data_bs$time),]
                        )
                      ))
                      # calculate average treatment effect
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
      if(.Platform$OS.type == "windows"){
        stopCluster(cl)
        registerDoSEQ()
      }
      # derive confidence intervals
      CI_EBS <- apply(res_EBS, 1, function(res_EBS_t){
        res_EBS_t <- res_EBS_t[is.finite(res_EBS_t) & !is.na(res_EBS_t)]
        if(length(res_EBS_t) == 0) return(list(rep(NA, 2), 0))
        return(list(c(max(quantile(res_EBS_t, 0.025), -1), 
                      min(quantile(res_EBS_t, 0.975), 1)), 
                    length(res_EBS_t)))
      })
      # determine quantile for confidence bands
      q_EBS_sup <- quantile(apply(abs(res_EBS - rowMeans(res_EBS, na.rm = TRUE)) / 
                                    sqrt(apply(res_EBS, 1, var, na.rm = TRUE)), 
                                  2, max, na.rm = TRUE), 0.95)
      EBS_time <- ifelse(any(!is.na(unlist(lapply(CI_EBS, `[[`, 1)))), 
                         (proc.time() - EBS_time)[3]*1e3, 
                         NA)
      
      
      # prepare IF/WBS
      Cox_time <- proc.time()
      if(is.null(m)){ # competing risks scenarios
        # fit cause-specific Cox models
        csc <- CSC(formula = list(Hist(time, event) ~ 
                                    A+z1+z2+z3+z4+z5+z6+z7+z8+z9+z10+z11+z12,
                                  Hist(time, event) ~ 
                                    A+z1+z2+z3+z4+z5+z6+z7+z8+z9+z10+z11+z12),
                   data[order(data$time),])
      }else{ # setting with staggered entry & type II censoring
        # fit Cox model
        cox <- coxph(Surv(time, event) ~ 
                       A+z1+z2+z3+z4+z5+z6+z7+z8+z9+z10+z11+z12,
                     data[order(data$time),])
      }
      Cox_time <- (proc.time() - Cox_time)[3]*1e3
      
      # create matrices to store average treatment effect, confidence intervals & bands
      res_IF_WBS <- list(rep(NA, length(t)),
                         matrix(NA, nrow=2, ncol=length(t)), 
                         matrix(NA, nrow=2, ncol=length(t)),
                         matrix(NA, nrow=8, ncol=length(t)), 
                         matrix(NA, nrow=12, ncol=length(t)))
      # perform IF/WBS
      IF_time <- NA
      WBS_time <- NA
      try({
        # determine number of observed (uncensored) events for multipliers
        event_num <- sum(data$event[order(data$time)] != 0)
        # compute weird bootstrap multipliers
        mult_Weird <- sapply(1:event_num, function(j){
          time <- sort(data$time)
          Y <- sum(time[data$event[order(data$time)] != 0][j] <= time)
          rbinom(BS_iter, Y, 1/Y) - 1
        })
        if(is.null(m)){ # competing risks scenarios
          # calculate average treatment effect, confidence intervals & bands
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
            bs_iter = BS_iter,
            G_IF = rnorm(n * BS_iter),
            G_Lin_init = rnorm(BS_iter * event_num),
            G_Bey_init = rpois(BS_iter * event_num, 1) - 1,
            G_Weird_init = c(t(mult_Weird))
          )
        }else{ # setting with staggered entry & type II censoring
          # calculate average treatment effect, confidence intervals & bands
          res_IF_WBS <- ATE_IF_WBS(
            Z = list(model.matrix(cox)),
            event = data$event[order(data$time)], 
            time = sort(data$time),
            t = t,
            beta = list(cox$coefficients),
            index_A = c(which(colnames(model.matrix(cox)) == 
                                paste0("A", levels(data$A)[2])), NA)[1],
            bs_iter = BS_iter,
            G_IF = rnorm(n * BS_iter),
            G_Lin_init = rnorm(BS_iter * event_num),
            G_Bey_init = rpois(BS_iter * event_num, 1) - 1,
            G_Weird_init = c(t(mult_Weird))
          )
        }
        IF_time <- unname(
          Cox_time + sum(clock$timer[clock$ticker %in% c("ATE","IF")]/1e6)
        )
        WBS_time <- unname(
          Cox_time + sum(clock$timer[clock$ticker %in% c("ATE","WBS")]/1e6)
        )
      }, silent = TRUE)
      
      
      # store confidence intervals
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
      # store confidence bands
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
      
      # store results
      res[[i]] <- list(CI = CI,
                       CB = CB,
                       event_distribution = event_prob,
                       treatment_probability = treat_prob,
                       EBS_valid = unlist(setNames(lapply(CI_EBS, `[[`, 2), 
                                                   paste0("t = ", t))),
                       times = c("EBS"=EBS_time, "IF"=IF_time, "WBS"=WBS_time))
      
      # store coverage probabilities
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
    # sample size
    n = n,
    # direction of treatment effect on type I events
    effect = ifelse(beta_0 == 0, "no", 
                    ifelse(beta_0 == 2, "yes", 
                           ifelse(beta_0 == -2, "adverse", NA))),
    # treatment probability
    treatProb = ifelse(alpha_0 == 0, "normal", 
                       ifelse(alpha_0 == -2, "low", 
                              ifelse(alpha_0 == 2, "high", NA))),
    # standard deviation of the standard normal covariates
    sdCovariates = ifelse(sd_cov == 1, "1", 
                          ifelse(sd_cov == 0.5, "0.5", 
                                 ifelse(sd_cov == 2, "2", NA))),
    # true average treatment effect
    ATE_true = unlist(setNames(ATE_true, paste0("t = ", t))),
    # extent of censoring
    censoring = ifelse(!is.null(m), "typeII", 
                       ifelse(cens == FALSE, "no", 
                              ifelse(cens_par == 200, "low", 
                                     ifelse(cens_par == 50, "high", NA)))),
    # confidence intervals & bands, event time distributions, treatment probabilities, 
    #   number of valid EBS samples, and computation times for each 
    #   iteration
    results = res,
    # number of iterations with errors
    errors = sum(is.na(res)),
    # average event time distribution
    mean_event_distribution = 
      apply(sapply(res, 
                   function(i){
                     if(length(i) > 1){
                       i[[3]]
                     }else{
                       matrix(NA, nrow=ifelse(is.null(m),3,2),ncol=length(t))
                     }
                   }, simplify = "array"), 
            c(1,2), mean, na.rm = TRUE),
    # average treatment probability
    mean_treatment_probability = 
      mean(sapply(res, 
                  function(i){
                    if(length(i) > 1){i[[4]]}else{NA}
                  }), na.rm=TRUE),
    # average number of valid EBS samples
    mean_EBS_valid = 
      rowMeans(sapply(res, 
                      function(i){
                        if(length(i) > 1){i[[5]]}else{rep(NA, length(t))}
                      }), na.rm = TRUE),
    # average computation times
    mean_time = 
      rowMeans(sapply(res, 
                      function(i){
                        if(length(i) > 1){i[[6]]}else{rep(NA, 3)}
                      }), na.rm = TRUE),
    # average confidence interval widths
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
    # average confidence band widths
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
    # confidence interval coverage probabilities
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
    # confidence band coverage probabilities
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

# create plots #################################################################

## confidence interval coverage ####

# function to plot the confidence interval coverage at a given time point
plot_coverage_CI_t <- function(scenario, effect, t, ylim){
  # input:
  ## scenario: simulation scenario
  ##   ('noCens', 'lowCens', 'highCens', 
  ##   'lowTreatProb', 'highTreatProb', 
  ##   'lowVarCov', 'highVarCov',
  ##   'typeII')
  ## effect: underlying ATE ('adverse', 'no', 'yes')
  ## t: time point
  ## ylim: common limit of the y axis
  
  # determine y axis breaks
  if(ylim[2]-ylim[1] <= 10){
    labels <- as.character(c(seq(0,50,10), 
                             seq(55,90,5), 
                             seq(91,100)))
  }else if(ylim[2]-ylim[1] <= 15){
    labels <- as.character(c(seq(0,50,10), 
                             seq(55,90,5), 
                             c(rbind(rep("", 5), 
                                     seq(92, 100, 2)))))
  }else{
    labels <- as.character(c(seq(0,50,10), 
                             seq(55,90,5), 
                             c(rbind(rep("", 2),rep("", 2),rep("", 2),rep("", 2), 
                                     c(95, 100)))))
  }
  # create plot
  ggplot(data = total_coverage_CI[
    total_coverage_CI$scenario == scenario & 
      total_coverage_CI$effect == effect & 
      total_coverage_CI$t == t & 
      total_coverage_CI$type != "WBS (calculated)",
  ], 
  aes(x = n, y = coverage)) +
    geom_hline(aes(yintercept = 95)) +
    geom_line(aes(color = type), linewidth=1.1) +
    geom_point(aes(color = type), size=2) +
    scale_color_manual("", 
                       values = c('orange','red','turquoise4','cyan3','chartreuse2')) +
    labs(title = paste0("t = ", t), 
         x = "n", y = ifelse(t == 
                               ifelse((scenario != "typeII") | (effect == "no"), 1,
                                      ifelse(effect == "adverse", 2, 0.5)), 
                             "coverage [%]", "")) +
    scale_x_continuous(breaks = c(50,75,100,200,300)) + 
    scale_y_continuous(limits = ylim, 
                       breaks = c(0,10,20,30,40,50,55,60,65,70,75,80,85,90:100),
                       labels = labels) +
    theme(text = element_text(size=15),
          plot.title = element_text(hjust=0.5),
          axis.title.x = element_text(vjust=-0.5),
          axis.text.x = element_text(angle=90, vjust=0.5),
          panel.grid.minor = element_line(colour=NA),
          panel.grid.major = element_line(colour='lightgrey'),
          panel.background = element_rect(fill=NA), 
          panel.border = element_rect(fill=NA),
          legend.key = element_rect(fill=NA),
          legend.background = element_rect(colour="black"))
}

# function to plot the confidence interval coverage
plot_coverage_CI <- function(scenario, effect){
  # input:
  ## scenario: simulation scenario
  ##   ('noCens', 'lowCens', 'highCens', 
  ##   'lowTreatProb', 'highTreatProb', 
  ##   'lowVarCov', 'highVarCov',
  ##   'typeII')
  ## effect: underlying ATE ('adverse', 'no', 'yes')
  
  # determine common limit of the y axis
  ylim <- total_coverage_CI[total_coverage_CI$scenario == scenario & 
                              total_coverage_CI$effect == effect & 
                              total_coverage_CI$type != "WBS (calculated)",
                            "coverage"]
  ylim <- c(floor(min(ylim)*10)/10, ceiling(max(ylim)*10)/10)
  ylim <- c(pmin(95, ylim[1]), pmax(95, ylim[2]))
  # determine time points considered in the respective scenario
  if(scenario != "typeII"){
    t <- c(1,3,5,7,9)
  }else if(effect == "adverse"){
    t <- c(2,4,6,8,10)
  }else if(effect == "no"){
    t <- c(1,2,3,4,5)
  }else if(effect == "yes"){
    t <- c(0.5,1,1.5,2,2.5)
  }
  # create plots of the confidence interval coverages for each time point
  p1 <- plot_coverage_CI_t(scenario, effect, t[1], ylim)
  p2 <- plot_coverage_CI_t(scenario, effect, t[2], ylim)
  p3 <- plot_coverage_CI_t(scenario, effect, t[3], ylim)
  p4 <- plot_coverage_CI_t(scenario, effect, t[4], ylim)
  p5 <- plot_coverage_CI_t(scenario, effect, t[5], ylim)
  # combine
  p1 + p2 + p3 + p4 + p5 + 
    plot_layout(nrow = 1, guides = "collect") + 
    plot_annotation(title = 
                      paste0("Confidence interval coverage (", 
                             ifelse(effect == "yes", "", paste0(effect, " ")), 
                             "treatment effect)")) & 
    theme(plot.title = element_text(hjust=0.5, size=18, face="bold")) &
    theme(legend.position = 'bottom')
  
}

## confidence band coverage ####

# function to plot the confidence band coverage
plot_coverage_CB <- function(scenario, effect){
  # input:
  ## scenario: simulation scenario
  ##   ('noCens', 'lowCens', 'highCens', 
  ##   'lowTreatProb', 'highTreatProb', 
  ##   'lowVarCov', 'highVarCov',
  ##   'typeII')
  ## effect: underlying ATE ('adverse', 'no', 'yes')
  
  # determine y axis breaks
  ylim <- total_coverage_CB[total_coverage_CB$scenario == scenario & 
                              total_coverage_CB$effect == effect & 
                              !(total_coverage_CB$type %in% 
                                  c("WBS - Lin et al. (calculated)",
                                    "WBS - Beyersmann et al. (calculated)",
                                    "WBS - Weird bootstrap (calculated)")),
                            "coverage"]
  ylim <- c(floor(min(ylim)*10)/10, ceiling(max(ylim)*10)/10)
  ylim <- c(pmin(95, ylim[1]), pmax(95, ylim[2]))
  if(ylim[2]-ylim[1] <= 10){
    labels <- as.character(c(seq(0,50,10), 
                             seq(55,90,5), 
                             seq(91,100)))
  }else if(ylim[2]-ylim[1] <= 15){
    labels <- as.character(c(seq(0,50,10), 
                             seq(55,90,5), 
                             c(rbind(rep("", 5), 
                                     seq(92, 100, 2)))))
  }else{
    labels <- as.character(c(seq(0,50,10), 
                             seq(55,90,5), 
                             c(rbind(rep("", 2),rep("", 2),rep("", 2),rep("", 2), 
                                     c(95, 100)))))
  }
  # create plot
  ggplot(data = total_coverage_CB[
    total_coverage_CB$scenario == scenario & 
      total_coverage_CB$effect == effect & 
      !(total_coverage_CB$type %in% 
          c("WBS - Lin et al. (calculated)",
            "WBS - Beyersmann et al. (calculated)",
            "WBS - Weird bootstrap (calculated)")),
  ], 
  aes(x = n, y = coverage)) +
    geom_hline(aes(yintercept = 95)) +
    geom_line(aes(color = type), linewidth=1.1) +
    geom_point(aes(color = type), size=2) +
    scale_color_manual("", 
                       values = c('orange','red','turquoise4','cyan3','chartreuse2')) +
    labs(title = paste0("Confidence band coverage (", 
                        ifelse(effect == "yes", "", paste0(effect, " ")), 
                        "treatment effect)"), 
         x = "n", y = "coverage [%]") +
    scale_x_continuous(breaks = c(50,75,100,200,300)) + 
    scale_y_continuous(limits = ylim, 
                       breaks = c(0,10,20,30,40,50,55,60,65,70,75,80,85,90:100),
                       labels = labels) +
    theme(text = element_text(size=15),
          plot.title = element_text(hjust=0.5, face="bold"),
          axis.title.x = element_text(vjust=-0.5),
          axis.text.x = element_text(angle=90, vjust=0.5),
          panel.grid.minor = element_line(colour=NA),
          panel.grid.major = element_line(colour='lightgrey'),
          panel.background = element_rect(fill=NA), 
          panel.border = element_rect(fill=NA),
          legend.key = element_rect(fill=NA),
          legend.background = element_rect(colour="black"),
          legend.position = "bottom")
}

## confidence interval/band width ####

# function to plot the confidence interval/band width at a given time point
plot_width_t <- function(scenario, effect, type, t){
  # input:
  ## scenario: simulation scenario
  ##   ('noCens', 'lowCens', 'highCens', 
  ##   'lowTreatProb', 'highTreatProb', 
  ##   'lowVarCov', 'highVarCov',
  ##   'typeII')
  ## effect: underlying ATE ('adverse', 'no', 'yes')
  ## type: confidence region type 
  ##   ('CI' = confidence interval, 'CB': confidence band)
  ## t: time point
  
  # determine time points considered in the respective scenario
  if(scenario != "typeII"){
    all_t <- c("1","3","5","7","9")
  }else if(effect == "adverse"){
    all_t <- c("2","4","6","8","10")
  }else if(effect == "no"){
    all_t <- c("1","2","3","4","5")
  }else if(effect == "yes"){
    all_t <- c("0.5","1","1.5","2","2.5")
  }
  # retrieve data on confidence interval/band widths
  data <- data.frame()
  if(type == "CI"){ # confidence intervals
    for(n in c(50, 75, 100, 200, 300)){
      data <- rbind(
        data,
        data.frame(
          n = n, 
          type = factor(rep(c("EBS",
                              "IF",
                              "WBS (calculated)",
                              "WBS - Lin et al.",
                              "WBS - Beyersmann et al.",
                              "WBS - Weird bootstrap"), 
                            5000), 
                        levels = c("EBS",
                                   "IF",
                                   "WBS - Lin et al.",
                                   "WBS - Beyersmann et al.",
                                   "WBS - Weird bootstrap")),
          width = c(unlist(sapply(eval(parse(
            text = paste0("res_", 
                          ifelse(effect == "no", "no", 
                                 ifelse(effect == "yes", "", "adv")), 
                          "ATE_", scenario, "_n", n, "$results")
          )), 
          function(i){
            if(length(i) > 1){
              i[[1]][seq(2,12,2),
                     which(as.character(t) == all_t)] - 
                i[[1]][seq(1,11,2),
                       which(as.character(t) == all_t)]
            }else{
              rep(NA, 6)
            }
          }
          )))
        )
      )
    }
  }else if(type == "CB"){ # confidence bands
    for(n in c(50, 75, 100, 200, 300)){
      data <- rbind(
        data,
        data.frame(
          n = n, 
          type = factor(rep(c("EBS", 
                              "IF",
                              "WBS - Lin et al. (calculated)",
                              "WBS - Lin et al.",
                              "WBS - Beyersmann et al. (calculated)",
                              "WBS - Beyersmann et al.",
                              "WBS - Weird bootstrap (calculated)",
                              "WBS - Weird bootstrap"), 
                            5000), 
                        levels = c("EBS", 
                                   "IF",
                                   "WBS - Lin et al.",
                                   "WBS - Beyersmann et al.",
                                   "WBS - Weird bootstrap")),
          width = c(unlist(sapply(eval(parse(
            text = paste0("res_", 
                          ifelse(effect == "no", "no", 
                                 ifelse(effect == "yes", "", "adv")), 
                          "ATE_", scenario, "_n", n, "$results")
          )), 
          function(i){
            if(length(i) > 1){
              i[[2]][seq(2,16,2),
                     which(as.character(t) == all_t)] - 
                i[[2]][seq(1,15,2),
                       which(as.character(t) == all_t)]
            }else{
              rep(NA, 8)
            }
          }
          )))
        )
      )
    }
  }
  data <- data[!is.na(data$type),]
  # identify mild & extreme outliers
  q1 <- aggregate(width ~ n + type, data, quantile, 0.25)
  q3 <- aggregate(width ~ n + type, data, quantile, 0.75)
  outliers <- data.frame(
    n = q1$n,
    type = q1$type,
    whisker_low = q1$width - 1.5 * (q3$width - q1$width),
    whisker_up = q3$width + 1.5 * (q3$width - q1$width),
    extreme_outlier_low = q1$width - 3 * (q3$width - q1$width),
    extreme_outlier_up = q3$width + 3 * (q3$width - q1$width)
  )
  data <- merge(data, outliers, by = c("n", "type"), all.x = TRUE)
  data$outlier <- ifelse(data$width < data$whisker_low | 
                           data$width > data$whisker_up, 
                         "mild", NA)
  data$outlier <- ifelse(data$width < data$extreme_outlier_low | 
                           data$width > data$extreme_outlier_up,
                         "extreme", data$outlier)
  data$outlier <- factor(data$outlier, levels = c("mild","extreme"))
  
  # create plot
  ggplot(data = rbind(data, 
                      c(125, rep(NA, 7)), 
                      c(150, rep(NA, 7)), 
                      c(175, rep(NA, 7)),
                      c(225, rep(NA, 7)),
                      c(250, rep(NA, 7)),
                      c(275, rep(NA, 7))), 
         aes(x = factor(n), y = width, fill = type)) + 
    stat_boxplot(geom="errorbar", position="dodge") +
    geom_boxplot(position="dodge", outlier.shape = NA) + 
    geom_point(aes(col = outlier, group = type), size = 0.5, 
               position = position_dodge(width=0.75), show.legend = FALSE) + 
    scale_fill_manual("", values = c('orange','red','turquoise4','cyan3','chartreuse2'),
                      breaks = c("EBS","IF","WBS - Lin et al.","WBS - Beyersmann et al.","WBS - Weird bootstrap")) +
    scale_color_manual("", values = c("black","gray50"), na.translate = FALSE) + 
    labs(title = paste0("Confidence ", ifelse(type == "CI", "interval", "band"), 
                        " width at t = ", t, 
                        " (", ifelse(effect == "yes", "", paste0(effect, " ")), 
                        "treatment effect)"), 
         x = "n", y = "width") +
    scale_x_discrete(breaks = c(50,75,100,200,300)) +
    theme(text = element_text(size=15),
          plot.title = element_text(hjust=0.5, face="bold"),
          axis.title.x = element_text(vjust=-0.5),
          axis.text.x = element_text(angle=90, vjust=0.5),
          panel.grid.minor = element_line(colour=NA),
          panel.grid.major = element_line(colour='lightgrey'),
          panel.background = element_rect(fill=NA), 
          panel.border = element_rect(fill=NA),
          legend.key = element_rect(fill=NA),
          legend.title = element_blank(),
          legend.background = element_rect(colour="black"),
          legend.position = "bottom")
  
}

## computation times ####

# function to plot the computation times
plot_computation_times <- function(scenario, effect){
  # input:
  ## scenario: simulation scenario
  ##   ('noCens', 'lowCens', 'highCens', 
  ##   'lowTreatProb', 'highTreatProb', 
  ##   'lowVarCov', 'highVarCov',
  ##   'typeII')
  ## effect: underlying ATE ('adverse', 'no', 'yes')
  
  # retrieve data on computation times
  data <- data.frame()
  for(n in c(50, 75, 100, 200, 300)){
    data <- rbind(
      data,
      data.frame(
        n = n, 
        type = rep(c("EBS","IF","WBS"), each = 5000),
        time = c(t(sapply(eval(parse(
          text = paste0("res_", 
                        ifelse(effect == "no", "no", 
                               ifelse(effect == "yes", "", "adv")), 
                        "ATE_", scenario, "_n", n, "$results")
        )),
        function(i){if(length(i) > 1){i[[6]]/1000}else{rep(NA, 3)}})))
      ))
  }
  data <- data[!is.na(data$type) & !is.na(data$time),]
  # identify mild & extreme outliers
  q1 <- aggregate(time ~ n + type, data, quantile, 0.25)
  q3 <- aggregate(time ~ n + type, data, quantile, 0.75)
  mean_val <- aggregate(time ~ n + type, data, mean)
  outliers <- data.frame(
    n = q1$n,
    type = q1$type,
    mean = mean_val$time,
    whisker_low = q1$time - 1.5 * (q3$time - q1$time),
    whisker_up = q3$time + 1.5 * (q3$time - q1$time),
    extreme_outlier_low = q1$time - 3 * (q3$time - q1$time),
    extreme_outlier_up = q3$time + 3 * (q3$time - q1$time)
  )
  data <- merge(data, outliers, by = c("n", "type"), all.x = TRUE)
  data$outlier <- ifelse(data$time < data$whisker_low | 
                           data$time > data$whisker_up, 
                         "mild", NA)
  data$outlier <- ifelse(data$time < data$extreme_outlier_low | 
                           data$time > data$extreme_outlier_up,
                         "extreme", data$outlier)
  data$outlier <- factor(data$outlier, levels = c("mild","extreme"))
  
  # create plot
  ggplot(data = rbind(data, 
                      c(125, rep(NA, 8)), 
                      c(150, rep(NA, 8)), 
                      c(175, rep(NA, 8)),
                      c(225, rep(NA, 8)),
                      c(250, rep(NA, 8)),
                      c(275, rep(NA, 8))), 
         aes(x = factor(n), y = time, fill = type)) +
    geom_bar(aes(y = mean), stat = "identity", position = position_dodge(width=0.9)) + 
    stat_boxplot(geom="errorbar", width = 0.5, position = position_dodge(width=0.9)) +
    geom_boxplot(width = 0.5, position = position_dodge(width=0.9), 
                 outlier.shape = NA, show.legend = FALSE) + 
    geom_point(aes(col = outlier, group = type), size = 1, 
               position = position_dodge(width=0.9), show.legend = FALSE) + 
    scale_fill_manual("", values = c('orange','red','turquoise'), 
                      na.translate = F) +
    scale_color_manual("", values = c("black","gray50"), na.translate = FALSE) + 
    guides(fill = guide_legend(override.aes = list(shape = NA))) +
    labs(title = paste0("Computation times (", 
                        ifelse(effect == "yes", "", paste0(effect, " ")), 
                        "treatment effect)"), 
         x = "n", y = "time [s]") +
    scale_x_discrete(breaks = c(50,75,100,200,300)) +
    theme(text = element_text(size=15),
          plot.title = element_text(hjust=0.5, face="bold"),
          axis.title.x = element_text(vjust=-0.5),
          axis.text.x = element_text(angle=90, vjust=0.5),
          panel.grid.minor = element_line(colour=NA),
          panel.grid.major = element_line(colour='lightgrey'),
          panel.background = element_rect(fill=NA), 
          panel.border = element_rect(fill=NA),
          legend.key = element_rect(fill=NA),
          legend.title = element_blank(),
          legend.background = element_rect(colour="black"),
          legend.position = "bottom")
  
}

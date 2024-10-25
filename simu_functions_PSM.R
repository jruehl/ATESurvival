# This script includes the functions to produce the results of the simulation study
# w.r.t. PS matching.

# generate data ################################################################

# function to generate data
generate_data <- function(n, 
                          alpha_0 = 0, beta_0, gamma_0 = 0, delta_0 = 0, 
                          sd_cov = 1, 
                          true_ATE = FALSE, 
                          cens = !true_ATE, cens_par = 200, m = NULL, 
                          method = "cs"){
  
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
  ## method: string denoting the data generating method 
  ##   (ozenne: latent failure time model, cs: cause-specific Cox model)
  
  
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
  # generate event times
  arg_T1 <- beta_0*A + colSums(rep(c(1,0,1,0,0,1),2)*log(2) * t(Z))
  if(is.null(m)){
    arg_T2 <- gamma_0*A + colSums(rep(c(-1,0,0,0,1,1),2)*log(2) * t(Z))
  }
  if(method  == "ozenne" | !is.null(m)){ # latent failure time model
    # generate type I event times
    T1 <- rweibull(n, shape = 2, scale = 10/sqrt(exp(arg_T1)))
    if(is.null(m)){ # competing risks scenarios
      # generate type II event times
      T2 <- rweibull(n, shape = 2, scale = 10/sqrt(exp(arg_T2)))
      # generate censoring times
      if(!cens){
        C <- Inf
      }else{
        arg_C <- delta_0*A + colSums(rep(c(-1,0,0,1,0,-1),2)*log(2) * t(Z))
        C <- rweibull(n, shape = 2, scale = sqrt(cens_par/exp(arg_C)))
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
        entry <- runif(n, 0, min(qweibull(m/n, 2, 10/sqrt(exp(arg_T1)))))
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
  }else{ # cause-specific model
    # generate waiting times
    time <- rweibull(n, shape = 2, scale = 10/sqrt(exp(arg_T1) + exp(arg_T2)))
    # generate cause
    cause <- rbinom(n, size = 1, prob = exp(arg_T2) / (exp(arg_T1)+exp(arg_T2))) + 1
    # generate censoring times
    if(!cens){
      C <- Inf
    }else{
      arg_C <- delta_0*A + colSums(rep(c(-1,0,0,1,0,-1),2)*log(2) * t(Z))
      C <- rweibull(n, shape = 2, scale = sqrt(cens_par/exp(arg_C)))
    }
    # save data
    data <- data.frame(time = pmin(time, C), 
                       event = ifelse(time <= C, cause, 0),
                       A = as.factor(A),
                       Z)
  }
  
  return(data)
  
}


# run simulations ##############################################################

# function to run naive EBS
EBS <- function(m, t, n, data, EBS_iter, i, seed, cores){
  
  # input:
  ## m: in case of type II censoring: number of events to observe before censoring
  ##   is imputed (else: NULL)
  ## t: time point(s) to evaluate
  ## n: sample size (of original data)
  ## data: if resample_method == 'naive': original data set, else: matched data
  ## EBS_iter: number of samples to use
  ## i: number of current iteration in the 'run' function (for reproducibility)
  ## seed: seed for reproducibility
  ## cores: number of CPU cores for parallel computations
  
  
  # output: matrix of bootstrap ATE for each time point (rows) & sample (columns)
  

  # perform EBS
  if(.Platform$OS.type == "windows"){
    cl <- makeCluster(cores)
    registerDoParallel(cl)
  }else{
    registerDoParallel(cores = cores)
  }
  if(is.null(m)){
    res_EBS <- 
      foreach(j = 1:EBS_iter, .combine = cbind, 
              .packages = c('riskRegression','ATESurvival'),
              .options.RNG = seed+i) %dorng% {
                tryCatch({
                  # draw from data with replacement
                  data_bs <- data[sample(1:n, n, replace = TRUE),]
                  data_bs <- data_bs[order(data_bs$time),]
                  # determine propensity scores
                  data_bs$PS <- predict(glm(A ~ z1+z2+z3+z4+z5+z6+z7+z8+z9+z10+z11+z12, data_bs, 
                                            family = "binomial"), 
                                        data_bs, type = "response")
                  # determine matching bias
                  PS_B_bs <- PS_bias(ID = data_bs$ID,
                                     A = as.numeric(data_bs$A)-1,
                                     event = data_bs$event,
                                     time = data_bs$time,
                                     beta = coefficients(CSC(
                                       list(Hist(time, event) ~ A + PS, 
                                            Hist(time, event) ~ A + PS), 
                                       data_bs, ties = "breslow")),
                                     PS = data_bs$PS,
                                     t = t)
                  # match data
                  data_matched_bs <- rbind(data_bs[,c("ID","time","event","A")],
                                           merge(merge(data_bs[!duplicated(data_bs$ID),c("ID","time","event","A")], 
                                                       matrix(PS_B_bs$matches, ncol=2, dimnames=list(NULL, c("ID_old", "ID"))), 
                                                       by="ID"),
                                                 matrix(data_bs$ID, ncol=1, dimnames=list(NULL, c("ID_old"))),
                                                 by = "ID_old")[,c("ID","time","event","A")])
                  data_matched_bs <- data_matched_bs[order(data_matched_bs$time, data_matched_bs$ID),]
                  # estimate ATE
                  cox_bs <- CSC(list(Hist(time, event) ~ A, 
                                     Hist(time, event) ~ A), 
                                data_matched_bs, ties = "breslow")
                  as.numeric(ATE(ID = data_matched_bs$ID,
                                 Z = list(cox_bs$models$`Cause 1`$x, 
                                          cox_bs$models$`Cause 2`$x),
                                 index_A = c(c(which(colnames(cox_bs$models$`Cause 1`$x) == 
                                                       paste0("A", levels(data_matched_bs$A)[2])), NA)[1],
                                             c(which(colnames(cox_bs$models$`Cause 2`$x) == 
                                                       paste0("A", levels(data_matched_bs$A)[2])), NA)[1]),
                                 event = data_matched_bs$event,
                                 time = data_matched_bs$time,
                                 beta = coefficients(cox_bs),
                                 t = t,
                                 IF = FALSE, WBS = FALSE)$ATE - PS_B_bs$bias)
                }, error = function(e){rep(NA, length(t))})
              }
  }else{
    res_EBS <- 
      foreach(j = 1:EBS_iter, .combine = cbind, 
              .packages = c('survival','ATESurvival'),
              .options.RNG = seed+i) %dorng% { 
                tryCatch({
                  # draw from data with replacement
                  data_bs <- data[sample(1:n, n, replace = TRUE),]
                  data_bs <- data_bs[order(data_bs$time),]
                  # determine propensity scores
                  data_bs$PS <- predict(glm(A ~ z1+z2+z3+z4+z5+z6+z7+z8+z9+z10+z11+z12, data_bs, 
                                            family = "binomial"), 
                                        data_bs, type = "response")
                  # determine matching bias
                  PS_B_bs <- PS_bias(ID = data_bs$ID,
                                     A = as.numeric(data_bs$A)-1,
                                     event = data_bs$event,
                                     time = data_bs$time,
                                     beta = list(coefficients(coxph(
                                       Surv(time, event) ~ A + PS, 
                                       data_bs, ties = "breslow"))),
                                     PS = data_bs$PS,
                                     t = t)
                  # match data
                  data_matched_bs <- rbind(data_bs[,c("ID","time","event","A")],
                                           merge(merge(data_bs[!duplicated(data_bs$ID),c("ID","time","event","A")], 
                                                       matrix(PS_B_bs$matches, ncol=2, dimnames=list(NULL, c("ID_old", "ID"))), 
                                                       by="ID"),
                                                 matrix(data_bs$ID, ncol=1, dimnames=list(NULL, c("ID_old"))),
                                                 by = "ID_old")[,c("ID","time","event","A")])
                  data_matched_bs <- data_matched_bs[order(data_matched_bs$time, data_matched_bs$ID),]
                  # estimate ATE
                  cox_bs <- coxph(Surv(time, event) ~ A, 
                                  data_matched_bs, ties = "breslow")
                  as.numeric(ATE(ID = data_matched_bs$ID,
                                 Z = list(model.matrix(cox_bs)),
                                 index_A = c(which(colnames(model.matrix(cox_bs)) == 
                                                     paste0("A", levels(data_matched_bs$A)[2])), NA)[1],
                                 event = data_matched_bs$event,
                                 time = data_matched_bs$time,
                                 beta = list(cox_bs$coefficients),
                                 t = t,
                                 IF = FALSE, WBS = FALSE)$ATE - PS_B_bs$bias)
                }, error = function(e){rep(NA, length(t))})
              }
  }
  if(.Platform$OS.type == "windows"){
    stopCluster(cl)
    registerDoSEQ()
  }
  return(res_EBS)
}

# function to organize parallelized outcome list
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

# function to run simulations and save summarized results
run <- function(n, 
                alpha_0 = 0, beta_0, gamma_0 = 0, delta_0 = 0, 
                sd_cov = 1, 
                cens = TRUE, cens_par = 200, m = NULL,
                method = "cs",
                t = NULL, 
                ATE_true, 
                run_EBS = TRUE,
                EBS_iter = 1e3, BS_iter = 1e3, 
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
  ## method: string denoting the data generating method 
  ##   (ozenne: latent failure time model, cs: cause-specific Cox model)
  ## t: time point(s) to evaluate
  ## ATE_true: true average treatment effect
  ## run_EBS: logical value indicating whether the naive EBS should be performed
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
    cov_CI_dr <-
    cov_CI_IF <- cov_CI_IF_cluster <- 
    cov_CI_WBS <- cov_CI_WBS_cluster <- 
    matrix(NA, nrow=iter, ncol=length(t))
  cov_CB_EBS <- 
    cov_CB_dr <-
    cov_CB_IF <- cov_CB_IF_cluster <-
    cov_CB_WBS <- cov_CB_WBS_cluster <- 
    rep(NA, iter)
  
  # prepare results list
  res <- list()
  set.seed(seed, kind = "Mersenne-Twister") # for reproducibility
  
  for(i in 1:iter){
    
    res[[i]] <- NA
    
    try({
      
      ## generate data ####
      if(is.null(m)){ # competing risks scenarios
        data <- generate_data(
          n = n, 
          alpha_0 = alpha_0, beta_0 = beta_0, 
          gamma_0 = gamma_0, delta_0 = delta_0, 
          sd_cov = sd_cov, 
          cens = cens, cens_par = cens_par,
          method = method
        )
        # determine transition frequencies at each time point by treatment group
        event_prob <- list(
          A0 = sapply(t, function(t){
            c(sum(data$event == 0 & data$time <= t & data$A == 0),
              sum(data$event == 1 & data$time <= t & data$A == 0),
              sum(data$event == 2 & data$time <= t & data$A == 0))/n
          }),
          A1 = sapply(t, function(t){
            c(sum(data$event == 0 & data$time <= t & data$A == 1),
              sum(data$event == 1 & data$time <= t & data$A == 1),
              sum(data$event == 2 & data$time <= t & data$A == 1))/n
          })
        )
        row.names(event_prob[[1]]) <- row.names(event_prob[[2]]) <- c("0","1","2")
      }else{ # setting with staggered entry & type II censoring
        data <- generate_data(n = n, beta_0 = beta_0, m = m)
        # determine transition frequencies at each time point by treatment group
        event_prob <- list(
          sapply(t, function(t){
            c(sum(data$event == 0 & data$time <= t & data$A == 0),
              sum(data$event == 1 & data$time <= t & data$A == 0))/n
          }),
          sapply(t, function(t){
            c(sum(data$event == 0 & data$time <= t & data$A == 1),
              sum(data$event == 1 & data$time <= t & data$A == 1))/n
          })
        )
        row.names(event_prob[[1]]) <- row.names(event_prob[[2]]) <- c("0","1")
      }
      colnames(event_prob[[1]]) <- colnames(event_prob[[2]]) <- paste0("t = ", t)
      # determine treatment frequency
      treat_prob <- sum(data$A == 1)/n

      # order data by time & assign unique patient identifiers
      data <- data[order(data$time),]
      data$ID <- 1:n
      ## determine propensity scores ####
      data$PS <- predict(glm(A ~ z1+z2+z3+z4+z5+z6+z7+z8+z9+z10+z11+z12, data, 
                             family = "binomial"), 
                         data, type = "response")
      
      ## determine matching bias ####
      bias_time <- proc.time()
      if(is.null(m)){
        PS_B <- PS_bias(ID = data$ID,
                        A = as.numeric(data$A)-1,
                        event = data$event,
                        time = data$time,
                        beta = coefficients(CSC(
                          list(Hist(time, event) ~ A + PS, 
                               Hist(time, event) ~ A + PS), 
                          data)),
                        PS = data$PS,
                        t = t)
      }else{
        PS_B <- PS_bias(ID = data$ID,
                        A = as.numeric(data$A)-1,
                        event = data$event,
                        time = data$time,
                        beta = list(coefficients(coxph(
                          Surv(time, event) ~ A + PS, 
                          data))),
                        PS = data$PS,
                        t = t)
      }
      ## match data ####
      data_matched <- merge(rbind(cbind(data[,c("ID","time","event","A")],
                                        "pair_ID"=1:n),
                                  cbind(merge(merge(data[!duplicated(data$ID),c("ID","time","event","A")], 
                                                    matrix(PS_B$matches, ncol=2, dimnames=list(NULL, c("ID_old", "ID"))), 
                                                    by="ID"),
                                              matrix(data$ID, ncol=1, dimnames=list(NULL, c("ID_old"))),
                                              by = "ID_old"),
                                        "pair_ID"=1:n)[,c("ID","time","event","A","pair_ID")]),
                            matrix(PS_B$n_i, ncol=2, dimnames=list(NULL, c("ID","n_i"))),
                            by="ID")
      data_matched$n_sampled <- data_matched$n_i + 1
      data_matched <- data_matched[order(data_matched$time, data_matched$ID),c("ID","time","event","A","n_sampled","pair_ID")]
      # estimate ATE
      if(is.null(m)){
        cox <- CSC(list(Hist(time, event) ~ A, 
                        Hist(time, event) ~ A), 
                   data_matched, ties = "breslow")
        ATE_m <- as.numeric(ATE(ID = data_matched$ID,
                                Z = list(cox$models$`Cause 1`$x, 
                                         cox$models$`Cause 2`$x),
                                index_A = c(c(which(colnames(cox$models$`Cause 1`$x) == 
                                                      paste0("A", levels(data_matched$A)[2])), NA)[1],
                                            c(which(colnames(cox$models$`Cause 2`$x) == 
                                                      paste0("A", levels(data_matched$A)[2])), NA)[1]),
                                event = data_matched$event,
                                time = data_matched$time,
                                beta = coefficients(cox),
                                t = t,
                                IF = FALSE, WBS = FALSE)$ATE - as.numeric(PS_B$bias))
      }else{
        cox <- coxph(Surv(time, event) ~ A, data_matched, ties = "breslow")
        ATE_m <- as.numeric(ATE(ID = data_matched$ID,
                                Z = list(model.matrix(cox)),
                                index_A = c(which(colnames(model.matrix(cox)) == 
                                                    paste0("A", levels(data_matched$A)[2])), NA)[1],
                                event = data_matched$event,
                                time = data_matched$time,
                                beta = list(cox$coefficients),
                                t = t,
                                IF = FALSE, WBS = FALSE)$ATE - as.numeric(PS_B$bias))
      }
      bias_time <- (proc.time() - bias_time)[3]*1e3
      
      ## perform EBS ####
      CI_EBS <- list(matrix(NA, nrow=2, ncol=length(t)), rep(NA, length(t)))
      CB_EBS <- matrix(NA, nrow=2, ncol=length(t))
      EBS_time <- NA
      if(run_EBS){
        EBS_time <- proc.time()
        res_EBS <- EBS(m=m, t=t, n=n, data=data, EBS_iter=EBS_iter, i=i, seed=seed, cores=cores)
        EBS_invalid <- apply(res_EBS, 1, function(t){sum(is.na(t))})
        res_EBS <- pmin(pmax(res_EBS, -1, na.rm=TRUE), 1, na.rm=TRUE)
        # derive confidence intervals
        CI_EBS <- list(rbind(apply(res_EBS, 1, quantile, 0.025, na.rm=TRUE),
                             apply(res_EBS, 1, quantile, 0.975, na.rm=TRUE)),
                       EBS_invalid)
        # determine quantile for confidence bands
        suppressWarnings(
          q_EBS_sup <- quantile(apply(abs(res_EBS - rowMeans(res_EBS, na.rm = TRUE)) / 
                                        apply(res_EBS, 1, sd, na.rm = TRUE), 
                                      2, FUN = max, na.rm = TRUE), 0.95)
        )
        CB_EBS <- rbind(pmax(ATE_m - q_EBS_sup * apply(res_EBS, 1, sd, na.rm = TRUE), -1), 
                        pmin(ATE_m + q_EBS_sup * apply(res_EBS, 1, sd, na.rm = TRUE), 1))
        EBS_time <- ifelse(any(!is.na(CI_EBS[[1]])), 
                           bias_time + (proc.time() - EBS_time)[3]*1e3, 
                           NA)
      } 


      ## perform double-resampling ####
      # create matrices to store double resampling results
      CI_dr <- CB_dr <- matrix(NA, nrow=2, ncol=length(t))
      # derive double-resampling confidence intervals & bands following Tang et al.
      dr_time <- NA
      try({
        # bootstrap treatment & multipliers
        dr_time <- proc.time()
        if(.Platform$OS.type == "windows"){
          cl <- makeCluster(cores)
          registerDoParallel(cl)
        }else{
          registerDoParallel(cores = cores)
        }
        # generate auxiliary quantities
        dr_quantities <- 
          foreach(j = 1:BS_iter, .combine='comb', .multicombine=TRUE, .init=list(list(),list(),list()),
                  .options.RNG = seed+i) %dorng% {
                    tryCatch({
                      A_dr <- rbinom(n, 1, data$PS)
                      PS_dr <- predict(glm(A_dr ~ z1+z2+z3+z4+z5+z6+z7+z8+z9+z10+z11+z12,
                                           data, family = "binomial"),
                                       data, type = "response")
                      u_dr <- rnorm(n)
                      list(A_dr, PS_dr, u_dr)
                    }, error = function(e){list(rep(NA, n), rep(NA, n), rep(NA, n))})
                  }
        if(.Platform$OS.type == "windows"){
          stopCluster(cl)
          registerDoSEQ()
        }
        # perform double resampling
        if(is.null(m)){
          res_dr <- PS_dr(ID = data$ID,
                          A = as.numeric(data$A)-1,
                          Z = as.matrix(data[,c("z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12")]),
                          event = data$event,
                          time = data$time,
                          beta = unlist(coefficients(cox)),
                          PS = data$PS,
                          n_i = PS_B$n_i[,2],
                          t = t,
                          A_BS = do.call(cbind, dr_quantities[[1]]), 
                          PS_BS = do.call(cbind, dr_quantities[[2]]),
                          u_BS = do.call(cbind, dr_quantities[[3]]))
        }else{
          res_dr <- PS_dr(ID = data$ID,
                          A = as.numeric(data$A)-1,
                          Z = as.matrix(data[,c("z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12")]),,
                          event = data$event,
                          time = data$time,
                          beta = cox$coefficients,
                          PS = data$PS,
                          n_i = PS_B$n_i[,2],
                          t = t,
                          A_BS = do.call(cbind, dr_quantities[[1]]), 
                          PS_BS = do.call(cbind, dr_quantities[[2]]),
                          u_BS = do.call(cbind, dr_quantities[[3]]))
        }
        # derive double resampling confidence intervals & bands
        CI_dr <- rbind(pmax(ATE_m + apply(res_dr, 2, quantile, 0.025, na.rm=TRUE), -1),
                       pmin(ATE_m + apply(res_dr, 2, quantile, 0.975, na.rm=TRUE), 1))
        # determine quantile for confidence bands
        suppressWarnings(
          q_dr_sup <- quantile(apply(abs(t(res_dr) - colMeans(res_dr, na.rm = TRUE)) / 
                                       apply(res_dr, 2, sd, na.rm = TRUE), 
                                     2, FUN = max, na.rm = TRUE), 0.95)
        )
        CB_dr <- rbind(pmax(ATE_m - q_dr_sup * apply(res_dr, 2, sd, na.rm = TRUE), -1), 
                       pmin(ATE_m + q_dr_sup * apply(res_dr, 2, sd, na.rm = TRUE), 1))
        dr_time <- bias_time + (proc.time() - dr_time)[3]*1e3
      }, silent = TRUE)

      
      ## perform (clustered) IF/WBS ####
      # create matrices to store (cluster-based) IF & WBS results
      CI_IF <- CI_IF_cluster <- 
        CI_WBS <- CI_WBS_cluster <- 
        matrix(NA, nrow=2, ncol=length(t))
      CB_IF <- CB_IF_cluster <- 
        CB_WBS <- CB_WBS_cluster <- 
        matrix(NA, nrow=2, ncol=length(t))
      # derive IF & WBS confidence intervals & bands as well as 
      # cluster-based confidence intervals & bands following Austin & Cafri
      IF_time <- NA
      IF_cluster_time <- NA
      WBS_time <- NA
      WBS_cluster_time <- NA  
      # perform IF/WBS
      try({
        
        # determine number of observed (uncensored) events for multipliers
        event_num <- sum(data$event != 0)
        
        if(is.null(m)){ # competing risks scenarios
          
          # calculate average treatment effect, considering unique (matched) observations
          ATE_independent <- ATE(ID = data_matched$ID,
                                 Z = list(cox$models$`Cause 1`$x, 
                                          cox$models$`Cause 2`$x),
                                 index_A = c(c(which(colnames(cox$models$`Cause 1`$x) == 
                                                       paste0("A", levels(data_matched$A)[2])), NA)[1],
                                             c(which(colnames(cox$models$`Cause 2`$x) == 
                                                       paste0("A", levels(data_matched$A)[2])), NA)[1]),
                                 event = data_matched$event,
                                 time = data_matched$time,
                                 beta = coefficients(cox),
                                 t = t,
                                 G_IF_init = rnorm(n * BS_iter), 
                                 G_WBS_init = rnorm(BS_iter * event_num),
                                 bs_iter = BS_iter)

          # derive confidence intervals & bands
          ## IF
          CI_IF <- rbind(pmax(ATE_independent$CI_IF[1,] - as.numeric(PS_B$bias), -1),
                         pmin(ATE_independent$CI_IF[2,] - as.numeric(PS_B$bias), 1))
          CB_IF <- rbind(pmax(ATE_independent$CB_IF[1,] - as.numeric(PS_B$bias), -1),
                         pmin(ATE_independent$CB_IF[2,] - as.numeric(PS_B$bias), 1))
          IF_time <- bias_time +
            sum(clock$timer[clock$ticker %in% c("ATE","IF")])/1e6
          ## WBS
          CI_WBS <- rbind(pmax(ATE_independent$CI_WBS[1,] - as.numeric(PS_B$bias), -1),
                          pmin(ATE_independent$CI_WBS[2,] - as.numeric(PS_B$bias), 1))
          CB_WBS <- rbind(pmax(ATE_independent$CB_WBS[1,] - as.numeric(PS_B$bias), -1),
                          pmin(ATE_independent$CB_WBS[2,] - as.numeric(PS_B$bias), 1))
          WBS_time <- bias_time +
          sum(clock$timer[clock$ticker %in% c("ATE","WBS")])/1e6
          
          ## clustered IF/WBS
          # save time for clustered CIs/CBs
          IF_cluster_time <-  
            bias_time + sum(clock$timer[clock$ticker %in% c("ATE","IF")])/1e6 - 
            clock$timer[clock$ticker == "IF_CI_CB"]/1e6
          WBS_cluster_time <- 
            bias_time + sum(clock$timer[clock$ticker %in% c("ATE","WBS")])/1e6 - 
            clock$timer[clock$ticker == "WBS_CI_CB"]/1e6
          # cluster by individual subjects
          beta_time <- proc.time()
          beta_cluster <- CSC(list(Hist(time, event) ~ A,
                                   Hist(time, event) ~ A),
                              data_matched, ties = "breslow", weight = 1/data_matched$n_sampled)
          beta_time <- (proc.time() - beta_time)[3]*1e3
          ATE_subjects <- ATE(ID = data_matched$ID,
                              Z = list(beta_cluster$models$`Cause 1`$x, 
                                       beta_cluster$models$`Cause 2`$x),
                              index_A = c(c(which(colnames(beta_cluster$models$`Cause 1`$x) == 
                                                    paste0("A", levels(data_matched$A)[2])), NA)[1],
                                          c(which(colnames(beta_cluster$models$`Cause 2`$x) == 
                                                    paste0("A", levels(data_matched$A)[2])), NA)[1]),
                              event = data_matched$event,
                              time = data_matched$time,
                              beta = coefficients(beta_cluster),
                              t = t,
                              cluster_init = data_matched$ID,
                              G_IF_init = rnorm(n * BS_iter), 
                              G_WBS_init = rnorm(BS_iter * event_num),
                              bs_iter = BS_iter)
          IF_cluster_time <- IF_cluster_time +
            beta_time + sum(clock$timer[clock$ticker %in% c("ATE","IF")])/1e6
          WBS_cluster_time <- WBS_cluster_time + 
            beta_time + sum(clock$timer[clock$ticker %in% c("ATE","WBS")])/1e6
          # cluster by matched pairs
          beta_time <- proc.time()
          beta_cluster <- CSC(list(Hist(time, event) ~ A,
                                   Hist(time, event) ~ A),
                              data_matched, ties = "breslow", weight = rep(1/2, 2*n))
          beta_time <- (proc.time() - beta_time)[3]*1e3
          ATE_pairs <- ATE(ID = data_matched$ID,
                           Z = list(beta_cluster$models$`Cause 1`$x, 
                                    beta_cluster$models$`Cause 2`$x),
                           index_A = c(c(which(colnames(beta_cluster$models$`Cause 1`$x) == 
                                                 paste0("A", levels(data_matched$A)[2])), NA)[1],
                                       c(which(colnames(beta_cluster$models$`Cause 2`$x) == 
                                                 paste0("A", levels(data_matched$A)[2])), NA)[1]),
                           event = data_matched$event,
                           time = data_matched$time,
                           beta = coefficients(beta_cluster),
                           t = t,
                           cluster_init = data_matched$pair_ID,
                           G_IF_init = rnorm(n * BS_iter), 
                           G_WBS_init = rnorm(BS_iter * event_num),
                           bs_iter = BS_iter)
          IF_cluster_time <- IF_cluster_time + 
            beta_time + sum(clock$timer[clock$ticker %in% c("ATE","IF")])/1e6
          WBS_cluster_time <- WBS_cluster_time + 
            beta_time + sum(clock$timer[clock$ticker %in% c("ATE","WBS")])/1e6          
          
        }else{ # setting with staggered entry & type II censoring

          # calculate average treatment effect, considering unique (matched) observations
          ATE_independent <- ATE(ID = data_matched$ID,
                                 Z = list(model.matrix(cox)),
                                 index_A = c(which(colnames(model.matrix(cox)) == 
                                                     paste0("A", levels(data_matched$A)[2])), NA)[1],
                                 event = data_matched$event,
                                 time = data_matched$time,
                                 beta = list(cox$coefficients),
                                 t = t,
                                 G_IF_init = rnorm(n * BS_iter), 
                                 G_WBS_init = rnorm(BS_iter * event_num),
                                 bs_iter = BS_iter)

          # derive confidence intervals & bands
          ## IF
          CI_IF <- rbind(pmax(ATE_independent$CI_IF[1,] - as.numeric(PS_B$bias), -1),
                         pmin(ATE_independent$CI_IF[2,] - as.numeric(PS_B$bias), 1))
          CB_IF <- rbind(pmax(ATE_independent$CB_IF[1,] - as.numeric(PS_B$bias), -1),
                         pmin(ATE_independent$CB_IF[2,] - as.numeric(PS_B$bias), 1))
          IF_time <- bias_time +
            sum(clock$timer[clock$ticker %in% c("ATE","IF")])/1e6
          ## WBS
          CI_WBS <- rbind(pmax(ATE_independent$CI_WBS[1,] - as.numeric(PS_B$bias), -1),
                          pmin(ATE_independent$CI_WBS[2,] - as.numeric(PS_B$bias), 1))
          CB_WBS <- rbind(pmax(ATE_independent$CB_WBS[1,] - as.numeric(PS_B$bias), -1),
                          pmin(ATE_independent$CB_WBS[2,] - as.numeric(PS_B$bias), 1))
          WBS_time <- bias_time +
            sum(clock$timer[clock$ticker %in% c("ATE","WBS")])/1e6

          ## clustered IF/WBS
          # save time for clustered CIs/CBs
          IF_cluster_time <- 
            bias_time + sum(clock$timer[clock$ticker %in% c("ATE","IF")])/1e6 - 
            clock$timer[clock$ticker == "IF_CI_CB"]/1e6
          WBS_cluster_time <- 
            bias_time + sum(clock$timer[clock$ticker %in% c("ATE","WBS")])/1e6 - 
            clock$timer[clock$ticker == "WBS_CI_CB"]/1e6
          # cluster by individual subjects
          beta_time <- proc.time()
          beta_cluster <- coxph(Surv(time, event) ~ A, data = data_matched, 
                                ties = "breslow", weight = 1/data_matched$n_sampled)
          beta_time <- (proc.time() - beta_time)[3]*1e3
          ATE_subjects <- ATE(ID = data_matched$ID,
                              Z = list(model.matrix(beta_cluster)),
                              index_A = c(which(colnames(model.matrix(beta_cluster)) == 
                                                  paste0("A", levels(data_matched$A)[2])), NA)[1],
                              event = data_matched$event,
                              time = data_matched$time,
                              beta = list(beta_cluster$coefficients),
                              t = t,
                              cluster_init = data_matched$ID,
                              G_IF_init = rnorm(n * BS_iter), 
                              G_WBS_init = rnorm(BS_iter * event_num),
                              bs_iter = BS_iter)
          IF_cluster_time <- IF_cluster_time + 
            beta_time + sum(clock$timer[clock$ticker %in% c("ATE","IF")])/1e6
          WBS_cluster_time <- WBS_cluster_time + 
            beta_time + sum(clock$timer[clock$ticker %in% c("ATE","WBS")])/1e6
          # cluster by matched pairs
          beta_time <- proc.time()
          beta_cluster <- coxph(Surv(time, event) ~ A, data = data_matched, 
                                ties = "breslow", weight = rep(1/2, 2*n))
          beta_time <- (proc.time() - beta_time)[3]*1e3
          ATE_pairs <- ATE(ID = data_matched$ID,
                           Z = list(model.matrix(beta_cluster)),
                           index_A = c(which(colnames(model.matrix(beta_cluster)) == 
                                               paste0("A", levels(data_matched$A)[2])), NA)[1],
                           event = data_matched$event,
                           time = data_matched$time,
                           beta = list(beta_cluster$coefficients),
                           t = t,
                           cluster_init = data_matched$pair_ID,
                           G_IF_init = rnorm(n * BS_iter), 
                           G_WBS_init = rnorm(BS_iter * event_num),
                           bs_iter = BS_iter)
          IF_cluster_time <- IF_cluster_time + 
            beta_time + sum(clock$timer[clock$ticker %in% c("ATE","WBS")])/1e6
          WBS_cluster_time <- WBS_cluster_time + 
            beta_time + sum(clock$timer[clock$ticker %in% c("ATE","WBS")])/1e6
        }

        IF_cluster_time_calc <- proc.time()
        # calculate combined standard error
        SE_IF_cluster <- sqrt((dim(data_matched)[1]^2 / 
                                 (dim(data_matched)[1] - 1) * 
                                 ATE_independent$SE_IF^2 +
                                 length(unique(data_matched$ID))^2 / 
                                 (length(unique(data_matched$ID)) - 1) * 
                                 ATE_subjects$SE_IF^2 + 
                                 (dim(data_matched)[1]/2)^2 / 
                                 (dim(data_matched)[1]/2 - 1) * 
                                 ATE_pairs$SE_IF^2) / 
                                dim(data_matched)[1])
        # derive confidence intervals & bands
        CI_IF_cluster <- t(cbind(
          pmax(ATE_m - SE_IF_cluster * qnorm(1-(1-0.95)/2), -1),
          pmin(ATE_m + SE_IF_cluster * qnorm(1-(1-0.95)/2), 1)))
        # determine quantile for confidence bands
        q_IF_cluster <- quantile(apply(
          dim(data_matched)[1]/(dim(data_matched)[1]-1) * ATE_independent$Un_std_IF + 
            length(unique(data_matched$ID))/(length(unique(data_matched$ID))-1) * ATE_subjects$Un_std_IF + 
            (dim(data_matched)[1]/2)/(dim(data_matched)[1]/2-1) * ATE_pairs$Un_std_IF, 
          1, max), 0.95)
        CB_IF_cluster <- t(cbind(
          pmax(ATE_m - SE_IF_cluster * q_IF_cluster, -1), 
          pmin(ATE_m + SE_IF_cluster * q_IF_cluster, 1)))
        IF_cluster_time <- IF_cluster_time + 
          (proc.time() - IF_cluster_time_calc)[3]*1e3
        WBS_cluster_time_calc <- proc.time()
        # calculate combined standard error
        SE_WBS_cluster <- sqrt((dim(data_matched)[1]^2 / 
                                  (dim(data_matched)[1] - 1) * 
                                  ATE_independent$SE_WBS^2 +
                                  length(unique(data_matched$ID))^2 / 
                                  (length(unique(data_matched$ID)) - 1) * 
                                  ATE_subjects$SE_WBS^2 + 
                                  (dim(data_matched)[1]/2)^2 / 
                                  (dim(data_matched)[1]/2 - 1) * 
                                  ATE_pairs$SE_WBS^2) / 
                                 dim(data_matched)[1])
        # derive confidence intervals & bands
        CI_WBS_cluster <- t(cbind(
          pmax(ATE_m - SE_WBS_cluster * qnorm(1-(1-0.95)/2), -1),
          pmin(ATE_m + SE_WBS_cluster * qnorm(1-(1-0.95)/2), 1)))
        # determine quantile for confidence bands
        q_WBS_cluster <- quantile(apply(
          dim(data_matched)[1]/(dim(data_matched)[1]-1) * ATE_independent$Un_std_WBS + 
            length(unique(data_matched$ID))/(length(unique(data_matched$ID))-1) * ATE_subjects$Un_std_WBS + 
            (dim(data_matched)[1]/2)/(dim(data_matched)[1]/2-1) * ATE_pairs$Un_std_WBS, 
          1, max), 0.95)
        CB_WBS_cluster <- t(cbind(
          pmax(ATE_m - SE_WBS_cluster * q_WBS_cluster, -1), 
          pmin(ATE_m + SE_WBS_cluster * q_WBS_cluster, 1)))
        WBS_cluster_time <- WBS_cluster_time + 
          (proc.time() - WBS_cluster_time_calc)[3]*1e3
        
      }, silent = TRUE)
      
      ## store results ####
      # store confidence intervals
      CI <- rbind(CI_EBS[[1]],
                  CI_dr,
                  CI_IF,
                  CI_IF_cluster,
                  CI_WBS,
                  CI_WBS_cluster)
      row.names(CI) <- c(
        "EBS (lower)", "EBS (upper)",
        "double resampling (lower)", "double resampling (upper)",
        "IF (lower)", "IF (upper)",
        "clustered IF (lower)", "clustered IF (upper)",
        "WBS (lower)", "WBS (upper)",
        "clustered WBS (lower)", "clustered WBS (upper)"
      )
      # store confidence bands
      CB <- rbind(CB_EBS,
                  CB_dr,
                  CB_IF,
                  CB_IF_cluster,
                  CB_WBS,
                  CB_WBS_cluster)
      row.names(CB) <- c(
        "EBS (lower)", "EBS (upper)",
        "double resampling (lower)", "double resampling (upper)",
        "IF (lower)", "IF (upper)",
        "clustered IF (lower)", "clustered IF (upper)",
        "WBS (lower)", "WBS (upper)",
        "clustered WBS (lower)", "clustered WBS (upper)"
      )
      colnames(CI) <- colnames(CB) <- paste0("t = ", t)
      
      # store results
      res[[i]] <- list(CI = CI,
                       CB = CB,
                       event_distribution = event_prob,
                       treatment_probability = treat_prob,
                       EBS_invalid = setNames(CI_EBS[[2]], paste0("t = ", t)),
                       times = c("EBS" = EBS_time,
                                 "double resampling" = dr_time,
                                 "IF" = IF_time,
                                 "clustered IF" = IF_cluster_time,
                                 "WBS" = WBS_time,
                                 "clustered WBS" = WBS_cluster_time))
      
      # store coverage probabilities
      cov_CI_EBS[i,] <- CI[1,] <= ATE_true & ATE_true <= CI[2,]
      cov_CI_dr[i,] <- CI[3,] <= ATE_true & ATE_true <= CI[4,]
      cov_CI_IF[i,] <- CI[5,] <= ATE_true & ATE_true <= CI[6,]
      cov_CI_IF_cluster[i,] <- CI[7,] <= ATE_true & ATE_true <= CI[8,]
      cov_CI_WBS[i,] <- CI[9,] <= ATE_true & ATE_true <= CI[10,]
      cov_CI_WBS_cluster[i,] <- CI[11,] <= ATE_true & ATE_true <= CI[12,]
      cov_CB_EBS[i] <- all(CB[1,] <= ATE_true & ATE_true <= CB[2,])
      cov_CB_dr[i] <- all(CB[3,] <= ATE_true & ATE_true <= CB[4,])
      cov_CB_IF[i] <- all(CB[5,] <= ATE_true & ATE_true <= CB[6,])
      cov_CB_IF_cluster[i] <- all(CB[7,] <= ATE_true & ATE_true <= CB[8,])
      cov_CB_WBS[i] <- all(CB[9,] <= ATE_true & ATE_true <= CB[10,])
      cov_CB_WBS_cluster[i] <- all(CB[11,] <= ATE_true & ATE_true <= CB[12,])
      
    })
    
    if(i%%(iter/20) == 0){print(paste0(100* i/iter, "%"))}
    
  }
  
  ## summarize & save results ####
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
    #   number of invalid EBS samples, and computation times for each iteration
    results = res,
    # number of iterations with errors
    errors = sum(is.na(res)),
    # average event time distribution
    mean_event_distribution = list(
      A0 = apply(sapply(res, function(i){
        if(length(i) > 1){
          i[[3]][[1]]
        }else{
          matrix(NA, nrow=ifelse(is.null(m),3,2),ncol=length(t))
        }
      }, simplify = "array"), 
      c(1,2), mean, na.rm = TRUE),
      A1 = apply(sapply(res, function(i){
        if(length(i) > 1){
          i[[3]][[1]]
        }else{
          matrix(NA, nrow=ifelse(is.null(m),3,2),ncol=length(t))
        }
      }, simplify = "array"), 
      c(1,2), mean, na.rm = TRUE)
    ),
    # average treatment probability
    mean_treatment_probability = 
      mean(sapply(res, 
                  function(i){
                    if(length(i) > 1){i[[4]]}else{NA}
                  }), na.rm=TRUE),
    # average number of invalid EBS samples
    mean_EBS_samples_invalid = 
      setNames(rowMeans(sapply(res, 
                               function(i){
                                 if(length(i) > 1){i[[5]]}else{rep(NA, length(t))}
                               }), na.rm = TRUE), 
               paste0(c(paste0("t = ", t)))),
    # number of invalid confidence regions
    CIs_invalid = matrix(rowSums(sapply(res, 
                                        function(i){
                                          if(length(i) > 1){
                                            is.na(i[[1]][seq(1,12,2),])
                                          }else{
                                            rep(TRUE, 6*length(t))
                                          }
                                        })), nrow=6,
                         dimnames=list(c("EBS",
                                         "double resampling",
                                         "IF",
                                         "clustered IF",
                                         "WBS",
                                         "clustered WBS"), 
                                       paste0("t = ", t))),
    # average computation times
    mean_time = 
      rowMeans(sapply(res, 
                      function(i){
                        if(length(i) > 1){i[[6]]}else{rep(NA, 6)}
                      }), na.rm = TRUE),
    # average confidence interval widths
    mean_width_CI = 
      matrix(rowMeans(sapply(res, 
                             function(i){
                               if(length(i) > 1){
                                 i[[1]][seq(2,12,2),] - i[[1]][seq(1,11,2),]
                               }else{
                                 rep(NA, 6*length(t))
                               }
                             }), na.rm = TRUE), 
             nrow=5, dimnames=list(c("EBS",
                                     "double resampling",
                                     "IF",
                                     "clustered IF",
                                     "WBS",
                                     "clustered WBS"), 
                                   paste0("t = ", t))),
    # average confidence band widths
    mean_width_CB = 
      matrix(rowMeans(sapply(res, 
                             function(i){
                               if(length(i) > 1){
                                 i[[2]][seq(2,12,2),] - i[[2]][seq(1,11,2),]
                               }else{
                                 rep(NA, 6*length(t))
                               }
                             }), na.rm = TRUE), 
             nrow=5, dimnames=list(c("EBS",
                                     "double resampling",
                                     "IF",
                                     "clustered IF",
                                     "WBS",
                                     "clustered WBS"), 
                                   paste0("t = ", t))),
    # confidence interval coverage probabilities
    coverage_CI = matrix(100*c(colMeans(cov_CI_EBS, na.rm=TRUE),
                               colMeans(cov_CI_dr, na.rm = TRUE), 
                               colMeans(cov_CI_IF, na.rm = TRUE),
                               colMeans(cov_CI_IF_cluster, na.rm = TRUE),
                               colMeans(cov_CI_WBS, na.rm = TRUE),
                               colMeans(cov_CI_WBS_cluster, na.rm = TRUE)), 
                         byrow=TRUE, nrow=6, 
                         dimnames=list(c("EBS",
                                         "double resampling",
                                         "IF",
                                         "clustered IF",
                                         "WBS",
                                         "clustered WBS"), 
                                       paste0("t = ", t))),
    # confidence band coverage probabilities
    coverage_CB = 100*c("EBS" = mean(cov_CB_EBS, na.rm=TRUE),
                        "double resampling" = mean(cov_CB_dr, na.rm = TRUE), 
                        "IF" = mean(cov_CB_IF, na.rm = TRUE),
                        "clustered IF" = mean(cov_CB_IF_cluster, na.rm = TRUE),
                        "WBS" = mean(cov_CB_WBS, na.rm = TRUE),
                        "clustered WBS" = mean(cov_CB_WBS_cluster, na.rm = TRUE))
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
      total_coverage_CI$t == t,], 
    aes(x = n, y = coverage)) +
    geom_hline(aes(yintercept = 95)) +
    geom_line(aes(color = type), linewidth=1.1) +
    geom_point(aes(color = type, fill = type, shape = type), size=2) +
    scale_color_manual("", 
                       values = c('orange','red','turquoise4','cyan3','chartreuse2')) +
    scale_fill_manual("", 
                      values = c('orange','red','turquoise4','cyan3','chartreuse2')) +
    scale_shape_manual("",
                       values = c(19,15,17,25,23)) +
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
                              total_coverage_CI$effect == effect,
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
                              total_coverage_CB$effect == effect,
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
      total_coverage_CB$effect == effect,
  ], 
  aes(x = n, y = coverage)) +
    geom_hline(aes(yintercept = 95)) +
    geom_line(aes(color = type), linewidth=1.1) +
    geom_point(aes(color = type, fill = type, shape = type), size=2) +
    scale_color_manual("", 
                       values = c('orange','red','turquoise4','cyan3','chartreuse2')) +
    scale_fill_manual("", 
                      values = c('orange','red','turquoise4','cyan3','chartreuse2')) +
    scale_shape_manual("",
                       values = c(19,15,17,25,23)) +
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
              i[[1]][seq(2,10,2),
                     which(as.character(t) == all_t)] - 
                i[[1]][seq(1,9,2),
                       which(as.character(t) == all_t)]
            }else{
              rep(NA, 5)
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
              i[[2]][seq(2,10,2),
                     which(as.character(t) == all_t)] - 
                i[[2]][seq(1,9,2),
                       which(as.character(t) == all_t)]
            }else{
              rep(NA, 5)
            }
          }
          )))
        )
      )
    }
  }
  # identify mild & extreme outliers
  mean_val <- aggregate(width ~ n + type, data, mean)
  q1 <- aggregate(width ~ n + type, data, quantile, 0.25)
  q3 <- aggregate(width ~ n + type, data, quantile, 0.75)
  sum <- merge(merge(mean_val, q1, by = c("n","type")), q3, by = c("n","type"))
  colnames(sum)[3:5] <- c("mean","q1","q3")
  sum$whisker_low <- sum$q1 - 1.5 * (sum$q3 - sum$q1)
  sum$whisker_up <- sum$q3 + 1.5 * (sum$q3 - sum$q1)
  
  outliers <- merge(data, sum[,c("n","type","q1","q3","whisker_low","whisker_up")], 
                    by = c("n","type"))
  outliers <- outliers[(outliers$width < outliers$whisker_low | 
                          outliers$width > outliers$whisker_up) & 
                         !is.na(outliers$width),]
  outliers$extreme_outlier_low <- outliers$q1 - 3 * (outliers$q3 - outliers$q1)
  outliers$extreme_outlier_up <- outliers$q3 + 3 * (outliers$q3 - outliers$q1)
  outliers$outlier <- ifelse(outliers$width < outliers$whisker_low | 
                               outliers$width > outliers$whisker_up, 
                             "mild", NA)
  outliers$outlier <- ifelse(outliers$width < outliers$extreme_outlier_low | 
                               outliers$width > outliers$extreme_outlier_up,
                             "extreme", outliers$outlier)
  outliers$outlier <- factor(outliers$outlier, levels = c("mild","extreme"))
  
  # create plot
  ggplot(data = sum, aes(x = factor(n), fill = type)) + 
    geom_errorbar(aes(ymin = whisker_low, ymax = whisker_up), position = "dodge") +
    geom_boxplot(aes(ymin = whisker_low, lower = q1, middle = mean, upper = q3, ymax = whisker_up),
                 stat = "identity", position = "dodge", outlier.shape = NA) + 
    geom_point(data = outliers, aes(x = factor(n), y = width, col = outlier, group = type), 
               size = 0.5, position = position_dodge(width=0.9), show.legend = FALSE) + 
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
  # identify mild & extreme outliers
  mean_val <- aggregate(time ~ n + type, data, mean)
  q1 <- aggregate(time ~ n + type, data, quantile, 0.25)
  q3 <- aggregate(time ~ n + type, data, quantile, 0.75)
  sum <- merge(merge(mean_val, q1, by = c("n","type")), q3, by = c("n","type"))
  colnames(sum)[3:5] <- c("mean","q1","q3")
  sum$whisker_low <- sum$q1 - 1.5 * (sum$q3 - sum$q1)
  sum$whisker_up <- sum$q3 + 1.5 * (sum$q3 - sum$q1)
  
  outliers <- merge(data, sum[,c("n","type","q1","q3","whisker_low","whisker_up")], 
                    by = c("n","type"))
  outliers <- outliers[(outliers$time < outliers$whisker_low | 
                          outliers$time > outliers$whisker_up) & 
                         !is.na(outliers$time),]
  outliers$extreme_outlier_low <- outliers$q1 - 3 * (outliers$q3 - outliers$q1)
  outliers$extreme_outlier_up <- outliers$q3 + 3 * (outliers$q3 - outliers$q1)
  outliers$outlier <- ifelse(outliers$time < outliers$whisker_low | 
                               outliers$time > outliers$whisker_up, 
                             "mild", NA)
  outliers$outlier <- ifelse(outliers$time < outliers$extreme_outlier_low | 
                               outliers$time > outliers$extreme_outlier_up,
                             "extreme", outliers$outlier)
  outliers$outlier <- factor(outliers$outlier, levels = c("mild","extreme"))
  
  # create plot
  ggplot(data = sum, aes(x = factor(n), fill = type)) + 
    geom_bar(aes(y = mean), stat = "identity", position = position_dodge(width=0.9)) + 
    geom_errorbar(aes(ymin = whisker_low, ymax = whisker_up), 
                  width = 0.5, position = position_dodge(width=0.9)) +
    geom_boxplot(aes(ymin = whisker_low, lower = q1, middle = mean, upper = q3, ymax = whisker_up),
                 stat = "identity", width = 0.5, position = position_dodge(width=0.9), 
                 outlier.shape = NA, show.legend = FALSE) + 
    geom_point(data = outliers, aes(x = factor(n), y = time, col = outlier, group = type), 
               size = 1, position = position_dodge(width=0.9), show.legend = FALSE) + 
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
          plot.title = element_text(hjust=0.5),
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

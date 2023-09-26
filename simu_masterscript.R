# This script reproduces the results of the simulation study.

# load functions
source("simu_functions.R")

# load packages
library(riskRegression)
library(prodlim)
library(survival)
library(doParallel)
library(rngtools)
library(doRNG)
library(ATESurvival)
library(ggplot2)
library(patchwork)


# find true average treatment effect ###########################################

# to save time, skip this section (l. 23-89) and use
# load("Results/ATE_true.Rda")

# create table of scenarios
ATE_true <- data.frame(m = c(rep("NULL", 9), rep("n/2", 3)),
                       sd_cov = rep(c(1, 0.5, 2, 1), each = 3), 
                       beta_0 = rep(c(-2, 0, 2), 4),
                       ATE_true = NA)
# create matrix to store average treatment effects for Figure 1
ATE_true_fig1 <- matrix(nrow = 2000, ncol = 101)

# determine true average treatment effect
set.seed(1234, kind = "Mersenne-Twister") # for reproducibility
for(i in 1:dim(ATE_true)[1]){
  if(ATE_true$beta_0[i] == 0){ # no treatment effect
    ATE_true$ATE_true[i] <- "0, 0, 0, 0, 0"
    next
  }
  m <- switch(ATE_true$m[i], 
              "NULL" = NULL,
              "n/2" = 50000)
  print(paste0("m = ", m, ", sd_cov = ", ATE_true$sd_cov[i], ", beta_0 = ", 
               ATE_true$beta_0[i], ":"))
  # create matrix to store average treatment effects
  ATEs <- matrix(nrow = 1000, ncol = 5)
  
  for(j in 1:1000){
    # generate data
    data <- generate_data(n = 100000, 
                          beta_0 = ATE_true$beta_0[i], 
                          sd_cov = ATE_true$sd_cov[i],
                          true_ATE = TRUE, m = m)
    if(is.null(m)){ # competing risks scenarios
      # calculate average treatment effect
      ATEs[j,] <- diff(predictRisk(prodlim(Hist(time, event) ~ A, data=data), 
                                   newdata=data.frame(A=0:1), times=c(1,3,5,7,9),
                                   cause=1))
      if(ATE_true$sd_cov[i] == 1){
        # calculate average treatment effect for Figure 1
        ATE_true_fig1[ifelse(ATE_true$beta_0[i] == -2, 0, 1000) + j,] <- 
          diff(predictRisk(prodlim(Hist(time, event) ~ A, data=data), 
                           newdata=data.frame(A=0:1), times=seq(0,10,0.1),
                           cause=1))
      }
    }else{ # setting with staggered entry & type II censoring
      if(ATE_true$beta_0[i] == -2){
        t <- c(2,4,6,8,10)
      }else if(ATE_true$beta_0[i] == 2){
        t <- c(0.5,1,1.5,2,2.5)
      }
      # calculate average treatment effect
      ATEs[j,] <- diff(predictRisk(survfit(Surv(time, event) ~ A, data=data), 
                                   newdata=data.frame(A=0:1), times=t))
    }
    if(j %% 100 == 0) print(paste0(j/10, "%"))
  }
  # compute median average treatment effect
  ATE_true$ATE_true[i] <- paste0(round(apply(ATEs, 2, median, na.rm = TRUE), 3), 
                                 collapse = ", ")
}

# save results
ATE_true_fig1 <- data.frame(t = rep(seq(0,10,0.1), 3),
                            ATE_true = c(apply(ATE_true_fig1[1:1000,], 2, median, 
                                               na.rm = TRUE),
                                         rep(0, 101), 
                                         apply(ATE_true_fig1[1001:2000,], 2, median, 
                                               na.rm = TRUE)),
                            beta_0 = c(rep("-2", 101), rep("0", 101), rep("2", 101)))
save(ATE_true, ATE_true_fig1, file = "Results/ATE_true.Rda")


# run simulations ##############################################################

# create table of scenarios
scenarios <- data.frame(order = 1:24,
                        scenario = rep(c("noCens", "lowCens", "highCens", 
                                         "lowTreatProb", "highTreatProb", 
                                         "lowVarCov", "highVarCov",
                                         "typeII"), 
                                       each = 3),
                        alpha_0 = c(rep(0, 3*3), rep(-2, 3), rep(2, 3), rep(0, 3*3)),
                        beta_0 = rep(c(-2, 0, 2), 8),
                        sd_cov = c(rep(1, 3*5), rep(0.5, 3), rep(2, 3), rep(1, 3*1)),
                        cens = c(rep(FALSE, 3), rep(TRUE, 3*7)),
                        cens_par = c(rep(200, 3*2), rep(50, 3), rep(200, 3*5)),
                        m = c(rep("NULL", 3*7), rep("n/2", 3)))
scenarios <- merge(scenarios, ATE_true, by = c("beta_0", "sd_cov", "m"))
scenarios <- scenarios[order(scenarios$order), 
                       c("scenario","alpha_0","beta_0","sd_cov","cens","cens_par","m","ATE_true")]
rownames(scenarios) <- NULL

for(i in 1:dim(scenarios)[1]){
  print(paste0("Scenario: ", 
               switch(scenarios$scenario[i],
                      "noCens" = "no censoring",
                      "lowCens" = "low censoring", 
                      "highCens" = "high censoring", 
                      "lowTreatProb" = "low treatment probability",
                      "highTreatProb" = "high treatment probability",
                      "lowVarCov" = "low variance of the covariates",
                      "highVarCov" = "high variance of the covariates",
                      "typeII" = "type II censoring with staggered entry")))
  for(n in c(50, 75, 100, 200, 300)){
    print(paste0("beta_0 = ", scenarios$beta_0[i], ", n = ", n, ":"))
    name <- paste0("res_",
                   ifelse(scenarios$beta_0[i] == 0, "no", 
                          ifelse(scenarios$beta_0[i] == 2, "", 
                                 ifelse(scenarios$beta_0[i] == -2, "adv", NA))), 
                   "ATE_", scenarios$scenario[i], "_n", n)
    # run simulation
    assign(name, 
           run(n = n, 
               alpha_0 = scenarios$alpha_0[i], 
               beta_0 = scenarios$beta_0[i], 
               sd_cov = scenarios$sd_cov[i], 
               cens = scenarios$cens[i], 
               cens_par = scenarios$cens_par[i],
               m = switch(scenarios$m[i], "NULL" = NULL, "n/2" = n/2),
               ATE_true = eval(parse(text = paste0("c(", scenarios$ATE_true[i], ")"))),
               iter = 5000))
    # save results
    save(eval(parse(text = name)), file = paste0("Results/", name, ".Rda"))
  }
}


# prepare outcomes #############################################################

# to save time, skip this section (l. 158-252) and use
# load("Results/total_coverages.Rda")

# summarize data for plots
total_coverage_CI <- data.frame()
total_coverage_CB <- data.frame()

for(i in 1:dim(scenarios)[1]){
  
  # determine time points considered in the respective scenario
  if(scenarios$scenario[i] != "typeII"){
    t <- c(1,3,5,7,9)
  }else if(scenarios$scenario[i] == "typeII"){
    if(scenarios$beta_0[i] == -2){
      t <- c(2,4,6,8,10)
    }else if(scenarios$beta_0[i] == 0){
      t <- c(1,2,3,4,5)
    }else if(scenarios$beta_0[i] == 2){
      t <- c(0.5,1,1.5,2,2.5)
    }
  }
  
  for(n in c(50, 75, 100, 200, 300)){
    
    # retrieve data on CI coverage
    total_coverage_CI <- rbind(
      total_coverage_CI, 
      data.frame(
        scenario = scenarios$scenario[i],
        effect = ifelse(scenarios$beta_0[i] == 0, "no", 
                        ifelse(scenarios$beta_0[i] == 2, "yes", 
                               ifelse(scenarios$beta_0[i] == -2, "adverse", NA))),
        n = n,
        t = rep(t, each=6),
        type = factor(rep(c("EBS",
                            "IF",
                            "WBS (calculated)",
                            "WBS - Lin et al.",
                            "WBS - Beyersmann et al.",
                            "WBS - Weird bootstrap"), 
                          5),
                      levels = c("EBS",
                                 "IF",
                                 "WBS (calculated)",
                                 "WBS - Lin et al.",
                                 "WBS - Beyersmann et al.",
                                 "WBS - Weird bootstrap")),
        coverage = c(eval(parse(
          text = paste0("res_", 
                        ifelse(scenarios$beta_0[i] == 0, "no", 
                               ifelse(scenarios$beta_0[i] == 2, "", 
                                      ifelse(scenarios$beta_0[i] == -2, "adv", NA))), 
                        "ATE_", scenarios$scenario[i], "_n", n, "[[15]]")
        )))
      )
    )
    
    # retrieve data on CB coverage
    total_coverage_CB <- rbind(
      total_coverage_CB, 
      data.frame(
        scenario = scenarios$scenario[i],
        effect = ifelse(scenarios$beta_0[i] == 0, "no", 
                        ifelse(scenarios$beta_0[i] == 2, "yes", 
                               ifelse(scenarios$beta_0[i] == -2, "adverse", NA))),
        n = n,
        type = factor(c("EBS", 
                        "IF",
                        "WBS - Lin et al. (calculated)",
                        "WBS - Lin et al.",
                        "WBS - Beyersmann et al. (calculated)",
                        "WBS - Beyersmann et al.",
                        "WBS - Weird bootstrap (calculated)",
                        "WBS - Weird bootstrap"),
                      levels = c("EBS",
                                 "IF",
                                 "WBS - Lin et al. (calculated)",
                                 "WBS - Lin et al.",
                                 "WBS - Beyersmann et al. (calculated)",
                                 "WBS - Beyersmann et al.",
                                 "WBS - Weird bootstrap (calculated)",
                                 "WBS - Weird bootstrap")),
        coverage = eval(parse(
          text = paste0("res_", 
                        ifelse(scenarios$beta_0[i] == 0, "no", 
                               ifelse(scenarios$beta_0[i] == 2, "", 
                                      ifelse(scenarios$beta_0[i] == -2, "adv", NA))), 
                        "ATE_", scenarios$scenario[i], "_n", n, "[[16]]")
        ))
      )
    )
    
  }
  
}
rownames(total_coverage_CB) <- NULL
# save results
save(total_coverage_CI, total_coverage_CB, file = "Results/total_coverages.Rda")


# create plots #################################################################

## Figure 1 ####
png("Results/figures/figure_1.png", width = 800, height = 500)
ggplot(data = ATE_true_fig1, aes(x = t, y = ATE_true)) + 
  geom_line(aes(linetype = beta_0), linewidth = 0.85) +
  geom_point(data = ATE_true_fig1[c(rep(0, 5), rep(101, 5), rep(2 * 101, 5)) + 
                                    rep(findInterval(c(1,3,5,7,9), seq(0,10,0.1)), 3),], 
             shape = 19, size = 1.7) + 
  scale_linetype_manual(name = "", values = c("solid","88","18"), 
                        labels = c(expression(paste(beta[0], " = -2")),
                                   expression(paste(beta[0], " = 0")),
                                   expression(paste(beta[0], " = 2")))) + 
  labs(x = "t", y = "ATE(t)") +
  scale_x_continuous(expand = expansion(mult = c(0.03,0)), limits = c(0, 10), 
                     breaks = c(1,3,5,7,9)) +
  theme(text = element_text(size = 15), 
        axis.title.x = element_text(margin = margin(t = 10, r = 10, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(colour = "gray95"),
        panel.background = element_blank(), 
        legend.key = element_rect(fill = NA),
        legend.key.width = unit(3,"line"),
        legend.background = element_rect(colour = "black"),
        legend.position = "bottom",
        legend.text = element_text(margin = margin(r = 2, unit = "line")))
dev.off()

## confidence interval coverages ####
for(scenario in c("noCens", "lowCens", "highCens", 
                  "lowTreatProb", "highTreatProb", 
                  "lowVarCov", "highVarCov", 
                  "typeII")){
  for(effect in c("adverse", "no", "yes")){
    if(scenario == "lowCens" & effect == "yes"){
      ## Figure 2 ####
      png("Results/figures/figure_2.png", width = 1200, height = 500)
    }else if(scenario == "noCens" & effect == "adverse"){
      ## Figure 3 ####
      png("Results/figures/figure_3.png", width = 1200, height = 500)
    }else{
      png(paste0("Results/figures/SuppInfo/coverage_CI_", 
                 ifelse(effect == "no", "no", 
                        ifelse(effect == "yes", "", "adv")), 
                 "ATE_", scenario, ".png"), 
          width = 1200, height = 500)
    }
    print(plot_coverage_CI(scenario, effect))
    dev.off()
  }
}

## confidence band coverages ####
for(scenario in c("noCens", "lowCens", "highCens", 
                  "lowTreatProb", "highTreatProb", 
                  "lowVarCov", "highVarCov",
                  "typeII")){
  for(effect in c("adverse", "no", "yes")){
    if(scenario == "highTreatProb" & effect == "yes"){
      ## Figure 4 ####
      png("Results/figures/figure_4.png", width = 800, height = 500)
    }else{
      png(paste0("Results/figures/SuppInfo/coverage_CB_", 
                 ifelse(effect == "no", "no", 
                        ifelse(effect == "yes", "", "adv")), 
                 "ATE_", scenario, ".png"), 
          width = 800, height = 500)
    }
    print(plot_coverage_CB(scenario, effect))
    dev.off()
  }
}

## Figure 5 ####
png("Results/figures/figure_5.png", width = 800, height = 500)
plot_width_t("lowCens", "yes", "CI", t = 5)
dev.off()

## Figure 6 ####
png("Results/figures/figure_6.png", width = 800, height = 500)
plot_computation_times("lowCens", "yes")
dev.off()

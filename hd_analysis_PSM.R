# This script reproduces the analysis of the hd data w.r.t. PS matching.

# load packages
library(doRNG)
library(doParallel)
library(riskRegression)
library(ATESurvival)
library(ggplot2)
library(cowplot)

# prepare data #################################################################

# prepare data
data(hd, package = "randomForestSRC")
hd$trtgiven <- factor(hd$trtgiven, levels = c("RT","CMT"))
hd$clinstg <- factor(hd$clinstg)
hd$medwidsi <- factor(hd$medwidsi, levels = c("N","S","L"))
# break ties
set.seed(1234)
hd$time[duplicated(hd$time)] <- 
  hd$time[duplicated(hd$time)] +
  rnorm(sum(duplicated(hd$time)), 0,0.001)
hd <- hd[order(hd$time),]
hd$ID <- 1:dim(hd)[1]
# determine propensity scores
hd$PS <- predict(glm(trtgiven ~ age+sex+clinstg+medwidsi+extranod, hd,
                     family = "binomial"), 
                 hd, type = "response")
  
# plot propensity scores
ggplot() +
  geom_histogram(data = hd[hd$trtgiven == "RT",], 
                 aes(x = PS, fill = trtgiven), color = "black",
                 bins = 80) + 
  geom_histogram(data = hd[hd$trtgiven == "CMT",],
                 aes(x = PS, y = -after_stat(count), fill = trtgiven), color = "black", 
                 bins = 80) + 
  scale_fill_manual(name = "", 
                    limits = c("RT", "CMT"),
                    labels = c("Radiation", "Radiation & chemotherapy"),
                    values=c("gray80","gray30")) +
  labs(x = "Propensity score", y = "Frequency") +
  scale_x_continuous(expand = expansion(mult = c(0.01,0.01)), limits = c(0,1),
                     breaks = c(0.0, 0.25, 0.50, 0.75, 1.0)) +
  scale_y_continuous(label = abs) +
  theme(text = element_text(size = 14),
        panel.grid.minor = element_line(colour=NA),
        panel.grid.major = element_line(colour='lightgrey'),
        panel.background = element_rect(fill=NA), 
        panel.border = element_rect(fill=NA),
        legend.key = element_rect(fill=NA),
        legend.background = element_rect(colour="black"),
        legend.position = "bottom")


# analyse data #################################################################

## Relapse ####

# determine matching bias
PS_B_r <- PS_bias(ID = hd$ID,
                  A = as.numeric(hd$trtgiven)-1,
                  event = hd$status,
                  time = hd$time,
                  beta = coefficients(CSC(
                    list(Hist(time, status) ~ trtgiven + PS, 
                         Hist(time, status) ~ trtgiven + PS), 
                    hd)),
                  PS = hd$PS,
                  t = 0:35, 
                  cause = 1)
# match data
hd_matched_r <- merge(rbind(cbind(hd[,c("ID","time","status","trtgiven")],
                                  "pair_ID"=1:dim(hd)[1]),
                            cbind(merge(merge(hd[!duplicated(hd$ID),c("ID","time","status","trtgiven")], 
                                              matrix(PS_B_r$matches, ncol=2, dimnames=list(NULL, c("ID_old", "ID"))), 
                                              by="ID"),
                                        matrix(hd$ID, ncol=1, dimnames=list(NULL, c("ID_old"))),
                                        by = "ID_old"),
                                  "pair_ID"=1:dim(hd)[1])[,c("ID","time","status","trtgiven","pair_ID")]),
                      matrix(PS_B_r$n_i, ncol=2, dimnames=list(NULL, c("ID","n_i"))),
                      by="ID")
hd_matched_r$n_sampled <- hd_matched_r$n_i + 1
hd_matched_r <- hd_matched_r[order(hd_matched_r$time, hd_matched_r$ID),c("ID","time","status","trtgiven","n_sampled","pair_ID")]
# estimate ATE
cox_r <- CSC(list(Hist(time, status) ~ trtgiven, 
                  Hist(time, status) ~ trtgiven), 
             hd_matched_r, ties = "breslow")
ATE_m_r <- as.numeric(ATE(ID = hd_matched_r$ID,
                          Z = list(cox_r$models$`Cause 1`$x, 
                                   cox_r$models$`Cause 2`$x),
                          index_A = c(1,1),
                          event = hd_matched_r$status,
                          time = hd_matched_r$time,
                          beta = coefficients(cox_r),
                          t = 0:35,
                          cause = 1,
                          IF = FALSE, WBS = FALSE)$ATE - as.numeric(PS_B_r$bias))

# perform EBS
cores <- detectCores() - 1
if(.Platform$OS.type == "windows"){
  cl <- makeCluster(cores)
  registerDoParallel(cl)
}else{
  registerDoParallel(cores = cores)
}
res_EBS_r <- 
  foreach(j = 1:1000, .combine = cbind, 
          .packages = c('riskRegression', 'ATESurvival'), 
          .options.RNG = 1234) %dorng% { 
            # draw from data with replacement
            hd_bs <- hd[sample(1:dim(hd)[1], dim(hd)[1], replace = TRUE),]
            hd_bs <- hd_bs[order(hd_bs$time),]
            # determine propensity scores
            hd_bs$PS <- predict(glm(trtgiven ~ age+sex+clinstg+medwidsi+extranod, hd_bs, 
                                    family = "binomial"), 
                                hd_bs, type = "response")
            # determine matching bias
            PS_B_bs <- PS_bias(ID = hd_bs$ID,
                               A = as.numeric(hd_bs$trtgiven)-1,
                               event = hd_bs$status,
                               time = hd_bs$time,
                               beta = coefficients(CSC(
                                 list(Hist(time, status) ~ trtgiven + PS, 
                                      Hist(time, status) ~ trtgiven + PS), 
                                 hd_bs, ties = "breslow")),
                               PS = hd_bs$PS,
                               t = 0:35,
                               cause = 1)
            # match data
            hd_matched_bs <- rbind(hd_bs[,c("ID","time","status","trtgiven")],
                                   merge(merge(hd_bs[!duplicated(hd_bs$ID),c("ID","time","status","trtgiven")], 
                                               matrix(PS_B_bs$matches, ncol=2, dimnames=list(NULL, c("ID_old", "ID"))), 
                                               by="ID"),
                                         matrix(hd_bs$ID, ncol=1, dimnames=list(NULL, c("ID_old"))),
                                         by = "ID_old")[,c("ID","time","status","trtgiven")])
            hd_matched_bs <- hd_matched_bs[order(hd_matched_bs$time, hd_matched_bs$ID),]
            # estimate ATE
            cox_bs <- CSC(list(Hist(time, status) ~ trtgiven, 
                               Hist(time, status) ~ trtgiven), 
                          hd_matched_bs, ties = "breslow")
            as.numeric(ATE(ID = hd_matched_bs$ID,
                           Z = list(cox_bs$models$`Cause 1`$x, 
                                    cox_bs$models$`Cause 2`$x),
                           index_A = c(1,1),
                           event = hd_matched_bs$status,
                           time = hd_matched_bs$time,
                           beta = coefficients(cox_bs),
                           t = 0:35,
                           cause = 1,
                           IF = FALSE, WBS = FALSE)$ATE - PS_B_bs$bias)
          }
if(.Platform$OS.type == "windows"){
  stopCluster(cl)
  registerDoSEQ()
}
# derive confidence intervals & bands
res_EBS_r <- pmin(pmax(res_EBS_r, -1, na.rm=TRUE), 1, na.rm=TRUE)
CI_EBS_r <- rbind(apply(res_EBS_r, 1, quantile, 0.025, na.rm=TRUE),
                  apply(res_EBS_r, 1, quantile, 0.975, na.rm=TRUE))
# determine quantile for confidence bands
q_EBS_sup_r <- quantile(apply(abs(res_EBS_r - rowMeans(res_EBS_r, na.rm = TRUE)) / 
                                apply(res_EBS_r, 1, sd, na.rm = TRUE), 
                              2, FUN = max, na.rm = TRUE), 0.95)
CB_EBS_r <- rbind(pmax(ATE_m_r - q_EBS_sup_r * apply(res_EBS_r, 1, sd, na.rm = TRUE), -1), 
                  pmin(ATE_m_r + q_EBS_sup_r * apply(res_EBS_r, 1, sd, na.rm = TRUE), 1))


# perform double-resampling
## function to organize parallelized outcome list
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}
if(.Platform$OS.type == "windows"){
  cl <- makeCluster(cores)
  registerDoParallel(cl)
}else{
  registerDoParallel(cores = cores)
}
# generate auxiliary quantities
dr_quantities <- 
  foreach(j = 1:1000, .combine='comb', .multicombine=TRUE, .init=list(list(),list(),list()),
          .options.RNG = 1234) %dorng% {
            A_dr <- rbinom(dim(hd)[1], 1, hd$PS)
            PS_dr <- predict(glm(A_dr ~ age+sex+clinstg+medwidsi+extranod,
                                 hd, family = "binomial"),
                             hd, type = "response")
            u_dr <- rnorm(dim(hd)[1])
            list(A_dr, PS_dr, u_dr)
          }
if(.Platform$OS.type == "windows"){
  stopCluster(cl)
  registerDoSEQ()
}
# perform double-resampling
res_dr_r <- PS_dr(ID = hd$ID,
                  A = as.numeric(hd$trtgiven)-1,
                  Z = model.matrix(glm(trtgiven ~ age+sex+clinstg+medwidsi+extranod,
                                       hd, family = "binomial"))[,-1],
                  event = hd$status,
                  time = hd$time,
                  beta = unlist(coefficients(cox_r)),
                  PS = hd$PS,
                  n_i = PS_B_r$n_i[,2],
                  t = 0:35,
                  A_BS = do.call(cbind, dr_quantities[[1]]), 
                  PS_BS = do.call(cbind, dr_quantities[[2]]),
                  u_BS = do.call(cbind, dr_quantities[[3]]),
                  cause = 1)
# derive confidence intervals & bands
CI_dr_r <- rbind(pmax(ATE_m_r + apply(res_dr_r, 2, quantile, 0.025, na.rm=TRUE), -1),
                 pmin(ATE_m_r + apply(res_dr_r, 2, quantile, 0.975, na.rm=TRUE), 1))
# determine quantile for confidence bands
q_dr_sup_r <- quantile(apply(abs(t(res_dr_r) - colMeans(res_dr_r, na.rm = TRUE)) / 
                               apply(res_dr_r, 2, sd, na.rm = TRUE), 
                             2, FUN = max, na.rm = TRUE), 0.95)
CB_dr_r <- rbind(pmax(ATE_m_r - q_dr_sup_r * apply(res_dr_r, 2, sd, na.rm = TRUE), -1), 
                 pmin(ATE_m_r + q_dr_sup_r * apply(res_dr_r, 2, sd, na.rm = TRUE), 1))                                  


# perform (clustered) IF/WBS
# determine number of observed (uncensored) events for multipliers
event_num <- sum(hd$status != 0)
# calculate average treatment effect, considering unique (matched) observations
ATE_independent_r <- ATE(ID = hd_matched_r$ID,
                         Z = list(cox_r$models$`Cause 1`$x, 
                                  cox_r$models$`Cause 2`$x),
                         index_A = c(1,1),
                         event = hd_matched_r$status,
                         time = hd_matched_r$time,
                         beta = coefficients(cox_r),
                         t = 0:35,
                         cause = 1,
                         G_IF_init = rnorm(dim(hd)[1] * 1000), 
                         G_WBS_init = rnorm(1000 * event_num),
                         bs_iter = 1000)
# derive confidence intervals & bands
# IF
CI_IF_r <- rbind(pmax(ATE_independent_r$CI_IF[1,] - as.numeric(PS_B_r$bias), -1),
                 pmin(ATE_independent_r$CI_IF[2,] - as.numeric(PS_B_r$bias), 1))
CB_IF_r <- rbind(pmax(ATE_independent_r$CB_IF[1,] - as.numeric(PS_B_r$bias), -1),
                 pmin(ATE_independent_r$CB_IF[2,] - as.numeric(PS_B_r$bias), 1))
# WBS
CI_WBS_r <- rbind(pmax(ATE_independent_r$CI_WBS[1,] - as.numeric(PS_B_r$bias), -1),
                  pmin(ATE_independent_r$CI_WBS[2,] - as.numeric(PS_B_r$bias), 1))
CB_WBS_r <- rbind(pmax(ATE_independent_r$CB_WBS[1,] - as.numeric(PS_B_r$bias), -1),
                  pmin(ATE_independent_r$CB_WBS[2,] - as.numeric(PS_B_r$bias), 1))
# clustered IF/WBS
# cluster by individual subjects
beta_cluster_r <- CSC(list(Hist(time, status) ~ trtgiven,
                           Hist(time, status) ~ trtgiven),
                      hd_matched_r, ties = "breslow", weight = 1/hd_matched_r$n_sampled)
ATE_subjects_r <- ATE(ID = hd_matched_r$ID,
                      Z = list(beta_cluster_r$models$`Cause 1`$x, 
                               beta_cluster_r$models$`Cause 2`$x),
                      index_A = c(1,1),
                      event = hd_matched_r$status,
                      time = hd_matched_r$time,
                      beta = coefficients(beta_cluster_r),
                      t = 0:35,
                      cluster_init = hd_matched_r$ID,
                      cause = 1,
                      G_IF_init = rnorm(dim(hd)[1] * 1000), 
                      G_WBS_init = rnorm(1000 * event_num),
                      bs_iter = 1000)
# cluster by matched pairs
beta_cluster_r <- CSC(list(Hist(time, status) ~ trtgiven,
                           Hist(time, status) ~ trtgiven),
                      hd_matched_r, ties = "breslow", weight = rep(1/2, 2*dim(hd)[1]))
ATE_pairs_r <- ATE(ID = hd_matched_r$ID,
                   Z = list(beta_cluster_r$models$`Cause 1`$x, 
                            beta_cluster_r$models$`Cause 2`$x),
                   index_A = c(1,1),
                   event = hd_matched_r$status,
                   time = hd_matched_r$time,
                   beta = coefficients(beta_cluster_r),
                   t = 0:35,
                   cluster_init = hd_matched_r$pair_ID,
                   cause = 1,
                   G_IF_init = rnorm(dim(hd)[1] * 1000), 
                   G_WBS_init = rnorm(1000 * event_num),
                   bs_iter = 1000)
# calculate combined standard error (IF)
SE_IF_cluster_r <- sqrt((dim(hd_matched_r)[1]^2 / 
                           (dim(hd_matched_r)[1] - 1) * 
                           ATE_independent_r$SE_IF^2 +
                           length(unique(hd_matched_r$ID))^2 / 
                           (length(unique(hd_matched_r$ID)) - 1) * 
                           ATE_subjects_r$SE_IF^2 + 
                           (dim(hd_matched_r)[1]/2)^2 / 
                           (dim(hd_matched_r)[1]/2 - 1) * 
                           ATE_pairs_r$SE_IF^2) / 
                          dim(hd_matched_r)[1])
# derive confidence intervals & bands
CI_IF_cluster_r <- t(cbind(
  pmax(ATE_m_r - SE_IF_cluster_r * qnorm(1-(1-0.95)/2), -1),
  pmin(ATE_m_r + SE_IF_cluster_r * qnorm(1-(1-0.95)/2), 1)))
# determine quantile for confidence bands
q_IF_cluster_r <- quantile(apply(
  dim(hd_matched_r)[1]/(dim(hd_matched_r)[1]-1) * ATE_independent_r$Un_std_IF + 
    length(unique(hd_matched_r$ID))/(length(unique(hd_matched_r$ID))-1) * ATE_subjects_r$Un_std_IF + 
    (dim(hd_matched_r)[1]/2)/(dim(hd_matched_r)[1]/2-1) * ATE_pairs_r$Un_std_IF, 
  1, max), 0.95)
CB_IF_cluster_r <- t(cbind(
  pmax(ATE_m_r - SE_IF_cluster_r * q_IF_cluster_r, -1), 
  pmin(ATE_m_r + SE_IF_cluster_r * q_IF_cluster_r, 1)))
# calculate combined standard error (WBS)
SE_WBS_cluster_r <- sqrt((dim(hd_matched_r)[1]^2 / 
                            (dim(hd_matched_r)[1] - 1) * 
                            ATE_independent_r$SE_WBS^2 +
                            length(unique(hd_matched_r$ID))^2 / 
                            (length(unique(hd_matched_r$ID)) - 1) * 
                            ATE_subjects_r$SE_WBS^2 + 
                            (dim(hd_matched_r)[1]/2)^2 / 
                            (dim(hd_matched_r)[1]/2 - 1) * 
                            ATE_pairs_r$SE_WBS^2) / 
                           dim(hd_matched_r)[1])
# derive confidence intervals & bands
CI_WBS_cluster_r <- t(cbind(
  pmax(ATE_m_r - SE_WBS_cluster_r * qnorm(1-(1-0.95)/2), -1),
  pmin(ATE_m_r + SE_WBS_cluster_r * qnorm(1-(1-0.95)/2), 1)))
# determine quantile for confidence bands
q_WBS_cluster_r <- quantile(apply(
  dim(hd_matched_r)[1]/(dim(hd_matched_r)[1]-1) * ATE_independent_r$Un_std_WBS + 
    length(unique(hd_matched_r$ID))/(length(unique(hd_matched_r$ID))-1) * ATE_subjects_r$Un_std_WBS + 
    (dim(hd_matched_r)[1]/2)/(dim(hd_matched_r)[1]/2-1) * ATE_pairs_r$Un_std_WBS, 
  1, max), 0.95)
CB_WBS_cluster_r <- t(cbind(
  pmax(ATE_m_r - SE_WBS_cluster_r * q_WBS_cluster_r, -1), 
  pmin(ATE_m_r + SE_WBS_cluster_r * q_WBS_cluster_r, 1)))

                                      
# store results
res_r <- data.frame(
  time = rep(0:35, 6),
  type = rep(c("EBS",
               "DR",
               "IF",
               "clustered IF",
               "WBS",
               "clustered WBS"), 
             each = 36),
  ATE_hat = rep(ATE_m_r, 6),
  CI_lower = c(CI_EBS_r[1,],
               CI_dr_r[1,],
               CI_IF_r[1,],
               CI_IF_cluster_r[1,],
               CI_WBS_r[1,],
               CI_WBS_cluster_r[1,]),
  CI_upper = c(CI_EBS_r[2,],
               CI_dr_r[2,],
               CI_IF_r[2,],
               CI_IF_cluster_r[2,],
               CI_WBS_r[2,],
               CI_WBS_cluster_r[2,]),
  CB_lower = c(CB_EBS_r[1,],
               CB_dr_r[1,],
               CB_IF_r[1,],
               CB_IF_cluster_r[1,],
               CB_WBS_r[1,],
               CB_WBS_cluster_r[1,]),
  CB_upper = c(CB_EBS_r[2,],
               CB_dr_r[2,],
               CB_IF_r[2,],
               CB_IF_cluster_r[2,],
               CB_WBS_r[2,],
               CB_WBS_cluster_r[2,])
)

# create plot
p1_r <- ggplot(data = res_r, aes(x = time)) +
  geom_line(aes(y = 100*ATE_hat), color = "black", linewidth=1.1) + 
  geom_line(aes(x = 0, y = 0, color = type), linewidth=1.1, lty = 1) + 
  geom_line(aes(y = 100*CI_lower, color = type), linewidth=1.1, lty = 5) +
  geom_line(aes(y = 100*CI_upper, color = type), linewidth=1.1, lty = 5) +
  scale_color_manual(name = "", 
                     values = c('orange','purple','red','darkred','turquoise4','darkblue'),
                     labels = c("EBS","DR","IF","clustered IF","WBS","clustered WBS")) +
  labs(title = paste0("Confidence intervals"), 
       x = "years since diagnosis", 
       y = "absolute risk difference [%]") +
  scale_x_continuous(breaks=seq(0,35,5), 
                     labels = c("0","","10","","20","","30","")) +
  scale_y_continuous(breaks = seq(-50,20,5), 
                     labels = c("-50","","-40","","-30","","-20","","-10","","0","","10","","20"), 
                     limits = c(-50,16.5)) +
  theme(text = element_text(size=15),
        plot.title = element_text(hjust=0.5, face="bold"),
        panel.grid.minor = element_line(colour=NA),
        panel.grid.major = element_line(colour='lightgrey'),
        panel.background = element_rect(fill=NA), 
        panel.border = element_rect(fill=NA),
        legend.key = element_rect(fill=NA),
        legend.background = element_rect(colour="black"),
        legend.position = "bottom")
p2_r <- ggplot(data = res_r, aes(x = time)) +
  geom_line(aes(y = 100*ATE_hat), color = "black", linewidth=1.1) + 
  geom_line(aes(x = 0, y = 0, color = type), linewidth=1.1, lty = 1) + 
  geom_line(aes(y = 100*CB_lower, color = type), linewidth=1.1, lty = 5) +
  geom_line(aes(y = 100*CB_upper, color = type), linewidth=1.1, lty = 5) +
  scale_color_manual(name = "", 
                     values = c('orange','purple','red','darkred','turquoise4','darkblue'),
                     labels = c("EBS","DR","IF","clustered IF","WBS","clustered WBS")) +
  labs(title = paste0("Confidence bands"), 
       x = "years since diagnosis", 
       y = "") + 
  scale_x_continuous(breaks=seq(0,35,5), 
                     labels = c("0","","10","","20","","30","")) +
  scale_y_continuous(breaks = seq(-50,20,5), 
                     labels = c("-50","","-40","","-30","","-20","","-10","","0","","10","","20"), 
                     limits = c(-50,16.5)) +
  theme(text = element_text(size=15),
        plot.title = element_text(hjust=0.5, face="bold"),
        panel.grid.minor = element_line(colour=NA),
        panel.grid.major = element_line(colour='lightgrey'),
        panel.background = element_rect(fill=NA), 
        panel.border = element_rect(fill=NA),
        legend.position = "none")
p_r <- plot_grid(p1_r + theme(legend.position = "none"), p2_r)
plot_grid(p_r, get_legend(p1_r), nrow=2, rel_heights = c(0.75, 0.2))

                                      
## Death ####

# determine matching bias
PS_B_d <- PS_bias(ID = hd$ID,
                  A = as.numeric(hd$trtgiven)-1,
                  event = hd$status,
                  time = hd$time,
                  beta = coefficients(CSC(
                    list(Hist(time, status) ~ trtgiven + PS, 
                         Hist(time, status) ~ trtgiven + PS), 
                    hd)),
                  PS = hd$PS,
                  t = 0:35, 
                  cause = 2)
# match data
hd_matched_d <- merge(rbind(cbind(hd[,c("ID","time","status","trtgiven")],
                                  "pair_ID"=1:dim(hd)[1]),
                            cbind(merge(merge(hd[!duplicated(hd$ID),c("ID","time","status","trtgiven")], 
                                              matrix(PS_B_d$matches, ncol=2, dimnames=list(NULL, c("ID_old", "ID"))), 
                                              by="ID"),
                                        matrix(hd$ID, ncol=1, dimnames=list(NULL, c("ID_old"))),
                                        by = "ID_old"),
                                  "pair_ID"=1:dim(hd)[1])[,c("ID","time","status","trtgiven","pair_ID")]),
                      matrix(PS_B_d$n_i, ncol=2, dimnames=list(NULL, c("ID","n_i"))),
                      by="ID")
hd_matched_d$n_sampled <- hd_matched_d$n_i + 1
hd_matched_d <- hd_matched_d[order(hd_matched_d$time, hd_matched_d$ID),c("ID","time","status","trtgiven","n_sampled","pair_ID")]
# estimate ATE
cox_d <- CSC(list(Hist(time, status) ~ trtgiven, 
                  Hist(time, status) ~ trtgiven), 
             hd_matched_d, ties = "breslow")
ATE_m_d <- as.numeric(ATE(ID = hd_matched_d$ID,
                          Z = list(cox_d$models$`Cause 1`$x, 
                                   cox_d$models$`Cause 2`$x),
                          index_A = c(1,1),
                          event = hd_matched_d$status,
                          time = hd_matched_d$time,
                          beta = coefficients(cox_d),
                          t = 0:35,
                          cause = 2,
                          IF = FALSE, WBS = FALSE)$ATE - as.numeric(PS_B_d$bias))

# perform EBS
cores <- detectCores() - 1
if(.Platform$OS.type == "windows"){
  cl <- makeCluster(cores)
  registerDoParallel(cl)
}else{
  registerDoParallel(cores = cores)
}
res_EBS_d <- 
  foreach(j = 1:1000, .combine = cbind, 
          .packages = c('riskRegression', 'ATESurvival'), 
          .options.RNG = 1234) %dorng% { 
            # draw from data with replacement
            hd_bs <- hd[sample(1:dim(hd)[1], dim(hd)[1], replace = TRUE),]
            hd_bs <- hd_bs[order(hd_bs$time),]
            # determine propensity scores
            hd_bs$PS <- predict(glm(trtgiven ~ age+sex+clinstg+medwidsi+extranod, hd_bs, 
                                    family = "binomial"), 
                                hd_bs, type = "response")
            # determine matching bias
            PS_B_bs <- PS_bias(ID = hd_bs$ID,
                               A = as.numeric(hd_bs$trtgiven)-1,
                               event = hd_bs$status,
                               time = hd_bs$time,
                               beta = coefficients(CSC(
                                 list(Hist(time, status) ~ trtgiven + PS, 
                                      Hist(time, status) ~ trtgiven + PS), 
                                 hd_bs, ties = "breslow")),
                               PS = hd_bs$PS,
                               t = 0:35,
                               cause = 2)
            # match data
            hd_matched_bs <- rbind(hd_bs[,c("ID","time","status","trtgiven")],
                                   merge(merge(hd_bs[!duplicated(hd_bs$ID),c("ID","time","status","trtgiven")], 
                                               matrix(PS_B_bs$matches, ncol=2, dimnames=list(NULL, c("ID_old", "ID"))), 
                                               by="ID"),
                                         matrix(hd_bs$ID, ncol=1, dimnames=list(NULL, c("ID_old"))),
                                         by = "ID_old")[,c("ID","time","status","trtgiven")])
            hd_matched_bs <- hd_matched_bs[order(hd_matched_bs$time, hd_matched_bs$ID),]
            # estimate ATE
            cox_bs <- CSC(list(Hist(time, status) ~ trtgiven, 
                               Hist(time, status) ~ trtgiven), 
                          hd_matched_bs, ties = "breslow")
            as.numeric(ATE(ID = hd_matched_bs$ID,
                           Z = list(cox_bs$models$`Cause 1`$x, 
                                    cox_bs$models$`Cause 2`$x),
                           index_A = c(1,1),
                           event = hd_matched_bs$status,
                           time = hd_matched_bs$time,
                           beta = coefficients(cox_bs),
                           t = 0:35,
                           cause = 2,
                           IF = FALSE, WBS = FALSE)$ATE - PS_B_bs$bias)
          }
if(.Platform$OS.type == "windows"){
  stopCluster(cl)
  registerDoSEQ()
}
# derive confidence intervals & bands
res_EBS_d <- pmin(pmax(res_EBS_d, -1, na.rm=TRUE), 1, na.rm=TRUE)
CI_EBS_d <- rbind(apply(res_EBS_d, 1, quantile, 0.025, na.rm=TRUE),
                  apply(res_EBS_d, 1, quantile, 0.975, na.rm=TRUE))
# determine quantile for confidence bands
q_EBS_sup_d <- quantile(apply(abs(res_EBS_d - rowMeans(res_EBS_d, na.rm = TRUE)) / 
                                apply(res_EBS_d, 1, sd, na.rm = TRUE), 
                              2, FUN = max, na.rm = TRUE), 0.95)
CB_EBS_d <- rbind(pmax(ATE_m_d - q_EBS_sup_d * apply(res_EBS_d, 1, sd, na.rm = TRUE), -1), 
                  pmin(ATE_m_d + q_EBS_sup_d * apply(res_EBS_d, 1, sd, na.rm = TRUE), 1))


# perform double-resampling
res_dr_d <- PS_dr(ID = hd$ID,
                  A = as.numeric(hd$trtgiven)-1,
                  Z = model.matrix(glm(trtgiven ~ age+sex+clinstg+medwidsi+extranod,
                                       hd, family = "binomial"))[,-1],
                  event = hd$status,
                  time = hd$time,
                  beta = unlist(coefficients(cox_d)),
                  PS = hd$PS,
                  n_i = PS_B_d$n_i[,2],
                  t = 0:35,
                  A_BS = do.call(cbind, dr_quantities[[1]]), 
                  PS_BS = do.call(cbind, dr_quantities[[2]]),
                  u_BS = do.call(cbind, dr_quantities[[3]]),
                  cause = 2)
# derive confidence intervals & bands
CI_dr_d <- rbind(pmax(ATE_m_d + apply(res_dr_d, 2, quantile, 0.025, na.rm=TRUE), -1),
                 pmin(ATE_m_d + apply(res_dr_d, 2, quantile, 0.975, na.rm=TRUE), 1))
# determine quantile for confidence bands
q_dr_sup_d <- quantile(apply(abs(t(res_dr_d) - colMeans(res_dr_d, na.rm = TRUE)) / 
                               apply(res_dr_d, 2, sd, na.rm = TRUE), 
                             2, FUN = max, na.rm = TRUE), 0.95)
CB_dr_d <- rbind(pmax(ATE_m_d - q_dr_sup_d * apply(res_dr_d, 2, sd, na.rm = TRUE), -1), 
                 pmin(ATE_m_d + q_dr_sup_d * apply(res_dr_d, 2, sd, na.rm = TRUE), 1))                                  


# perform (clustered) IF/WBS
# calculate average treatment effect, considering unique (matched) observations
ATE_independent_d <- ATE(ID = hd_matched_d$ID,
                         Z = list(cox_d$models$`Cause 1`$x, 
                                  cox_d$models$`Cause 2`$x),
                         index_A = c(1,1),
                         event = hd_matched_d$status,
                         time = hd_matched_d$time,
                         beta = coefficients(cox_d),
                         t = 0:35,
                         cause = 2,
                         G_IF_init = rnorm(dim(hd)[1] * 1000), 
                         G_WBS_init = rnorm(1000 * event_num),
                         bs_iter = 1000)
# derive confidence intervals & bands
# IF
CI_IF_d <- rbind(pmax(ATE_independent_d$CI_IF[1,] - as.numeric(PS_B_d$bias), -1),
                 pmin(ATE_independent_d$CI_IF[2,] - as.numeric(PS_B_d$bias), 1))
CB_IF_d <- rbind(pmax(ATE_independent_d$CB_IF[1,] - as.numeric(PS_B_d$bias), -1),
                 pmin(ATE_independent_d$CB_IF[2,] - as.numeric(PS_B_d$bias), 1))
# WBS
CI_WBS_d <- rbind(pmax(ATE_independent_d$CI_WBS[1,] - as.numeric(PS_B_d$bias), -1),
                  pmin(ATE_independent_d$CI_WBS[2,] - as.numeric(PS_B_d$bias), 1))
CB_WBS_d <- rbind(pmax(ATE_independent_d$CB_WBS[1,] - as.numeric(PS_B_d$bias), -1),
                  pmin(ATE_independent_d$CB_WBS[2,] - as.numeric(PS_B_d$bias), 1))
# clustered IF/WBS
# cluster by individual subjects
beta_cluster_d <- CSC(list(Hist(time, status) ~ trtgiven,
                           Hist(time, status) ~ trtgiven),
                      hd_matched_d, ties = "breslow", weight = 1/hd_matched_d$n_sampled)
ATE_subjects_d <- ATE(ID = hd_matched_d$ID,
                      Z = list(beta_cluster_d$models$`Cause 1`$x, 
                               beta_cluster_d$models$`Cause 2`$x),
                      index_A = c(1,1),
                      event = hd_matched_d$status,
                      time = hd_matched_d$time,
                      beta = coefficients(beta_cluster_d),
                      t = 0:35,
                      cluster_init = hd_matched_d$ID,
                      cause = 2,
                      G_IF_init = rnorm(dim(hd)[1] * 1000), 
                      G_WBS_init = rnorm(1000 * event_num),
                      bs_iter = 1000)
# cluster by matched pairs
beta_cluster_d <- CSC(list(Hist(time, status) ~ trtgiven,
                           Hist(time, status) ~ trtgiven),
                      hd_matched_d, ties = "breslow", weight = rep(1/2, 2*dim(hd)[1]))
ATE_pairs_d <- ATE(ID = hd_matched_d$ID,
                   Z = list(beta_cluster_d$models$`Cause 1`$x, 
                            beta_cluster_d$models$`Cause 2`$x),
                   index_A = c(1,1),
                   event = hd_matched_d$status,
                   time = hd_matched_d$time,
                   beta = coefficients(beta_cluster_d),
                   t = 0:35,
                   cluster_init = hd_matched_d$pair_ID,
                   cause = 2,
                   G_IF_init = rnorm(dim(hd)[1] * 1000), 
                   G_WBS_init = rnorm(1000 * event_num),
                   bs_iter = 1000)
# calculate combined standard error (IF)
SE_IF_cluster_d <- sqrt((dim(hd_matched_d)[1]^2 / 
                           (dim(hd_matched_d)[1] - 1) * 
                           ATE_independent_d$SE_IF^2 +
                           length(unique(hd_matched_d$ID))^2 / 
                           (length(unique(hd_matched_d$ID)) - 1) * 
                           ATE_subjects_d$SE_IF^2 + 
                           (dim(hd_matched_d)[1]/2)^2 / 
                           (dim(hd_matched_d)[1]/2 - 1) * 
                           ATE_pairs_d$SE_IF^2) / 
                          dim(hd_matched_d)[1])
# derive confidence intervals & bands
CI_IF_cluster_d <- t(cbind(
  pmax(ATE_m_d - SE_IF_cluster_d * qnorm(1-(1-0.95)/2), -1),
  pmin(ATE_m_d + SE_IF_cluster_d * qnorm(1-(1-0.95)/2), 1)))
# determine quantile for confidence bands
q_IF_cluster_d <- quantile(apply(
  dim(hd_matched_d)[1]/(dim(hd_matched_d)[1]-1) * ATE_independent_d$Un_std_IF + 
    length(unique(hd_matched_d$ID))/(length(unique(hd_matched_d$ID))-1) * ATE_subjects_d$Un_std_IF + 
    (dim(hd_matched_d)[1]/2)/(dim(hd_matched_d)[1]/2-1) * ATE_pairs_d$Un_std_IF, 
  1, max), 0.95)
CB_IF_cluster_d <- t(cbind(
  pmax(ATE_m_d - SE_IF_cluster_d * q_IF_cluster_d, -1), 
  pmin(ATE_m_d + SE_IF_cluster_d * q_IF_cluster_d, 1)))
# calculate combined standard error (WBS)
SE_WBS_cluster_d <- sqrt((dim(hd_matched_d)[1]^2 / 
                            (dim(hd_matched_d)[1] - 1) * 
                            ATE_independent_d$SE_WBS^2 +
                            length(unique(hd_matched_d$ID))^2 / 
                            (length(unique(hd_matched_d$ID)) - 1) * 
                            ATE_subjects_d$SE_WBS^2 + 
                            (dim(hd_matched_d)[1]/2)^2 / 
                            (dim(hd_matched_d)[1]/2 - 1) * 
                            ATE_pairs_d$SE_WBS^2) / 
                           dim(hd_matched_d)[1])
# derive confidence intervals & bands
CI_WBS_cluster_d <- t(cbind(
  pmax(ATE_m_d - SE_WBS_cluster_d * qnorm(1-(1-0.95)/2), -1),
  pmin(ATE_m_d + SE_WBS_cluster_d * qnorm(1-(1-0.95)/2), 1)))
# determine quantile for confidence bands
q_WBS_cluster_d <- quantile(apply(
  dim(hd_matched_d)[1]/(dim(hd_matched_d)[1]-1) * ATE_independent_d$Un_std_WBS + 
    length(unique(hd_matched_d$ID))/(length(unique(hd_matched_d$ID))-1) * ATE_subjects_d$Un_std_WBS + 
    (dim(hd_matched_d)[1]/2)/(dim(hd_matched_d)[1]/2-1) * ATE_pairs_d$Un_std_WBS, 
  1, max), 0.95)
CB_WBS_cluster_d <- t(cbind(
  pmax(ATE_m_d - SE_WBS_cluster_d * q_WBS_cluster_d, -1), 
  pmin(ATE_m_d + SE_WBS_cluster_d * q_WBS_cluster_d, 1)))

                                      
# store results
res_d <- data.frame(
  time = rep(0:35, 6),
  type = rep(c("EBS",
               "DR",
               "IF",
               "clustered IF",
               "WBS",
               "clustered WBS"), 
             each = 36),
  ATE_hat = rep(ATE_m_d, 6),
  CI_lower = c(CI_EBS_d[1,],
               CI_dr_d[1,],
               CI_IF_d[1,],
               CI_IF_cluster_d[1,],
               CI_WBS_d[1,],
               CI_WBS_cluster_d[1,]),
  CI_upper = c(CI_EBS_d[2,],
               CI_dr_d[2,],
               CI_IF_d[2,],
               CI_IF_cluster_d[2,],
               CI_WBS_d[2,],
               CI_WBS_cluster_d[2,]),
  CB_lower = c(CB_EBS_d[1,],
               CB_dr_d[1,],
               CB_IF_d[1,],
               CB_IF_cluster_d[1,],
               CB_WBS_d[1,],
               CB_WBS_cluster_d[1,]),
  CB_upper = c(CB_EBS_d[2,],
               CB_dr_d[2,],
               CB_IF_d[2,],
               CB_IF_cluster_d[2,],
               CB_WBS_d[2,],
               CB_WBS_cluster_d[2,])
)
                                      
# create plot
p1_d <- ggplot(data = res_d, aes(x = time)) +
  geom_line(aes(y = 100*ATE_hat), color = "black", linewidth=1.1) + 
  geom_line(aes(x = 0, y = 0, color = type), linewidth=1.1, lty = 1) + 
  geom_line(aes(y = 100*CI_lower, color = type), linewidth=1.1, lty = 5) +
  geom_line(aes(y = 100*CI_upper, color = type), linewidth=1.1, lty = 5) +
  scale_color_manual(name = "", 
                     values = c('orange','purple','red','darkred','turquoise4','darkblue'),
                     labels = c("EBS","DR","IF","clustered IF","WBS","clustered WBS")) +
  labs(title = paste0("Confidence intervals"), 
       x = "years since diagnosis", 
       y = "absolute risk difference [%]") +
  scale_x_continuous(breaks=seq(0,35,5), 
                     labels = c("0","","10","","20","","30","")) +
  scale_y_continuous(breaks = seq(-30,45,5), 
                     labels = c("-30","","-20","","-10","","0","","10","","20","","30","","40",""), 
                     limits = c(-30,41.5)) +
  theme(text = element_text(size=15),
        plot.title = element_text(hjust=0.5, face="bold"),
        panel.grid.minor = element_line(colour=NA),
        panel.grid.major = element_line(colour='lightgrey'),
        panel.background = element_rect(fill=NA), 
        panel.border = element_rect(fill=NA),
        legend.key = element_rect(fill=NA),
        legend.background = element_rect(colour="black"),
        legend.position = "bottom")
p2_d <- ggplot(data = res_d, aes(x = time)) +
  geom_line(aes(y = 100*ATE_hat), color = "black", linewidth=1.1) + 
  geom_line(aes(x = 0, y = 0, color = type), linewidth=1.1, lty = 1) + 
  geom_line(aes(y = 100*CB_lower, color = type), linewidth=1.1, lty = 5) +
  geom_line(aes(y = 100*CB_upper, color = type), linewidth=1.1, lty = 5) +
  scale_color_manual(name = "", 
                     values = c('orange','purple','red','darkred','turquoise4','darkblue'),
                     labels = c("EBS","DR","IF","clustered IF","WBS","clustered WBS")) +
  labs(title = paste0("Confidence bands"), 
       x = "years since diagnosis", 
       y = "") + 
  scale_x_continuous(breaks=seq(0,35,5), 
                     labels = c("0","","10","","20","","30","")) +
  scale_y_continuous(breaks = seq(-30,45,5), 
                     labels = c("-30","","-20","","-10","","0","","10","","20","","30","","40",""), 
                     limits = c(-30,41.5)) +
  theme(text = element_text(size=15),
        plot.title = element_text(hjust=0.5, face="bold"),
        panel.grid.minor = element_line(colour=NA),
        panel.grid.major = element_line(colour='lightgrey'),
        panel.background = element_rect(fill=NA), 
        panel.border = element_rect(fill=NA),
        legend.position = "none")
p_r <- plot_grid(p1_r + theme(legend.position = "none"), p2_r)
plot_grid(p_r, get_legend(p1_r), nrow=2, rel_heights = c(0.75, 0.2))

# This script reproduces the analysis of the hd data.

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


## Table 2 ####
n_treat <- table(hd$trtgiven)
descrip <- rbind(paste0(round(aggregate(age ~ trtgiven, hd, mean)$age, 2), " (",
                        round(aggregate(age ~ trtgiven, hd, sd)$age, 2), ")"), 
                 paste0(table(hd$sex, hd$trtgiven)["M",], " (", 
                        round(100*table(hd$sex, hd$trtgiven)["M",] / n_treat, 2), 
                        "%)"),
                 paste0(table(hd$clinstg, hd$trtgiven)["1",], " (", 
                        round(100*table(hd$clinstg, hd$trtgiven)["1",] / n_treat, 2), 
                        "%)"),
                 matrix(paste0(table(hd$medwidsi, hd$trtgiven), " (",
                               round(100 * table(hd$medwidsi, hd$trtgiven) / 
                                       matrix(n_treat, nrow = 3, ncol = 2, byrow = TRUE), 
                                     2), 
                               "%)"), ncol = 2),
                 paste0(table(hd$extranod, hd$trtgiven)["Y",], " (", 
                        round(100*table(hd$extranod, hd$trtgiven)["Y",] / n_treat, 2), 
                        "%)"))
colnames(descrip) <- paste0(c("Radiation alone", "Radiation & chemotherapy"), 
                            " (n = ", n_treat, ")")
rownames(descrip) <- c("age, mean (sd)", 
                       "sex: male", 
                       "lymphoma stage: I", 
                       "no mediastinum involvement",
                       "small mediastinum involvement",
                       "large mediasstinum involvement",
                       "extranodal disease")
descrip


# analyse data #################################################################

# perform EBS
cores <- detectCores() - 1
if(.Platform$OS.type == "windows"){
  cl <- makeCluster(cores)
  registerDoParallel(cl)
}else{
  registerDoParallel(cores = cores)
}
res_EBS <- 
  foreach(j = 1:1000, .combine = cbind, 
          .packages = c('riskRegression', 'ATESurvival'), 
          .options.RNG = 1234) %dorng% { 
            # draw from data with replacement
            hd_bs <- hd[sample(1:dim(hd)[1], dim(hd)[1], replace = TRUE),]
            # fit cause-specific Cox models
            invisible(capture.output(
              csc_bs <- CSC(
                formula = 
                  list(Hist(time, status) ~ trtgiven+age+sex+clinstg+medwidsi+extranod,
                       Hist(time, status) ~ trtgiven+age+sex+clinstg+medwidsi+extranod),
                hd_bs[order(hd_bs$time),])
            ))
            # calculate average treatment effect
            EBS(
              Z = list(csc_bs$models$`Cause 1`$x,
                       csc_bs$models$`Cause 2`$x),
              event = hd_bs$status[order(hd_bs$time)],
              time = sort(hd_bs$time),
              t = 0:35,
              beta = coefficients(csc_bs),
              index_A = c(1,1),
              cause = 2
            )
          }
if(.Platform$OS.type == "windows"){
  stopCluster(cl)
  registerDoSEQ()
}
# determine quantile for confidence bands
q_EBS_sup <- quantile(apply(abs(res_EBS - rowMeans(res_EBS, na.rm = TRUE)) / 
                              sqrt(apply(res_EBS, 1, var, na.rm = TRUE)), 
                            2, max, na.rm = TRUE), 0.95)

# fit cause-specific Cox models
csc <- CSC(formula = list(Hist(time, status) ~ trtgiven+age+sex+clinstg+medwidsi+extranod,
                          Hist(time, status) ~ trtgiven+age+sex+clinstg+medwidsi+extranod),
           hd[order(hd$time),])

# perform IF/WBS
set.seed(1234, kind = "Mersenne-Twister") # for reproducibility
# determine number of observed (uncensored) events for multipliers
event_num <- sum(hd$status[order(hd$time)] != 0)
# compute weird bootstrap multipliers
mult_Weird <- sapply(1:event_num, function(j){
  time <- sort(hd$time)
  Y <- sum(time[hd$status[order(hd$time)] != 0][j] <= time)
  rbinom(1e4, Y, 1/Y) - 1
})
# calculate average treatment effect, confidence intervals & bands
res_IF_WBS <- ATE_IF_WBS(
  Z = list(csc$models$`Cause 1`$x,
           csc$models$`Cause 2`$x),
  event = hd$status[order(hd$time)], 
  time = sort(hd$time),
  t = 0:35,
  beta = coefficients(csc),
  index_A = c(1,1),
  cause = 2,
  G_IF = rnorm(dim(hd)[1] * 1e4),
  G_Lin_init = rnorm(1e4 * event_num),
  G_Bey_init = rpois(1e4 * event_num, 1) - 1,
  G_Weird_init = c(t(mult_Weird))
)

# store results
res <- data.frame(
  time = rep(0:35, 5),
  type = rep(c("EBS",
               "IF",
               "WBS - Lin et al.",
               "WBS - Beyersmann et al.",
               "WBS - Weird bootstrap"), 
             each = 36),
  ATE_hat = rep(res_IF_WBS[[1]], 5),
  CI_lower = c(apply(res_EBS, 1, 
                     function(res_EBS_t){max(quantile(res_EBS_t, 0.025), -1)}),
               res_IF_WBS[[2]][1,],
               res_IF_WBS[[4]][3,],
               res_IF_WBS[[4]][5,],
               res_IF_WBS[[4]][7,]),
  CI_upper = c(apply(res_EBS, 1, 
                     function(res_EBS_t){min(quantile(res_EBS_t, 0.975), 1)}),
               res_IF_WBS[[2]][2,],
               res_IF_WBS[[4]][4,],
               res_IF_WBS[[4]][6,],
               res_IF_WBS[[4]][8,]),
  CB_lower = c(as.numeric(res_IF_WBS[[1]]) - 
                 q_EBS_sup * sqrt(apply(res_EBS, 1, var, na.rm = TRUE)),
               res_IF_WBS[[3]][1,],
               res_IF_WBS[[5]][3,],
               res_IF_WBS[[5]][7,],
               res_IF_WBS[[5]][11,]),
  CB_upper = c(as.numeric(res_IF_WBS[[1]]) + 
                 q_EBS_sup * sqrt(apply(res_EBS, 1, var, na.rm = TRUE)),
               res_IF_WBS[[3]][2,],
               res_IF_WBS[[5]][4,],
               res_IF_WBS[[5]][8,],
               res_IF_WBS[[5]][12,])
)

## Figure 7 ####

# create plot
p1 <- ggplot(data = res, aes(x = time)) +
  geom_line(aes(y = 100*ATE_hat), color = "black", linewidth=1.1) + 
  geom_line(aes(x = 0, y = 0, color = type), linewidth=1.1, lty = 1) + 
  geom_line(aes(y = 100*CI_lower, color = type), linewidth=1.1, lty = 5) +
  geom_line(aes(y = 100*CI_upper, color = type), linewidth=1.1, lty = 5) +
  scale_color_manual(name = "", 
                     values = c('orange','red','turquoise4','cyan3','chartreuse2'),
                     labels = c('EBS',
                                'IF',
                                'WBS - Lin et al.',
                                'WBS - Beyersmann et al.',
                                'WBS - Weird bootstrap')) +
  labs(title = paste0("Confidence intervals"), 
       x = "years since diagnosis", 
       y = "absolute risk difference [%]") +
  scale_x_continuous(breaks=seq(0,35,5), 
                     labels = c("0","","10","","20","","30","")) +
  theme(text = element_text(size=15),
        plot.title = element_text(hjust=0.5, face="bold"),
        panel.grid.minor = element_line(colour=NA),
        panel.grid.major = element_line(colour='lightgrey'),
        panel.background = element_rect(fill=NA), 
        panel.border = element_rect(fill=NA),
        legend.key = element_rect(fill=NA),
        legend.background = element_rect(colour="black"),
        legend.position = "bottom")
p2 <- ggplot(data = res, aes(x = time)) +
  geom_line(aes(y = 100*ATE_hat), color = "black", linewidth=1.1) + 
  geom_line(aes(x = 0, y = 0, color = type), linewidth=1.1, lty = 1) + 
  geom_line(aes(y = 100*CB_lower, color = type), linewidth=1.1, lty = 5) +
  geom_line(aes(y = 100*CB_upper, color = type), linewidth=1.1, lty = 5) +
  scale_color_manual(name = "", 
                     values = c('orange','red','turquoise4','cyan3','chartreuse2'), 
                     #drop = FALSE,
                     labels = c('EBS',
                                'IF',
                                'WBS - Lin et al.',
                                'WBS - Beyersmann et al.',
                                'WBS - Weird bootstrap')) +
  labs(title = paste0("Confidence bands"), 
       x = "years since diagnosis", 
       y = "") + 
  scale_x_continuous(breaks=seq(0,35,5), 
                     labels = c("0","","10","","20","","30","")) +
  theme(text = element_text(size=15),
        plot.title = element_text(hjust=0.5, face="bold"),
        panel.grid.minor = element_line(colour=NA),
        panel.grid.major = element_line(colour='lightgrey'),
        panel.background = element_rect(fill=NA), 
        panel.border = element_rect(fill=NA),
        legend.position = "none")
p <- plot_grid(p1 + theme(legend.position = "none"), p2)
png("Results/figures/figure_7.png", width = 1200, height = 500)
plot_grid(p, get_legend(p1), nrow=2, rel_heights = c(0.75, 0.2))
dev.off()

# This script includes the functions for the analysis of the hd data.

# prepare data
data(hd, package = "randomForestSRC")
hd$trtgiven <- factor(hd$trtgiven, levels = c("RT","CMT"))
hd$clinstg <- factor(hd$clinstg)
hd$medwidsi <- factor(hd$medwidsi, levels = c("N","S","L"))

# prepare IF/WBS: fit (cause-specific) Cox models
library(riskRegression)
library(ATESurvival)
library(doRNG)
csc <- CSC(formula = list(Hist(time, status) ~ 
                            trtgiven+age+sex+clinstg+medwidsi+extranod,
                          Hist(time, status) ~ 
                            trtgiven+age+sex+clinstg+medwidsi+extranod),
           hd[order(hd$time),])

# perform IF/WBS
set.seed(1234)
res_IF_WBS <- ATE_IF_WBS(
  Z = list(csc$models$`Cause 1`$x,
           csc$models$`Cause 2`$x),
  event = hd$status[order(hd$time)], 
  time = sort(hd$time),
  t = 0:35,
  beta = coefficients(csc),
  index_A = c(1,1),
  cause = 2
)

# perform EBS
cl <- makeCluster(7)
registerDoParallel(cl)
res_EBS <- 
  foreach(j = 1:1000, .combine = cbind, 
          .packages = c('riskRegression', 'ATESurvival'), 
          .options.RNG = 1234) %dorng% { 
            hd_bs <- hd[sample(1:dim(hd)[1], dim(hd)[1], replace = TRUE),]
            invisible(capture.output(
              csc_bs <- CSC(
                formula = 
                  list(Hist(time, status) ~ 
                         trtgiven+age+sex+clinstg+medwidsi+extranod,
                       Hist(time, status) ~ 
                         trtgiven+age+sex+clinstg+medwidsi+extranod),
                hd_bs[order(hd_bs$time),])
            ))
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
stopCluster(cl)
q_EBS_sup <- quantile(apply(abs(res_EBS - rowMeans(res_EBS, na.rm = TRUE)) / 
                              sqrt(apply(res_EBS, 1, var, na.rm = TRUE)), 
                            2, max, na.rm = TRUE), 0.95)

# save results
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

# plot results: confidence intervals
library(ggplot2)
ggplot(data = res, aes(x = time)) +
  geom_line(aes(y = 100*ATE_hat)) + 
  geom_line(aes(y = 100*CI_lower, color = type), lty = 5) +
  geom_line(aes(y = 100*CI_upper, color = type), lty = 5) +
  scale_color_manual(name = "", 
                     values = c('orange',
                                'red',
                                'turquoise4',
                                'cyan3',
                                'chartreuse2'),
                     labels = c('EBS',
                                'IF',
                                'WBS - Lin et al.',
                                'WBS - Beyersmann et al.',
                                'WBS - Weird bootstrap')) +
  scale_x_continuous(breaks=seq(0,35,5), 
                     labels = c("0","","10","","20","","30","")) +
  labs(title="Confidence Intervals", 
       x="years since diagnosis", y="absolute risk difference [%]") + 
  theme(text = element_text(size=15),
        plot.title = element_text(hjust=0.5, face="bold"),
        panel.grid.minor = element_line(colour=NA),
        panel.grid.major = element_line(colour='lightgrey'),
        panel.background = element_rect(fill=NA), 
        panel.border = element_rect(fill=NA),
        legend.key = element_rect(fill=NA),
        legend.background = element_rect(colour="black"),
        legend.position = "bottom")

# plot results: confidence bands
ggplot(data = res, aes(x = time)) +
  geom_line(aes(y = 100*ATE_hat)) + 
  geom_line(aes(y = 100*CB_lower, color = type), lty = 5) +
  geom_line(aes(y = 100*CB_upper, color = type), lty = 5) +
  scale_color_manual(name = "", 
                     values = c('orange',
                                'red',
                                'turquoise4',
                                'cyan3',
                                'chartreuse2'), 
                     labels = c('EBS',
                                'IF',
                                'WBS - Lin et al.',
                                'WBS - Beyersmann et al.',
                                'WBS - Weird bootstrap')) +
  scale_x_continuous(breaks=seq(0,35,5), 
                     labels = c("0","","10","","20","","30","")) +
  labs(title="Confidence Bands", 
         x="years since diagnosis", y="absolute risk difference [%]") + 
  theme(text = element_text(size=15),
        plot.title = element_text(hjust=0.5, face="bold"),
        panel.grid.minor = element_line(colour=NA),
        panel.grid.major = element_line(colour='lightgrey'),
        panel.background = element_rect(fill=NA), 
        panel.border = element_rect(fill=NA),
        legend.key = element_rect(fill=NA),
        legend.background = element_rect(colour="black"),
        legend.position = "bottom")

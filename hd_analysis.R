# This script reproduces the analysis of the hd data.

# load packages
library(survival)
library(splines)
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


## Table 3 ####
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

# break ties
set.seed(1234)
hd$time[duplicated(hd$time)] <- 
  hd$time[duplicated(hd$time)] +
  rnorm(sum(duplicated(hd$time)), 0,0.001)

## check PH assumption ####

# fit cause-specific Cox models for all causes
csc <- CSC(formula = list(Hist(time, status) ~ trtgiven+age+sex+clinstg+medwidsi+extranod,
                          Hist(time, status) ~ trtgiven+age+sex+clinstg+medwidsi+extranod),
           hd[order(hd$time),])
PH_test_relapse <- cox.zph(csc$models[[1]])
PH_test_death <- cox.zph(csc$models[[2]])

# function to plot outcome of PH Test
zph_plot <- function(cause, var){
  # input:
  ## cause: event type ('relapse', 'death')
  ## var: covariate ('trtgiven', 'age', 'sex', 'clinstg', 'medwidsi', 'extranod')
  
  # rename variable for plot title
  var_name <- switch(var, 
                     trtgiven = "Treatment",
                     age = "Age",
                     sex = "Sex", 
                     clinstg = "Lymphoma stage",
                     medwidsi = "Mediastinum involvement",
                     extranod = "Extranodal disease")
  
  # save Beta(t)
  data <- data.frame(time = eval(parse(text = paste0("PH_test_", cause, "$time"))),
                     x = eval(parse(text = paste0("PH_test_", cause, "$x"))),
                     y = eval(parse(text = paste0("PH_test_", cause, "$y[,'", var,"']"))))
  
  # determine x axis ticks
  xval <- signif(suppressWarnings(approx(data$x, 
                                         data$time, 
                                         seq(min(data$x), 
                                             max(data$x), 
                                             length=17)[2*(1:8)])$y), 2)
  
  # calculate y hat, ylow, yup
  pred <- data.frame(x = seq(from=min(data$x), to=max(data$x), length=40))
  lmat <- ns(c(pred$x, data$x), df=4, intercept=TRUE)
  pmat <- lmat[1:40,]
  xmat <- lmat[-(1:40),]
  qmat <- qr(xmat)
  bk <- backsolve(qmat$qr[1:4, 1:4], diag(4))
  seval <- ((pmat%*% (bk %*% t(bk))) * pmat) %*% rep(1, 4)
  index <- which(dimnames(eval(parse(text = paste0("PH_test_", cause, "$y"))))[[2]] == var)
  variance <- eval(parse(text = paste0("PH_test_", cause, "$var[", index, ",", index, "]")))
  pred$yhat <- (pmat %*% qr.coef(qmat, data$y))[,1]
  pred$ylow <- (pred$yhat - 2*sqrt(variance*seval))[,1]
  pred$yup <- (pred$yhat + 2*sqrt(variance*seval))[,1]
  
  # plot results
  ggplot() + 
    geom_point(data=data, aes(x=x, y=y)) + 
    geom_line(data=pred, aes(x=x, y=yhat)) +
    geom_line(data=pred, aes(x=x, y=ylow), linetype = "dashed") + 
    geom_line(data=pred, aes(x=x, y=yup), linetype = "dashed") +
    labs(title = paste0("Variable: ", var_name),
         subtitle = paste0("Schoenfeld Individual Test: p = ", 
                           format(round(
                             eval(parse(text = paste0("PH_test_", cause, "$table['", var, "','p']"))), 
                             4), nsmall = 4)),
         x = ifelse(var %in% c("medwidsi","extranod"), "Time [years]", ""), 
         y = ifelse(var %in% c("trtgiven","sex","medwidsi"), "Beta(t)", "")) +
    scale_x_continuous(breaks = suppressWarnings(approx(data$time, data$x, xval)$y), 
                       labels = sapply(xval, format)) +
    theme(text = element_text(size=14),
          plot.title = element_text(hjust=0.5, face="bold"),
          plot.subtitle = element_text(hjust=0.5),
          panel.grid.minor = element_line(colour=NA),
          panel.grid.major = element_line(colour='lightgrey'),
          panel.background = element_rect(fill=NA), 
          panel.border = element_rect(fill=NA),
          legend.key = element_rect(fill=NA),
          legend.background = element_rect(colour="black"),
          legend.position = "bottom")
  
}

p_r1 <- zph_plot("relapse","trtgiven")
p_r2 <- zph_plot("relapse","age")
p_r3 <- zph_plot("relapse","sex")
p_r4 <- zph_plot("relapse","clinstg")
p_r5 <- zph_plot("relapse","medwidsi")
p_r6 <- zph_plot("relapse","extranod")
p <- plot_grid(p_r1, p_r2, p_r3, p_r4, p_r5, p_r6, ncol = 2)
t <- ggdraw() + draw_label(paste0("Global Schoenfeld test: p = ", 
                                  format(round(PH_test_relapse$table["GLOBAL","p"], 4), nsmall=4)), 
                           fontface='bold', size=18)
png("Results/PH_relapse.png", width = 1200, height = 1500)
plot_grid(t, p, ncol = 1, rel_heights = c(0.1, 1))
dev.off()

p_d1 <- zph_plot("death","trtgiven")
p_d2 <- zph_plot("death","age")
p_d3 <- zph_plot("death","sex")
p_d4 <- zph_plot("death","clinstg")
p_d5 <- zph_plot("death","medwidsi")
p_d6 <- zph_plot("death","extranod")
p <- plot_grid(p_d1, p_d2, p_d3, p_d4, p_d5, p_d6, ncol = 2)
t <- ggdraw() + draw_label(paste0("Global Schoenfeld test: p = ", 
                                  format(round(PH_test_death$table["GLOBAL","p"], 4), nsmall=4)), 
                           fontface='bold', size=18)
png("Results/PH_death.png", width = 1200, height = 1500)
plot_grid(t, p, ncol = 1, rel_heights = c(0.1, 1))
dev.off()


# analyse data #################################################################

## Relapse ####

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
            # fit cause-specific Cox models for all causes
            invisible(capture.output(
              csc_bs <- CSC(
                formula = 
                  list(Hist(time, status) ~ trtgiven+age+sex+clinstg+medwidsi+extranod,
                       Hist(time, status) ~ trtgiven+age+sex+clinstg+medwidsi+extranod),
                hd_bs[order(hd_bs$time),])
            ))
            # calculate average treatment effect
            ATE(ID = order(hd_bs$time),
                Z = list(csc_bs$models$`Cause 1`$x,
                         csc_bs$models$`Cause 2`$x),
                index_A = c(1,1),
                event = hd_bs$status[order(hd_bs$time)],
                time = sort(hd_bs$time),
                beta = coefficients(csc_bs),
                t = 0:35,
                cause = 1, # cause 1 indicates relapse
                IF = FALSE, WBS = FALSE)$ATE
          }
if(.Platform$OS.type == "windows"){
  stopCluster(cl)
  registerDoSEQ()
}
# determine quantile for confidence bands
q_EBS_sup_r <- quantile(apply(abs(res_EBS_r - rowMeans(res_EBS_r, na.rm = TRUE)) / 
                                sqrt(apply(res_EBS_r, 1, var, na.rm = TRUE)), 
                              2, max, na.rm = TRUE), 0.95)

# perform IF/WBS
set.seed(1234, kind = "Mersenne-Twister") # for reproducibility
# determine number of observed (uncensored) events for multipliers
event_num <- sum(hd$status[order(hd$time)] != 0)
# compute weird bootstrap multipliers
mult_Weird <- sapply(1:event_num, function(j){
  time <- sort(hd$time)
  Y <- sum(time[hd$status[order(hd$time)] != 0][j] <= time)
  rbinom(1e3, Y, 1/Y) - 1
})
# calculate average treatment effect, confidence intervals & bands
## IF & WBS (Lin multipliers)
res_IF_WBS_Lin_r <- ATE(ID = order(hd$time),
                        Z = list(csc$models$`Cause 1`$x,
                                 csc$models$`Cause 2`$x),
                        index_A = c(1,1),
                        event = hd$status[order(hd$time)],
                        time = sort(hd$time),
                        beta = coefficients(csc),
                        t = 0:35,
                        cause = 1, # cause 1 indicates relapse
                        G_IF_init = rnorm(dim(hd)[1] * 1e3),
                        G_WBS_init = rnorm(1e3 * event_num),
                        bs_iter = 1e3)
## WBS (Beyersmann multipliers)
res_WBS_Bey_r <- ATE(ID = order(hd$time),
                     Z = list(csc$models$`Cause 1`$x,
                              csc$models$`Cause 2`$x),
                     index_A = c(1,1),
                     event = hd$status[order(hd$time)],
                     time = sort(hd$time),
                     beta = coefficients(csc),
                     t = 0:35,
                     cause = 1, # cause 1 indicates relapse
                     IF = FALSE,
                     G_WBS_init = rpois(1e3 * event_num, 1) - 1, 
                     bs_iter = 1e3)
## WBS (weird bootstrap multipliers)
res_WBS_Weird_r <- ATE(ID = order(hd$time),
                       Z = list(csc$models$`Cause 1`$x,
                                csc$models$`Cause 2`$x),
                       index_A = c(1,1),
                       event = hd$status[order(hd$time)],
                       time = sort(hd$time),
                       beta = coefficients(csc),
                       t = 0:35,
                       cause = 1, # cause 1 indicates relapse
                       IF = FALSE,
                       G_WBS_init = c(t(mult_Weird)),
                       bs_iter = 1e3)
# store results
res_r <- data.frame(
  time = rep(0:35, 5),
  type = rep(c("EBS",
               "IF",
               "WBS - Lin et al.",
               "WBS - Beyersmann et al.",
               "WBS - Weird bootstrap"), 
             each = 36),
  ATE_hat = rep(res_IF_WBS_Lin_r[[1]], 5),
  CI_lower = c(apply(res_EBS_r, 1, 
                     function(res_EBS_t){max(quantile(res_EBS_t, 0.025), -1)}),
               res_IF_WBS_Lin_r[[4]][1,],
               res_IF_WBS_Lin_r[[8]][1,],
               res_WBS_Bey_r[[8]][1,],
               res_WBS_Weird_r[[8]][1,]),
  CI_upper = c(apply(res_EBS_r, 1, 
                     function(res_EBS_t){min(quantile(res_EBS_t, 0.975), 1)}),
               res_IF_WBS_Lin_r[[4]][2,],
               res_IF_WBS_Lin_r[[8]][2,],
               res_WBS_Bey_r[[8]][2,],
               res_WBS_Weird_r[[8]][2,]),
  CB_lower = c(as.numeric(res_IF_WBS_Lin_r[[1]]) - 
                 q_EBS_sup_r * sqrt(apply(res_EBS_r, 1, var, na.rm = TRUE)),
               res_IF_WBS_Lin_r[[5]][1,],
               res_IF_WBS_Lin_r[[9]][1,],
               res_WBS_Bey_r[[9]][1,],
               res_WBS_Weird_r[[9]][1,]),
  CB_upper = c(as.numeric(res_IF_WBS_Lin_r[[1]]) + 
                 q_EBS_sup_r * sqrt(apply(res_EBS_r, 1, var, na.rm = TRUE)),
               res_IF_WBS_Lin_r[[5]][2,],
               res_IF_WBS_Lin_r[[9]][2,],
               res_WBS_Bey_r[[9]][2,],
               res_WBS_Weird_r[[9]][2,])
)

# create plot
p1_r <- ggplot(data = res_r, aes(x = time)) +
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
p2_r <- ggplot(data = res_r, aes(x = time)) +
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
p_r <- plot_grid(p1_r + theme(legend.position = "none"), p2_r)
png("Results/ATE_relapse.png", width = 1200, height = 500)
plot_grid(p_r, get_legend(p1_r), nrow=2, rel_heights = c(0.75, 0.2))
dev.off()

## Death ####

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
            # fit cause-specific Cox models for all causes
            invisible(capture.output(
              csc_bs <- CSC(
                formula = 
                  list(Hist(time, status) ~ trtgiven+age+sex+clinstg+medwidsi+extranod,
                       Hist(time, status) ~ trtgiven+age+sex+clinstg+medwidsi+extranod),
                hd_bs[order(hd_bs$time),])
            ))
            # calculate average treatment effect
            ATE(ID = order(hd_bs$time),
                Z = list(csc_bs$models$`Cause 1`$x,
                         csc_bs$models$`Cause 2`$x),
                index_A = c(1,1),
                event = hd_bs$status[order(hd_bs$time)],
                time = sort(hd_bs$time),
                beta = coefficients(csc_bs),
                t = 0:35,
                cause = 2, # cause 2 indicates death
                IF = FALSE, WBS = FALSE)$ATE
          }
if(.Platform$OS.type == "windows"){
  stopCluster(cl)
  registerDoSEQ()
}
# determine quantile for confidence bands
q_EBS_sup_d <- quantile(apply(abs(res_EBS_d - rowMeans(res_EBS_d, na.rm = TRUE)) / 
                                sqrt(apply(res_EBS_d, 1, var, na.rm = TRUE)), 
                              2, max, na.rm = TRUE), 0.95)

# perform IF/WBS
set.seed(1234, kind = "Mersenne-Twister") # for reproducibility
# determine number of observed (uncensored) events for multipliers
event_num <- sum(hd$status[order(hd$time)] != 0)
# compute weird bootstrap multipliers
mult_Weird <- sapply(1:event_num, function(j){
  time <- sort(hd$time)
  Y <- sum(time[hd$status[order(hd$time)] != 0][j] <= time)
  rbinom(1e3, Y, 1/Y) - 1
})
# calculate average treatment effect, confidence intervals & bands
## IF & WBS (Lin multipliers)
res_IF_WBS_Lin_d <- ATE(ID = order(hd$time),
                        Z = list(csc$models$`Cause 1`$x,
                                 csc$models$`Cause 2`$x),
                        index_A = c(1,1),
                        event = hd$status[order(hd$time)],
                        time = sort(hd$time),
                        beta = coefficients(csc),
                        t = 0:35,
                        cause = 2, # cause 2 indicates death
                        G_IF_init = rnorm(dim(hd)[1] * 1e3),
                        G_WBS_init = rnorm(1e3 * event_num),
                        bs_iter = 1e3)
## WBS (Beyersmann multipliers)
res_WBS_Bey_d <- ATE(ID = order(hd$time),
                     Z = list(csc$models$`Cause 1`$x,
                              csc$models$`Cause 2`$x),
                     index_A = c(1,1),
                     event = hd$status[order(hd$time)],
                     time = sort(hd$time),
                     beta = coefficients(csc),
                     t = 0:35,
                     cause = 2, # cause 2 indicates death
                     IF = FALSE,
                     G_WBS_init = rpois(1e3 * event_num, 1) - 1, 
                     bs_iter = 1e3)
## WBS (weird bootstrap multipliers)
res_WBS_Weird_d <- ATE(ID = order(hd$time),
                       Z = list(csc$models$`Cause 1`$x,
                                csc$models$`Cause 2`$x),
                       index_A = c(1,1),
                       event = hd$status[order(hd$time)],
                       time = sort(hd$time),
                       beta = coefficients(csc),
                       t = 0:35,
                       cause = 2, # cause 2 indicates death
                       IF = FALSE,
                       G_WBS_init = c(t(mult_Weird)),
                       bs_iter = 1e3)
# store results
res_d <- data.frame(
  time = rep(0:35, 5),
  type = rep(c("EBS",
               "IF",
               "WBS - Lin et al.",
               "WBS - Beyersmann et al.",
               "WBS - Weird bootstrap"), 
             each = 36),
  ATE_hat = rep(res_IF_WBS_Lin_d[[1]], 5),
  CI_lower = c(apply(res_EBS_d, 1, 
                     function(res_EBS_t){max(quantile(res_EBS_t, 0.025), -1)}),
               res_IF_WBS_Lin_d[[4]][1,],
               res_IF_WBS_Lin_d[[8]][1,],
               res_WBS_Bey_d[[8]][1,],
               res_WBS_Weird_d[[8]][1,]),
  CI_upper = c(apply(res_EBS_d, 1, 
                     function(res_EBS_t){min(quantile(res_EBS_t, 0.975), 1)}),
               res_IF_WBS_Lin_d[[4]][2,],
               res_IF_WBS_Lin_d[[8]][2,],
               res_WBS_Bey_d[[8]][2,],
               res_WBS_Weird_d[[8]][2,]),
  CB_lower = c(as.numeric(res_IF_WBS_Lin_d[[1]]) - 
                 q_EBS_sup_d * sqrt(apply(res_EBS_d, 1, var, na.rm = TRUE)),
               res_IF_WBS_Lin_d[[5]][1,],
               res_IF_WBS_Lin_d[[9]][1,],
               res_WBS_Bey_d[[9]][1,],
               res_WBS_Weird_d[[9]][1,]),
  CB_upper = c(as.numeric(res_IF_WBS_Lin_d[[1]]) + 
                 q_EBS_sup_d * sqrt(apply(res_EBS_d, 1, var, na.rm = TRUE)),
               res_IF_WBS_Lin_d[[5]][2,],
               res_IF_WBS_Lin_d[[9]][2,],
               res_WBS_Bey_d[[9]][2,],
               res_WBS_Weird_d[[9]][2,])
)

### Figure 8 ####

# create plot
p1_d <- ggplot(data = res_d, aes(x = time)) +
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
p2_d <- ggplot(data = res_d, aes(x = time)) +
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
p_d <- plot_grid(p1_d + theme(legend.position = "none"), p2_d)
png("Results/figure_7.png", width = 1200, height = 500)
plot_grid(p_d, get_legend(p1_d), nrow=2, rel_heights = c(0.75, 0.2))
dev.off()

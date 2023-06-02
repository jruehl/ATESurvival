\name{EBS}
\alias{EBS}
\title{Compute the ATE for data where the observations were drawn with replacement}
\usage{
EBS(Z, event, time, t, beta, index_A, factor_init = NULL, trunc_time_init = NULL, cause = 1)
}
\arguments{
\item{Z}{A list with covariate matrices for each cause (with rows representing subjects, ordered by their event times)}
\item{event}{A vector of event indicators, ordered by the subjects' event times (\code{0} indicating censoring)}
\item{time}{A (sorted) vector of event times}
\item{t}{A (sorted) vector of time points at which the data should be evaluated}
\item{beta}{A list with Cox regression coefficients for each cause}
\item{index_A}{A list with column indexes at which the treatment variable is stored in the covariate matrices for each cause, respectively}
\item{factor_init}{A vector of (multinomial) bootstrap weights, adding up to the sample size (default is \code{NULL}, indicating that the data have already been resampled)}
\item{trunc_time_init}{A vector of left-truncation times, ordered by the subjects' event times (default is \code{NULL}, indicating that there is no left-truncation)}
\item{cause}{An integer denoting the cause of interest (as specified in the \code{'event'} vector and corresponding to the list position in \code{'Z'}/\code{'beta'}/\code{'index_A'})}
}
\value{
A vector of the values of ATE at the specified time points.
}
\description{
Compute the ATE for time-to-event data with competing risks that have been drawn with replacement (i.e. in the context of bootstrapping). 
}
\examples{
\dontrun{
data(hd, package = "randomForestSRC")
hd$trtgiven <- factor(hd$trtgiven, levels = c("RT","CMT"))
hd$clinstage <- factor(hd$clinstg)
hd$medwidsi <- factor(hd$medwidsi, levels = c("N","S","L"))

library(doParallel)
registerDoParallel(cores=detectCores())
res <- 
  foreach(j = 1:1000, .combine = cbind, 
          .packages = c('riskRegression', 'ATESurvival')) %dopar% {
            hd_bs <- hd[sample(1:dim(hd)[1], dim(hd)[1], replace = TRUE),]
            csc_bs <- CSC(
              formula = 
                list(Hist(time, status) ~ trtgiven+age+sex+clinstage+medwidsi+extranod,
                     Hist(time, status) ~ trtgiven+age+sex+clinstage+medwidsi+extranod),
              hd_bs[order(hd_bs$time),])
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
rbind(apply(res, 1, function(res_EBS_t){max(quantile(res_EBS_t, 0.025), -1)}),
      apply(res, 1, function(res_EBS_t){min(quantile(res_EBS_t, 0.975), 1)}))
}
}
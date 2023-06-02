\name{ATE_IF_WBS}
\alias{ATE_IF_WBS}
\title{Compute the ATE and the associated confidence intervals/bands using influence functions and/or the wild bootstrap}
\usage{
ATE_IF_WBS(Z, event, time, t, beta, index_A, trunc_time_init = NULL, cause = 1, IF = TRUE, WBS = TRUE, bs_iter = 1e4, level = 0.95, arms = FALSE)
}
\arguments{
\item{Z}{A list with covariate matrices for each cause (with rows representing subjects, ordered by their event times)}
\item{event}{A vector of event indicators, ordered by the subjects' event times (\code{0} indicating censoring)}
\item{time}{A (sorted) vector of event times}
\item{t}{A (sorted) vector of time points at which the data should be evaluated}
\item{beta}{A list with Cox regression coefficients for each cause}
\item{index_A}{A list with column indexes at which the treatment variable is stored in the covariate matrices for each cause, respectively}
\item{trunc_time_init}{A vector of left-truncation times, ordered by the subjects' event times (default is \code{NULL}, indicating that there is no left-truncation)}
\item{cause}{An integer denoting the cause of interest (as specified in the \code{'event'} vector and corresponding to the list position in \code{'Z'}/\code{'beta'}/\code{'index_A'})}
\item{IF}{(logical) If \code{TRUE}, confidence intervals/bands are determined based on influence functions}
\item{WBS}{(logical) If \code{TRUE}, confidence intervals/bands are determined based on the wild bootstrap}
\item{bs_iter}{An integer denoting the number of bootstrap samples to use}
\item{level}{A double denoting the significance level for the confidence intervals/bands}
\item{arms}{(logical) If \code{TRUE}, the wild bootstrap confidence intervals/bands for each treatment arm are computed in addition to the intervals/bands for their difference}
}
\value{
A list including the ATE at the specified time points as first entry.
If \code{IF = TRUE}, the second element is a matrix with the lower (first row) and upper (second row) limit of the confidence interval derived using influence functions, and the third element is a matrix with the lower (first row) and upper (second row) limit of the time-simultaneous confidence band derived using influence functions.
If \code{WBS = TRUE}, the fourth element is a matrix with the limits of the confidence intervals derived using the wild bootstrap. If \code{arms = FALSE}, the first and second row include the lower and upper limit considering the standard deviation according to the asymptotic distribution of the stochastic process for the average treatment effect ('calculated' SD), the third and fourth row include the lower and upper limit considering standard normal multipliers, the fifth and sixth row include the lower and upper limit considering centered Poisson multipliers, and the seventh and eighth row include the lower and upper limit considering weird bootstrap multipliers. If \code{arms = TRUE}, the first and second, third & fourth, and fifth & sixth row include the lower and upper limit considering the calculated SD for A = 1, A = 0, and for their difference, respectively. This order continues in the subsequent rows, with rows 7-12 referring to standard normal multipliers, rows 13-18 covering centered Poisson multipliers, and rows 19-24 being associated with weird bootstrap multipliers.
Besides, if \code{WBS = TRUE}, the fifth element is a matrix with the limits of the time-simultaneous confidence bands derived using the wild bootstrap. If \code{arms = FALSE}, the first and second row include the lower and upper limit considering the calculated SD but quantiles based on standard normal multipliers, the third and fourth row include the lower and upper limit considering standard normal multipliers to compute both the SD and quantiles, the fifth and sixth row include the lower and upper limit considering the calculated SD but quantiles based on centered Poisson multipliers, the seventh and eighth row include the lower and upper limit considering centered Poisson multipliers to compute both the SD and quantiles, the nineth and tenth row include the lower and upper limit considering the calculated SD but quantiles based on weird boostrap multipliers, and the eleventh and twelveth row include the lower and upper limit considering weird bootstrap multipliers to compute both the SD and quantiles. If \code{arms = TRUE}, the first and second row include the lower and upper limit considering the calculated SD but quantiles based on standard normal multipliers for A = 1, the third and fourth row include the lower and upper limit considering standard normal multipliers to compute both the SD and quantiles for A = 1, the fifth and sixth row include the lower and upper limit considering the calculated SD but quantiles based on standard normal multipliers for A = 0, the seventh and eighth row include the lower and upper limit considering standard normal multipliers to compute both the SD and quantiles for A = 0, the nineth and tenth row include the lower and upper limit considering the calculated SD but quantiles based on standard normal multipliers for the difference between A = 1 and A = 0, the eleventh and twelveth row include the lower and upper limit considering standard normal multipliers to compute both the SD and quantiles for the difference between A = 1 and A = 0. This order continues in the subsequent rows, with rows 13-24 referring to centered Poisson multipliers and rows 25-36 covering weird bootstrap multipliers.
}
\description{
Compute the ATE and the associated confidence intervals/bands for time-to-event data with competing risks, using influence functions and/or the wild bootstrap. 
}
\examples{
\dontrun{
data(hd, package = "randomForestSRC")
hd$trtgiven <- factor(hd$trtgiven, levels = c("RT","CMT"))
hd$clinstage <- factor(hd$clinstg)
hd$medwidsi <- factor(hd$medwidsi, levels = c("N","S","L"))

library(riskRegression)
csc <- CSC(formula = list(Hist(time, status) ~ trtgiven+age+sex+clinstage+medwidsi+extranod,
                          Hist(time, status) ~ trtgiven+age+sex+clinstage+medwidsi+extranod),
           hd[order(hd$time),])

res <- ATE_IF_WBS(
  Z = list(csc$models$`Cause 1`$x,
           csc$models$`Cause 2`$x),
  event = hd$status[order(hd$time)], 
  time = sort(hd$time),
  t = 0:35,
  beta = coefficients(csc),
  index_A = c(1,1),
  cause = 2
)
}
}
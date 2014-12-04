# Implementing Method Proposed on November 6th

# November 7, 2014

#rm(list=ls())
set.seed(5)
source("functions.R")
library(waveslim)

update_tau <- function(nvec, theta_p, theta_m, tstat, dist_ic, ...){
  # Input:
  # rlength: T, time
  # nvec: new vector recieved at time T
  # theta - 2 x (T-2) matrix containing old theta estimates
  rlength <- length(theta_p) +1 
  theta_p <- append(theta_p, 0)
  theta_m <- append(theta_m, 0)
  nvec_null <- tstat(nvec, dist_ic, ...)
  #print(nvec_null)
  theta_m[rlength] <- theta_m[(rlength-1)] + theta_p[(rlength-1)]
  theta_p <- theta_p + nvec_null
  return(cbind(theta_p,theta_m))
}
# theta_p
x <- seq(.1,.5,length.out=5)
#theta_m
y <- 1:5
#nvec
t_obs <- rnorm(50)

theta <- update_tau(nvec=t_obs, theta_p=x, theta_m=y, tstat=wave.den, dist_ic="pnorm")
theta
find_max_dif(theta)


Find_RL_Fast(num.samp=30, dist="norm", params=list(mean=0, sd=1), tstat=wave.den, UCL=c(.1,.2))
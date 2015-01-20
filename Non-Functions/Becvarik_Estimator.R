# Rachel Becvarik Estimator

# wave.bec

# rm(list=ls())
library(waveslim)
source("functions.R")
x <- rnorm(100)
y <- rnorm(100)
wave.bec(x,y)





x <- rnorm(64)
y <- rnorm(64,0,2)
tran <- wave.bec(x,y)
tran

x <- rnorm(64,0,2)
tran <- wave.bec(x,"qnorm")
tran

x <- rnorm(382,0,2)
tran <- wave.bec(x,"qnorm")
tran


x <- rnorm(382,0,2)
tran <- wave.bec(x,"qnorm")
tran2 <- wave.energy(x,"pnorm")
tran
tran2
x <- matrix(rnorm(382*1000,0,2), nrow=382, ncol=1000)
becvarik_big <- mean(apply(x, 2, function(x) wave.bec(x, "qnorm")))
energy_curve_big <- mean(apply(x, 2, function(x) wave.energy(x, "pnorm")))


x1 <- rnorm(64,0,2)
tran <- wave.bec(x1,"qnorm")
tran2 <- wave.energy(x1,"pnorm")
tran
tran2
x1 <- matrix(rnorm(64*1000,0,2), nrow=64, ncol=1000)
becvarik_small <- mean(apply(x1, 2, function(x) wave.bec(x, "qnorm")))
energy_curve_small <- mean(apply(x1, 2, function(x) wave.energy(x, "pnorm")))


becvarik_small
energy_curve_small
becvarik_big
energy_curve_big


Find_IC_RL_Slow(num.samp=64, tstat=wave.bec, UCL=80)
# Monday, January 12, 2015

# GOAL: Get ARL's for the (200) for the Becvarik Estimator

#Note: wave.bec had reduce=2

rm(list=ls())
source("functions.R")
set.seed(502)
library(waveslim)
library(doMC)

######## Finding 200 ARL for Slow Method, Under N(0,1)

UCL <- seq(0,4.3, length.out=100)
Method_Name <- "Slow Method, N(0,1), Bec Est"
Cores <- 6

registerDoMC(Cores)
ARL_Curve_SA <- times(Cores) %dopar% ARL_Proc(UCL=UCL, time=8*60*60, 
                               method=Find_IC_RL_Slow, 
                               tstat=wave.bec, 
                               dist="qnorm",
                               params=list(mean=0, sd=1))
save.image(file="ARL_Baseline_1.RData")
q(save="no")


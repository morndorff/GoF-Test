# Monday, January 12, 2015

# GOAL: Get ARL's for the (200) for the Becvarik Estimator

#Note: wave.bec had reduce=2

rm(list=ls())
source("functions.R")
set.seed(502)
library(waveslim)

ARL_Matrices <- vector(mode="list")
ARL_Functions <- vector(mode="list")
UCL_200 <- vector(mode="list")
######## Finding 200 ARL for Slow Method, Under N(0,1)

UCL <- seq(0,4.3, length.out=100)
Method_Name <- "Slow Method, N(0,1), Bec Est"
ARL_Curve_SA <- ARL_Proc(UCL=UCL, time=12*60*60, 
                         method=Find_IC_RL_Slow, 
                         tstat=wave.bec, 
                         dist="qnorm",
                         params=list(mean=0, sd=1))

### Automatic Finding ##
ARL_Matrices[[Method_Name]] <- ARL_Curve_SA[[1]]
ARL_Functions[[Method_Name]] <- approxfun(as.numeric(ARL_Matrices[[Method_Name]][1,]), UCL)
UCL_200[[Method_Name]] <- ARL_Functions[[Method_Name]](200)
save.image(file="ARL_Baseline_2.RData")
#q(save="no")


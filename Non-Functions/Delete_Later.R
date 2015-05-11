# Yet Another Test Script

rm(list=ls())
source("functions.R")
set.seed(50)
library(waveslim)

ARL_Matrices <- vector(mode="list")
ARL_Functions <- vector(mode="list")
UCL_200 <- vector(mode="list")
######## Finding 200 ARL for Windowed Method, Under N(0,1)

UCL <- seq(0,6, length.out=100)
Method_Name <- "Slow Method, N(0,1), WSize=15"
ARL_Curve_SA <- ARL_Proc(UCL=UCL, time=60, 
                         method=Find_IC_RL_Slow, 
                         tstat=wave.bec, 
                         dist="qnorm",
                         params=list(mean=0, sd=1),
                         num.samp=64)

# While I'm at it, lets fix all the possible problems
# Fixed IC_RL_Slow first
# TOfix: IC_RL_Fast second


UCL <- seq(0,5, length.out=100)
Method_Name <- "Fast Method, N(0,1), WSize=15"
ARL_Curve_SA <- ARL_Proc(UCL=UCL, time=8, 
                         method=Find_IC_RL_Fast, 
                         tstat=wave.bec, 
                         dist="qnorm",
                         params=list(mean=0, sd=1),
                         num.samp=64)


UCL <- seq(0,4, length.out=100)
Method_Name <- "Fast Method, N(0,1), WSize=15"
ARL_Curve_SA <- ARL_Proc(UCL=UCL, time=8, 
                         method=Find_IC_RL_Windowed, 
                         tstat=wave.bec, 
                         dist="qnorm",
                         params=list(mean=0, sd=1),
                         num.samp=64)






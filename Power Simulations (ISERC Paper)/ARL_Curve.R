# ARL Simulations for ISERC 2015 Paper

# rm(list=ls())
library(waveslim)
set.seed(50)
source("functions.R")

# Initial Code: December 31, 2015


UCL <- seq(.01,.09, length.out=100)
ARL_Curve_SA <- ARL_Proc(UCL=UCL, time=60*60*8, 
                         method=Find_IC_RL_Slow, 
                         tstat=wave.energy, 
                         dist="norm",
                         params=list(mean=0, sd=1))
ARL_Curve_SA_Mat <- ARL_Curve_SA[[1]] # Sample Average ARL
ARL <- as.numeric(ARL_Curve_SA_Mat[1,])
f <- approxfun(ARL,UCL)
ARL_200 <- f(200)

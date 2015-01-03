# Some ARL simulation code:

# rm(list=ls())
source("functions.R")
set.seed(50)
library(waveslim)

ARL_Matrices <- vector(mode="list")
ARL_Functions <- vector(mode="list")
ARL_200 <- vector(mode="list")
######## Finding 200 ARL for Windowed Method, Under N(0,1)

UCL <- seq(.01,.08, length.out=100)
Method_Name <- "Windowed Method, N(0,1), WSize=15"
ARL_Curve_SA <- ARL_Proc(UCL=UCL, time=60*60, 
                         method=Find_IC_RL_Windowed, 
                         tstat=wave.energy, 
                         dist="norm",
                         params=list(mean=0, sd=1),
                         WSize=15)

### Automatic Finding ##
ARL_Matrices[[Method_Name]] <- ARL_Curve_SA[[1]]
ARL_Functions[[Method_Name]] <- approxfun(as.numeric(ARL_Matrices[[Method_Name]][1,]))
ARL_200[[Method_Name]] <- ARL_Functions[[Method_Name]](200)
save.image(file="ARL_Baseline.RData")




  
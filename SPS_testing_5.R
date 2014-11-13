# SPS function testing, take 5

# November 11, 2014
# rm(list=ls())
source("functions.R")
set.seed(5)
UCL <- seq(.1, .4, length.out=10)

test <- Find_IC_RL_Fast(num.samp=30, dist="norm", 
                     params=list(mean=0, sd=1), 
                     tstat=wave.den, 
                     UCL=UCL,
                     detail=TRUE)
plot(test[[1]])

test2 <- ARL_Proc(UCL=c(.25, .75, 1), method=Find_IC_RL_Fast)
test2[[1]]

test3 <- Sim_CP_Process_Iter(num.samp=30, dist_one="norm", param_one=list(mean=0, sd=1),
dist_two="norm", param_two=list(mean=0, sd=2), cp=50)

UCL <- seq(.1, 16, length.out=10)

test3 <- Find_CP_RL_Fast(num.samp=30, dist_one="norm", param_one=list(mean=0, sd=1), 
                         dist_two="norm", param_two=list(mean=0, sd=2), cp=100,
                         tstat=wave.den, 
                         UCL=UCL,
                         detail=TRUE)
plot(test3[[1]])

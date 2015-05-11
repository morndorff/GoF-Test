# Testing Out of Control ARL Code
# Friday, January 16, 2015

# rm(list=ls())
source("functions.R")
library(wavethresh)
library(waveslim)

set.seed(1000)

Find_ARL_OOC(tstat=wave.bec,UCL=4.2, method=Find_CP_RL_Slow, dist_one="qnorm", dist_two="qnorm", param_two=list(mean=0, sd=2))



# Find_ARL_OOC <- function(num.samp=32, dist_one="norm", param_one=list(mean=0, sd=1),
#                          dist_two="norm", param_two=list(mean=0, sd=2), cp=1,
#                          tstat=wave.energy, UCL, time = 30, method=Find_CP_RL_Fast, ...)
  
Find_CP_RL_Slow
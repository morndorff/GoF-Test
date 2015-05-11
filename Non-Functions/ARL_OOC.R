# Finding ARL's for various out of control distributions
# For details, see Becvarik, page 36
# In Control distribution is N(0,1)

# rm(list=ls())
load("ARL_Baseline.RData")
source("functions.R")
set.seed(50)

##### Out of Control Distribution ####

# N(0,2)
test <- Find_CP_RL_Windowed(UCL=.06, num.samp=64)
plot(test[[1]])
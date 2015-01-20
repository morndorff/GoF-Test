# SPS Testing, Take 6

# November 14, 2014
# Find RLs for 'fast' process
rm(list=ls())
source("functions.R")
library(waveslim)
set.seed(5)

Plots <- function (test, weight="Biased") {
  plot(test[[1]], main=paste("Plot of Test Statistic h(t), Weight=",weight), type="l") # Weight = Biased Towards Extreme Obs?
  hist(test[[7]], main="Tau Chosen (Raw)") # Which Tau is Chosen?
  hist(test[[8]], main="Tau Chosen (Chosen/Time)", bin=20) # Tau Chosen Divided by Time (Scaled 0-1)
  d <- density(test[[8]], adjust=.5)
  plot(d, , "gaussian", xlim=c(0,1), main="Tau/T Chosen at Each Step") #Kernel Density Estimate
  quantile(test[[8]], probs=seq(0,1,.1)) # Quantiles of Tau/T
}
# Biased Weighting Scheme
UCL <- seq(.1, .3, length.out=10)
set.seed(5)
test <- Find_IC_RL_Fast(num.samp=30, dist="norm", 
                        params=list(mean=2,sd=1), # mean=0, sd=1), 
                        tstat=wave.energy, 
                        UCL=UCL,
                        detail=TRUE, weight=TRUE)


Plots(test, weight="Biased")
# Average Weighting Scheme
set.seed(5)
UCL <- seq(.01, .06, length.out=10)

test <- Find_IC_RL_Fast(num.samp=30, dist="norm", 
                        params=list(mean=0, sd=1), 
                        tstat=wave.energy, 
                        UCL=UCL,
                        detail=TRUE, weight=FALSE)
Plots(test, weight="Average")



# 
# Testing and Creating Functions for Statistical Process Control
rm(list=ls())
source("functions.R")
library(waveslim)
set.seed(5)
Process <- Sim_IC_Process(num.samp=20, run.length=40, dist="norm", param=list(mean=0, sd=1))

Track_Stat(Process, mean, doplot=TRUE)
testfun <- function(z){
  x <- runif(1)
  if(x > .5) y <- c(1)
  if(x < .5) y <- c(1, 1)
  y
}
testproc <- Sim_CP_Process(num.samp=20, run.length=40, 
                           dist_one="norm", param_one=list(mean=0, sd=1), 
                           dist_two="norm", param_two=list(mean=0, sd=2), bpoint=20)
a <- Track_Stat_Over(testproc, tstat=wave.den)
a

testproc <- Sim_CP_Process(num.samp=50, run.length=40,
                           dist_one="norm", param_one=list(mean=0, sd=1), 
                           dist_two="norm", param_two=list(mean=0, sd=2), bpoint=20)
a <- Track_Stat_Over(testproc, tstat=wave.den)
a

a2 <- Estimate_Tau(proc=testproc, tstat=wave.den, dist_ic=pnorm, mean=0, sd=1)
a2
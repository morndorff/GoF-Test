# November 2, 2014

# Statistical Process Control, Take 2
# rm(list=ls())
source("functions.R")
set.seed(5)
library(waveslim)

icproc <- Sim_IC_Process(num.samp=20, run.length=200, dist="norm", param=list(mean=0, sd=1))

# For a particular t, output the test statistic
Track_Stat_Over2(icproc, tstat=wave.den, dist_ic="pnorm", mean=0, sd=1)
Track_Stat_Over2(icproc[,1], tstat=wave.den, dist_ic="pnorm", mean=0, sd=1)

Test_Stat <- vector(mode="numeric", length=dim(icproc)[2])
for(i in 1:dim(icproc)[2]){
  Test_Stat[i] <- Process_Stat(proc=icproc[,1:i], tstat=wave.den, dist_ic="pnorm", mean=0, sd=1)
}


plot(Test_Stat)

oocproc <- Sim_CP_Process(num.samp=50, run.length=200,
                           dist_one="norm", param_one=list(mean=0, sd=1), 
                           dist_two="norm", param_two=list(mean=0, sd=2), bpoint=100)
Test_Stat <- vector(mode="numeric", length=dim(oocproc)[2])
for(i in 1:dim(oocproc)[2]){
  Test_Stat[i] <- Process_Stat(proc=oocproc[,1:i], tstat=wave.den, dist_ic="pnorm", mean=0, sd=1)
}

plot(Test_Stat)


x <- 1:dim(icproc)[2]
plot(x,Test_Stat)
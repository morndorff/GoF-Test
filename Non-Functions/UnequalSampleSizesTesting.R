# December 5th, 2014

# Unequal Sample Sizes Assessment
# rm(list=ls())
source("functions.R")
library(waveslim)
set.seed(50)

lenproc <- 200
Reps <- 60
Results <- matrix(, ncol=(lenproc-1), nrow=Reps) #rows= replicates
for(i in 1:Reps){
  Proc <- Sim_IC_Process(num.samp=32, run.length=lenproc, dist="norm", params=list(mean=0, sd=1))
  for(j in 1:(lenproc-1)){
    Results[i, j] <- wave.den(Proc[, 1:j], Proc[, (j+1):lenproc])
  }
}

AvgVals <- colMeans(Results)
pdf('WDif.pdf')
WDif <- plot(AvgVals, main="Average Values for Wavelet Difference Test")
dev.off


lenproc <- 200
Reps <- 60
Results <- matrix(, ncol=(lenproc-1), nrow=Reps) #rows= replicates
for(i in 1:Reps){
  Proc <- Sim_IC_Process(num.samp=32, run.length=lenproc, dist="norm", params=list(mean=0, sd=1))
  for(j in 1:(lenproc-1)){
    Results[i, j] <- wave.curve(Proc[, 1:j], Proc[, (j+1):lenproc])
  }
}

AvgVals <- colMeans(Results)
pdf('WCurve.pdf')
WCurve <- plot(AvgVals, main="Average Values for Wavelet Power Curve Test")
dev.off



lenproc <- 200
Reps <- 60
Results <- matrix(, ncol=(lenproc-1), nrow=Reps) #rows= replicates
for(i in 1:Reps){
  Proc <- Sim_IC_Process(num.samp=32, run.length=lenproc, dist="norm", params=list(mean=0, sd=1))
  for(j in 1:(lenproc-1)){
    Results[i, j] <- ks.res.simp(Proc[, 1:j], Proc[, (j+1):lenproc])
  }
}

AvgVals <- colMeans(Results)
pdf('KS.pdf')
Ks <- plot(AvgVals, main="Average Values for Kolmogorov-Smirnov Test")
dev.off

# Creating Multicore Functions

# rm(list=ls())
set.seed(500)
source("functions.R")

library(doMC)
registerDoMC(2)
foreach(i=1:3) %dopar% sqrt(i)

x <- iris[which(iris[,5] != "setosa"), c(1,5)]
trials <- 10000
ptime <- system.time({
  r <- foreach(icount(trials), .combine=cbind) %dopar% {
    ind <- sample(100, 100, replace=TRUE)
    result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
    coefficients(result1)
    }
  })[3]
ptime


stime <- system.time({
 r <- foreach(icount(trials), .combine=cbind) %do% {
 ind <- sample(100, 100, replace=TRUE)
 result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
 coefficients(result1)
 }
 })[3]
stime

UCL <- 3
ptm <- proc.time()
x <- times(2) %dopar% ARL_Proc(UCL=UCL, time=15, 
                  method=Find_IC_RL_Slow, 
                  tstat=wave.bec, 
                  dist="qnorm",
                  params=list(mean=0, sd=1))
etime <- proc.time()- ptm
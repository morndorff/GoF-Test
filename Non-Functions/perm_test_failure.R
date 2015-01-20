# Making a figure for explaining permutation test failure 
# 
rm(list=ls())
source("functions.R")
set.seed(11)
x <- rnorm(100,0,sqrt(3))
y <- rt(100,3)
perm.test(x,y,f=Max_Quan_TS, doplot=TRUE)
Max_Quan_TS(x,y,do.plot=TRUE)
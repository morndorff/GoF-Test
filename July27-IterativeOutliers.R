# iterative outlier function
# July 27, 2014

# rm(list=ls())

set.seed(5)
source("functions.R")
x <- rnorm(100)
y <- rnorm(100)

myts.out(x,y)
myts(x,y)

perm.test.out(x,y,f=myts.out)

x <- rnorm(100,0,sqrt(3))
y <-rt(100,3)

# Testing (Official)

x <- c(1,2,3,4,5,6,7,8,9,100)
y <- c(1,2,3,4,5,6,7,8,9,10)
new.perm.test.out(x,y,f=myts.out)
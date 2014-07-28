# Bug Testing the new.perm.test function for 
# two samples. Will be obsolete later.

# July 27, 2014

# rm(list=ls())

set.seed(5)
source("functions.R")

x <- rnorm(100)
y <- rnorm(100)

myts(x,y)
new.perm.test(x,y, f=myts)

perm.test(x,y, f=myts)

# seems to work now

x <- rnorm(100)
y <- rnorm(101)

myts(x,y)
new.perm.test(x,y, f=myts)

perm.test(x,y, f=myts)
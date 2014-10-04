rm(list=ls())
source("functions.R")
source("functions_testing.R")
set.seed(5)

# Checking speed
library(microbenchmark)
x <- rnorm(20)
y <- rnorm(20)
microbenchmark(perm.test(x,y,f=wave.energy), perm.test2(x,y,f=wave.energy), times=15)

Rprof(filename="test.out")
for(i in 1:900) wave.energy(x,y)
Rprof(NULL)
summaryRprof(filename="test.out")
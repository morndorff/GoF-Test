# Testing the wavelet functions 
# September 24, 2014
source("functions.R")
library(waveslim)
# source("functions_wavelets.R")

# NOTE: changed function to use ks.test output
# More accurate than previous way.

# wave.den appears to be working correctly
# permutation test also appears to be working correctly

# Changing functions to my schema
set.seed(5)
x <- rnorm(100)
y <- rnorm(100)
testvec <- wave.den(x,y)
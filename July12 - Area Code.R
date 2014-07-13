# July 12, 2014
# Goals: Make diag.plot better
# Goals: Get that damn trapezoid code working

source("functions.R")
set.seed(5)
x <- rnorm(100)

myts.par(x, qt, 3)

# Updated myts.par to use quantiles based on null hypothesisized distribution

# Now I want to make myts.par compatible with perm.test if I can
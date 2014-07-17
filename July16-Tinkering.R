# Tinkering with the new perm.test function
source("functions.R")
x <- rnorm(100)

new.perm.test(x, qnorm, f = myts.par)

new.perm.test(x, "qnorm", f = myts.par)

new.perm.test(x, "qnorm", f = ks.res)
# Doesn't work for now

y <- rnorm(100)

perm.test(x, y, f = ks.res)

# Current problem: I've managed to screw up perm.test for ks.res 

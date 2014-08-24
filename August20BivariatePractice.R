# August 20, 2014
# Doing some testing using the new Bivariate PIT test
# Can be archived later

source("functions.R")
x <- rnorm(100)
y <- runif(100,0,1)
a <- Bi_Var_PIT_ks(x,y)
#a1 <- c(unlist(a))
#a1 <- sort(a1)

# This will be EXACTLY uniform, by calculation
# Want to look at U_x and U_y

ks.test(a[[1]],punif,0,1)
ks.test(a[[2]],punif,0,1)


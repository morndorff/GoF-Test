# Making Bivariate Probability Integral Transform Test

# rm(list=ls())
library(ggplot2)
set.seed(5)


# First, lets make a univariate PIT:

# Generate Sample:
x <- rt(100,3)
U <- pt(x,3)
qplot(x=1,y=U, alpha=.5, size=5)

# That was easy enough.

Make_Inter_QFun <- function(x,y,interp=4){
  z <-c(x,y)
  lenz <- length(z)
  z1 <- seq(1/(lenz+1),lenz/(lenz+1), length.out=lenz)
  q1 <- quantile(z, probs = z1, type = interp)
  q_inter <- approxfun(z1, q1, yleft = min(q1), yright = max(q1))
  return(q_inter)
}

Make_Inter_CDF <- function(x,y,interp=4){
  z <-c(x,y)
  z <- sort(z)
  lenz <- length(z)
  z1 <- seq(1/(lenz+1),lenz/(lenz+1), length.out=lenz)
  #q1 <- quantile(z, probs = z1, type = interp)
  #q_inter <- approxfun(z1, q1, yleft = min(q1), yright = max(q1))
  cdf_inter <- approxfun(z,z1, yleft=0, yright=1)
  return(cdf_inter)
}
test <- Make_Inter_CDF(x,y)
U_x <- test(x)
qplot(y=U_x,x=1)

Bi_Var_PIT <- function(x,y){
  #com_QFun <- Make_Inter_QFun(x,y)
  inter_CDF <- Make_Inter_CDF(x,y)
  U_x <- inter_CDF(x)
  U_y <- inter_CDF(y)
  return(list(U_x,U_y))
}

a <- Bi_Var_PIT(x,y)
a1 <- c(unlist(a))
a1 <- sort(a1)

# This will be EXACTLY uniform, by calculation
# Want to look at U_x and U_y

ks.test(a[[1]],punif,0,1)
ks.test(a[[2]],punif,0,1)

# of course 
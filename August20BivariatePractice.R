# August 20, 2014
# Doing some testing using the new Bivariate PIT test

source("functions.R")
x <- rnorm(100)
y <- runif(100,0,1)
a <- Bi_Var_PIT(x,y)
#a1 <- c(unlist(a))
#a1 <- sort(a1)

# This will be EXACTLY uniform, by calculation
# Want to look at U_x and U_y

ks.test(a[[1]],punif,0,1)
ks.test(a[[2]],punif,0,1)

Bi_Var_PIT_ks <- function(x,y, alpha){
  a <- Bi_Var_PIT(x,y)
  p1 <- ks.test(a[[1]],punif,0,1)$p.value
  p2 <- ks.test(a[[2]],punif,0,1)$p.value
  FWER <- alpha/2
  pvals <- c(p1,p2)
  if(p1 < FWER & p2 < FWER){
    result <- "REJECT NULL"
  }else{
    results <- "FAIL TO REJECT NULL"
  } 
}

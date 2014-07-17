# Some updates, since this is getting complicated quick
# Note: I'm definitely breaking some old code
# new.perm.test(testx,y=paste("q",dist[[i]],sep=""),param[[i]],f=myts)
# This works
# now to make it work with power.res

# new.perm.test(testx,"qunif",list(min=0,max=2),f=myts.par)
# Correct way to use new.perm.test

# Thoughts since last night:

# Perhaps I should have two seperate permutation tests, one to use for 
# power simulation studies, and one to use for general audience stuff

# Bencharmking studies
set.seed(5)
testx <- runif(8,min=0,max=1)

library(microbenchmark)
compare <- microbenchmark(perm.test(testx,qunif,0,2,f=myts),new.perm.test(testx,"qunif",list(min=0,max=2),f=myts),times=200)
compare
# Unit: milliseconds
# expr
# perm.test(testx, qunif, 0, 2, f = myts)
# new.perm.test(testx, "qunif", list(min = 0, max = 2), f = myts)
# min       lq   median       uq      max neval
# 480.0719 483.8697 491.0205 499.8999 704.7425   200
# 278.4365 281.3906 285.1507 289.1433 363.8789   200

Rprof("out.out")
for (i in 1:200) pos = new.perm.test(testx,"qunif",list(min=0,max=2),f=myts)
Rprof(NULL)
summaryRprof("out.out")

Rprof("oldperm")
for (i in 1:200) pos = perm.test(testx,qunif,0,2,f=myts)
Rprof(NULL)
summaryRprof("oldperm")

# Seeing if a larger vector changes anything
testx2 <- runif(100,min=0,max=2)
compare <- microbenchmark(perm.test(testx2,qunif,0,2,f=myts),new.perm.test(testx2,"qunif",list(min=0,max=2),f=myts),times=100)
compare
# Unit: milliseconds
# expr
# perm.test(testx, qunif, 0, 2, f = myts)
# new.perm.test(testx, "qunif", list(min = 0, max = 2), f = myts)
# min       lq   median       uq      max neval
# 480.0719 483.8697 491.0205 499.8999 704.7425   200
# 278.4365 281.3906 285.1507 289.1433 363.8789   200

Rprof("out.out")
for (i in 1:100) pos = new.perm.test(testx2,"qunif",list(min=0,max=2),f=myts)
Rprof(NULL)
summaryRprof("out.out")

Rprof("oldperm")
for (i in 1:100) pos = perm.test(testx2,qunif,0,2,f=myts)
Rprof(NULL)
summaryRprof("oldperm")


testx3 <- runif(100,min=0,max=2)
testy3 <- runif(100,min=0,max=1)
compare <- microbenchmark(perm.test(testx3,testy3,f=myts),new.perm.test(testx3,testy3,f=myts),times=100)
compare

# as expected, these have near identical run times
# Difference appears only in handling one sample arguments

# Switching topics, back to the ever popular 'take a function', get a string
# topic

funtoli <- function(fun){
  a <- list(fun)
  names(a) <- as.character(substitute(fun))
  return(a)
}

# Doing some experimentation with ellipsis 
funex <- function(fun,...){
  a <- list(...)
  return(a)
}
d <- funex(qnorm,0,1)
d

funalt <- function(fun,funopts){
  a <- funopts
  return(a)
}
d <- funex(qnorm, min=0,max=1)
d

do.call(runif,c(list(n=100),d))
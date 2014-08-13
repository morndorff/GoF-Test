# August 13 - Iterative Outliers
# The final battle

# testing for bugs in outlier code
source("functions.R")
set.seed(5)
testx <- 1:60
testy <- 1:60
perm.test.out.iter(testx,testy, f=myts.out)
# Number of outliers expected: 0

testx <- c(1:59,100)
testy <- c(1:60)
perm.test.out.iter(testx,testy, f=myts.out)
# Number of outliers expected: 1

testx <- c(1:59, 100)
testy <- c(1:61)
perm.test.out.iter(testx,testy, f=myts.out)
# Number of outliers expected: 1

testx <- c(1:58, 100, 200)
testy <- c(1:60)
perm.test.out.iter(testx,testy, f=myts.out)
# Number of outliers expected: 2

testx <- c(1:58, 100, 200)
testy <- c(1:61)
perm.test.out.iter(testx,testy, f=myts.out)
# Number of outliers expected: 2

testx <- rt(100,3)
testy <- rnorm(100,0,sqrt(3))
perm.test.out.iter(testx,testy, f=myts.out)
# Number of outliers expected: 2

myts(testx,testy,do.plot=TRUE)

# Delete Later
# Want to use replace new.perm.test with perm.test, but there's work
# to be done

perm.test(testx,testy,f=myts)
new.perm.test(testx,testy,f=myts)

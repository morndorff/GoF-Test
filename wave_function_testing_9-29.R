# Wavelet Function Testing, September 29, 2014
# rm(list=ls())
library(waveslim)
source("functions.R")
set.seed(5)
x <- rnorm(50)
y <- rnorm(50)
#### Using different wavelet basis --------
wave.den(x, y)
wave.den(x, y, wf="la8")
perm.test(x,y, f=wave.den, fops=list(wf="la8"))
perm.test(x,y, f=wave.den, fops=list(wf="haar"))

y <- rnorm(50,0,2)
perm.test(x,y, f=wave.den, fops=list(wf="la8"))
perm.test(x,y, f=wave.den, fops=list(wf="haar"))

x <- matrix(rnorm(200), 10, 20)
y <- matrix(rnorm(200,0,2), 10, 20)
haar <- power.res.twosamp(x,y, f="wave.den", fops=list(wf="la8"))
la8 <- power.res.twosamp(x,y, f="wave.den", fops=list(wf="haar"))
haar > la8

# seems like quite a bit of variability depending on which wavelet basis you use
cor(haar, la8)

haar_curve <- power.res.twosamp(x,y, f="wave.energy", fops=list(wf="la8"))
la8_curve <- power.res.twosamp(x,y, f="wave.energy", fops=list(wf="haar"))

cor(haar_curve, la8_curve) #yikes

x <- matrix(rnorm(400), 10, 40)
y <- matrix(rnorm(400,0,2), 10, 40)
haar <- power.res.twosamp(x,y, f="wave.den", fops=list(wf="la8"))
la8 <- power.res.twosamp(x,y, f="wave.den", fops=list(wf="haar"))
haar > la8

# seems like quite a bit of variability depending on which wavelet basis you use
cor(haar, la8)

haar_curve <- power.res.twosamp(x,y, f="wave.energy", fops=list(wf="la8"))
la8_curve <- power.res.twosamp(x,y, f="wave.energy", fops=list(wf="haar"))

cor(haar_curve, la8_curve) #yikes

haar_curve_sum <- power.res.twosamp(x,y, f="wave.energy", fops=list(wf="la8", opt="sum"))
la8_curve_sum<- power.res.twosamp(x,y, f="wave.energy", fops=list(wf="haar", opt="sum"))





  
####---- Reconstruction Code ------


# Lets apply this to our data
x <- rnorm(50)
F.x <- ecdf(x)
n <- 2^5
z <- seq(min(x), max(x), length.out=n)
F.dwt <- dwt(F.x(z), wf="haar", n.levels=log(n, 2))
F.dwt_4 <- dwt(F.x(z), wf="haar", n.levels=4)
# Reconstruction
wd.F <- wd(F.x(z), 
           filter.number=1, 
           family="DaubExPhase", 
           bc="periodic")
wd.thresh <- threshold(wd.F, 
                  levels=3:(nlevelsWT(wd.F)- 1),
                  #type="hard", 
                  policy="mannum", 
                  by.level=FALSE, 
                  value=30, 
                  dev=var, 
                  boundary=FALSE,     
                  verbose = TRUE, 
                  return.threshold=F)
F.thresh<- wr(wd.thresh)
plot(F.x)
mtext("N=4, Threshold=7.32032")
lines(z, F.thresh, col= "violetred" , lwd=2,type="l")


# Changing 'n' (number of sample points) in wave.den
# Number is currently 32, regardless of sample size. Seems like it should be smallest
# diadic number
# Can use log2

floor(log2(33))


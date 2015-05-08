# Testing Eubank Tests
source("functions.R")
x <- rnorm(50)
fourier_test(x,pnorm, lenj=200)
cvm.res(x,pnorm)
smooth_test(x,pnorm,M=5)
smooth_test2(x,pnorm,lambda=2)




x <- rnorm(50)
x_test <- fourier_test(x,pnorm)
x_test
x_test1 <- fourier_test(x,pnorm, lenj=1)
x_test1
x_test2 <- cvm.res(x,pnorm)
x_test2




perm.test(x, pnorm,f="fourier_test", num.perm=2000)
perm.test(x, pnorm, distops=list(mean=0, sd=2),f="fourier_test", num.perm=2000)
perm.test(x, pnorm, distops=list(mean=0, sd=2),f="cvm.res", num.perm=2000)

x_mat <- matrix(rnorm(5000), nrow=100, ncol=50)

cvm <- apply(x_mat, 1, cvm.res, "pnorm")
fourier <- apply(x_mat, 1, fourier_test, y="pnorm", lenj=500)

plot(cvm, type="l")
lines(fourier, col="red")

plot(cvm-fourier, type="l")
mean(cvm-fourier)

cvm <- apply(x_mat, 1, cvm.res, "pnorm")
fourier <- apply(x_mat, 1, fourier_test, y="pnorm", lenj=200)

plot(cvm, type="l")
lines(fourier, col="red")

plot(cvm-fourier, type="l")
mean(cvm-fourier)

# x <- c(1,2,3,4,2,3,5,2)
# fx <- fft(x)
# spec.pgram(x)
# 
# abs(fx)

x <- 1:4
fx <- fft(x)
N <- length(x)

# Veryifying Parseval's Identity
sum_x2 <- sum(x^2)
sum_x2
sum_ft <- sum(abs(fx)^2)/N
sum_ft

# absolute value of a complex number is sqrt(Re^2 + Im^2)

# Amplitude of Frequency Components is:
amp <- abs(fx)/ N

# For a flat function, these had better be 0
x <- rep(5,20)
fx <- fft(x)
# And they are

x <- rnorm(30)
fx <- fft(x)
amp <- abs(fx)/ N

# Test Fan Power
set.seed(5)
x <- matrix(rnorm(100*50,0,1.3), ncol=100, nrow=50)
pow_fan <- power.res.onesamp(x,pnorm, f=Fan_Test, num.perm=1000) 
pow_fan_1 <- power.res.onesamp(x,pnorm, f=Fan_Test, num.perm=1000, fops=list(m=1))
pow_ks <- power.res.onesamp(x,pnorm, f=ks.res.simp, num.perm=1000)

pow_fan_50 <- sum(pow_fan < .05)/length(pow_fan)
pow_ks_50 <- sum(pow_ks < .05)/length(pow_ks)
pow_fan1_50 <- sum(pow_fan_1 < .05)/length(pow_fan_1)

pow_fan_50
pow_ks_50
pow_fan1_50

x <- matrix(rnorm(100*50,0.2,1), ncol=100, nrow=50)

pow_fan <- power.res.onesamp(x,pnorm, f=Fan_Test, num.perm=1000)
pow_fan_1 <- power.res.onesamp(x,pnorm, f=Fan_Test, num.perm=1000, fops=list(m=1))
pow_ks <- power.res.onesamp(x,pnorm, f=ks.res.simp, num.perm=1000)

pow_fan_50 <- sum(pow_fan < .05)/length(pow_fan)
pow_fan1_50 <- sum(pow_fan_1 < .05)/length(pow_fan_1)
pow_ks_50 <- sum(pow_ks < .05)/length(pow_ks)

pow_fan_50
pow_fan1_50
pow_ks_50


x <- matrix(rnorm(100*200,0,1.3), ncol=100, nrow=200)

pow_fan <- power.res.onesamp(x,pnorm, f=Fan_Test, num.perm=1000) 
pow_fan_1 <- power.res.onesamp(x,pnorm, f=Fan_Test, num.perm=1000, fops=list(m=1))
pow_ks <- power.res.onesamp(x,pnorm, f=ks.res.simp, num.perm=1000)

pow_fan_200_sd <- sum(pow_fan < .05)/length(pow_fan)
pow_fan1_200_sd <- sum(pow_fan_1 < .05)/length(pow_fan_1)
pow_ks_200_sd <- sum(pow_ks < .05)/length(pow_ks)

pow_fan_200_sd
pow_fan1_200_sd
pow_ks_200_sd


x <- matrix(rnorm(100*200,0.2,1), ncol=100, nrow=200)

pow_fan <- power.res.onesamp(x,pnorm, f=Fan_Test, num.perm=1000)
pow_fan_1 <- power.res.onesamp(x,pnorm, f=Fan_Test, num.perm=1000, fops=list(m=1))
pow_ks <- power.res.onesamp(x,pnorm, f=ks.res.simp, num.perm=1000)

pow_fan_200 <- sum(pow_fan < .05)/length(pow_fan)
pow_fan1_200 <- sum(pow_fan_1 < .05)/length(pow_fan_1)
pow_ks_200 <- sum(pow_ks < .05)/length(pow_ks)

pow_fan_200
pow_fan1_200
pow_ks_200
# Fourier Series
# Functions
plot.frequency.spectrum <- function(X.k, xlimits=c(0,length(X.k))) {
  plot.data  <- cbind(0:(length(X.k)-1), Mod(X.k))
  
  # TODO: why this scaling is necessary?
  plot.data[2:length(X.k),2] <- 2*plot.data[2:length(X.k),2] 
  
  plot(plot.data, t="h", lwd=2, main="", 
       xlab="Frequency (Hz)", ylab="Strength", 
       xlim=xlimits, ylim=c(0,max(Mod(plot.data[,2]))))
}

# Plot the i-th harmonic
# Xk: the frequencies computed by the FFt
#  i: which harmonic
# ts: the sampling time points
# acq.freq: the acquisition rate
plot.harmonic <- function(Xk, i, ts, acq.freq, color="red") {
  Xk.h <- rep(0,length(Xk))
  Xk.h[i+1] <- Xk[i+1] # i-th harmonic
  harmonic.trajectory <- get.trajectory(Xk.h, ts, acq.freq=acq.freq)
  points(ts, harmonic.trajectory, type="l", col=color)
}





acq.freq <- 100                    # data acquisition (sample) frequency (Hz)
time     <- 6                      # measuring time interval (seconds)
ts       <- seq(0,time-1/acq.freq,1/acq.freq) # vector of sampling time-points (s) 
f.0 <- 1/time

dc.component <- 1
component.freqs <- c(3,7,10)        # frequency of signal components (Hz)
component.delay <- c(0,0,0)         # delay of signal components (radians)
component.strength <- c(1.5,.5,.75) # strength of signal components

f   <- function(t,w) { 
  dc.component + 
    sum( component.strength * sin(component.freqs*w*t + component.delay)) 
}

plot.fourier(f,f.0,ts=ts)

w <- 2*pi*f.0
trajectory <- sapply(ts, function(t) f(t,w))
head(trajectory,n=30)

X.k <- fft(trajectory)                   # find all harmonics with fft()
plot.frequency.spectrum(X.k, xlimits=c(0,20))

### Other Stuff

x <- rnorm(20)
fx <- fft(x)

sum(abs(fx)^2)/length(fx)
sum(x^2)


sum_cn2 <- sum(abs(fx)^2)
sum_cn <- sum((abs(fx)/(pi)^2))
sum_cn
sum_cn2
sum_x2 <- sum(x^2)
sum_x2
# Computing real coefficients from complex

an <- Re(fx)
bn <- Im(fx)

abs_an <- abs(an)
abs_bn <- abs(bn)


# sum(an) ~~ sum(x^2)

x <- 1:4
fx <- fft(x)

an <- Re(fx)
an
bn <- Im(fx)
bn

abs_an <- abs(an)
abs_an
abs_bn <- abs(bn)
abs_bn

# IMPORTANT
sum_x2 <- sum(x^2)
sum_x2
# IMPORTANT
sum_ft <- sum(abs(fx)^2)/4

# DWT Stuff

x <- 1:16
library(waveslim)
dwt(x)
sum(unlist(dwt(x))^2)
sum(x^2)

x <- seq(0, 1 length.out=16)
fx <- fft(x)
sum(abs(fx)^2)/length(fx)
sum(x^2)

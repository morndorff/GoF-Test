# December 5th, 2014

# Finding/Fixing Problems with wave.den

# Suspicion: May be very bad when sample size is very large
Results <- rep(0,500)
for(i in 1:500){
  x <- rnorm(1000)
  y <- rnorm(1000)
  Results[i] <- wave.den(x,y)
}
mean(Results)
var(Results)

Results2 <- rep(0,500)
for(i in 1:500){
  x1 <- rnorm(1000)
  y1 <- rnorm(1000, 0, 1.1)
  Results2[i] <- wave.den(x1,y1)
}
mean(Results2)
var(Results2)

ResultsC <- rep(0,500)
for(i in 1:500){
  x <- rnorm(1000)
  y <- rnorm(1000)
  ResultsC[i] <- wave.energy(x,y)
}
mean(ResultsC)
var(ResultsC)

Results2C <- rep(0,500)
for(i in 1:500){
  x1 <- rnorm(1000)
  y1 <- rnorm(1000, 0, 1.1)
  Results2C[i] <- wave.energy(x1,y1)
}
mean(Results2C)
var(Results2C)

# TAKEAWAY: Much lower variance


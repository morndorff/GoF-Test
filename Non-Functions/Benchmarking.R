# Quick Benchmarking
x <- rnorm(50)
x1 <- rnorm(500)
x2 <- rnorm(10000)
x3 <- rnorm(30000)
library(microbenchmark)
library(waveslim)
microbenchmark(
  wave.bec(x,"qnorm"),
  wave.energy(x,"pnorm")
)

microbenchmark(
  wave.bec(x1,"qnorm"),
  wave.energy(x1,"pnorm")
)

microbenchmark(
  wave.bec(x2,"qnorm"),
  wave.energy(x2,"pnorm")
)
# September 22, 2014
# Updating the Area Approximation Test
# TODO: Change density estimate to CDF using cumsum
# TODO: Split into two variants: one based on order statistics, one not.
# TODO: Bug Test

# rm(list=ls())
source("functions.R")
set.seed(5)
x <- rnorm(50)
x_pdf <- density(x)

y <- cumsum(x_pdf$y)
x_cdf <- x_pdf
x_cdf$y <- y/max(y)
plot(x_cdf)
a_cs <- Kernel_CDF_Estimate(x, opt="cumsum")
a_int <- Kernel_CDF_Estimate(x, opt="integrate")
a_cs(2)
a_int(2)

# TODO: Split
# Split is achieved. New TODO: Make this function efficient. 

library(microbenchmark)
x <- rnorm(100)
y <- rnorm(100)
microbenchmark(Quad_OS_Area_TS(x,y), Quad_Quan_Area_TS(x,y), neval=2000)

# Profiling
Rprof("file.out")
for(i in 1:1000) Quad_Quan_Area_TS(x,y)
Rprof(NULL)
summaryRprof("file.out")

# Important do.call info
test <- list(maxval=TRUE)
print(do.call(Quad_Quan_Area_TS, c(list(x,y),test)))
print(do.call(Quad_Quan_Area_TS, list(x,y,maxval=TRUE)))
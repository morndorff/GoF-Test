source("functions.R")

set.seed(1234)
x <- rnorm(45)
y <- rnorm(145)
sort(x)
wave.bec.test(x,y)

y <- "qnorm"
wave.bec.test(x,"qnorm")

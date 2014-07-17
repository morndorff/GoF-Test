# June 3, 2014

# Developing suitable functions and exploring ways of computing the Cramer von Mises
# test for one and two samples

# Results: For one sample tests, use of: cvm.test(x, 'pnorm', mean='2', sd='1') is
# recommended library(goftest)

# For two sample tests, use of: CramerVonMisesTwoSamples(x,y) is recommended
# library(CDFt)

# Unsolved: Getting asymtotic CVM tests.






library(dgof)

x <- rpois(100, 1)


# for discrete distributions
cvm.test(x, ecdf(rpois(100, 1)))

library(sos)

library(goftest)

# for continuous distributions
a <- cvm.test(x, "pnorm", mean = 2, sd = 1)
# has class htest (good)

# Note: library(goftest) seems to have what we want in the one sample case

# The two sample version is out of date, but that shouldn't matter...

library(CDFt)
x <- rnorm(100)
y <- rnorm(100)
CramerVonMisesTwoSamples(x, y)

# Calculates just the statistics (works for me)!

# Below is some code investigation
sortx <- sort(x)
sorty <- sort(y)
M <- length(sortx)
N <- length(sorty)
a <- data.frame(val = sortx, rang = seq(M), ens = rep(1, M))
b <- data.frame(val = sorty, rang = seq(N), ens = rep(2, N))
head(a)
head(b)
d <- rbind(a, b)
head(d)
d <- d[order(d$val), ]
head(d)
d <- data.frame(d, rangTot = seq(M + N))
head(d)

dtfM = d[which(d$ens == 1), ]
head(dtfM)

dtfN = d[which(d$ens == 2), ]

somN = sum((dtfN$rang - dtfN$rangTot)^2)
somM = sum((dtfM$rang - dtfM$rangTot)^2)
U = N * somN + M * somM
CvM = ((U/(N * M))/(N + M)) - ((4 * M * N - 1)/(6 * (M + N))) 

# May 19 Forging ahead with some power studies This uses some code from May18.R, in
# particular the 'finished' Power Code Cleanup from Earlier
set.seed(5)
source("functions.R")
liResults <- list(NULL)
liData200 <- list(NULL)
##################################################### N(0,3) vs t(3) Notes that these distributions have different variances
z <- 250  #Number of Samples
lens <- c(10, 20, 50, 100, 200)  #Number of Draws in each sample
num.dif.draws <- length(lens)  #How many different draws?
num.test <- 4  #Number of Different tests
results <- matrix(, nrow = num.dif.draws, ncol = num.test)  #Will record rejection rate later
for (i in 1:num.dif.draws) {
    # Generating Data
    x <- matrix(rnorm(z * lens[i], 0, 3), nrow = z, ncol = lens[i])  #X ~ N(0,3)
    y <- matrix(rt(z * lens[i], 3), nrow = z, ncol = lens[i])  #Y ~ t(3)
    # OBC
    pvals <- power.res(x, y, myts)
    results[i, 1] <- sum(pvals < 0.05)/z  #rejection Rate for OBC 
    # OBC
    pvals <- power.res(x, y, ks.res.simp)
    results[i, 2] <- sum(pvals < 0.05)/z  #rejection Rate for KS
    # Max(KS,OBC)
    pvals <- power.res(x, y, myts.max.simp)
    results[i, 3] <- sum(pvals < 0.05)/z  #rejection Rate for max of OBC, KS
    # Max(KS,OBC) using range adjustment
    pvals <- power.res(x, y, myts.max.range.simp)
    results[i, 4] <- sum(pvals < 0.05)/z  #rejection Rate for max of OBC, KS, range adjusted
}

rownames(results) <- lens
colnames(results) <- c("max_q", "perm_KS", "max_q_ks", "max_q_ks_adj")
liResults[[1]] <- results
names(liResults)[1] <- c("N(0,3) vs t(3)")
liData200[1] <- list(list(x, y))
x <- liData200[[1]][[1]]
y <- liData200[[1]][[2]]
perm.test(x[2, ], y[2, ], myts)[1]
diag.plot(x[2, ], y[2, ])  #p-val ~.002
perm.test(x[3, ], y[3, ], myts)[1]
diag.plot(x[3, ], y[3, ])  #p-val ~0
perm.test(x[4, ], y[4, ], myts)[1]
diag.plot(x[4, ], y[4, ])  #p-val ~0

################################################ 

################################################# N(0,sqrt(3)) vs t(3) Notes that these distributions have identical variances
ptm <- proc.time()
z <- 250  #Number of Samples
lens <- c(10, 20, 50, 100, 200)  #Number of Draws in each sample
num.dif.draws <- length(lens)  #How many different draws?
num.test <- 4  #Number of Different tests
results <- matrix(, nrow = num.dif.draws, ncol = num.test)  #Will record rejection rate later
for (i in 1:num.dif.draws) {
    # Generating Data
    x <- matrix(rnorm(z * lens[i], 0, sqrt(3)), nrow = z, ncol = lens[i])  #X ~ N(0,sqrt(3))
    y <- matrix(rt(z * lens[i], 3), nrow = z, ncol = lens[i])  #Y ~ t(3)
    # OBC
    pvals <- power.res(x, y, myts)
    results[i, 1] <- sum(pvals < 0.05)/z  #rejection Rate for OBC 
    # OBC
    pvals <- power.res(x, y, ks.res.simp)
    results[i, 2] <- sum(pvals < 0.05)/z  #rejection Rate for KS
    # Max(KS,OBC)
    pvals <- power.res(x, y, myts.max.simp)
    results[i, 3] <- sum(pvals < 0.05)/z  #rejection Rate for max of OBC, KS
    # Max(KS,OBC) using range adjustment
    pvals <- power.res(x, y, myts.max.range.simp)
    results[i, 4] <- sum(pvals < 0.05)/z  #rejection Rate for max of OBC, KS, range adjusted
}
rownames(results) <- lens
colnames(results) <- c("max_q", "perm_KS", "max_q_ks", "max_q_ks_adj")
liResults[[2]] <- results
names(liResults)[2] <- c("N(0,sqrt(3)) vs t(3)")
liData200[2] <- list(list(x, y))
proc.time() - ptm
# Diagnostic plots
x <- liData200[[2]][[1]]
y <- liData200[[2]][[2]]
diag.plot(x[2, ], y[2, ])  #p-val .062
diag.plot(x[3, ], y[3, ])  #p-val 1
diag.plot(x[4, ], y[4, ])  #p-val .074

############################################################ N(0,2) vs exp(1/2) Notes that these distributions have identical variances
############################################################ (different means)
ptm <- proc.time()
z <- 250  #Number of Samples
lens <- c(10, 20, 50, 100, 200)  #Number of Draws in each sample
num.dif.draws <- length(lens)  #How many different draws?
num.test <- 4  #Number of Different tests
results <- matrix(, nrow = num.dif.draws, ncol = num.test)  #Will record rejection rate later
for (i in 1:num.dif.draws) {
    # Generating Data
    x <- matrix(rnorm(z * lens[i], 0, 2), nrow = z, ncol = lens[i])  #X ~ N(0,2)
    y <- matrix(rexp(z * lens[i], 1/2), nrow = z, ncol = lens[i])  #Y ~ exp(1/2)
    # OBC
    pvals <- power.res(x, y, myts)
    results[i, 1] <- sum(pvals < 0.05)/z  #rejection Rate for OBC 
    # OBC
    pvals <- power.res(x, y, ks.res.simp)
    results[i, 2] <- sum(pvals < 0.05)/z  #rejection Rate for KS
    # Max(KS,OBC)
    pvals <- power.res(x, y, myts.max.simp)
    results[i, 3] <- sum(pvals < 0.05)/z  #rejection Rate for max of OBC, KS
    # Max(KS,OBC) using range adjustment
    pvals <- power.res(x, y, myts.max.range.simp)
    results[i, 4] <- sum(pvals < 0.05)/z  #rejection Rate for max of OBC, KS, range adjusted
}
rownames(results) <- lens
colnames(results) <- c("max_q", "perm_KS", "max_q_ks", "max_q_ks_adj")
liResults[[3]] <- results
names(liResults)[3] <- c("N(0,2) vs exp(1/2)")
liData200[3] <- list(list(x, y))



proc.time() - ptm
# Diagnostic plots
x <- liData200[[3]][[1]]
y <- liData200[[3]][[2]]
# perm.test(x[2,],y[2,],myts)[1]
diag.plot(x[2, ], y[2, ])
# perm.test(x[3,],y[3,],myts)[1]
diag.plot(x[3, ], y[3, ])
# perm.test(x[4,],y[4,],myts)[1]
diag.plot(x[4, ], y[4, ])
######################################### N(0,1) vs N(0,2) Notes that these distributions have different variances
ptm <- proc.time()
z <- 250  #Number of Samples
lens <- c(10, 20, 50, 100, 200)  #Number of Draws in each sample
num.dif.draws <- length(lens)  #How many different draws?
num.test <- 4  #Number of Different tests
results <- matrix(, nrow = num.dif.draws, ncol = num.test)  #Will record rejection rate later
for (i in 1:num.dif.draws) {
    # Generating Data
    x <- matrix(rnorm(z * lens[i], 0, 1), nrow = z, ncol = lens[i])  #X ~ N(0,1)
    y <- matrix(rnorm(z * lens[i], 0, 2), nrow = z, ncol = lens[i])  #Y ~ N(0,3)
    # OBC
    pvals <- power.res(x, y, myts)
    results[i, 1] <- sum(pvals < 0.05)/z  #rejection Rate for OBC 
    # OBC
    pvals <- power.res(x, y, ks.res.simp)
    results[i, 2] <- sum(pvals < 0.05)/z  #rejection Rate for KS
    # Max(KS,OBC)
    pvals <- power.res(x, y, myts.max.simp)
    results[i, 3] <- sum(pvals < 0.05)/z  #rejection Rate for max of OBC, KS
    # Max(KS,OBC) using range adjustment
    pvals <- power.res(x, y, myts.max.range.simp)
    results[i, 4] <- sum(pvals < 0.05)/z  #rejection Rate for max of OBC, KS, range adjusted
}
rownames(results) <- lens
colnames(results) <- c("max_q", "perm_KS", "max_q_ks", "max_q_ks_adj")
liResults[[4]] <- results
names(liResults)[4] <- c("N(0,1) vs N(0,2)")
liData200[4] <- list(list(x, y))
proc.time() - ptm
# Diagnostic plots
x <- liData200[[4]][[1]]
y <- liData200[[4]][[2]]
perm.test(x[2, ], y[2, ], myts)[1]
diag.plot(x[2, ], y[2, ])  #p-val ~.002
perm.test(x[3, ], y[3, ], myts)[1]
diag.plot(x[3, ], y[3, ])  #p-val ~0
perm.test(x[4, ], y[4, ], myts)[1]
diag.plot(x[4, ], y[4, ])  #p-val ~0
############################################ Cleaning up some earlier Naming mistakes
names(liResults)[2:3] <- c("N(0,sqrt(3)) vs t(3)", "N(0,2) vs exp(1/2)")
########################################### N(3,3) vs Gamma(3,1) Notes that these distributions have different variances
ptm <- proc.time()
z <- 250  #Number of Samples
lens <- c(10, 20, 50, 100, 200)  #Number of Draws in each sample
num.dif.draws <- length(lens)  #How many different draws?
num.test <- 4  #Number of Different tests
results <- matrix(, nrow = num.dif.draws, ncol = num.test)  #Will record rejection rate later
for (i in 1:num.dif.draws) {
    # Generating Data
    x <- matrix(rnorm(z * lens[i], 0, sqrt(3)), nrow = z, ncol = lens[i])  #X ~ N(3,3)
    y <- matrix(rgamma(z * lens[i], 3, shape = 1), nrow = z, ncol = lens[i])  #Y ~ Gamma(alpha=3,beta=1) (shape=1,rate=1)
    # OBC
    pvals <- power.res(x, y, myts)
    results[i, 1] <- sum(pvals < 0.05)/z  #rejection Rate for OBC 
    # OBC
    pvals <- power.res(x, y, ks.res.simp)
    results[i, 2] <- sum(pvals < 0.05)/z  #rejection Rate for KS
    # Max(KS,OBC)
    pvals <- power.res(x, y, myts.max.simp)
    results[i, 3] <- sum(pvals < 0.05)/z  #rejection Rate for max of OBC, KS
    # Max(KS,OBC) using range adjustment
    pvals <- power.res(x, y, myts.max.range.simp)
    results[i, 4] <- sum(pvals < 0.05)/z  #rejection Rate for max of OBC, KS, range adjusted
}
rownames(results) <- lens
colnames(results) <- c("max_q", "perm_KS", "max_q_ks", "max_q_ks_adj")
liResults[[5]] <- results
names(liResults)[5] <- c("N(3,sqrt(3)) vs Gamma(3,1)")
liData200[5] <- list(list(x, y))
proc.time() - ptm
# Diagnostic plots
x <- liData200[[5]][[1]]
y <- liData200[[5]][[2]]
perm.test(x[2, ], y[2, ], myts)[1]
diag.plot(x[2, ], y[2, ])  #p-val ~.002
perm.test(x[3, ], y[3, ], myts)[1]
diag.plot(x[3, ], y[3, ])  #p-val ~0
perm.test(x[4, ], y[4, ], myts)[1]
diag.plot(x[4, ], y[4, ])  #p-val ~0
hist(rgamma(3000, 3, 1))
hist(rnorm(3000, 3, sqrt(3)))
################################ 

save.image("~/Dropbox/Research/GoF Test/May20.RData")
savehistory("~/Dropbox/Research/GoF Test/May20.Rhistory") 

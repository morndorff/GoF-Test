# A script illustrating problems with the one sample
# test statistic for myts

source("functions.R")
set.seed(5)

tx <- rt(100,3)
normx <- rnorm(100,0,sqrt(3))

new.perm.test(tx,qnorm,list(mean=0,sd=3) ,f=myts)

new.perm.test(normx,qnorm,list(mean=0,sd=3), f=myts)

new.perm.test(normx,qnorm,list(mean=0,sd=3800), f=myts)
# If you run the above, you'll notice that still occasionally,
# you'll get a non-zero p-value (kind of odd)
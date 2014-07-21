# A script illustrating problems with the one sample
# test statistic for myts

rm(list=ls())
library(ggplot2)
source("functions.R")
set.seed(5)

tx <- rt(100,3)
normx <- rnorm(100,0,sqrt(3))

new.perm.test(tx,qnorm,list(mean=0,sd=3) ,f=myts)

new.perm.test(normx,qnorm,list(mean=0,sd=3), f=myts)

new.perm.test(normx,qnorm,list(mean=0,sd=3800), f=myts)
# If you run the above, you'll notice that still occasionally,
# you'll get a non-zero p-value (kind of odd)

normx <- rnorm(50,0,3)

new.perm.test(normx, qt, list(df=3), f=myts)

normx1 <- rnorm(100,0,3)

new.perm.test(normx1, qt, list(df=3), f=myts)

normx1 <- rnorm(200,0,3)

new.perm.test(normx1, qt, list(df=3), f=myts)

normx1 <- rnorm(500,0,3)

new.perm.test(normx1, qt, list(df=3), f=myts)

# Note that p-value increases as number of samples increases
# major problem

# An illustratio nof why this happens
myts(rt(500,3), qt,3,do.plot=TRUE)
myts(rt(500,3), qt,3,do.plot=TRUE)
myts(rt(500,3), qt,3,do.plot=TRUE)
myts(rt(500,3), qt,3,do.plot=TRUE)

# Distribution of my ts under the null that F0=t(3)
matx_t_3 <- matrix(rt(1000000,3),nrow=2000,ncol=500)
dist_t_3 <- apply(matx_t_3,MARGIN=1, FUN=myts, y=qt,df=3)
plot_t_3 <- hist(dist_t_3,breaks=20, main="Distribution of myts under F0=t(3)")
plot(plot_t_3)

qplot(dist_t_3)

# distribution of myts when F=N(0,3) and the null is F0=t(3)
matx_n_03_t_3 <- matrix(rnorm(1000000,0,3),nrow=2000,ncol=500)
dist_n_03_t_3 <- apply(matx_n_03_t_3,MARGIN=1, FUN=myts, y=qt,df=3)
plot_n_03_t_3 <- hist(dist_n_03_t_3,breaks=20, main="Distribution of myts under F0=t(3)")
plot(plot_n_03_t_3)

qplot(dist_n_03_t_3)

x <- c(dist_t_3,dist_n_03_t_3)
y <- rep(c("a","b"),each=length(dist_t_3))
overlay_data <- data.frame(x,y)

quantile(dist_n_03_t_3,probs=.95)
quantile(dist_t_3,probs=.95)
mean(dist_t_3)
mean(dist_n_03_t_3)

ggplot(data=overlay_data,aes(x=x)) +
  geom_density(data=subset(overlay_data, y=='a'), fill="red", alpha=.2) + 
  geom_density(data=subset(overlay_data, y=='b'), fill="blue", alpha=.2) +
  coord_cartesian(xlim=c(0,5))

# ggplot(overlay_data, aes(x=dist_t_3, fill=yy)) + geom_histogram(alpha=0.2, position="identity")
# 
# + theme_bw()



overlay_plot
overlay_plot + geom_histogram(aes(x=dist_n_03_t_3))
# The distribution of myts (obviously) depends on the 
# underlying null distribution
# Many of the distributions are pretty interesting
matx <- matrix(runif(1000000,min=0,max=2),nrow=2000,ncol=500)
dist <- apply(matx,MARGIN=1, FUN=myts, y=qunif,min=0, max=2)
plot_u_0_2 <- hist(dist,breaks=20, main="Distribution of myts under F0=U(0,2)")
plot(plot_u_0_2)
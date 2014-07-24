# A script illustrating problems with the one sample
# test statistic for myts

rm(list=ls())
library(ggplot2)
source("functions.R")
set.seed(5)

tx <- rt(100,3)
normx <- rnorm(100,0,sqrt(3))

ptest_t3_n03 <- new.perm.test(tx,qnorm,list(mean=0,sd=3) ,f=myts)
ptest_t3_n03

ptest_n0rt3_n03 <- new.perm.test(normx,qnorm,list(mean=0,sd=3), f=myts)
ptest_n0rt3_n03

ptest_n0rt3_n02800 <- new.perm.test(normx,qnorm,list(mean=0,sd=3800), f=myts)
ptest_n0rt3_n02800
# If you run the above, you'll notice that still occasionally,
# you'll get a non-zero p-value (kind of odd)

normx <- rnorm(50,0,3)

pteset_n03_rt3_50 <- new.perm.test(normx, qt, list(df=3), f=myts)
pteset_n03_rt3_50

normx1 <- rnorm(100,0,3)

ptest_n03_t3 <- new.perm.test(normx1, qt, list(df=3), f=myts)
ptest_n03_t3_100

normx1 <- rnorm(200,0,3)

ptest_n03_t3 <- new.perm.test(normx1, qt, list(df=3), f=myts)
ptest_n03_t3_100
normx1 <- rnorm(500,0,3)

ptest_n03_t3_500 <- new.perm.test(normx1, qt, list(df=3), f=myts)
ptest_n03_t3_500

# Note that p-value increases as number of samples increases
# major problem

# An illustratio nof why this happens
myts(rt(500,3), qt,3,do.plot=TRUE)
myts(rt(500,3), qt,3,do.plot=TRUE)
myts(rt(500,3), qt,3,do.plot=TRUE)
myts(rt(500,3), qt,3,do.plot=TRUE)

myts(rt(20000,3), qt,3,do.plot=TRUE)
myts(rnorm(20000,0,3), qt,3,do.plot=TRUE)

myts(rt(50000,3), qt,3)
myts(rnorm(50000,0,3), qt,3)

myts(rt(100000,3), qt,3) # Finally starts to win?
myts(rnorm(100000,0,3), qt,3)


graphfunc <- function(sampsize){
  
  # Distribution of my ts under the null that F0=t(3)
  matx_t_3 <- matrix(rt(sampsize*2000,3),nrow=2000,ncol=sampsize)
  dist_t_3 <- apply(matx_t_3,MARGIN=1, FUN=myts, y=qt,df=3)
  #plot_t_3 <- hist(dist_t_3,breaks=20, main="Distribution of myts under F0=t(3)")
  #plot(plot_t_3)
  
  #qplot(dist_t_3, title="Distribution of myts under F0=t(3)")
  
  # distribution of myts when F=N(0,3) and the null is F0=t(3)
  matx_n_03_t_3 <- matrix(rnorm(sampsize*2000,0,3),nrow=2000,ncol=sampsize)
  dist_n_03_t_3 <- apply(matx_n_03_t_3,MARGIN=1, FUN=myts, y=qt,df=3)
  # plot data
  x <- c(dist_t_3,dist_n_03_t_3)
  y <- rep(c("a","b"),each=length(dist_t_3))
  overlay_data <- data.frame(x,y)
  
  alt_lim <- as.numeric(quantile(dist_n_03_t_3,probs=.95))
  null_lim <- as.numeric(quantile(dist_t_3,probs=.95))
  plotxlim <- max(alt_lim,null_lim)
  
  limdata <- data.frame(c(alt_lim,null_lim),c("alt_lim","null_lim"))
  
  mean(dist_t_3)
  mean(dist_n_03_t_3)
  
  twoplottitle <- paste("Distribution of max quantile distance under Null (red) and Normal (Blue). Sample Size = ",sampsize)
  
  twoplot <- ggplot(data=overlay_data,aes(x=x)) +
    geom_density(data=subset(overlay_data, y=='a'), fill="red", alpha=.2) + 
    geom_density(data=subset(overlay_data, y=='b'), fill="blue", alpha=.2) +
    coord_cartesian(xlim=c(0,plotxlim)) +
    labs(title=twoplottitle) + 
    xlab("Max Quantile Distance") +
    geom_segment(aes(x=null_lim, y=0, xend=null_lim, yend =.2), color="red")
  return(twoplot)
}

plot200 <- graphfunc(200)
plot200
plot500 <- graphfunc(500)
plot500
plot2000 <- graphfunc(2000)
plot2000
plot10000 <- graphfunc(10000)
plot10000













overlay_plot
overlay_plot + geom_histogram(aes(x=dist_n_03_t_3))
# The distribution of myts (obviously) depends on the 
# underlying null distribution
# Many of the distributions are pretty interesting
matx <- matrix(runif(1000000,min=0,max=2),nrow=2000,ncol=500)
dist <- apply(matx,MARGIN=1, FUN=myts, y=qunif,min=0, max=2)
plot_u_0_2 <- hist(dist,breaks=20, main="Distribution of myts under F0=U(0,2)")
plot(plot_u_0_2)
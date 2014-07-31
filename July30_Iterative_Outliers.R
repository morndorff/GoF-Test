# Iterative Outliers
# Take 2
# July 30, 2014
# Goals: Add counter to outlier code
source("functions.R")
set.seed(5)
# Did small vector for exact --------
# x <- c(1,2,3,4,5,6,7,8,9,10)
# y <- c(1,2,3,4,5,6,7,8,9,100)
# 
# new.perm.test.out(x,y,f=myts.out) #(1)
# # expected output is: one outlier removed, 
# # p-value high
# new.perm.test(x,y,f=myts) #(2)
# 
# new.perm.test(x,y[1:9],f=myts) #(3)
# # as of right now, there is a difference b/w (1) and (3)
# new.perm.test.out(x,y[1:9],f=myts.out) #(4)
# 
# new.perm.test(x[1:9],y[1:9],f=myts) #(5)
# 
# # (1) and (5) are the same
# # (3) and (4) are also the same
# 
# # current outlier test deletes both entries of the offending
# # observation
# Larger vector
# Larger vector --------
x <- c(1:15)
y <- c(1:14,100)
new.perm.test(x,y,f=myts)
new.perm.test.out(x,y,f=myts)

# Current task: idientify x or y in a non-dumb way ----------
# temporarily set new.perm.test.out to return quan.mat
quanmat <- new.perm.test.out(x,y,f=myts.out)
quantable <- table(quanmat)
str(quantable)
attributes(quantable)
typeof(quantable)
sample(quantable, size=4, replace=TRUE)
# testing out schemes
x <- 1:10
y <- 11:20
names(x) <- rep("x",length(x))
names(y) <- rep("y",length(y))
newvec <- c(x,y)
newsam <- sample(newvec,size=20,replace=TRUE)
# above scheme seems promising
tabsam <- table(newsam)
which.max(tabsam)
which(x == max(tabsam)) # also includes the fact that it is from x

# slightly different scheme:
x <- 1:10
y <- 1:10
names(x) <- rep("x",length(x))
names(y) <- rep("y",length(y))
newvec <- c(x,y)
newsam <- sample(newvec,size=20,replace=TRUE)
lisam <- as.list(newsam)
# above scheme seems promising
tabsam <- table(newsam)
which.max(tabsam)
which(x == max(tabsam))
tablisam <- table(lisam)

newmat <- matrix(c(names(x),x),nrow=2,ncol=length(x), byrow=TRUE)
as.list(newmat)

test <- lapply(x,function(x) list(x,"x"))
test
str(test)
samtest <- sample(test, size=15, replace=TRUE)
table(samtest[][[]][1])
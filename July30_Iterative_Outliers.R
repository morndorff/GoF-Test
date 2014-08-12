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


######### Implementing new scheme based on matching:
x <- 1:10
y <- 1:10
z <- c(x,y)
z_seq <- seq_along(z)
seq_sam <- sample(z_seq, size=20, replace=TRUE)
tabsam <- table(seq_sam)
maxtab <- which.max(seq_sam)
pos <- as.numeric(names(seq_sam))
z_sam <- z[seq_sam]


from <- 


str(maxtab)
which(x == max(tabsam))


# Figuring out why perm.test.match.out isn't working

# x <- runif(100,0,1)
# y <- runif(100,0,1)
# lenx <- length(x)
# leny <- length(y)
# z <- c(x,y)
# z_seq <- seq_along(z)
# lenz <- length(z)
# 
# i <- 1
# seq_sam <- sample(z_seq, size = lenz, replace= FALSE)
# sam_x <- z[seq_sam[1:lenx]]
# sam_y <- z[seq_sam[(lenx+1):lenz]]
# sam_x <- sort(sam_x)
# sam_y <- sort(sam_y)
# res <- f(sam_x, sam_y)
# ts.random[i] <- res[[1]]
# ord <- res[[2]]
# x_out <- sam_x[ord]
# y_out <- sam_y[ord]
# x_out_ind <- match(x_out,z)
# y_out_ind <- match(y_out,z)
# quan.mat[i, ] <- c(x_out_ind,y_out_ind)



# New Testing
set.seed(5)
x <- runif(100,0,1)
y <- runif(100,0,1)
perm.test.match.out(x,y,f=myts.out)
perm.test.out(x,y,f=myts.out)

x <- rt(100,3)
y <- rnorm(100,0,sqrt(3))
perm.test.match.out(x,y,f=myts.out)
perm.test.out(x,y,f=myts.out)

x <- c(1:20)
y <- c(1:19,200)
perm.test.match.out(x,y,f=myts.out)
perm.test.out(x,y,f=myts.out)

# Creating a toy function to do what we want

iter.test <- function(x){
  if(exists("i", where = -1, inherits=FALSE) == FALSE){
    i <- 0
  }
    x <- x/2
  while(x > 1 ){
    x <- iter.test(x)
  }
  i <<- i+1
  x
  return(x)
}


# This works!
iter.test2 <- function(x, count=0){
  x <- x/2
  count <- count+1
  while(x > 1 ){
    stor <- do.call(iter.test2,list(x=x, count=count))
    x <- stor[[1]]
    count <- stor[[2]]
  }
  return(list(x,count))
}
  
    

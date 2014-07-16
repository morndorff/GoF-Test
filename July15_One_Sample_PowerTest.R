# One Sample Power Code

# Also taking the opportunity to generate data better than before
# Generate all the data first, store it in a list, then apply the methods
# to the stored data
# Avoids messes when I want to add a method later and use the same data

# rm(list=ls())

set.seed(5)
source("functions.R")

# Initializing Power Simulation Parameters
z <- 10  #Number of samples
lens <- c(4, 5, 6, 7, 8)  #Number of draws in each sample
num.samps <- length(lens)  #How many different draws?


# Distributions to test
dist <- list("norm",
             "t",
             "unif")
param <- list(list(mean=0,sd=1),
              list(df=3),
              list(min=0,max=1))
num.dist.test <- length(dist) #Number of different distributions to test

#Methods to use
liMethods <- list(myts, 
                  ks.res.simp)#,myts.max.simp,myts.max.range.simp,myts.out,myts.par)
num.methods <- length(liMethods)


# 
# Generating Data
# 

liData <- vector(mode="list", length=num.dist.test)
# results <- matrix(, nrow = num.samps ncol = num.test)  #Will record rejection rate later

for (i in 1:num.dist.test){
  samps <- vector(mode="list", length=num.samps)
  for (j in 1:num.samps) {
    # Generating Data
    #x<- matrix(distf[[i]](z * lens[j], param[[i]]), nrow = z, ncol = lens[j])
    x<- matrix(make_sample(n = z * lens[j], dist=dist[[i]], param=param[[i]]), nrow = z, ncol = lens[j])
    
    samps[[j]] <- x
  }
  liData[[i]] <- samps
}

# Naming
for(i in 1:num.dist.test){
  names(liData)[i] <-paste(dist[[i]], "(", paste(param[[i]], collapse=","), ")", sep="")
}
# Doing the tests


liPVals <- vector(mode="list", length=num.dist.test)
# Number of distributions to test

for (i in 1:num.dist.test) {
  pval_samp <- vector(mode="list",length=num.samps)
    # Number of sample sizes (10, 20, 50, 100)
  for (j in 1:num.samps) {
        # Looking up Data to Use
    x <- liData[[i]][[j]]   
    pval_methods <- vector(mode="list",length=num.methods)
    
    # Number of methods to test with
    for (k in 1:num.methods){
      # Evaluating using the correct method
      pvals <- power.res.onesamp(x, y=paste("q",dist[[i]],sep=""), param[[i]], f=liMethods[[k]], g=new.perm.test)
      pval_methods[[k]] <- pvals
      }
    pval_samp[[j]] <- pval_methods
  }
  liPVals[[i]] <- pval_samp
}
# Making table of results from the p-value Data




# results <- matrix(, nrow = num.dif.draws, ncol = num.test)  #Will record rejection rate later



rownames(results) <- lens
colnames(results) <- c("max_q", "perm_KS", "max_q_ks", "max_q_ks_adj", "q_out")
liResults[[1]] <- results
names(liResults)[1] <- c("N(0,3) vs t(3)")
liData200[1] <- list(list(x, y))

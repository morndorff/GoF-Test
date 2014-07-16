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
liMethods <- list(myts=myts, 
                  ks.res.simp=ks.res.simp)#,myts.max.simp,myts.max.range.simp,myts.out,myts.par)
num.methods <- length(liMethods)

# P-Value Cutoff to assess power 
cutoff <- .05


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

# Naming the PVals function
for(i in 1:num.dist.test){
  names(liPVals)[i] <-paste(dist[[i]], "(", paste(param[[i]], collapse=","), ")", sep="")
    for (j in 1:num.samps) {
      names(liPVals[[i]])[j] <- paste("Samp Size =", lens[j] ,  sep=" ")
      for (k in 1:num.methods){
        names(liPVals[[i]][[j]])[k] <- names(liMethods[k])
      }
    }
}

# Making table of results from the p-value Data

liResults <- vector(mode="list", length=num.dist.test)

# Number of distributions to test

for (i in 1:num.dist.test) {
  res_samp <- vector(mode="list",length=num.samps)
  # Number of sample sizes (10, 20, 50, 100)
  for (j in 1:num.samps) {
    # Looking up Pvals to evaluate
    x <- liPVals[[i]][[j]]   
    res_methods <- vector(mode="list",length=num.methods)
    
    # Number of methods to test with
    for (k in 1:num.methods){
      # Evaluating using the correct method
      pvals <- liPVals[[i]][[j]][[k]]
      res_methods[[k]] <- sum(pvals>.05)/length(pvals)
    }
    res_samp[[j]] <- res_methods
  }
  liResults[[i]] <- res_samp
}

# Naming the results list
for(i in 1:num.dist.test){
  names(liResults)[i] <-paste(dist[[i]], "(", paste(param[[i]], collapse=","), ")", sep="")
  for (j in 1:num.samps) {
    names(liResults[[i]])[j] <- paste("Samp Size =", lens[j] ,  sep=" ")
    for (k in 1:num.methods){
      names(liResults[[i]][[j]])[k] <- names(liMethods[k])
    }
  }
}

# Making into a more human readable set of tables
liTables <- vector(mode="list", length=num.dist.test)

for (i in 1:num.dist.test) {
  liTables[[i]] <- matrix(,ncol=num.methods,nrow=num.samps)
  for (j in 1:num.samps){
    liTables[[i]][j,] <- as.numeric(liResults[[i]][[j]])
    
  }
  names(liTables)[i] <-paste(dist[[i]], "(", paste(param[[i]], collapse=","), ")", sep="")
  rownames(liTables[[i]]) <- paste("Samp Size =", lens ,  sep=" ")
  colnames(liTables[[i]]) <- paste(names(liMethods))
}

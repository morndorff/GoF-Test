# August 13, 2014

# Creating two sample power code.


# Also taking the opportunity to generate data better than before 
# Generate all the  data first, store it in a list, 
# then apply the methods to the stored data 
#Avoids messes when I want to add a method later and use the same data

# rm(list=ls())

set.seed(5)
source("functions.R")

# Initializing Power Simulation Parameters
z <- 14  #Number of samples
lens <- c(4, 5, 6, 7, 8)  #Number of draws in each sample
num.samps <- length(lens)  #How many different draws?

# Distribution 1
dist_one <- list("norm", "norm", "norm", "norm", "norm", "unif", "unif", "unif")
param_one <- list(list(mean = 0, sd = 3), 
              list(mean = 0, sd = sqrt(3)), 
              list(mean = 0, sd = 2), 
              list(mean = 0, sd = 1), 
              list(mean = 0, sd = sqrt(3)), 
              list(min = 0.5, max = 1.5), 
              list(min = 0, max = 1), 
              list(min = 0, max = 1))

num.dist.test <- length(dist_one)  #Number of different distributions to test

# Distribution 2
dist_two <- list("t", "t", "exp", "norm", "gamma", "unif", "beta", "beta")
param_two <- list(list(df = 3), 
                   list(df = 3), 
                   list(rate = 0.5), 
                   list(mean = 0, sd = 2), 
                   list(shape = 1, rate = 3), 
                   list(min = 0, max = 2), 
                   list(shape1 = 0.5, shape2 = 0.5), 
                   list(shape1 = 0.2, shape2 = 2))

dist <-list(dist_one, dist_two)
param <- list(param_one, param_two)
# Methods to use
liMethods <- list(myts = myts, 
                  ks.res.simp = ks.res.simp,
                  myts.max.simp=myts.max.simp, 
                  myts.par=myts.par)

num.methods <- length(liMethods)

method_param <- list(NULL, 
                     NULL, 
                     NULL, 
                     NULL)

# P-Value Cutoff to assess power
cutoff <- 0.05

# Generating Data

liData <- vector(mode = "list", length = num.dist.test)
# results <- matrix(, nrow = num.samps ncol = num.test) #Will record rejection rate
# later
for (i in 1:num.dist.test) {
  liDist <- vector(mode="list", length=2)
  for (j in 1:2){
    samps <- vector(mode = "list", length = num.samps)
    for (k in 1:num.samps) {
      x <- matrix(make_sample(n = z * lens[k], dist = dist[[j]][[i]], param = param[[j]][[i]]), 
                  nrow = z, ncol = lens[k])
      
      samps[[k]] <- x
    }
    liDist[[j]] <- samps
  }
  liData[[i]] <- liDist
}
# Naming
for (i in 1:num.dist.test) {
  for (j in 1:2){
    names(liData[[i]])[j] <- paste(dist[[j]][[i]], "(", paste(param[[j]][[i]], collapse = ","), ")", 
                                 sep = "")
  }
}


# Doing tests and generating P-Values
liPVals <- vector(mode = "list", length = num.dist.test)
for (i in 1:num.dist.test) {
  pval_samp <- vector(mode = "list", length = num.samps)
  # Number of sample sizes (10, 20, 50, 100)
  for (j in 1:num.samps) {
    # Looking up Data to Use
    x <- liData[[i]][[1]][[j]]
    y <- liData[[i]][[2]][[j]]
    pval_methods <- vector(mode = "list", length = num.methods)
    
    # Number of methods to test with
    for (k in 1:num.methods) {
      # Evaluating using the correct method if(liMethods[[k]]=='myts.out'){ pvals <-
      # power.res.onesamp(x, y=paste('q',null_dist[[i]],sep=''), null_param[[i]],
      # f=liMethods[[k]], g=new.perm.test.out) }
      pvals <- power.res.twosamp(x, y, 
                                 f = names(liMethods[k]), 
                                 fops = method_param[[k]], 
                                 g = perm.test)
      pval_methods[[k]] <- pvals
    }
    pval_samp[[j]] <- pval_methods
  }
  liPVals[[i]] <- pval_samp
}
# Naming the PVals function
for (i in 1:num.dist.test) {
  names(liPVals)[i] <- paste(dist[[i]], "(", paste(param[[i]], collapse = ","), ")", " vs ",
                             null_dist[[i]], " (", paste(null_param[[i]], collapse = ","), ")", 
                             sep = "")
  for (j in 1:num.samps) {
    names(liPVals[[i]])[j] <- paste("Samp Size =", lens[j], sep = " ")
    for (k in 1:num.methods) {
      names(liPVals[[i]][[j]])[k] <- names(liMethods[k])
    }
  }
}


# Making table of results from the p-value Data
liResults <- vector(mode = "list", length = num.dist.test)
for (i in 1:num.dist.test) {
  res_samp <- vector(mode = "list", length = num.samps)
  # Number of sample sizes (10, 20, 50, 100)
  for (j in 1:num.samps) {
    # Looking up Pvals to evaluate
    x <- liPVals[[i]][[j]]
    res_methods <- vector(mode = "list", length = num.methods)
    
    # Number of methods to test with
    for (k in 1:num.methods) {
      # Evaluating using the correct method
      pvals <- liPVals[[i]][[j]][[k]]
      res_methods[[k]] <- sum(pvals < cutoff)/length(pvals)
    }
    res_samp[[j]] <- res_methods
  }
  liResults[[i]] <- res_samp
}
# Naming the results list
for (i in 1:num.dist.test) {
  names(liResults)[i] <- paste(dist[[i]], "(", paste(param[[i]], collapse = ","), ")", " vs ",
                               null_dist[[i]], " (", paste(null_param[[i]], collapse = ","), ")", 
                               sep = "")
  for (j in 1:num.samps) {
    names(liResults[[i]])[j] <- paste("Samp Size =", lens[j], sep = " ")
    for (k in 1:num.methods) {
      names(liResults[[i]][[j]])[k] <- names(liMethods[k])
    }
  }
}




# Making into a more human readable set of tables
liTables <- vector(mode = "list", length = num.dist.test)
for (i in 1:num.dist.test) {
  liTables[[i]] <- matrix(, ncol = num.methods, nrow = num.samps)
  for (j in 1:num.samps) {
    liTables[[i]][j, ] <- as.numeric(liResults[[i]][[j]])
    
  }
  names(liTables)[i] <- paste(dist[[i]], "(", paste(param[[i]], collapse = ","), ")", " vs ",
                              null_dist[[i]], " (", paste(null_param[[i]], collapse = ","), ")", 
                              sep = "")
  rownames(liTables[[i]]) <- paste("Samp Size =", lens, sep = " ")
  colnames(liTables[[i]]) <- paste(names(liMethods))
} 
# save.image("~/Dropbox/Research/GoF Test/Power Simulations/One Sample 4 Tests - July 17/July17OneSamplePower.RData")
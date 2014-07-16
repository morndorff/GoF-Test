# One Sample Power Code

# Also taking the opportunity to generate data better than before
# Generate all the data first, store it in a list, then apply the methods
# to the stored data
# Avoids messes when I want to add a method later and use the same data

set.seed(5)
source("functions.R")

# Initializing Power Simulation Parameters
z <- 1000  #Number of samples
lens <- c(10, 20, 50, 100, 200)  #Number of draws in each sample
num.samps <- length(lens)  #How many different draws?
num.dist.test <- 3 #Number of different distributions to test


# Distributions to test
dist <- list(rnorm,rt,runif)
param <- list(c(3,1),3,c(0,1))

#Methods to use
methods <- list(myts,ks.res.simp,myts.max.simp,myts.max.range.simp,myts.out,myts.par)
num.methods <- length(methods)

liData <- vector(mode="list", length=num.dist.test)
# results <- matrix(, nrow = num.samps ncol = num.test)  #Will record rejection rate later

for (i in 1:num.dist.test){
  samps <- vector(mode="list", length=num.samps)
  for (j in 1:num.samps) {
    # Generating Data
    x<- matrix(dist[[1]](z * lens[j], param[[i]]), nrow = z, ncol = lens[j])
    samps[[j]] <- x
  }
  liData[[i]] <- samps
}

names(liData)[1] <- c("N(3,1)")
names(liData)[2] <- c("t(3)")
names(liData)[3] <- c("U(0,1)")

# Doing the tests

# results <- matrix(, nrow = num.dif.draws, ncol = num.test)  #Will record rejection rate later

liPVals <- vector(mode="list", length=num.dist.test)
# Number of distributions to test

for (i in 1:num.dist.test) {
  
  # Number of sample sizes (10, 20, 50, 100)
  for (j in 1:num.samps) {
    
    # Looking up Data to Use
    x <- liData[[i]][[j]]   
    pval_samp <- vector(mode="list",length=num.samps)
    
    # Number of methods to test with
    for (k in 1:num.methods){
      pval_methods <- vector(mode="list",length=num.methods)
      # Evaluating using the correct method
      
      pvals <- power.res.onesamp(x, dist[[k]],param[[k]], methods[[k]], g=new.perm.test)
      
      pval_methods[[k]] <- pvals
      }
    pval_samp[[j]] <- pval_methods
  }
  liPvals[[i]] <- pval_samp
}
# Making table of results from the p-value Data







rownames(results) <- lens
colnames(results) <- c("max_q", "perm_KS", "max_q_ks", "max_q_ks_adj", "q_out")
liResults[[1]] <- results
names(liResults)[1] <- c("N(0,3) vs t(3)")
liData200[1] <- list(list(x, y))

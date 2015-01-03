# Finding ARL's for various out of control distributions
# For details, see Becvarik, page 36
# In Control distribution is N(0,1)

# rm(list=ls())
load("ARL_Baseline.RData")
source("functions.R")
set.seed(50)

Find_ARL_OC

##### Out of Control Distribution ####

# N(0,2)

Find_ARL_OOC <- function(num.samp=32, dist_one="norm", param_one=list(mean=0, sd=1),
                         dist_two="norm", param_two=list(mean=0, sd=2), cp=1,
                         tstat=wave.energy, UCL, time = 30, method=Find_CP_RL_Fast, ...){
  ptm <- proc.time()
  RLs <- vector(mode="list", length=0)
  e_time <- 0
  len_UCL <- length(UCL)
  # Calculating ARLs
  while(e_time < time){
    RLs_det <- method(num.samp=num.samp, 
                      dist_one=dist_one, param_one=param_one,
                      dist_two=dist_two, param_two=param_two,
                      tstat=tstat, UCL=UCL, cp=cp, ...)
    RLs[[length(RLs)+1]] <- RLs_det[[3]] # Append list of RL's
    howlong <- proc.time()-ptm
    e_time <- howlong["elapsed"]
  }
  
  matRL <- matrix(unlist(RLs), ncol=len_UCL, byrow=TRUE)# Making ARL Matrix
  ARL <- matrix(,nrow=2, ncol=len_UCL) 
  ARL[1,] <- colMeans(matRL)
  ARL[2,] <- apply(matRL, 2, sd)
  colnames(ARL) <- as.character(round(UCL,3))
  rownames(ARL) <- c("mean", "sd")
  print(ARL)
  return(list(ARL, RLs, e_time))
}





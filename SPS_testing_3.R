# SPS Testing 3

# Making Function to determine ARL for a set of given UCL's

Find_RL <- function(num.samp, dist, params, tstat, UCL){
  # num.samp - Number of Samples
  # dist - Incontrol distribution
  # params - Incontrol distribution parameters
  # tstat - Test Statistic for the process
  # UCL - Vector containing Upper Control Limits
  Proc <- NULL # Initializing
  count <- 1
  track_stat <- NULL
  tstat_proc <- 0
  while(tstat_proc < max(UCL)){
  Proc <- Sim_IC_Process_Iter(num.samp=num.samp, dist="norm", param=list(mean=0, sd=1), proc=Proc)
  tstat_proc <- Process_Stat(proc=Proc, tstat=wave.den, dist_ic="pnorm", mean=0, sd=1)
  track_stat <- append(track_stat, tstat_proc)
  count <- count +1
  }
  RL <- sapply(UCL, function(x) min(which(track_stat >= x)))
  res <- list("h(t)"=track_stat, "Length of Process"=count, "RL for Corresponding UCL"=RL, "UCLS"=UCL)
  return(res)
}

test <- Find_RL(num.samp=30, dist="norm", params=list(mean=0, sd=1), tstat=wave.den, UCL=c(.1,.2,.3))
  
ARL_Proc <- function(UCL){
  ptm <- proc.time()
  RLs <- vector(mode="list", length=0)
  e_time <- 0
  while(e_time < 2*60*60){
    RLs_det <- Find_RL(num.samp=30, dist="norm", params=list(mean=0, sd=1), tstat=wave.den, UCL=UCL)
    RLs[[length(RLs)+1]] <- RLs_det[[3]] # Append list of RL's
    howlong <- proc.time()-ptm
    e_time <- howlong["elapsed"]
  }
  return(list(RLs, e_time))
}

test2 <- ARL_Proc(UCL=seq(.1,.4, length.out=7))

numRuns <- length(test2[[1]])
matRL <- matrix(unlist(test2[[1]]),nrow=numRuns,ncol=7, byrow=TRUE)

ARL <- colMeans(matRL)
RL_sd <- apply(matRL,2,sd)
ARL <- rbind(ARL,RL_sd)
colnames(ARL) <- c(seq(.1,.4,length.out=7))

ARL
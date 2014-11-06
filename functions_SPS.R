# SPS Functions
Sim_IC_Process <- function(num.samp, run.length, dist, params){
  # Simulates a process of run length N with samples of size k
  # Inputs:
  # num.samp: The size of the samples. 1 corresponds to one data point per run unit
  # run.length: Length of the run
  # dist: what distribution do the data come from
  # params: what parameters to call from the dist function
  total_samps <- num.samp*run.length
  Proc <- matrix(make_sample(total_samps,dist,params),nrow=num.samp,ncol=run.length)
}

Sim_IC_Process_Iter <- function(proc=NULL, num.samp, dist, params){
  if(is.null(proc)){
    proc <- matrix(make_sample(num.samp, dist, params), nrow=num.samp, ncol=1)
    return(proc)
  }
  new_samp <- make_sample(num.samp, dist, params)
  proc <- cbind(proc, new_samp)
}



Sim_CP_Process <- function(num.samp,run.length,dist_one, param_one, 
                           dist_two, param_two, bpoint){
 samp_one <- (num.samp*bpoint)
 samp_two <- (num.samp*(run.length-bpoint))
 Proc_One <- matrix(make_sample(samp_one, dist_one, param_one),
                    nrow=num.samp, ncol=bpoint)
 Proc_Two <- matrix(make_sample(samp_two, dist_two, param_two),
                    nrow=num.samp, ncol=run.length-bpoint)
 CP_Proc <- cbind(Proc_One, Proc_Two)
 attr(CP_Proc, "bp") <- bpoint
 return(CP_Proc)
}

Process_Stat <- function(proc, tstat, dist_ic, ..., doplot=FALSE, detail=FALSE){
  # Tracks the value of a statistic for a process
  # Inputs:
  # proc: A matrix containing the process
  # stat: A two-sample test statistic
  lenproc <- dim(proc)[2] # run length of the process
  if(is.null(lenproc)) {
    ts <- NULL
    ts[1] <- 0
    ts[2] <- tstat(proc, dist_ic, ...)
    ts[3] <- ts[2]
    STAT <- ts[3]
    return(STAT)
  }# if process length is 1, this is needed
  
  ts <- matrix(nrow=3,ncol=lenproc)
  tau <- 0
  ts[1, (tau + 1)] <- 0
  ts[2, (tau + 1)] <- tstat(proc[, (tau+1):lenproc], dist_ic, ...)
  ts[3, (tau + 1)] <- ts[2, (tau+1)] - ts[1, (tau+1)]
  if(lenproc!=1) {
  for(tau in 1:(lenproc-1)){
    ts[1, (tau+1)] <- tstat(proc[, 1:tau], dist_ic, ...) #g(H0, tau-)
    
    ts[2, (tau+1)] <- tstat(proc[, (tau+1):lenproc], dist_ic, ...) #g(H0, tau+)
    
    ts[3, (tau+1)] <- ts[2, (tau+1)] - ts[1, (tau+1)] #g(H0, tau+) - g(H0, tau-)
  }
  }
  STAT <- max(ts[3,])
  
  if(detail){
    deta <- list("Statistic"=STAT,"All Stats"=ts)
    return(deta)
  }
  STAT
}

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

JKnife_Est <- function(data, tstat, ...){
  jack.est <- vector(mode="numeric", length=length(data))
  for(i in seq_along(data)){
    jack.est[i] <- tstat(x[-i], ...)
  }
  jack.mean <- mean(jack.est)
}

Old_TSO <- function(proc, tstat, dist_ic, ..., doplot=FALSE, detail=FALSE){
  # Tracks the value of a statistic for a process
  # Inputs:
  # proc: A matrix containing the process
  # stat: A two-sample test statistic
  lenproc <- dim(proc)[2] # run length of the process
  if(is.null(lenproc)) lenproc <- 1
  
  ts <- matrix(nrow=3,ncol=lenproc)
  for(i in 1:(lenproc-1)){
    ts[1, i] <- tstat(proc[, 1:i], dist_ic, ...)
    ts[2, i] <- tstat(proc[, (i+1):(lenproc-1)], dist_ic, ...)
    ts[3, i] <- ts[2,i] - ts[1,i]
  }
  STAT <- max(ts[3, i])
  if(detail){
    namets <- vector(mode="character", length=(lenproc-1))
    for(i in 1:(lenproc-1)){
      namets[i] <- paste(1,":",i, " vs " , i+1, ":", lenproc, sep="")
    }
    names(ts) <- namets
  }
  STAT
}
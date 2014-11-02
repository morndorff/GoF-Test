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

Track_Stat <- function(proc,tstat,doplot=FALSE){
  # Tracks the value of a statistic for a process
  # Inputs:
  # proc: A matrix containing the process
  # stat: the statistic tracked

  x <- apply(proc,2,tstat)
  if(doplot){
    plot(seq(1:dim(proc)[2]),x)
    title("Statistic on Each Sample")
  }
  x
}

Track_Stat_Over <- function(proc, tstat, dist_ic, ..., doplot=FALSE, detail=FALSE){
  # Tracks the value of a statistic for a process
  # Inputs:
  # proc: A matrix containing the process
  # stat: A two-sample test statistic
  
  
  
  
  
  lenproc <- dim(proc)[2] # run length of the process
  
  ts <- matrix(nrow=2,ncol=lenproc)
  for(i in 1:(lenproc-1)){
    ts[1, i] <- tstat(proc[, 1:i], dist_ic, ...)
    ts[2, i] <- tstat(proc[, 
  }
  
  
  
  large <- which.max(ts)
  if(detail){
    namets <- vector(mode="character", length=(lenproc-1))
    for(i in 1:(lenproc-1)){
      namets[i] <- paste(1,":",i, " vs " , i+1, ":", lenproc, sep="")
    }
    names(ts) <- namets
  }
  out <- list("Estimated Tau"=large, "Test Stat"=ts)
}

Estimate_Tau <- function(proc, tstat, dist_ic, ...){
  tau_stat <- NULL
  tau_stat[1] <- tstat(proc[,1], y=dist_ic, ...)
  tau_est <- Track_Stat_Over(proc, tstat)
  tau_stat <- append(tau_stat, tau_est[[2]])
  which.max(tau_stat)-1
}
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

Track_Stat <- function(proc,stat,doplot=FALSE){
  # Tracks the value of a statistic for a process
  # Inputs:
  # proc: A matrix containing the process
  # stat: the statistic tracked
  track <- vector(mode="list",length=dim(proc)[2])
  x <- apply(proc,2,stat)
  if(doplot){
    plot(seq(1:dim(proc)[1]),x)
    title("Statistic on Each Sample")
  }
}

Find_UCL <- function(Proc, stat, ARL=200){
  
}
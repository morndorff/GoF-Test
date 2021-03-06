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
  # Simulate a Incontrol Process Iteratively
  if(is.null(proc)){
    proc <- matrix(make_sample(num.samp, dist, params), nrow=num.samp, ncol=1)
    return(proc)
  }
  new_samp <- make_sample(num.samp, dist, params)
  proc <- cbind(proc, new_samp)
}

Sim_CP_Process_Iter <- function(proc=NULL, num.samp, cp, dist_one, param_one,
                    dist_two, param_two){
  # Args:
  # dist_one: "norm" "unif"
  # param_one: list(mean=0, sd=1)
  # cp: When does the distribution change.  If it is one, it does so immediately
  if(is.null(proc)){
    if(cp==1){
      proc <- matrix(make_sample(n=num.samp, dist=dist_two, params=param_two), nrow=num.samp, ncol=1)
      return(proc)     
    }
    proc <- matrix(make_sample(n=num.samp, dist=dist_one, params=param_one), nrow=num.samp, ncol=1)
    return(proc)
  }
  lenproc <- dim(proc)[2]
  if(lenproc < cp){
    new_samp <- make_sample(n=num.samp, dist=dist_one, params=param_one)
    proc <- cbind(proc, new_samp)
    return(proc)
  }
  if(lenproc >= cp){
    new_samp <- make_sample(n=num.samp, dist=dist_two, params=param_two)
    proc <- cbind(proc, new_samp)
    return(proc)
  }
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

Process_Stat <- function(proc, tstat, dist_ic, params, doplot=FALSE, detail=FALSE){
  # Tracks the value of a statistic for a process
  # Inputs:
  # proc: A matrix containing the process
  # tstat: A test statistic
  # dist_ic: The incontrol distribution - for wavelet and ks, should be "pnorm" or similar
  # params: parameters of the in control distribution
  lenproc <- dim(proc)[2] # run length of the process
  
  # If the proc is initially a vector, its dim value will be 0, so we do this
  if(is.null(lenproc)) {
    ts <- NULL
    ts[1] <- 0
    ts[2] <- do.call(tstat, c(list(proc), list(dist_ic), params))
    #ts[2] <- tstat(proc, dist_ic, ...)
    ts[3] <- ts[2]
    STAT <- ts[3]
    return(STAT)
  }
  
  ts <- matrix(nrow=3,ncol=lenproc)
  tau <- 0
  ts[1, (tau + 1)] <- 0
  ts[2, (tau + 1)] <- do.call(tstat, c(list(proc[, (tau+1):lenproc]), list(dist_ic), params))
  # ts[2, (tau + 1)] <- tstat(proc[, (tau+1):lenproc], dist_ic, ...)
  ts[3, (tau + 1)] <- ts[2, (tau+1)] - ts[1, (tau+1)]
  if(lenproc!=1) {
  for(tau in 1:(lenproc-1)){
    ts[1, (tau+1)] <- do.call(tstat, c(list(proc[, 1:tau]), list(dist_ic), params)) #g(H0, tau-)
    #ts[1, (tau+1)] <- tstat(proc[, 1:tau], dist_ic, ...) #g(H0, tau-)
    ts[2, (tau+1)] <- do.call(tstat, c(list(proc[, (tau+1):lenproc]), list(dist_ic), params)) #g(H0, tau+)
    #ts[2, (tau+1)] <- tstat(proc[, (tau+1):lenproc], dist_ic, ...) #g(H0, tau+)
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

Find_IC_RL_Slow <- function(num.samp=32, 
                            dist="norm", params=list(mean=0, sd=1), 
                            tstat=wave.energy, UCL){
  # num.samp - Number of Samples
  # dist - Incontrol distribution
  # params - Incontrol distribution parameters
  # tstat - Test Statistic for the process
  # UCL - Vector containing Upper Control Limits
  Proc <- NULL # Initializing
  count <- 1
  track_stat <- NULL
  tstat_proc <- 0
  dist_ic <- paste("p", dist, sep="")
  while(tstat_proc < max(UCL)){
    Proc <- Sim_IC_Process_Iter(proc=Proc, num.samp=num.samp, dist=dist, 
                                param=params) # good
    tstat_proc <- Process_Stat(proc=Proc, tstat=tstat, 
                               dist_ic=dist_ic, params=params)
    track_stat <- append(track_stat, tstat_proc)
    count <- count +1
  }
  RL <- sapply(UCL, function(x) min(which(track_stat >= x)))
  res <- list("h(t)"=track_stat, "Length of Process"=count, "RL for Corresponding UCL"=RL, "UCLS"=UCL)
  return(res)
}

Find_CP_RL_Slow <- function(num.samp=32, dist_one="norm", param_one=list(mean=0, sd=1),
                            dist_two="norm", param_two=list(mean=0, sd=2), cp=1,
                            tstat=wave.energy, UCL){
  # num.samp - Number of Samples
  # dist - Incontrol distribution
  # params - Incontrol distribution parameters
  # tstat - Test Statistic for the process
  # UCL - Vector containing Upper Control Limits
  Proc <- NULL # Initializing
  count <- 1
  track_stat <- NULL
  tstat_proc <- 0
  dist_ic <- paste("p", dist_one, sep="")
  while(tstat_proc < max(UCL)){
    Proc <- Sim_CP_Process_Iter(proc=Proc, num.samp=num.samp, cp=cp, 
                                dist_one=dist_one, param_one=param_one,
                                dist_two=dist_two, param_two=param_two)# good
    tstat_proc <- Process_Stat(proc=Proc, tstat=tstat, 
                               dist_ic=dist_ic, params=param_one)
    track_stat <- append(track_stat, tstat_proc)
    count <- count +1
  }
  RL <- sapply(UCL, function(x) min(which(track_stat >= x)))
  res <- list("h(t)"=track_stat, "Length of Process"=count, "RL for Corresponding UCL"=RL, "UCLS"=UCL)
  return(res)
}

Find_IC_RL_Fast <- function(num.samp=32, 
                            dist="norm", params=list(mean=0, sd=1), 
                            tstat=wave.energy,
                            UCL, detail=FALSE, weight=FALSE){
  # num.samp - numeric, # size of RS from each dist
  # dist - Incontrol Distribution - prefer string
  # params - Incontrol Distribution Parameters -- MUST BE LIST
  # tstat - Test Statistic to use, i.e. ks.res.simp, wave.den
  # UCL - vector of potential upper control limits -- max controls ends of loop
  # detail - Gives more output
  # weight - Weighting scheme used. Currently if weight=FALSE, then uses sample mean
  Proc <- NULL
  h_t <- 0
  h_t_new <- 0
  g_t_p <- NULL
  g_t_m <- NULL
  which_tau <- NULL
  which_tau_weight <- NULL
  dist_ic <- paste("p", dist, sep="")
  while(h_t_new < max(UCL)){
    Proc <- Sim_IC_Process_Iter(proc=Proc, num.samp=num.samp, dist=dist, 
                                param=params) # New Realization @ Time T
    Time <- dim(Proc)[2] # Can probably just be replaced by count var
    g_t <- update_g_t(nvec=Proc[, Time], 
                      theta_p = g_t_p, 
                      theta_m = g_t_m, 
                      tstat = tstat, 
                      dist_ic = dist_ic,
                      params = params) # get new estimates for g+ and g-
    g_t_p <- g_t[, 1]
    g_t_m <- g_t[, 2]
    # Scaling
    if(weight==TRUE){
      pos_scale <- seq(1/Time, 1, length.out=Time)
      neg_scale <- rev(pos_scale)
    }
    if(weight==FALSE){
      pos_scale <- rev(1 / 1:Time)
      neg_scale <- 1 / 1:Time
    }
    h_g_p <- g_t_p * pos_scale 
    h_g_m <- g_t_m * neg_scale
    which_tau <- append(which_tau, which.max(h_g_p - h_g_m)) # Diagnostics for Tau selection
    which_tau_weight <- append(which_tau_weight, tail(which_tau/Time,1))
    h_t <- append(h_t, max(h_g_p - h_g_m)) # For finding RLs
    h_t_new <- tail(h_t, 1) # When does loop stop?
  }
  RL <- sapply(UCL, function(x) min(which(h_t >= x))) # Finding RLs to Corresponding UCL
  lenproc <- dim(Proc)[2]
  if(detail){
    res <- list("h_t"=h_t, 
                "Length of Process"=lenproc, 
                "RL for Corresponding UCL"=RL, 
                "UCLs"=UCL,
                "hgp"=h_g_p,
                "hmp"=h_g_m,
                "Tau Chosen at each Time Period"=which_tau,
                "Tau Chosen Divided by Time" = which_tau_weight)
    return(res)
  }
  res <- list("h_t"=h_t, "Length of Process"=lenproc, "RL for Corresponding UCL"=RL, "UCLs"=UCL)
  return(res)
}

Find_CP_RL_Fast <- function(num.samp=32, dist_one="norm", param_one=list(mean=0, sd=1),
                            dist_two="norm", param_two=list(mean=0, sd=2), cp=1,
                            tstat=wave.energy, UCL, detail=FALSE, weight=TRUE){
  Proc <- NULL
  h_t <- 0
  h_t_new <- 0
  g_t_p <- NULL
  g_t_m <- NULL
  which_tau <- NULL
  which_tau_weight <- NULL
  pdist_one <- paste("p", dist_one, sep="")
  while(h_t_new < max(UCL)){
    Proc <- Sim_CP_Process_Iter(proc=Proc, num.samp=num.samp, cp=cp, 
                                dist_one=dist_one, param_one=param_one,
                                dist_two=dist_two, param_two=param_two) # New Realization @ Time T
    Time <- dim(Proc)[2] # Can probably just be replaced by count var
    g_t <- update_g_t(nvec=Proc[, Time], 
                      theta_p = g_t_p, 
                      theta_m = g_t_m, 
                      tstat = tstat, 
                      dist_ic = pdist_one,
                      params=param_one) # get new estimates for g+ and g-
    g_t_p <- g_t[, 1]
    g_t_m <- g_t[, 2]
    if(weight==TRUE){
      pos_scale <- seq(1/Time, 1, length.out=Time)
      neg_scale <- rev(pos_scale)
    }
    if(weight==FALSE){
      pos_scale <- rev(1 / 1:Time)
      neg_scale <- 1 / 1:Time
    }
    h_g_p <- g_t_p * pos_scale
    h_g_m <- g_t_m * neg_scale
    which_tau <- append(which_tau, which.max(h_g_p - h_g_m))
    which_tau_weight <- append(which_tau_weight, tail(which_tau/Time,1))
    h_t <- append(h_t, max(h_g_p - h_g_m))
    #h_t <- append(h_t, find_max_dif(theta)) # New estimates for h(t)
    h_t_new <- tail(h_t, 1) # For finding RLs
  }
  RL <- sapply(UCL, function(x) min(which(h_t >= x))) # Finding RLs to Corresponding UCL
  lenproc <- dim(Proc)[2]
  if(detail){
    res <- list("h_t"=h_t, 
                "Length of Process"=lenproc, 
                "RL for Corresponding UCL"=RL, 
                "UCLs"=UCL,
                "hgp"=h_g_p,
                "hmp"=h_g_m,
                "Tau Chosen at each Time Period"=which_tau,
                "Tau Chosen Divided by Time" = which_tau_weight)
    return(res)
  }
  res <- list("h_t"=h_t, "Length of Process"=lenproc, "RL for Corresponding UCL"=RL, "UCLs"=UCL)
  return(res)
}

update_g_t <- function(nvec, theta_p, theta_m, tstat, dist_ic, params, ...){
  # Updates 'Fast' Algorithm for changepoint process
  # Input:
  # rlength: T, time
  # nvec: new vector recieved at time T
  # theta_p tau+ estimates
  # theta_m tau- estimates
  # output: 2 column matrix of updated parameters 
  rlength <- length(theta_p) +1 
  theta_p <- append(theta_p, 0)
  theta_m <- append(theta_m, 0)# 
  nvec_null <- do.call(tstat, c(list(nvec), list(dist_ic), params))
  #nvec_null <- tstat(nvec, dist_ic, ...)
  if(rlength==1){
    # Special case when T=1
    # TODO: Averaging
    # Generate random sample from null distribution here, test against null
    rdist_ic <- dist.conv.str(dist_ic, "r") # Convert from pnorm to rnorm
    y <- get(rdist_ic, mode="function", envir=parent.frame())
    rvec <- y(length(nvec), ...)
    theta_m <- tstat(rvec, dist_ic, ...)
    theta_p <- nvec_null
    return(cbind(theta_p,theta_m))
  }
  theta_m[rlength] <- theta_m[(rlength-1)] + theta_p[(rlength-1)] # Will Break if T=1
  theta_p <- theta_p + nvec_null
  return(cbind(theta_p,theta_m))
}

# Temporary function
ARL_Proc <- function(UCL, num.samp = 30,
                     time=60, 
                     method=Find_IC_RL_Fast, 
                     tstat=wave.energy,
                     dist="norm",
                     params=list(mean=0, sd=1), 
                     ...){
  ptm <- proc.time()
  RLs <- vector(mode="list", length=0)
  e_time <- 0
  len_UCL <- length(UCL)
  # Calculating ARLs
  while(e_time < time){
    RLs_det <- method(num.samp=num.samp, 
                      dist=dist, 
                      params=params,
                      tstat=tstat, 
                      UCL=UCL, ...)
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

find_max_dif <- function(theta){
  # NOT Absolute Values
  h_t <- max(theta[,1]-theta[,2])
  return(h_t)
}
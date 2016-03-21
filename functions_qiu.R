# Qiu 2009 Functions

get_y <- function(x, boundaries, sum){
  int <- findInterval(x, boundaries)
  tab <- tabulate(int + 1, nbins=length(boundaries) + 1) # m*f
  emp_probs <- tab/length(x)  #f0
  return(emp_probs)
}


get_exact <- function(num_bps){
  exact_probs <- rep(1/num_bps, num_bps)
  return(exact_probs)
}

get_g <- function(data, quan, past_cols=0){
  #pastcols must be b/w 0 and dim(data)[2] - 1
  num_col <- dim(data)[2]
  num_quan <- length(quan)
  
  if(is.null(num_col) || num_col==1){
    if(past_cols!=0) stop("Too many columns on vector!")
    data_subset <- data
  }else{
    data_subset <- data[, (num_col - past_cols):num_col] 
  }
  int <- findInterval(data_subset, quan)
  tab <- tabulate(int + 1, nbins=num_quan + 1) # b/c needs positive, see ?tabulate for details
}

qiu_ARL <- function(ic_data=rnorm(500), kP=.1, num_bps=4, control_limit, m=5, exact=FALSE, s=.01) {
  #ic_data: a sample of in control data
  # kP: allowance parameter (Qiu says ~.05 is good)
  # num_bps: number of break points.  (Poorly named)
  # control_limit: upper control limit (manual input)
  # m: batch size
  # exact: Calculate based on N(0,1) quantiles. Will have to extend later.
  
  # Checked 3/2: Sobs and Sexp are calculated correctly
  
  ic_data_length <- length(ic_data)
  
  if(is.null(num_bps)){
    num_bps <- floor(sqrt(ic_data_length))
  }
  
  if(exact){
    stop("Not Working Now")
    ic_probs <- get_exact(num_bps + 1)
    boundaries <- get_exact_boundaries(num_bps + 1)
    
  }else{
    boundaries <- quantile(ic_data, probs=seq(1/(num_bps + 1), (num_bps)/(num_bps+1), length.out=num_bps))
    ic_probs <- get_y(ic_data, boundaries) #f0
  }

  data <- NULL

  # Initializing S_exp and S_obs
  S_obs <- matrix(0, nrow=(num_bps + 1), ncol=1)
  S_exp <- matrix(0, nrow=(num_bps + 1), ncol=1)
  u <- 0
  mf0 <- m * ic_probs
  i <- 0
  num_bins <- num_bps + 1
  
  while(tail(u, 1) < control_limit) {
    i <- i + 1
    data <- cbind(data, rnorm(m))
    g_n <- get_g(data = data[, i], quan=boundaries, past_cols=0) #g(n)
    g_n <- g_n + rnorm(num_bins, 0, s)
    # print(g_n)
    # C1 <- t((S_obs[, i] - S_exp[, i]) + (g_n -  mf0)) # check tranposition
    # C2 <- diag(1 / (S_exp[, i] + mf0))
    # C3 <- t(C1)
    C_n <- sum(((S_obs[, i] - S_exp[, i]) + (g_n -  mf0))^2 / (S_exp[, i] + mf0))
    
    
    # C_n <- C1 %*% C2 %*% C3
    # C_n <- as.vector(C_n)
    
    
    # print(paste("Cn", C_n))
    # print(paste("Ct", Ct))
    
    if(C_n <= kP){
      S_o_new = numeric(num_bps + 1)
      S_e_new = numeric(num_bps + 1)
      print(paste("Cn < kP at time", i)) 
    }else{
      S_o_new <- (S_obs[, i] + g_n) * ((C_n - kP) / C_n)
      S_e_new <- (S_exp[, i] + mf0) * ((C_n - kP) / C_n)
    }
    S_obs <- cbind(S_obs, S_o_new) 
    S_exp <- cbind(S_exp, S_e_new) 
    if(all(S_o_new==0)){
      u <- append(u, 0)
    }else{
      U1 <- t(S_obs[, (i + 1)] - S_exp[, (i + 1)])
      U2 <- diag( 1 / (S_exp[, (i + 1)]))
      U3 <- t(U1)
    }
    u <- append(u, U1 %*% U2 %*% U3)
    #print(u)
    if(i ==103){
      #print(data)
      # print(g_n)
      print(paste("Cn =" , C_n))
      #print(S_obs)
      #print(S_exp)
    }
  }
  return(list("uP"=u, "Time OOC"=i))
}

qiu_Phase_II <- function(ic_data=rnorm(500), kP=.1, num_bps=10, control_limit=20, m=5, exact=TRUE, 
                         tau=0, 
                         ICdist="qnorm", IC_dist_ops=NULL,
                         OOCdist="qnorm", OOC_dist_ops=NULL){
  #ic_data: a sample of in control data
  # kP: allowance parameter (Qiu says ~.05 is good)
  # num_bps: number of break points.  (Poorly named)
  # control_limit: upper control limit (manual input)
  # m: batch size
  # exact: Calculate based on N(0,1) quantiles. Will have to extend later.
  # IC_dist_ops= NULL or list(mean=100, sd=2) or similar
  # ICdist/OOCdist= "qnorm" (for now, must be character
  
  # Checked 3/2: Sobs and Sexp are calculated correctly
  
  ic_data_length <- length(ic_data)
  
  if(is.null(num_bps)){
    num_bps <- floor(sqrt(ic_data_length))
  }
  
  if(exact){
    ic_probs <- get_exact(num_bps + 1)
    
    boundaries <- do.call(ICdist, c(list(seq(1/num_bps, (num_bps-1)/num_bps, length.out=num_bps)),
                                    IC_dist_ops))
  }else{
    boundaries <- quantile(ic_data, probs=seq(1/(num_bps + 1), (num_bps)/(num_bps+1), length.out=num_bps))
    ic_probs <- get_y(ic_data, boundaries) #f0
  }
  
  IC_gen <- dist.conv.str(ICdist, "r")
  OOC_gen <- dist.conv.str(ICdist, "r")
  rIC <- get(IC_gen, mode = "function", envir = parent.frame())
  rOOC <- get(OOC_gen, mode = "function", envir = parent.frame())
  # print(ic_emp_probs)
  # print(sum(ic_emp_probs))
  # boundaries <- quantile(ic_data, probs=seq(1/num_bps, (num_bps-1)/num_bps, length.out=num_bps))
  # print(boundaries)
  # int <- findInterval(ic_data, quan)
  m <- 5 # data size
  
  #data <- matrix(rnorm(m), nrow=m, ncol=1)
  data <- NULL
  #n <- 10 # current time
  
  # g1
  #g <- get_g(data=data, quan=boundaries, past_cols=0)
  
  # Initializing S_exp and S_obs
  S_obs <- matrix(0, nrow=(num_bps + 1), ncol=1)
  S_exp <- matrix(0, nrow=(num_bps + 1), ncol=1)
  u <- 0
  mf0 <- m * ic_probs
  #print(mf0)
  i <- 0
  print(boundaries)
  while(tail(u, 1) < control_limit) {
    i <- i + 1
    if(i <= tau){
      data <- cbind(data, do.call(rIC, c(list(m), IC_dist_ops) ) )
      #print(data)
    }else{
      data <- cbind(data, do.call(rOOC, c(list(m), OOC_dist_ops)) )
    }
    g_n <- get_g(data = data[, i], quan=boundaries, past_cols=0) #g(n)
    g_n <- g_n + rnorm(num_bins, 0, s)
    C1 <- t((S_obs[, i] - S_exp[, i]) + (g_n -  mf0)) # check tranposition
    #print(g_n)
    #print(mf0)
    #print(g_n-mf0)
    C2 <- diag(1 / (S_exp[, i] + mf0))
    # print(C2)
    C3 <- t(C1)
    C_n <- C1 %*% C2 %*% C3
    C_n <- as.vector(C_n)
    #print(C_n)
    if(C_n <= kP){
      S_o_new = numeric(num_bps+1)
      S_e_new = numeric(num_bps +1)
      
    }else{
      S_o_new <- (S_obs[, i] + g_n) * ((C_n - kP) / C_n)
      S_e_new <- (S_exp[, i] + mf0) * ((C_n - kP) / C_n)
    }
    S_obs <- cbind(S_obs, S_o_new) 
    #print(S_obs)
    S_exp <- cbind(S_exp, S_e_new) 
    #print(S_exp)
    U1 <- t(S_obs[, (i + 1)] - S_exp[, (i + 1)])
    U2 <- diag( 1 / (S_exp[, (i + 1)]))
    U3 <- t(U1)
    u <- append(u, U1 %*% U2 %*% U3)
    
  }
  # Diagnostics
  # print(S_obs)
  # return(list(data, boundaries))
  return(list("uP"=u, "Time OOC"=i))
}

qiu_L_ARL <- function(ic_data=rnorm(500), kL=1, num_bps=5, control_limit, m=5, exact=FALSE,
                      additive_constant=.001, s=.01) {
  #ic_data: a sample of in control data
  # kP: allowance parameter (Qiu says ~.05 is good)
  # num_bps: number of break points.  (Poorly named)
  # control_limit: upper control limit (manual input)
  # m: batch size
  # exact: Calculate based on N(0,1) quantiles. Will have to extend later.
  
  # Checked 3/2: Sobs and Sexp are calculated correctly
  
  ic_data_length <- length(ic_data)
  
  if(is.null(num_bps)){
    num_bps <- floor(sqrt(ic_data_length))
  }
  
  if(exact){
    ic_probs <- get_exact(num_bps + 1)
    boundaries <- get_exact_boundaries(num_bps + 1)
    
  }else{
    # Choose boundaries based on quantiles
    boundaries <- quantile(ic_data, probs=seq(1/(num_bps + 1), (num_bps)/(num_bps+1), length.out=num_bps))
    ic_probs <- get_y(ic_data, boundaries) #f0
  }


  data <- NULL
  # Initializing S_exp and S_obs
  S_obs <- matrix(0, nrow=(num_bps + 1), ncol=1)
  S_exp <- matrix(0, nrow=(num_bps + 1), ncol=1)
  u <- 0
  mf0 <- m * ic_probs
  #print(mf0)
  i <- 0

  while(tail(u, 1) < control_limit) {
    i <- i + 1

    data <- cbind(data, rnorm(m))
    g_n <- get_g(data = data[, i], quan=boundaries, past_cols=0) #g(n)
    if(all(S_obs[, i]==0)){
      S_obs[, i] <- S_obs[, i] + additive_constant # avoids probs with log(0)
    }
    #print(mf0)
    C1 <- 2 * t(S_obs[, i] + g_n)
    C2 <- log( (S_obs[, i] + g_n) / (S_exp[, i] + mf0))
    #print(S_obs[, i])
    #print(S_exp[, i])
    #print(S_obs[, i] + g_n)
    #print(S_exp[, i] + mf0)
    #print(C1)
    #print(C2)
    C_n <- C1 %*% C2
    C_n <- as.vector(C_n)
    
    #print(C_n)
    if(C_n <= kL){
      S_o_new = numeric(num_bps + 1)
      S_e_new = numeric(num_bps + 1)
    }else{
      S_o_new <- (S_obs[, i] + g_n) * ((C_n - kL) / C_n)
      S_e_new <- (S_exp[, i] + mf0) * ((C_n - kL) / C_n)
    }
    S_obs <- cbind(S_obs, S_o_new) 
    S_exp <- cbind(S_exp, S_e_new) 
    if(all(S_o_new==0)){
      u <- append(u, 0)
    }else{
      U1 <- t(S_obs[, (i + 1)])
      U2 <- log(S_obs[, (i+1)] / mf0)
      print(U1)
      print(U2)
      u_new <- 2 * U1 %*% U2
      u <- append(u, u_new)
    }
  }
  # Diagnostics
  # print(S_obs)
  # return(list(data, boundaries))
  return(list("uP"=u, "Time OOC"=i))
}




qiu_KS_ARL <- function(ic_data=rnorm(500), kK=.02, control_limit, m=5,
                       ICdist="rnorm", IC_dist_ops=NULL,
                       bootstrap_samples=3000, keep_data=FALSE) {
  #ic_data: a sample of in control data
  # kP: allowance parameter (Qiu says ~.05 is good)
  # num_bps: number of break points.  (Poorly named)
  # control_limit: upper control limit (manual input) (hK)
  # m: batch size
  # no exact option available here
  
  # Checked 3/2: Sobs and Sexp are calculated correctly
  
  fhat_ic <- ecdf(ic_data)

  # Calculate d0 via bootstrapping
  
  j <- 1:m
  D_n <- numeric(bootstrap_samples)
  for(i in 1:bootstrap_samples){
    f0 <- fhat_ic(sort.int(sample(ic_data, replace=T, size=m)))
    d_n1 <- f0 - (j-1)/m 
    d_n2 <- (j/m) - f0
    D <- max(d_n1, d_n2)
    D_n[i] <- D
  }
  d0=mean(D_n)
  
  # Initializing Variables  
  data <- NULL
  u <- 0
  i <- 0
  j <- 1:m
  
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  while(tail(u, 1) < control_limit) {
    i <- i + 1
    if(keep_data){ 
    data <- cbind(data, do.call(rIC, c(list(m), IC_dist_ops) ) )
    
    f0_n <- fhat_ic(sort.int(data[, i]))
    d_n1 <- f0_n - (j-1)/m 
    d_n2 <- (j/m) - f0_n
    D_n <- max(d_n1, d_n2)
    u_nK <- max(0,
                u[i] + (D_n - d0) - kK)
    u <- append(u, u_nK)
    }else{
    data <- do.call(rIC, c(list(m), IC_dist_ops)) 
    f0_n <- fhat_ic(sort.int(data))
    d_n1 <- f0_n - (j-1)/m 
    d_n2 <- (j/m) - f0_n
    D_n <- max(d_n1, d_n2)
    u_nK <- max(0,
                u[i] + (D_n - d0) - kK)
    u <- append(u, u_nK)
    }
  }
  return(list("uP"=u, "Time OOC"=i))
}


qiu_KS_PhaseII <- function(ic_data=rnorm(500), kK=.02, control_limit=20, m=5, exact=TRUE, 
                           tau=3, 
                           ICdist="rnorm", IC_dist_ops=NULL,
                           OOCdist="rnorm", OOC_dist_ops=NULL,
                           bootstrap_samples=3000,
                           keep_data=FALSE){
  
  fhat_ic <- ecdf(ic_data)

  # Calculate d0 via bootstrapping
  j <- 1:m
  D_n <- numeric(bootstrap_samples)
  for(i in 1:bootstrap_samples){
    f0 <- fhat_ic(sort.int(sample(ic_data, replace=T, size=m)))
    d_n1 <- f0 - (j-1)/m 
    d_n2 <- (j/m) - f0
    D <- max(d_n1, d_n2)
    D_n[i] <- D
  }
  d0=mean(D_n)
  
  # Initializing Variables  
  data <- NULL
  u <- 0
  i <- 0
  j <- 1:m
  
  rOOC <- get(OOCdist, mode = "function", envir = parent.frame())
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  while(tail(u, 1) < control_limit) {
    i <- i + 1
    if(keep_data){
    if(i <= tau){
      data <- cbind(data, do.call(rIC, c(list(m), IC_dist_ops) ) )

    }else{
      data <- cbind(data, do.call(rOOC, c(list(m), OOC_dist_ops)) )
    }
    f0_n <- fhat_ic(sort.int(data[, i]))
    d_n1 <- f0_n - (j-1)/m 
    d_n2 <- (j/m) - f0_n
    D_n <- max(d_n1, d_n2)
    u_nK <- max(0,
                u[i] + (D_n - d0) - kK)
    u <- append(u, u_nK)
    }else{
     if(i <= tau){
      data <-  do.call(rIC, c(list(m), IC_dist_ops) ) 

    }else{
      data <- do.call(rOOC, c(list(m), OOC_dist_ops)) 
    }
    f0_n <- fhat_ic(sort.int(data))
    d_n1 <- f0_n - (j-1)/m 
    d_n2 <- (j/m) - f0_n
    D_n <- max(d_n1, d_n2)
    u_nK <- max(0,
                u[i] + (D_n - d0) - kK)
    u <- append(u, u_nK) 
    }
  }
  return(list("uP"=u, "Time OOC"=i))
}


qiu_CVM_ARL <- function(ic_data, kK=.02, control_limit, m, bootstrap_samples=1000){
  
  
  j <- 1:m
  D_n <- numeric(bootstrap_samples)
  for(i in 1:bootstrap_samples){
    D <- CvMTwoSamp.res(ic_data, sample(ic_data, replace=T, size=m))
    D_n[i] <- D
  }
  d0=mean(D_n)
  
  # Initializing Variables  
  data <- NULL
  u <- 0
  i <- 0
  j <- 1:m
  
  while(tail(u, 1) < control_limit) {
    i <- i + 1
    
    data <- cbind(data, rnorm(m))
    
    D_n <- CvMTwoSamp.res(data[, i], ic_data)
    u_nK <- max(0,
                u[i] + (D_n - d0) - kK)
    u <- append(u, u_nK)
  }
  # Diagnostics
  # print(S_obs)
  # return(list(data, boundaries))
  return(list("uP"=u, "Time OOC"=i))
  
}


qiu_CVM_PhaseII <- function(ic_data=rnorm(500), kK=.02,  control_limit=20, m=5, exact=TRUE, 
                            tau=3, 
                            ICdist="rnorm", IC_dist_ops=NULL,
                            OOCdist="rnorm", OOC_dist_ops=NULL,
                            bootstrap_samples=1000){
  
  j <- 1:m
  D_n <- numeric(bootstrap_samples)
  for(i in 1:bootstrap_samples){
    D <- CvMTwoSamp.res(ic_data, sample(ic_data, replace=T, size=m))
    D_n[i] <- D
  }
  d0=mean(D_n)
  print(d0)
  print(sd(D_n)*1.96/sqrt(length(D_n)))
  
  # Initializing Variables  
  data <- NULL
  u <- 0
  i <- 0
  j <- 1:m 
  
  # Random variable generation functions
  rOOC <- get(OOCdist, mode = "function", envir = parent.frame())
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  while(tail(u, 1) < control_limit) {
    i <- i + 1
    if(i <= tau){
      data <- cbind(data, do.call(rIC, c(list(m), IC_dist_ops) ) )
      #print(data)
    }else{
      data <- cbind(data, do.call(rOOC, c(list(m), OOC_dist_ops)) )
    }
    D_n <- CvMTwoSamp.res(data[, i], ic_data)
    u_nK <- max(0,
                u[i] + (D_n - d0) - kK)
    u <- append(u, u_nK)
  }
  return(list("uP"=u, "Time OOC"=i))
}


R_EWMA_PhaseII <- function(lambda=.05, 
                           control_limit=.129375, 
                           m=5, 
                           tau=3, 
                           ICdist="rnorm", IC_dist_ops=NULL,
                           OOCdist="rnorm", OOC_dist_ops=NULL){
  
  
  
  # Initializing Variables  
  data <- NULL
  v <- 0 # Initialize v at 0
  i <- 0
  j <- 1:m
  
  rOOC <- get(OOCdist, mode = "function", envir = parent.frame())
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  while(tail(v, 1) < control_limit) {
    i <- i + 1
    if(i <= tau){
      data <- cbind(data, do.call(rIC, c(list(m), IC_dist_ops) ) )
      #print(data)
    }else{
      data <- cbind(data, do.call(rOOC, c(list(m), OOC_dist_ops)) )
    }
    
    x_bar <- mean(data[, i])
    v_N <- lambda * x_bar + (1-lambda) * tail(v, 1)
    v <- append(v, v_N)
  }
  return(list("v_EWMA"=v, "Time OOC"=i))
}

R_EWMA_Find_CL <- function(lambda=.05, 
                           control_limit=.12, 
                           m=5, 
                           ICdist="rnorm", 
                           IC_dist_ops=NULL){
  # Initializing Variables  
  data <- NULL
  v <- 0 # Initialize v at 0
  i <- 0
  j <- 1:m
  
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  while(tail(v, 1) < control_limit) {
    i <- i + 1
    data <- cbind(data, do.call(rIC, c(list(m), IC_dist_ops) ) )
    
    x_bar <- mean(data[, i])
    v_N <- lambda * x_bar + (1-lambda) * tail(v, 1)
    v <- append(v, v_N)
  }
  return(list("v_EWMA"=v, "Time OOC"=i))
}

Find_ARL <- function(arl=200, lcl=0, ucl=.15, N_max=15, tol=1, sd_tol=3,
                          ICdist="rnorm", IC_dist_ops=NULL,
                          f=R_EWMA_Find_CL, N2_min=300){
  # This can be adapted for a number of problems later
  # Works for any function where the output has "Time OOC"
  
  N <- 1
  current_lcl <- lcl
  current_ucl <- ucl
  arl_track <- arl + tol * 1.5 # will get overwritten later, jf 1st while iter
  
  while(N < N_max){
    
    # Is the arl found on the previous iteration near 200?
    # If yes, then output the list
    if(abs(arl_track - arl) < tol ){
      
      return(list("Calculated Control Limit"=new_cl, 
                  "Number of Iterations"=N,
                  "Calculated ARL"=arl_track,
                  "ARL_SD"=sd_arl))
    }
    
    
    new_cl <- mean(c(current_lcl, current_ucl)) # New control limit via bisection
    
    # Calculating control limit based on new control limit    
    # Two triggers: run for at least 500 iterations
    # and s.dev within sd_tol
    sd_arl <- sd_tol + 1
    arl_track <- NULL
    N2 <- 0
    while(sd_arl > sd_tol){
      new_arl <- f(control_limit=new_cl, ICdist = ICdist, IC_dist_ops = IC_dist_ops)
      
      new_arl <- new_arl[["Time OOC"]]
      
      arl_track <- append(arl_track, new_arl)
      sd_arl <- sd(arl_track)/sqrt(length(arl_track))
      if(is.na(sd_arl)) sd_arl <- sd_tol + 1
      N2 <- N2 + 1
      if(N2 < N2_min) sd_arl <- sd_tol + 1 # don't stop until N2min
    }
    # output mean of arl_track
    arl_track <- mean(arl_track) # f(new_cl)
    print(paste("New Control Limit", new_cl, "has ARL of:", arl_track))
    print(paste("This took", N2, "iterations"))
    # Create new estimates
    if(arl_track < arl){
      current_lcl <- new_cl
    }else{
      current_ucl <- new_cl
    }
    
    N <- N + 1
  }
  stop("Did Not Converge")
}


EWMA_KS_Find_CL <- function(ic_data=rnorm(500),
                            lambda=.05, 
                            control_limit=.12, 
                            m=5, 
                            ICdist="rnorm", 
                            IC_dist_ops=NULL,
                            bootstrap_samples=3000,
                            keep_data=FALSE){
  
  fhat_ic <- ecdf(ic_data)

  # Calculate d0 via bootstrapping
  j <- 1:m
  D_n <- numeric(bootstrap_samples)
  for(i in 1:bootstrap_samples){
    f0 <- fhat_ic(sort.int(sample(ic_data, replace=T, size=m)))
    d_n1 <- f0 - (j-1)/m 
    d_n2 <- (j/m) - f0
    D <- max(d_n1, d_n2)
    D_n[i] <- D
  }
  d0=mean(D_n)
  
  # Initializing Variables  
  data <- NULL
  u<- 0 
  i <- 0
  j <- 1:m
  
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  while(tail(u, 1) < control_limit) {
    i <- i + 1
    if(keep_data){
    data <- cbind(data, do.call(rIC, c(list(m), IC_dist_ops) ) )
    f0_n <- fhat_ic(sort.int(data[, i]))
    d_n1 <- f0_n - (j-1)/m 
    d_n2 <- (j/m) - f0_n
    D_n <- max(d_n1, d_n2)
    u_nK <- lambda*(D_n - d0) + (1 - lambda)*tail(u, 1)
    u <- append(u, u_nK)
    }else{
    data <- do.call(rIC, c(list(m), IC_dist_ops))
    f0_n <- fhat_ic(sort.int(data))
    d_n1 <- f0_n - (j-1)/m
    d_n2 <- (j/m) - f0_n
    D_n <- max(d_n1, d_n2)
    u_nK <- lambda*(D_n - d0) + (1 - lambda)*tail(u, 1)
    u <- append(u, u_nK)
    
    }   
  }
  return(list("u_EWMA_KS"=u, "Time OOC"=i))
}

EWMA_KS_PhaseII <- function(ic_data=rnorm(500),
                            lambda=.05, 
                           control_limit=0.033125, 
                           m=5, 
                           tau=3, 
                           ICdist="rnorm", IC_dist_ops=NULL,
                           OOCdist="rnorm", OOC_dist_ops=NULL,
                           bootstrap_samples=3000,
                           keep_data=FALSE){
  fhat_ic <- ecdf(ic_data)

  # Calculate d0 via bootstrapping
  j <- 1:m
  D_n <- numeric(bootstrap_samples)
  for(i in 1:bootstrap_samples){
    f0 <- fhat_ic(sort(sample(ic_data, replace=T, size=m)))
    d_n1 <- f0 - (j-1)/m 
    d_n2 <- (j/m) - f0
    D <- max(d_n1, d_n2)
    D_n[i] <- D
  }
  d0=mean(D_n)
  
  # Initializing Variables  
  data <- NULL
  u <- 0 # Initialize u at 0
  i <- 0
  j <- 1:m
  
  rOOC <- get(OOCdist, mode = "function", envir = parent.frame())
  rIC <- get(ICdist, mode = "function", envir = parent.frame())
  
  while(tail(u, 1) < control_limit) {
    i <- i + 1
    if(keep_data){
      if(i <= tau){
        data <- cbind(data, do.call(rIC, c(list(m), IC_dist_ops) ) )
      }else{
        data <- cbind(data, do.call(rOOC, c(list(m), OOC_dist_ops)) )
      }
      f0_n <- fhat_ic(sort.int(data[, i]))
      d_n1 <- f0_n - (j-1)/m 
      d_n2 <- (j/m) - f0_n
      D_n <- max(d_n1, d_n2)
      u_nK <- lambda*(D_n - d0) + (1 - lambda)*tail(u, 1)
      u <- append(u, u_nK)
    }else{
      if(i <= tau){
        data <- do.call(rIC, c(list(m), IC_dist_ops))
      }else{
        data <- do.call(rOOC, c(list(m), OOC_dist_ops))
      }
      f0_n <- fhat_ic(sort.int(data))
      d_n1 <- f0_n - (j-1)/m 
      d_n2 <- (j/m) - f0_n
      D_n <- max(d_n1, d_n2)
      u_nK <- lambda*(D_n - d0) + (1 - lambda)*tail(u, 1)
      u <- append(u, u_nK)
    }
  }
  return(list("u_EWMA"=u, "Time OOC"=i))
}

rmixnorm <- function(N, u1, u2, sd1, sd2){
  # Does some mixture distributions
  components <- sample(1:2,prob=c(.5, .5),size=N,replace=TRUE)
  mus <- c(u1, u2)
  sds <- c(sd1, sd2)
  samples <- rnorm(n=N,mean=mus[components],sd=sds[components]) 
  return(samples)
}


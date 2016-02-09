# Testing Uniformity
A_Uniform <- function(z, k=1){
  if(z < 0 || z > 1 || k<0 ) stop("z must be between 0 and 1, k>0")
  F_z <- 1-(1-z)^k
  F_z
}
rA_Uniform <- function(n, k=1){
  if(k<0 ) stop("z must be between 0 and 1, k>0")
  # Inverse Transform method
  y <- runif(n)
  z = 1-(1-y)^(1/k)
}


B_Uniform <- function(z,k=1){
  if(z < 0 || z > 1 || k<0 ) stop("z must be between 0 and 1, k>0")
  if(length(k) > 1) stop("Haven't vectorized appropriately")
  if(length(z) > 1){
    F_z <- vector(mode="numeric", length=length(z))
    for(i in seq_along(z)){
      if(z[i] <= .5){
        F_z[i] <- 2^(k-1)*z[i]^k
      }
      if(z[i] > .5){
        F_z[i] <- 1 - 2^(k - 1) * (1 - z[i])^k
      }
    }
    return(F_z)
  }
  
  if(z<=.5){
    F_z <- 2^(k - 1) * z^k
    return(F_z)
  }
  if(z >.5){
    F_z <- 1 - 2^(k - 1) * (1 - z)^k
    return(F_z)
  }
}

rB_Uniform <- function(n, k=1){
  if(k < 0) stop("z must be between 0 and 1, k>0")
  # Inverse Transform method
  y <- runif(n)
  z <- vector(mode="numeric", length=n)
  for(i in seq_along(z)){
    if(y[i] <= .5){
      z[i] <- (y[i]/(2^(k-1)))^(1/k)
    }
    if(y[i] > .5){
      z[i] <- 1 - ( (1-y[i]) / (2^(k - 1)) )^(1/k)
    }
  }
  return(z)
}




rC_Uniform <- function(n, k=1){
  if(k<0 ) stop("z must be between 0 and 1, k>0")
  # Inverse Transform method
  y <- runif(n)
  z <- vector(mode="numeric", length=n)
  for(i in seq_along(z)){
    if(y[i] <= .5){
      z[i] <- .5 - ((.5-y[i])/(2^(k-1)))^(1/k)
    }
    if(y[i] > .5){
      z[i] <- .5 + ((y[i]-.5) / (2^(k-1)))^(1/k)
    }
  }
  return(z)
}







C_Uniform <- function(z,k=1){
  if(z < 0 || z > 1 || k<0 ) stop("z must be between 0 and 1, k>0")
  if(length(k)>1) stop("Haven't vectorized appropriately")
  if(length(z) > 1){
    F_z <- vector(mode="numeric", length=length(z))
    for(i in seq_along(z)){
      if( z[i] <=.5)
        F_z[i] <- .5 - 2^(k - 1)*(.5 - z[i])^k
      if(z[i] > .5){
        F_z[i] <- .5 + 2^(k - 1)*(z[i] - .5)^k
      }
    }
    return(F_z)
  }
  if(z <= .5){
    F_z <- .5 - 2^(k-1)*(.5 - z)^k
    return(F_z)
  }
  if(z > .5){
    F_z <- .5 + 2^(k-1)*(z - .5)^k
    return(F_z)
  }
}
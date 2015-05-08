# Developing One Sample Test based on LaRiccia

fourier_test <- function(x, y, ..., doplot=FALSE){
  # Args: x - vector
  # y : pnorm or similar
  
  # Is equivalent to CvM
  
  # Checking Inputs
  if(is.numeric(y)) stop("One Sample Test Only")
  if (is.list(y)) 
    y <- names(y)
  if (is.function(y)) 
    funname <- as.character(substitute(y))
  if (is.character(y)) 
    funname <- y
  y <- get(funname, mode = "function", envir = parent.frame())
  if (!is.function(y)) 
    stop("'y' must be a function or a string naming a valid function")
  lenx <- length(x)
  x_unit <- y(x, ...)
  # x_unit
  # x_u_ecdf <- ecdf(x_unit)
  a_jn <- vector(mode="numeric", length=lenx)
  for(j in 1:lenx){
    a_jn[j] <- (sqrt(2) / lenx) * sum(cos(j * pi * x_unit))
  }
  
  C_n <- a_jn^2 * (1 / ((1:lenx) * pi)^2 )
  if(doplot){
    plot(C_n, type="l", xlab="J Index", ylab="Component Value")
    plot(cumsum(C_n), type="l", xlab="J Index", ylab="Cum Sum")
  }
  C_n <- lenx * sum(C_n)
  C_n
}

smooth_test <- function(x, y, ..., M=5, doplot=FALSE){
  # Args: x - vector
  # y : pnorm or similar
  
  # Checking Inputs
  if(is.numeric(y)) stop("One Sample Test Only")
  if (is.list(y)) 
    y <- names(y)
  if (is.function(y)) 
    funname <- as.character(substitute(y))
  if (is.character(y)) 
    funname <- y
  y <- get(funname, mode = "function", envir = parent.frame())
  if (!is.function(y)) 
    stop("'y' must be a function or a string naming a valid function")
  lenx <- length(x)
  x_unit <- y(x, ...)
  # x_unit
  # x_u_ecdf <- ecdf(x_unit)
  a_jn <- vector(mode="numeric", length=M)
  for(j in 1:M){
    a_jn[j] <- (sqrt(2) / lenx) * sum(cos(j * pi * x_unit))
  }
  T_nm <- lenx * sum(a_jn^2)
  T_nm
}

smooth_test2 <- function(x, y, ..., lambda=5, doplot=FALSE){
  # Args: x - vector
  # y : pnorm or similar
  # lambda: smoothing parameter > 0
  # Checking Inputs
  if(is.numeric(y)) stop("One Sample Test Only")
  if (is.list(y)) 
    y <- names(y)
  if (is.function(y)) 
    funname <- as.character(substitute(y))
  if (is.character(y)) 
    funname <- y
  y <- get(funname, mode = "function", envir = parent.frame())
  if (!is.function(y)) 
    stop("'y' must be a function or a string naming a valid function")
  lenx <- length(x)
  x_unit <- y(x, ...)
  # x_unit
  # x_u_ecdf <- ecdf(x_unit)
  a_jn <- vector(mode="numeric", length=lenx)
  for(j in 1:lenx){
    a_jn[j] <- (sqrt(2) / lenx) * sum(cos(j * pi * x_unit))
  }
  seq_n <- 1:lenx
  S_lambda <- lenx * sum(a_jn^2 / (1 + (lambda * seq_n^2))^2)
  S_lambda
}

fft_test <- function(x,y, ...){
  if(is.numeric(y)) stop("One Sample Test Only")
  if (is.list(y)) 
    y <- names(y)
  if (is.function(y)) 
    funname <- as.character(substitute(y))
  if (is.character(y)) 
    funname <- y
  y <- get(funname, mode = "function", envir = parent.frame())
  if (!is.function(y)) 
    stop("'y' must be a function or a string naming a valid function")
  lenx <- length(x)
  x_unit <- y(x, ...)
  fft(x_unit)
}

Fan_Test <- function(x, y, ...){
  if(is.numeric(y)) stop("One Sample Test Only")
  if (is.list(y)) 
    y <- names(y)
  if (is.function(y)) 
    funname <- as.character(substitute(y))
  if (is.character(y)) 
    funname <- y
  y <- get(funname, mode = "function", envir = parent.frame())
  if (!is.function(y)) 
    stop("'y' must be a function or a string naming a valid function")
  lenx <- length(x)
  x_unit <- y(x, ...)
  
  # Creating Sub Intervals
  nbin <- ceiling(sqrt(lenx)) #ad hoc. Will work decently well for samples around 100
  bp <- seq(0, 1, length.out=nbin) # Bins are equally spaced
  int <- findInterval(x_unit, bp) # Binning
  counts <- as.numeric(table(int))   # Gives counts in each interval
  print(counts)
  Y_j <- 2 * (sqrt(counts)- sqrt(lenx/nbin)) # Square root transform
  print(Y_j)
  # Now we take orthogonal transform
  X_lambda <- fft(Y_j)
  X_lambda
}
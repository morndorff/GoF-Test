# Test statistic functions
Max_Quan_TS <- function(x, y, ..., interp = 4, do.plot = FALSE) {
  # Computes maximum difference of quantiles Comments: OBC Using Linear Interpolation
  # of the ECDF Recall: ECDF range is [1/n,1]. This is probably not realistic Args: x:
  # A vector of observations for a R.V. (must be numeric) y: Either (1) Another vector
  # of observations (two sample) (2) A quantile function such as qnorm (one sample)
  # interp: method of interpolation used. For more details, see ?quantile do.plot:
  # Creates a plot illustrating the statistic Returns: The value of the statistic
  
  # Error Handling
  if (is.numeric(x) != TRUE) 
    stop("x must be numeric")
  
  x <- sort(x)
  lenx <- length(x)
  # Two Sample Test
  if (is.numeric(y)) {
    leny <- length(y)
    y <- sort(y)
    # x1 and y1 are the quantile values
    x1 <- seq(1/lenx, 1, 1/lenx)
    y1 <- seq(1/leny, 1, 1/leny)
    if (lenx == leny) {
      # In the case of equal sample sizes take max of abs. val
      z <- max(abs(y - x))
      if (do.plot == TRUE) 
        plot.ts.2sam(x, y, x1, y1, lenx, leny)
    } else if (lenx > leny) {
      # If there are unequal sample sizes, then interpolation is necessary
      q1 <- quantile(x, probs = x1, type = interp)
      q_inter <- approxfun(x1, q1, yleft = min(q1), yright = max(q1))
      z <- max(abs(q_inter(y1) - y))
      if (do.plot == TRUE) 
        plot.ts.2sam(x, y, x1, y1, lenx, leny)
    } else {
      # So length y>x
      q2 <- quantile(y, probs = y1, type = interp)
      q_inter <- approxfun(y1, q2, yleft = min(q2), yright = max(q2))
      z <- max(abs(q_inter(x1) - x))
      if (do.plot == TRUE) 
        plot.ts.2sam(x, y, x1, y1, lenx, leny)
    }
    return(z)
  }
  
  # One Sample
  if (is.list(y)) 
    y <- names(y)
  if (is.function(y)) 
    funname <- as.character(substitute(y))
  if (is.character(y)) 
    funname <- y
  y <- get(funname, mode = "function", envir = parent.frame())
  if (!is.function(y)) 
    stop("'y' must be numeric or a function or a string naming a valid function")
  z <- max(abs(x - y(seq(1/(lenx + 1), lenx/(lenx + 1), length.out = lenx), ...)))  #Note: quantiles up for debate
  if (do.plot == TRUE) {
    plot.ts.1sam(x, y, ..., funname = funname, lenx = lenx)
  }
  return(z)
}

Trap_Quan_Area_TS <- function(x, y, ..., interp = 4, do.plot = FALSE, size = 0.25) {
  # Computes area based on trapezoid areas Args: x: A vector of observations for a
  # R.V. (must be numeric) y: Either (1) Another vector of observations (two sample)
  # (2) A quantile function such as qnorm (one sample) size: controls height of
  # trapezoids Returns: The value of the statistic
  
  # Error Handling
  if (is.numeric(x) != TRUE) 
    stop("x must be numeric")
  
  x <- sort(x)
  lenx <- length(x)
  
  # One Sample
  if (is.function(y)) 
    funname <- as.character(substitute(y))
  if (is.character(y)) 
    funname <- y
  y <- get(funname, mode = "function", envir = parent.frame())
  if (!is.function(y)) 
    stop("'y' must be numeric or a function or a string naming a valid function")
  # Quantiles
  x1 <- seq(1/(lenx + 1), lenx/(lenx + 1), length.out = lenx)
  # Parameter which determines height of trapezoids
  delta <- size * (1/lenx)
  # Creating quantile function q1 <- quantile(x, probs = x1, type = interp) q_inter <-
  # approxfun(x1, q1, yleft = min(q1), yright = max(q1)) Length of horizontal segments
  x_p_d <- x + (y(x1 + delta, ...) - y(x1, ...))
  x_m_d <- x + (y(x1 - delta, ...) - y(x1, ...))
  y_p_d <- y(x1 + delta, ...)
  y_m_d <- y(x1 + delta, ...)
  # sum of trapezoid areas
  z <- sum((abs(x_p_d - y_p_d) + abs(x_m_d - y_m_d)) * delta)
  if (do.plot == TRUE) {
    plot.ts.1sam(x, y, ..., funname = funname, lenx = lenx)
    
  }
  return(z)
}
Com_KS_Max_Quan_TS <- function(x, y, ...) {
  obc.stat <- Max_Quan_TS(x, y)
  ks.stat <- ks.res.simp(x, y)
  if (obc.stat > ks.stat) {
    z <- list(obc.stat, c("OBC"))
  } else if (ks.stat > obc.stat) {
    z <- list(ks.stat, c("KS"))
  } else z <- list(c(NULL), c("NEITHER"))
  
  return(z)
}
Com_KS_Max_Quan_TS_Simp <- function(x, y, ...) {
  obc.stat <- myts(x, y, ...)
  ks.stat <- ks.res.simp(x, y, ...)
  if (obc.stat > ks.stat) {
    z <- obc.stat
  } else if (ks.stat > obc.stat) {
    z <- ks.stat
  } else z <- 1000  #######FIX THIS LATER
  
  return(z)
}
Com_KS_Max_Quan_Range_TS <- function(x, y) {
  # Returns value of the statistic and whether or not it from OBC
  com <- c(x, y)
  range <- max(com) - min(com)
  obc.stat <- Max_Quan_TS(x, y)/range
  ks.stat <- ks.res(x, y)
  if (obc.stat > ks.stat) {
    z <- list(obc.stat, c("OBC"))
  } else if (ks.stat > obc.stat) {
    z <- list(ks.stat, c("KS"))
  } else z <- list(c(NULL), c("NEITHER"))
  
  return(z)
}
Com_KS_Max_Quan_Range_TS_Simp <- function(x, y) {
  # Returns value of the statistic and whether or not it from OBC
  com <- c(x, y)
  range <- max(com) - min(com)
  obc.stat <- Max_Quan_TS(x, y)/range
  ks.stat <- ks.res.simp(x, y)
  z <- max(obc.stat, ks.stat)
  return(z)
}

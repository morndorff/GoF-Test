# New Function, takes as input two samples and a delta value?
source("functions.R")
myts.quad <- function(x, y, ..., interp = 4, do.plot = FALSE, size = .25) {
  # Computes area based on trapezoid areas
  # Args: 
  # x: A vector of observations for a R.V. (must be numeric)
  #
  # y: Either 
  # (1) Another vector of observations (two sample) 
  # (2) A quantile function such as qnorm (one sample)
  #
  # interp: method of interpolation used. For more details, see ?quantile 
  # do.plot: Creates a plot illustrating the statistic 
  #
  # Returns: The value of the statistic
  
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
  # quantiles
  x1 <- seq(1/(lenx+1),lenx/(lenx+1),length.out=lenx)
  # Parameter which determines height of trapezoids
  delta <- size * (1 / lenx)
  # Creating quantile function
  q1 <- quantile(x, probs = x1, type = interp)
  q_inter <- approxfun(x1, q1, yleft = min(q1), yright = max(q1))
  # Length of horizontal segments
  x_p_d <- q_inter(x1 + delta)
  x_m_d <- q_inter(x1 - delta)
  y_p_d <- y(x1 + delta)
  y_m_d <- y(x1 + delta)
  # sum of trapezoid areas
  z <- sum((abs(x_p_d - y_p_del) + abs(x_m_d - y_m_del)) * delta)
     return(z)
}

q1 <- quantile(x, probs = x1, type = interp)
q_inter <- approxfun(x1, q1, yleft = min(q1), yright = max(q1))
z <- max(abs(q_inter(y1) - y))



# Calculating the area of a quadrilateral
quad.area <- function(x1, x2, y1, y2){
  t1 <- tri.area(x1, x2, y1)
  t2 <- tri.area(x2, y1, y2)
  area <- t1 + t2
  return(area)
}

tri.area <- function(x, y, z) {
  area <- 0.5 * abs((x[1] - z[1]) * (y[2] - x[1]) - (x[1] - y[1]) * (z[2] - x[2]))
  return(area)
}


# Testing

x <- c(1,1)
y <- c(2,2)
z <- c(3,1)
z1 <- c(4,2)
tri.area(x,y,z)
quad.area(x,y,z,z1)


# Testing the trapezoid code
x <- rnorm(100)
y <- rnorm(100)
myts.par(x,y)
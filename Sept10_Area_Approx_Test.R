# Testing the Area Approximation Function
# September 9th, 2014


Density_Estimate <- function(x,interp=4){
  # Take in data (x)
  # Output: Function 
  lenx <- length(x)
  p_x <- seq(1/lenx, 1, 1/lenx)
  q1 <- quantile(x, probs = p_x, type = interp)
  dens_est <- approxfun(p_x, q1, yleft = min(q1), yright = max(q1))
}

Kernel_Density_Estimate <- function(x){
  pdf <- density(x)
  f <- approxfun(pdf$x, pdf$y, yleft=0, yright=0)
  function(lim) {
    cdf <- integrate(f, -Inf, lim)
  }
}


Delta_Calc <- function(x, y, lenx, leny){
  x_min <- min(abs(diff(x)))
  y_min <- min(abs(diff(y)))
  delta <- min(x_min, y_min)
}

# Let's do one sample test first

Quad_Quan_Area_TS <- function(x, y, ..., interp = 4, do.plot = FALSE, size = 0.25) {
  # Computes area based on trapezoid areas 
  # Args: x: A vector of observations for a R.V. (must be numeric) 
  # y: Either (1) Another vector of observations (two sample)
  #           (2) A quantile function such as qnorm (one sample) 
  # size: controls height of trapezoids 
  # Returns: The value of the statistic
  
  # Error Handling
  if (is.numeric(x) != TRUE) 
    stop("x must be numeric")

  x <- sort(x)
  lenx <- length(x)
  y <- sort(y)
  leny <- length(y)
  
  # Placeholder for Delta
  delta <- Delta_Calc(x,y)
  
  
  # Two Sample Test
  
  # Density estimate for x and y samples
  x_density <- Kernel_Density_Estimate(x)
  y_density <- Kernel_Density_Estimate(y)
  
  # Calculating Test Statistic
  x_minus_delta <- x - delta
  x_plus_delta <- x + delta
  
  #p_x_minus_delta <- x_density(x_minus_delta)$value
  p_x_minus_delta <- sapply(x_minus_delta, function(x) x_density(x)$value)

  #p_x_plus_delta <- x_density(x_plus_delta)$value
  p_x_plus_delta <- sapply(x_plus_delta, function(x) x_density(x)$value)
  
  y_minus_delta <- y - delta
  y_plus_delta <- y + delta
  
  #p_y_minus_delta <- y_density(y_minus_delta)$value
  p_y_minus_delta <- sapply(y_minus_delta, function(y) y_density(y)$value)
  
  #p_y_plus_delta <- y_density(y_plus_delta)$values
  p_y_plus_delta <- sapply(y_minus_delta, function(y) y_density(y)$value)
  
  x1 <- as.list(as.data.frame(rbind(x_minus_delta, p_x_minus_delta)))
  x2 <- as.list(as.data.frame(rbind(x_plus_delta, p_x_plus_delta)))
  y1 <- as.list(as.data.frame(rbind(y_minus_delta, p_y_minus_delta)))
  y2 <- as.list(as.data.frame(rbind(y_plus_delta, p_y_plus_delta)))
  
  coords <- rbind(x1,x2,y1,y2)
  coords <- as.list(as.data.frame(coords))
  areas <- lapply(coords, function(x) quad.area(x[[1]],x[[2]],x[[3]],x[[4]]))
  total_area <- sum(unlist(areas))
}

x <- rnorm(100)
y <- rnorm(100)
results <- Quad_Quan_Area_TS(x,y)
results

x <- rnorm(100)
y <- runif(100,0,1)
results <- Quad_Quan_Area_TS(x,y)
results

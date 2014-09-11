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
  # Handling String/Function Inputs
  if (is.function(y)) 
    funname <- as.character(substitute(y))
  if (is.character(y)) 
    funname <- y
  y <- get(funname, mode = "function", envir = parent.frame())
  if (!is.function(y)) 
    stop("'y' must be numeric or a function or a string naming a valid function")
  
  x <- sort(x)
  lenx <- length(x)
  y <- sort(y)
  leny <- length(y)
  
  # Two Sample Test
  
  # Density estimate for x and y samples
  x_density <- Density_Estimate(x)
  y_density <- Density_Estimate(y)
  
  # Calculating Test Statistic
  x_minus_delta <- x - delta
  x_plus_delta <- x + delta
  
  p_x_minus_delta <- x_density(x_minus_delta)
  p_x_plus_delta <- x_density(x_plus_delta)
  
  y_minus_delta <- y - delta
  y_plus_delta <- y + delta
  
  p_y_minus_delta <- y_density(y_minus_delta)
  p_y_plus_delta <- y_density(y_plus_delta)
  
  x1 <- as.list(as.data.frame(rbind(x_minus_delta, p_x_minus_delta)))
  x2 <- as.list(as.data.frame(rbind(x_plus_delta, p_x_plus_delta)))
  y1 <- as.list(as.data.frame(rbind(y_minus_delta, p_y_minus_delta)))
  y2 <- as.list(as.data.frame(rbind(y_plus_delta, p_y_plus_delta)))
  
  coords <- rbind(x1,x2,y1,y2)
  coords <- as.list(as.data.frame(coords))
  areas <- lapply(coords, function(x) quad.area(x[[1]],x[[2]],x[[3]],x[[4]]))
  total_area <- sum(areas)
}
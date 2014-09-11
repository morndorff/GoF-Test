# Area Approximation Pt. 2
# Delta around x-axis

# GOALS: Prototype area approximation script
# GOALS: Find good delta?
# Goals: Find mean of kth order statistics, k in 1,...n

# Finding mean of kth order statistic
x <- rt(7, 3)
sortx <- sort(x)

samp_size <- 7
num_sims <- 10000

x <- matrix(rt(samp_size * num_sims, 3), nrow=samp_size, ncol=num_sims)
sortx <- apply(x, 2, sort)
x_means <- rowMeans(sortx)
x_means

# Compare this to:
x_seq <- seq(1/(samp_size+1),samp_size/(samp_size+1),length.out=7)
x_quantiles <- qt(x_seq,df=3)
x_quantiles

samp_size <- 50
num_sims <- 10000

x <- matrix(rt(samp_size * num_sims, 3), nrow=samp_size, ncol=num_sims)
sortx <- apply(x, 2, sort)
x_means <- rowMeans(sortx)
x_means

# Compare this to:
x_seq <- seq(1/(samp_size+1),samp_size/(samp_size+1),length.out=samp_size)
x_quantiles <- qt(x_seq,df=3)
x_quantiles

# Less error as n increases





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




# Calculating Delta


# Calculating Test Statistic Given Delta
x_minus_delta <- x - delta
x_plus_delta <- x + delta

p_x_minus_delta <- y(x_minus_delta)
p_x_plus_delta <- y(x_plus_delta)

# Calculating Mean0 of the Order Statistics











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
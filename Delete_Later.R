# Delete Later
# rm(list=ls())
set.seed(5)


# Adding Plot to Area Approximation Test
quad.area <- function(x1, x2, y1, y2) {
  t1 <- tri.area(x1, x2, y1)
  t2 <- tri.area(x2, y1, y2)
  area <- t1 + t2
  return(area)
}
tri.area <- function(x, y, z) {
  area <- 0.5 * abs((x[1] - z[1]) * (y[2] - x[1]) - (x[1] - y[1]) * (z[2] - x[2]))
  return(area)
}

Kernel_CDF_Estimate <- function(x){
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


Quad_Quan_Area_TS <- function(x, y, ..., interp = 4, do.plot = FALSE, size = 0.25, opt=FALSE) {
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
  x_density <- Kernel_CDF_Estimate(x)
  y_density <- Kernel_CDF_Estimate(y)
  
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
  p_y_plus_delta <- sapply(y_plus_delta, function(y) y_density(y)$value)
  
  x1 <- as.list(as.data.frame(rbind(x_minus_delta, p_x_minus_delta)))
  x2 <- as.list(as.data.frame(rbind(x_plus_delta, p_x_plus_delta)))
  y1 <- as.list(as.data.frame(rbind(y_minus_delta, p_y_minus_delta)))
  y2 <- as.list(as.data.frame(rbind(y_plus_delta, p_y_plus_delta)))
  
  coords <- rbind(x1,x2,y2,y1)
  coords <- as.list(as.data.frame(coords))
  if(opt=="coords"){
    return(coords)
  }
  areas <- lapply(coords, function(x) quad.area(x[[1]],x[[2]],x[[3]],x[[4]]))
  total_area <- sum(unlist(areas))
}


x <- rnorm(3)
y <- rnorm(3)
x <- sort(x)
y <- sort(y)
lenx <- length(x)
leny <- length(y)
x_density <- Kernel_CDF_Estimate(x)
y_density <- Kernel_CDF_Estimate(y)
a <- Quad_Quan_Area_TS(x,y)

x_points <- seq(min(x)-sd(x), max(x)+sd(x),length.out=500)
x_out <- sapply(x_points, function(x) x_density(x)$value)

y_points <- seq(min(y)-sd(y), max(y)+sd(y),length.out=500)
y_out <- sapply(y_points, function(x) y_density(x)$value)

plot(x_points,x_out, xlim = c(min(x[1]-sd(x), y[1]-sd(y)), max(x[lenx]+sd(x), y[lenx]+sd(y))), ylab = "Probs", xlab = "Data", 
     main = "ECDF and Points", type="l")  #plot 1st sample points
coords <- Quad_Quan_Area_TS(x,y,opt="coords")
coords[[1]]
 lines(y_points,y_out,col="red")
 points(-1.67,0.133)
 points(-0.84,0.40)
 points(-1.017,.097)
 points(-.188,.362)

plotx <- sapply(coords[[1]],function(x) x[1])  
ploty <- sapply(coords[[1]],function(x) x[2])  
polygon(plotx,ploty, col="blue")





#points(y, y1, col = "red")  #put 2nd sample on graph






plot.ts.2sam <- function(x, y, x1, y1, lenx, leny) {
  if (lenx == leny) {
    par(mfrow = c(2, 1))
    # Making Quantile Graph: Note: Because Lengths are the same, don't need
    # iterpolation. (use type=1,for quantiles) For graphical purposes only
    q1 <- quantile(x, probs = x1, type = 1)
    q_inter <- approxfun(x1, q1, yleft = min(q1), yright = max(q1), method = "constant")
    plot(q_inter, ylim = c(min(y[1], x[1]) - 1, max(x[lenx], y[leny]) + 1), main = "Interpolated Quartile Function", 
         xlab = "Probs", ylab = "Data")  #quartile plot for z2
    points(y1, y, col = "red")
    points(x1, q1)
    # Making ECDF Graph
    f1 <- ecdf(x)
    plot(f1, xlim = c(min(x[1], y[1]), max(x[lenx], y[lenx])), ylab = "Probs", xlab = "Data", 
         main = "ECDF and Points")  #plot 1st sample points
    points(y, y1, col = "red")  #put 2nd sample on graph
  } else if (lenx > leny) {
    q1 <- quantile(x, probs = x1, type = 4)
    q_inter <- approxfun(x1, q1, yleft = min(q1), yright = max(q1))
    z <- max(abs(q_inter(y1) - y))
    f <- approxfun(x, x1, yleft = 0, yright = 1, ties = max)  #ties=max makes sure cdf jumps to .4
    par(mfrow = c(2, 1))
    plot(q_inter, ylim = c(min(y[1], x[1]) - 1, max(x[lenx], y[leny]) + 1), main = "Interpolated Quartile Function", 
         xlab = "Probs", ylab = "Data")  #quartile plot for z2
    points(y1, y, col = "red")
    plot(f, ylim = c(0, 1), xlim = c(min(y[1], x[1]) - 1, max(x[lenx], y[leny]) + 
                                       1), main = "Linearly Interpolated ECDF", xlab = "Data", ylab = "Probs")  #plot of ECDF 
    # Want the limits to include the points, so above w/ maxes and mins is necessary
    lines(ecdf(x))  #draws in original ECDF
    points(y, y1, col = "red")
  } else {
    # So length y>x
    q2 <- quantile(y, probs = y1, type = 4)
    q_inter <- approxfun(y1, q2, yleft = min(q2), yright = max(q2))
    z <- max(abs(q_inter(x1) - x))
    f <- approxfun(y, y1, yleft = 0, yright = 1, ties = max)  #ties=max makes sure cdf jumps to .4
    par(mfrow = c(2, 1))
    plot(q_inter, ylim = c(min(y[1], x[1]) - 1, max(x[lenx], y[leny]) + 1), main = "Interpolated Quartile Function", 
         xlab = "Probs", ylab = "Data")  #quartile plot for z2
    points(x1, x, col = "red")
    plot(f, ylim = c(0, 1), xlim = c(min(y[1], x[1]) - 1, max(x[lenx], y[leny]) + 
                                       1), main = "Linearly Interpolated ECDF", xlab = "Data", ylab = "Probs")  #plot of ECDF 
    # Want the limits to include the points, so above w/ maxes and mins is necessary
    lines(ecdf(y))  #draws in original ECDF
    points(x, x1, col = "red")
  }
}
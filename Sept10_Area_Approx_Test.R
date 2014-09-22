# Testing the Area Approximation Function
# September 9th, 2014
rm(list=ls())
source("functions.R")
# Goal: Make Plot Function

set.seed(5)
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
lines(y_points,y_out,col="red")

coords <- Quad_Quan_Area_TS(x,y,opt="coords")
coords[[1]]
# points(-1.67,0.133)
# points(-0.84,0.40)
# points(-1.017,.097)
# points(-.188,.362)
plot_polygons <- function(liCoords){
plotx <- sapply(liCoords,function(x) x[1])  
ploty <- sapply(liCoords,function(x) x[2])  
my_plot <- polygon(plotx,ploty, col="blue")
my_plot
}

sapply(coords, plot_polygons)

source("functions.R")
set.seed(5)
x <- rnorm(5)
y <- rnorm(5)
Quad_Quan_Area_TS(x,y,do.plot=TRUE)


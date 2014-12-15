wave.den <- function (x, y, ..., doplot=F, wf="haar")  #n=2^5
{
  # Inputs:
  # Y must either be numeric or "p"dist
  F.x <- ecdf(x)
  
  # Two Sample 
  if(is.numeric(y)){
  F.y <- ecdf(y)
  ml <- min(length(x),length(y)) 
  n <- 2^floor(log2(ml))
  z <- seq(range(x, y)[1], range(x, y)[2], length.out=n)
  F.dwt <- dwt(F.x(z) - F.y(z), wf=wf, n.levels=log(n, 2))
  } else{
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
    
    
    n <- 2^floor(log2(length(x)))
    z <- seq(min(x),max(x),length=n)
    F.dwt <- dwt(F.x(z) - y(z,...), wf=wf, n.levels=log(n,2))
  }
  return(F.dwt)
  
  oc <- unlist(F.dwt)
  oc <- max(abs(oc))
  #test_ks <- max(abs(F.x(z)-F.y(z)))

  if(doplot)
  {
    plot(z, F.x(z), type="l", 
         ylim=range(1.1, F.x(z), F.y(z), F.x(z) - F.y(z)))
    lines(z, F.y(z), col=3)
    lines(z, F.x(z) - F.y(z), col=4)
    abline(h=0, lty=3)
    abline(h=1, lty=3)
    segments(z[1], 0, z[1], oc, lwd=2)
    #segments(z[3], 0, z[3], ks, lwd=2)
    title(paste("oc = ", round(oc, 2),sep=""))
  }
  oc
}


wave.energy <- function (x, y, ...,
                         #n=2^5, 
                         doplot=F, 
                         opt="max", 
                         wf="haar",
                         square=FALSE,
                         norm=TRUE)
{
  # Args:
  # y must be numeric of a "p" distribution function
  if(opt !="max" & opt != "sum") stop ("invalid option value")
  # Get cdfs:
  F.x <- ecdf(x)
  # Two Sample
  if(is.numeric(y)){
  F.y <- ecdf(y)
  
  ml <- min(length(x),length(y))
  n <- 2^floor(log2(ml))
  z <- seq(range(x, y)[1], range(x, y)[2], length=n)
  
  F.x.dwt <- dwt(F.x(z), wf=wf, n.levels=log(n, 2))
  F.y.dwt <- dwt(F.y(z), wf=wf, n.levels=log(n, 2))
  x.dwt <- unlist(F.x.dwt)
  y.dwt <- unlist(F.y.dwt)
  }else{
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
    
    n <- 2^floor(log2(length(x)))
    z <- seq(min(x), max(x), length=n)
    F.x.dwt <- dwt(F.x(z), wf=wf, n.levels=log(n, 2))
    F.y.dwt <- dwt(y(z,...), wf=wf, n.levels=log(n, 2))
    x.dwt <- unlist(F.x.dwt)
    y.dwt <- unlist(F.y.dwt)
  }
  
  if(square==FALSE){
  x.dwt <- sort(abs(x.dwt), decreasing = TRUE)
  y.dwt <- sort(abs(y.dwt), decreasing = TRUE)
  }
  # What if we square the coefficients?
  if(square){
  x.dwt <- sort(x.dwt^2, decreasing=TRUE)
  y.dwt <- sort(y.dwt^2, decreasing=TRUE)
  }

  xx <- cumsum(x.dwt)
  yy <- cumsum(y.dwt)
  max_x <- tail(xx, 1)
  max_y <- tail(yy, 1)
  # return(max(yy))
  if(norm==FALSE){
    if(max_x > max_y){
      xx <- xx / max_x
      yy <- yy / max_x
    } else{
      xx <- xx / max_y
      yy <- yy / max_y
    }
  }
  if(norm==TRUE){  
  xx <- xx / max_x
  yy <- yy / max_y
  }
  if(doplot){
    plot(xx, type="l", ylim=range(xx, yy, xx - yy))
    lines(yy, col=2)
    lines(xx - yy, col=3)
    abline(h=0, lty=3)
  }
  
  # Returning test satistic values
  if(opt=="max"){
    ts <- max(abs(xx-yy))
    names(ts) <- "max"
    return(ts)
  }
  if(opt=="sum"){
    ts <- sum(abs(xx-yy))
    names(ts) <- "sum"
    return(ts)
  }
}


wave.den.perm <- function (x=rnorm(10), y=rnorm(10), n=2^5, p=100, doplot=F, doplot1=F) 
{
  z <- c(x, y)
  n.x <- length(x)
  n.z <- length(z)
  ts.vec <- matrix(,nrow=p, ncol=2)
  
  for(i in 1:p)
  {
    perm <- sample(z, n.z, replace=F)
    #a <- wave.den(perm[1:n.x], perm[(n.x + 1):n.z])
    #b <- ks.res.simp(perm[1:n.x], perm[(n.x + 1):n.z])
    #n=n, 
    #doplot=doplot) 
    a <- wave.energy(perm[1:n.x], perm[(n.x + 1):n.z], n=n, 
                     doplot=doplot) 
    ts.vec[i, ] <- a
  }
  
  colnames(ts.vec) <- c("oc", "ks")
  rownames(ts.vec) <- rep("", p)
  
  d.oc <- density(ts.vec[, 1])
  d.ks <- density(ts.vec[, 2])
  
  pp <- ceiling(0.95 * p)
  oc.crit <- sort(ts.vec[, 1])[pp]
  ks.crit <- sort(ts.vec[, 2])[pp]
  # Plots
  a <- wave.den(x, y, n=n, doplot=doplot)
  
  if(doplot1)
  {
    plot(d.oc, xlim=range(d.oc$x, d.ks$x), 
         ylim=range(d.oc$y, d.ks$y),main="")
    lines(d.ks, col=3)
    abline(v=oc.crit, lty=3)
    abline(v=ks.crit, col=3, lty=3)
    abline(v=a[1])
    abline(v=a[2], col=3)
  }
  
  oc.dec <- ifelse(a[1] < oc.crit, 0, 1)
  ks.dec <- ifelse(a[2] < ks.crit, 0, 1)
  c(oc.dec, ks.dec)
}

wave.den.power <- function (n=2^5, p=500, doplot=F, m=100, doplot1=F) 
{
  b <- numeric(0)
  power.vec <- matrix(,nrow=m, ncol=2)
  for(i in 1:m)
  {
    x <- rnorm(100)
    #		y <- rt(100, df=20)
    y <- rnorm(100, sd=1)
    a <- wave.den.perm(x=x, y=y, p=p)#, n=n, p=p, doplot=doplot, doplot1=doplot1)
    power.vec[i, ] <- a 
  }
  c(oc=mean(power.vec[, 1]), ks=mean(power.vec[, 2]), frac=mean(power.vec[, 1]) / mean(power.vec[, 2]))
}
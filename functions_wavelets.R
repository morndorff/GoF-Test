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
                         opt="sum", 
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

wave.bec <- function(x,y, ..., interp = 4, doplot=F, wf="haar", reduce=2)
{
  library(wavethresh)
  library(waveslim)
  x <- sort(x)
  lenx <- length(x)
  # Two Sample Test
  if (is.numeric(y)) {
    leny <- length(y)
    y <- sort(y)

    num_quan <- min(2^floor(log(lenx/reduce,2)), 2^floor(log(leny/reduce,2)))
    prob <- seq((1-.5)/num_quan, (num_quan-.5)/num_quan, length.out=num_quan)

    # Because # of quantiles is < data points, need to interpolate
    qx <- quantile(x, probs = prob, type = interp)
    qy <- quantile(y, probs = prob, type = interp)
    q_inter_x <- approxfun(prob, qx, yleft = min(qx), yright = max(qx))
    q_inter_y <- approxfun(prob, qy, yleft = min(qy), yright = max(qy))
    quan_dif <- q_inter_y(prob) - q_inter_x(prob)
    
    test_wave <- dwt(quan_dif,wf="haar")
    w_coef <- unlist(test_wave)
    # Take abs. value of wavelet coefficients, sort in decreasing order
    w_coef <- sort(abs(w_coef), decreasing=TRUE) 
    # Sum of them
    sum_coef <- sum(w_coef)
    # Get number of coefficients to keep, based on 90% thresholding
    num_coef <- which(cumsum(abs(w_coef))/sum_coef<=.9)
    th_w_coef<- w_coef[num_coef]
    STAT <- sum(th_w_coef^2)
    return(STAT)
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
  # Assuring diadic lengths
  #nx_p2 <- 2^floor(log(lenx,2))
  
  num_quan <- 2^floor(log(lenx,2))
  prob <- seq((1-.5)/num_quan, (num_quan-.5)/num_quan, length.out=num_quan)
  
  
  # Estimated Quantiles:
  qx <- quantile(x, probs = prob, type = interp)
  q_inter_x <- approxfun(prob, qx, yleft = min(qx), yright = max(qx))
  x_quans <- q_inter_x(prob)
  # True Quantiles
  true_quan <- y(prob, ...)
  quan_dif <- x_quans - true_quan
  
  test_wave <- dwt(quan_dif,wf="haar")
  w_coef <- unlist(test_wave)
  # Take abs. value of wavelet coefficients, sort in decreasing order
  w_coef <- sort(abs(w_coef), decreasing=TRUE) 
  # Sum of them
  sum_coef <- sum(w_coef)
  # Get number of coefficients to keep, based on 90% thresholding
  num_coef <- which(cumsum(abs(w_coef))/sum_coef<=.9)
  th_w_coef<- w_coef[num_coef]
  STAT <- sum(th_w_coef^2)
  return(STAT)
}

# 
# wave.bec2 <- function(x,y, ..., interp = 4, doplot=F, wf="haar", reduce=2)
# {
#   library(wavethresh)
#   library(waveslim)
#   x <- sort(x)
#   lenx <- length(x)
#   # Two Sample Test
#   if (is.numeric(y)) {
#     leny <- length(y)
#     y <- sort(y)
#     
#     num_quan <- min(2^floor(log(lenx/reduce,2)), 2^floor(log(leny/reduce,2)))
#     prob <- seq((1-.5)/num_quan, (num_quan-.5)/num_quan, length.out=num_quan)
#     
#     # Because # of quantiles is < data points, need to interpolate
#     qx <- quantile(x, probs = prob, type = interp)
#     qy <- quantile(y, probs = prob, type = interp)
#     q_inter_x <- approxfun(prob, qx, yleft = min(qx), yright = max(qx))
#     q_inter_y <- approxfun(prob, qy, yleft = min(qy), yright = max(qy))
#     quan_dif <- q_inter_y(prob) - q_inter_x(prob)
#     print(quan_dif)
#     n <- length(quan_dif)
#     wave_tran <- wd(quan_dif, family="DaubExPhase", filter.number=1) #haar wavelet
#     scoef <- sort(abs(wave_tran$D), decreasing=TRUE)
#     sum_coef <- sum(abs(scoef))
#     num_coef <- min(which(cumsum(abs(scoef))/sum_coef>.9))
#     wave_th <- threshold(wave_tran, mannum= num_coef)
#     STAT <- sum(wave_th$D^2)
#     return(STAT)
#   }
#   
#   # One Sample
#   if (is.list(y)) 
#     y <- names(y)
#   if (is.function(y)) 
#     funname <- as.character(substitute(y))
#   if (is.character(y)) 
#     funname <- y
#   y <- get(funname, mode = "function", envir = parent.frame())
#   if (!is.function(y)) 
#     stop("'y' must be numeric or a function or a string naming a valid function")
#   # Assuring diadic lengths
#   #nx_p2 <- 2^floor(log(lenx,2))
#   
#   num_quan <- 2^floor(log(lenx,2))
#   prob <- seq((1-.5)/num_quan, (num_quan-.5)/num_quan, length.out=num_quan)
#   
#   
#   # Estimated Quantiles:
#   qx <- quantile(x, probs = prob, type = interp)
#   q_inter_x <- approxfun(prob, qx, yleft = min(qx), yright = max(qx))
#   x_quans <- q_inter_x(prob)
#   # True Quantiles
#   true_quan <- y(prob, ...)
#   
#   quan_dif <- x_quans - true_quan
#   print(quan_dif)
#   
#   n <- length(quan_dif)
#   wave_tran <- wd(quan_dif, family="DaubExPhase", filter.number=1) #haar wavelet
#   test_wave <- dwt(quan_dif,wf="haar")
#   return(list(test_wave,wave_tran))
#   
#   print(quan_dif)
#   scal_coef <- tail(wave_tran$C,1) # Scaling Coefficient
#   print(scal_coef)
#   # Take abs. value of wavelet coefficients, sort in decreasing order
#   w_coef <- sort(abs(c(scal_coef,wave_tran$D)), decreasing=TRUE) 
#   # Sum of them
#   sum_coef <- sum(w_coef)
#   # Get number of coefficients to keep, based on 90% thresholding
#   num_coef <- min(which(cumsum(abs(w_coef))/sum_coef>.9))
#   # Threshold according
#   wave_th <- threshold(wave_tran, policy="mannum",value= num_coef)
#   return(wave_th)
#   # Calculate square of remaining wavelet coefficients
#   th_scal_coef <- tail(wave_th$C,1)^2
#   th_w_coef <- c(th_scal_coef, wave_th$D)
#   print(th_w_coef)
#   STAT <- sum(th_w_coef^2)
#   return(STAT)
# }

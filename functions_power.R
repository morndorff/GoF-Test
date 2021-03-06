# Functions for Power Simulation Studies (ks, cramer von mises, etc.)

# KS-Test
ks.res <- function(x, y, ..., alternative = c("two.sided", "less", "greater"), exact = NULL) {
    alternative <- match.arg(alternative)
    DNAME <- deparse(substitute(x))
    x <- x[!is.na(x)]  #removing missing entries
    n <- length(x)
    # need more than one
    if (n < 1L) 
        stop("not enough 'x' data")
    
    PVAL <- NULL
    
    if (is.numeric(y)) {
        # Good case, y is a numeric vector
        DNAME <- paste(DNAME, "and", deparse(substitute(y)))
        y <- y[!is.na(y)]
        n.x <- as.double(n)
        n.y <- length(y)
        if (n.y < 1L) 
            stop("not enough 'y' data")
        # Checking to see how big the vectors are
        if (is.null(exact)) 
            exact <- (n.x * n.y < 10000)
        METHOD <- "Two-sample Kolmogorov-Smirnov test"
        TIES <- FALSE
        n <- n.x * n.y/(n.x + n.y)
        w <- c(x, y)
        z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))  #This order function is nice!
        if (length(unique(w)) < (n.x + n.y)) {
            if (exact) {
                warning("cannot compute exact p-value with ties")
                exact <- FALSE
            } else warning("p-value will be approximate in the presence of ties")
            z <- z[c(which(diff(sort(w)) != 0), n.x + n.y)]
            TIES <- TRUE
        }
        STATISTIC <- switch(alternative, two.sided = max(abs(z)), greater = max(z), 
            less = -min(z))
        nm_alternative <- switch(alternative, two.sided = "two-sided", less = "the CDF of x lies below that of y", 
            greater = "the CDF of x lies above that of y")
        # if (exact && (alternative == 'two.sided') && !TIES) #NOTE: DID NOT CALL P-VALUE
        # PVAL <- 1 - .Call(C_pSmirnov2x, STATISTIC, n.x, n.y)
    } else {
        if (is.list(y)) 
            y <- names(y)
        if (is.character(y)) 
            y <- get(y, mode = "function", envir = parent.frame())
        
        if (!is.function(y)) 
            stop("'y' must be numeric or a function or a string naming a valid function")
        METHOD <- "One-sample Kolmogorov-Smirnov test"
        TIES <- FALSE
        if (length(unique(x)) < n) {
            # Another ties warning
            warning("ties should not be present for the Kolmogorov-Smirnov test")
            TIES <- TRUE
        }
        
        nm_alternative <- switch(alternative, two.sided = "two-sided", less = "the CDF of x lies below the null hypothesis", 
            greater = "the CDF of x lies above the null hypothesis")
    }
    names(STATISTIC) <- switch(alternative, two.sided = "D", greater = "D^+", less = "D^-")
    
    
    RVAL <- list(statistic = STATISTIC, alternative = nm_alternative, method = METHOD, 
        data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
}

ad.res <- function(x,y, ...){
  if(is.numeric(y)) stop("Haven't Done Two Sample Yet")
  
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
  x <- sort(x)
  n <- length(x)
  if (n < 8) 
    stop("sample size must be greater than 7")
  logp1 <- pnorm((x - mean(x))/sd(x), log.p = TRUE)
  logp2 <- pnorm(-(x - mean(x))/sd(x), log.p = TRUE)
  h <- (2 * seq(1:n) - 1) * (logp1 + rev(logp2))
  STAT <- -n - mean(h)
  return(STAT)
}

# K test (Just outputting statistic)
ks.res.simp <- function(x, y, ..., alternative = c("two.sided", "less", "greater"), 
    exact = NULL) {
    alternative <- match.arg(alternative)
    DNAME <- deparse(substitute(x))
    x <- x[!is.na(x)]
    n <- length(x)
    if (n < 1L) 
        stop("not enough 'x' data")
    PVAL <- NULL
    if (is.numeric(y)) {
        DNAME <- paste(DNAME, "and", deparse(substitute(y)))
        y <- y[!is.na(y)]
        n.x <- as.double(n)
        n.y <- length(y)
        if (n.y < 1L) 
            stop("not enough 'y' data")
        if (is.null(exact)) 
            exact <- (n.x * n.y < 10000)
        METHOD <- "Two-sample Kolmogorov-Smirnov test"
        TIES <- FALSE
        n <- n.x * n.y/(n.x + n.y)
        w <- c(x, y)
        z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))
        if (length(unique(w)) < (n.x + n.y)) {
            if (exact) {
                warning("cannot compute exact p-value with ties")
                exact <- FALSE
            } else warning("p-value will be approximate in the presence of ties")
            z <- z[c(which(diff(sort(w)) != 0), n.x + n.y)]
            TIES <- TRUE
        }
        STATISTIC <- switch(alternative, two.sided = max(abs(z)), greater = max(z), 
            less = -min(z))
        nm_alternative <- switch(alternative, two.sided = "two-sided", less = "the CDF of x lies below that of y", 
            greater = "the CDF of x lies above that of y")
    } else {
        if (is.list(y)) 
            y <- names(y)
        if (is.character(y)) {
            y <- dist.conv.str(y, type = "p")
            y <- get(y, mode = "function", envir = parent.frame())
        }
        if (!is.function(y)) 
            stop("'y' must be numeric or a function or a string naming a valid function")
        METHOD <- "One-sample Kolmogorov-Smirnov test"
        TIES <- FALSE
        if (length(unique(x)) < n) {
            warning("ties should not be present for the Kolmogorov-Smirnov test")
            TIES <- TRUE
        }
        if (is.null(exact)) 
            exact <- (n < 100) && !TIES
        x <- y(sort(x), ...) - (0:(n - 1))/n
        STATISTIC <- switch(alternative, two.sided = max(c(x, 1/n - x)), greater = max(1/n - 
            x), less = max(x))
        nm_alternative <- switch(alternative, two.sided = "two-sided", less = "the CDF of x lies below the null hypothesis", 
            greater = "the CDF of x lies above the null hypothesis")
    }
    names(STATISTIC) <- switch(alternative, two.sided = "D", greater = "D^+", less = "D^-")
    
    RVAL <- list(statistic = STATISTIC, alternative = nm_alternative, method = METHOD, 
        data.name = DNAME)
    class(RVAL) <- "htest"
    return(STATISTIC)
}

CvMTwoSamp.res <- function(S1, S2) {
    sortx = sort(S1)
    lenx = length(sortx)
    sorty = sort(S2)
    leny = length(sorty)
    a = data.frame(val = sortx, rang = seq(lenx), ens = rep(1, lenx))
    b = data.frame(val = sorty, rang = seq(leny), ens = rep(2, leny))
    d = rbind(a, b)
    d = d[order(d$val), ]
    d = data.frame(d, rangTot = seq(lenx + leny))
    dtfM = d[which(d$ens == 1), ]
    dtfN = d[which(d$ens == 2), ]
    somN = sum((dtfN$rang - dtfN$rangTot)^2)
    somM = sum((dtfM$rang - dtfM$rangTot)^2)
    U = leny * somN + lenx * somM
    CvM = ((U/(leny * lenx))/(leny + lenx)) - ((4 * lenx * leny - 1)/(6 * (lenx + leny)))
    return(CvM)
} 

perm.test2 <- function(x, y, distops = NULL, f, fops = NULL, num.perm = 2001, diag = FALSE, 
                       exact = FALSE, out=FALSE, do.plot=FALSE, ...) {
  #Args: 
  # x: numeric vector
  # y: numeric vector or quantile distribution (qgamma, qnorm)
  # distops: list describing quantile parameters (e.g. for qunif, list(min=0,max=2)
  # f: function outputting a test statistic
  # fops: if the test statistic has options, put them here. (e.g. for myts.out, size=.2)
  # num.perm: number of permutations to assess p-values
  # 
  # Output:
  # list containing observed test statistic, 
  #if (is.null(distops)==FALSE){
  #  if(is.list(distops)==FALSE) stop("distops must be a list")
  #}
  library(doMC)
  library(foreach)
  registerDoMC(2)
  if(! (is.null(distops) || is.list(distops))){stop("distops should be NULL or a list")}
  if(! (is.null(fops) || is.list(fops))){stop("fops should be NULL or a list")}
  
  if (out==TRUE){
    res_out <- perm.test.out(x,y,distops,f,fops,num.perm,diag,exact)
    return(res_out)
  }
  
  
  # Handling function inputs for y
  if (is.function(y)) y <- as.character(substitute(y))
  if (is.character(y)) y <- chartoli(y)
  
  # Calculating observed test statistic
  # One Sample
  if (is.list(y)) {
    if (length(fops) == 0) fops <- NULL
    if (length(distops) == 0) distops <- NULL
    ts.obs <- do.call(f, c(list(x), list(names(y)), distops, fops))
  }
  # Two sample
  if (is.numeric(y)){
    if (length(fops)== 0) fops <- NULL
    ts.obs <- do.call(f, c(list(x,y), fops))
  }
  
  lenx <- length(x)
  ts.random <- vector(mode = "numeric", length = num.perm)
  # Two sample
  if (is.numeric(y)) {
    z <- c(x, y)
    lenz <- length(z)
    if (lenz < 11 & exact == TRUE) {
      require(gtools)
      
      all.perm <- permutations(n = lenz, r = lenz, v = z, repeats.allowed = FALSE, 
                               set = FALSE)
      all.permx <- all.perm[, 1:lenx]
      all.permy <- all.perm[, (lenx + 1):lenz]
      exact.perm <- dim(all.perm)[1]
      for (i in 1:exact.perm) {
        ts.random[i] <- f(all.permx[i, ], all.permy[i, ])
      }
      p.val <- sum(abs(ts.random) >= abs(ts.obs))/exact.perm
      c.val <- quantile(ts.random, probs = 0.95)
    } else {
#       for (i in 1:num.perm) {
#         z1 <- sample(z, size = lenz, replace = FALSE)
#         a <- z1[1:lenx]
#         b <- z1[(lenx + 1):lenz]
#         ts.random[i] <- do.call(f, c(list(a,b), fops))
#       }
#       foreach(i=1:num.perm) %do% {
#         z1 <- sample(z, size = lenz, replace = FALSE)
#         a <- z1[1:lenx]
#         b <- z1[(lenx + 1):lenz]
#         ts.random[i] <- do.call(f, c(list(a,b), fops))
#       }
      z_mat <- replicate(num.perm,sample(z,size=lenz, replace=FALSE))
      ts.random <- apply(z_mat,2, function(x) do.call(f, c( list(x[1:lenx], x[(lenx +1):lenz], fops ))))
      # Slower
      
      
      p.val <- sum(abs(ts.random) >= abs(ts.obs)) / num.perm
      c.val <- quantile(ts.random, probs = 0.95)
    }
    if (do.plot == TRUE){
      hplot <- hist(ts.random, prob=TRUE)
      hplot
      segments(ts.obs, 0, x1=ts.obs, y1=max(hplot$density))
    }
    if (diag == TRUE) {
      return(list("p-value" = p.val, "95% crit val" = c.val, "Obs. TS" = ts.obs, 
                  "ts.dist" = ts.random))
    } else {
      return(list("p-value" = p.val, "95% crit val" = c.val, "Obs. TS" = ts.obs))
    }
  }
  # One Sample
  if (is.list(y)) {
    ry <- dist.conv(funname = names(y), type = "r")
    fy <- get(names(y), mode = "function", envir = parent.frame())
    for (i in 1:num.perm) {
      z <- do.call(ry, c(list(lenx), distops))  #(lenx,...)
      ts.random[i] <- do.call(f, c(list(z), list(names(y)), distops, fops))
    }
    p.val <- sum(abs(ts.random) >= abs(ts.obs)) / num.perm
    c.val <- quantile(ts.random, probs = 0.95)
    if (do.plot == TRUE){
      hplot <- hist(ts.random, prob=TRUE, main="Histogram of Permuted Test Statistic (line=Observed TS)")
      hplot
      segments(ts.obs,0,x1=ts.obs,y1=max(hplot$density))
    }
    if (diag == TRUE) {
      # 1st value of output is p value, 2nd is 95% critical value, 3rd is the actual test
      # statistic
      return(list("p-value" = p.val, "95% crit val" = c.val, "Obs. TS" = ts.obs, 
                  "ts.dist" = ts.random))
    } else {
      return(list("p-value" = p.val, "95% crit val" = c.val, "Obs. TS" = ts.obs))
    }
  }
}

cvm.res <- function(x, y, ..., test=FALSE){
  # x: numeric vector
  # y: numeric vector or prob distribution e.g. "pnorm" 
  if (is.numeric(x) != TRUE) 
    stop("x must be numeric")
  x <- sort(x)
  lenx <- length(x)
  # Two Sample Test
  if (is.numeric(y)) {
    sortx <- x
    lenx = length(sortx)
    sorty = sort(y)
    leny = length(sorty)
    a = data.frame(val = sortx, rang = seq(lenx), ens = rep(1, lenx))
    b = data.frame(val = sorty, rang = seq(leny), ens = rep(2, leny))
    d = rbind(a, b)
    d = d[order(d$val), ]
    d = data.frame(d, rangTot = seq(lenx + leny))
    dtfM = d[which(d$ens == 1), ]
    dtfN = d[which(d$ens == 2), ]
    somN = sum((dtfN$rang - dtfN$rangTot)^2)
    somM = sum((dtfM$rang - dtfM$rangTot)^2)
    U = leny * somN + lenx * somM
    CvM = ((U/(leny * lenx))/(leny + lenx)) - ((4 * lenx * leny - 1)/(6 * (lenx + leny)))
    return(CvM)
  }
  # One Sample Test  
  if (is.list(y)) 
    y <- names(y)
  if (is.function(y)) 
    funname <- as.character(substitute(y))
  if (is.character(y)) 
    funname <- y
  y <- get(funname, mode = "function", envir = parent.frame())
  if (!is.function(y)) 
    stop("'y' must be numeric or a function or a string naming a valid function")
  F_x <- y(x, ...)
  i <- 1:lenx
  STAT <- (1/(12*lenx)) + sum( (F_x -(1:lenx - .5)/ lenx)^2)
  if(test){
    # print(F_x)
    # print(F_x - (1:lenx -.5)/lenx)
    # return(F_x - (1:lenx -.5)/lenx)
    STAT <- (1/(12*lenx)) + sum((F_x - (1:lenx - .5)/lenx)^2)
    return(STAT)
  }
  STAT
}
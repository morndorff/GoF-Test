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
    STATISTIC <- switch(alternative, two.sided = max(abs(z)), greater = max(z), less = -min(z))
    nm_alternative <- switch(alternative, two.sided = "two-sided", less = "the CDF of x lies below that of y", 
                             greater = "the CDF of x lies above that of y")
    # if (exact && (alternative == 'two.sided') && !TIES) #NOTE: DID NOT CALL P-VALUE PVAL <- 1 -
    # .Call(C_pSmirnov2x, STATISTIC, n.x, n.y)
  } else {
    if (is.list(y)) y <- names(y)
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
  
  
  RVAL <- list(statistic = STATISTIC, alternative = nm_alternative, method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}

# K test (Just outputting statistic)
ks.res.simp <- function (x, y, ..., 
                         alternative = c("two.sided", "less", "greater"), 
                         exact = NULL) 
{
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
      }
      else warning("p-value will be approximate in the presence of ties")
      z <- z[c(which(diff(sort(w)) != 0), n.x + n.y)]
      TIES <- TRUE
    }
    STATISTIC <- switch(alternative, two.sided = max(abs(z)), 
                        greater = max(z), less = -min(z))
    nm_alternative <- switch(alternative, two.sided = "two-sided", 
                             less = "the CDF of x lies below that of y", greater = "the CDF of x lies above that of y")
  }
  else {
    if (is.list(y)) y <- names(y)
    if (is.character(y)) {
      y <- dist.conv.str(y,type="p")
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
    STATISTIC <- switch(alternative, 
                        two.sided = max(c(x, 1/n - x)), 
                        greater = max(1/n - x), 
                        less = max(x))
    nm_alternative <- switch(alternative, two.sided = "two-sided", 
                             less = "the CDF of x lies below the null hypothesis", 
                             greater = "the CDF of x lies above the null hypothesis")
  }
  names(STATISTIC) <- switch(alternative, two.sided = "D", 
                             greater = "D^+", less = "D^-")

  RVAL <- list(statistic = STATISTIC, alternative = nm_alternative, 
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  return(STATISTIC)
}

CvMTwoSamp.res <- function (S1, S2) {
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
  CvM = ((U/(leny * lenx))/(leny + lenx)) - ((4 * lenx * leny - 1)/(6 * (lenx + 
                                                                           leny)))
  return(CvM)
}
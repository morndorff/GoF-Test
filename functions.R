# Useful Functions
source("functions_out.R")
source("plot_functions.R")
source("functions_power.R")
source("functions_string.R")

myts <- function(x, y, ..., interp = 4, do.plot = FALSE) {
    # Computes maximum difference of quantiles 
    # Comments:
    # OBC Using Linear Interpolation of the ECDF Recall: ECDF range is [1/n,1]. This is probably not
    # realistic
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
    leny <- length(y)
    # Two Sample Test
    if (is.numeric(y)) {
        y <- sort(y)
        # x1 and y1 are the quantile values
        x1 <- seq(1/lenx, 1, 1/lenx)
        y1 <- seq(1/leny, 1, 1/leny)
        if (lenx == leny) {
            # In the case of equal sample sizes take max of abs. val
            z <- max(abs(y - x))
            if (do.plot == TRUE) 
                plot.ts.2sam(x, y, x1, y1, lenx, leny)
        } else if (lenx > leny) {
            # If there are unequal sample sizes, then interpolation is necessary
            q1 <- quantile(x, probs = x1, type = interp)
            q_inter <- approxfun(x1, q1, yleft = min(q1), yright = max(q1))
            z <- max(abs(q_inter(y1) - y))
            if (do.plot == TRUE) 
                plot.ts.2sam(x, y, x1, y1, lenx, leny)
        } else {
            # So length y>x
            q2 <- quantile(y, probs = y1, type = interp)
            q_inter <- approxfun(y1, q2, yleft = min(q2), yright = max(q2))
            z <- max(abs(q_inter(x1) - x))
            if (do.plot == TRUE) 
                plot.ts.2sam(x, y, x1, y1, lenx, leny)
        }
        return(z)
    }
    
    # One Sample
    if (is.function(y)) 
        funname <- as.character(substitute(y))
    if (is.character(y)) 
        funname <- y
    y <- get(funname, mode = "function", envir = parent.frame())
    if (!is.function(y)) 
        stop("'y' must be numeric or a function or a string naming a valid function")
    z <- max(x - y(seq(1/(lenx + 1), lenx/(lenx + 1), length.out = lenx), ...))  #Note: quantiles up for debate
    if (do.plot == TRUE) {
        plot.ts.1sam(x, y, ..., funname = funname, lenx = lenx)
    }
    return(z)
}
power.res <- function(x, y, f, g = perm.test, boot = FALSE) {
  # A nice function for power calculations Takes as input a two matrices. Each row a different draw of
  # col's from the row dist.
    dim.x <- dim(x)
    fun <- f
    pv <- NULL
    if (boot) {
        for (i in 1:dim.x[1]) {
            a <- boot.test(x[i, ], y[i, ], f)
            pv[i] <- a[[1]]
        }
        z <- pv
        return(z)
    }
    for (i in 1:dim.x[1]) {
        a <- g(x[i, ], y[i, ], f = fun)
        pv[i] <- a[[1]]
    }
    z <- pv
    return(z)
}
power.res.exact <- function(x, dist, ..., f, boot = FALSE) {
  # Does test of distribution based on EXACT quantiles of hypothesized distribution
  
    dim.x <- dim(x)
    pv <- NULL
    if (boot) {
        for (i in 1:dim.x[1]) {
            a <- boot.test(x[i, ], dist, ..., f)
            pv[i] <- a[[1]]
        }
        z <- pv
        return(z)
    }
    for (i in 1:dim.x[1]) {
        a <- perm.test(x[i, ], dist, ..., f = f)
        pv[i] <- a[[1]]
    }
    z <- pv
    return(z)
}
myts.par <- function(x, y, ..., interp = 4, do.plot = FALSE, size = .25) {
  # Computes area based on trapezoid areas
  # Args: 
  # x: A vector of observations for a R.V. (must be numeric)
  #
  # y: Either 
  # (1) Another vector of observations (two sample) 
  # (2) A quantile function such as qnorm (one sample)
  #
  # size: controls height of trapezoids
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
  # Quantiles
  x1 <- seq(1/(lenx+1),lenx/(lenx+1),length.out=lenx)
  # Parameter which determines height of trapezoids
  delta <- size * (1 / lenx)
  # Creating quantile function
  # q1 <- quantile(x, probs = x1, type = interp)
  # q_inter <- approxfun(x1, q1, yleft = min(q1), yright = max(q1))
  # Length of horizontal segments
  x_p_d <- x + ( y(x1+delta, ...) - y(x1, ...))
  x_m_d <- x + ( y(x1-delta, ...) - y(x1, ...))
  y_p_d <- y(x1 + delta, ...)
  y_m_d <- y(x1 + delta, ...)
  # sum of trapezoid areas
  z <- sum((abs(x_p_d - y_p_d) + abs(x_m_d - y_m_d)) * delta)
  if(do.plot==TRUE){
    plot.ts.1sam(x, y, ..., funname = funname, lenx = lenx)
    
  }
  return(z)
}
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
myts.max <- function(x, y) {
    obc.stat <- myts(x, y)
    ks.stat <- ks.res(x, y)$statistic
    if (obc.stat > ks.stat) {
        z <- list(obc.stat, c("OBC"))
    } else if (ks.stat > obc.stat) {
        z <- list(ks.stat, c("KS"))
    } else z <- list(c(NULL), c("NEITHER"))
    
    return(z)
}
myts.max.simp <- function(x, y) {
    obc.stat <- myts(x, y)
    ks.stat <- ks.res(x, y)$statistic
    if (obc.stat > ks.stat) {
        z <- obc.stat
    } else if (ks.stat > obc.stat) {
        z <- ks.stat
    } else z <- 1000  #######FIX THIS LATER
    
    return(z)
}
myts.max.range <- function(x, y) {
    # Returns value of the statistic and whether or not it from OBC
    com <- c(x, y)
    range <- max(com) - min(com)
    obc.stat <- myts(x, y)/range
    ks.stat <- ks.res(x, y)$statistic
    if (obc.stat > ks.stat) {
        z <- list(obc.stat, c("OBC"))
    } else if (ks.stat > obc.stat) {
        z <- list(ks.stat, c("KS"))
    } else z <- list(c(NULL), c("NEITHER"))
    
    return(z)
}
myts.max.range.simp <- function(x, y) {
    # Returns value of the statistic and whether or not it from OBC
    com <- c(x, y)
    range <- max(com) - min(com)
    obc.stat <- myts(x, y)/range
    ks.stat <- ks.res(x, y)$statistic
    z <- max(obc.stat, ks.stat)
    return(z)
}
perm.test <- function(x, y, ..., f, num.perm = 2001, diag = FALSE, exact = TRUE) {
  # Runs a permutation test for a given function 
  #(Note: Input function should JUST output p-value) Takes
  # as input f, a function that returns a statistic 
  # Computes the results of a permutation test for a
  # test statistic 
  # Args: x: A vector of observations for a R.V. (must be numeric) y: Either (1) Another
  # vector of observations (two sample) (2) A quantile function such as qnorm (one sample) f: a function
  # which outputs a test statistic. MUST output only a numeric test statistic num.perm: Number of
  # permutations. To avoid messiness, is best if not a multiple of 5 diag: if TRUE, also outputs the
  # resulting values of the test statistics from each permutation exact: if TRUE, then for small samples
  # <11 the function calculates every possible permutation 
  #Returns: A list containing: 
  # [1]: The p-value
  # [2]: The critical value 
  # [3]: The test statistic 
  # [4](if diag=TRUE): The values of the test statistic
  # in each permutation   
    require(gtools)
    lenx <- length(x)
    # This code snippet allows us to take in a function argument such as qgamma
    if (is.character(y)) 
        y <- get(y, mode = "function", envir = parent.frame())
    if (is.function(y)) 
        y <- y(seq(1/(lenx + 1), lenx/(lenx + 1), length.out = lenx), ...)  #Note: quantiles up for debate
    leny <- length(y)
    
    # First step, combine into one dataset
    z <- c(x, y)
    lenz <- length(z)
    # Calculate TS for the ACTUAL data
    ts.obs <- f(x, y)
    ts.random <- c(NULL)
    ### 
    if (lenz < 11 & exact == TRUE) {
        all.perm <- permutations(n = lenz, r = lenz, v = z, repeats.allowed = FALSE, set = FALSE)
        all.permx <- all.perm[, 1:lenx]
        all.permy <- all.perm[, (lenx + 1):lenz]
        exact.perm <- dim(all.perm)[1]
        for (i in 1:exact.perm) {
            ts.random[i] <- f(all.permx[i, ], all.permy[i, ])
        }
        p.val <- sum(abs(ts.random) >= abs(ts.obs))/exact.perm
        c.val <- quantile(ts.random, probs = 0.95)
    } else {
        for (i in 1:num.perm) {
            z1 <- sample(z, size = lenz, replace = FALSE)
            a <- z1[1:lenx]
            b <- z1[(lenx + 1):lenz]
            ts.random[i] <- f(a, b)
        }
        p.val <- sum(abs(ts.random) >= abs(ts.obs))/num.perm
        c.val <- quantile(ts.random, probs = 0.95)
    }
    
    # 1st value of output is p value, 2nd is 95% critical value, 3rd is the actual test statistic
    if (diag == TRUE) {
        return(list(`p-value` = p.val, `95% crit val` = c.val, `Obs. TS` = ts.obs, ts.dist = ts.random))
    } else {
        return(list(`p-value` = p.val, `95% crit val` = c.val, `Obs. TS` = ts.obs))
    }
}
boot.test <- function(x, y, f, num.perm = 1000, diag = FALSE, exact = TRUE) {
  # Runs a BOOTSTRAP test for a given function (Note: Input function should JUST output p-value) Takes
  # as input f, a function that returns a statistic
    require(gtools)
    lenx <- length(x)
    leny <- length(y)
    # First Step, combine into one dataset
    z <- c(x, y)
    lenz <- length(z)
    # Calculate TS for the ACTUAL data
    ts.obs <- f(x, y)
    ts.random <- c(NULL)
    ### 
    if (lenz < 10 & exact == TRUE) {
        all.perm <- permutations(n = lenz, r = lenz, v = z, repeats.allowed = FALSE, set = FALSE)
        all.permx <- all.perm[, 1:lenx]
        all.permy <- all.perm[, (lenx + 1):lenz]
        exact.perm <- dim(all.perm)[1]
        for (i in 1:exact.perm) {
            ts.random[i] <- f(all.permx[i, ], all.permy[i, ])
        }
        p.val <- sum(abs(ts.random) >= abs(ts.obs))/exact.perm
        c.val <- quantile(ts.random, probs = 0.95)
    } else {
        for (i in 1:num.perm) {
            z1 <- sample(z, size = lenz, replace = TRUE)
            a <- z1[1:lenx]
            b <- z1[(lenx + 1):lenz]
            ts.random[i] <- f(a, b)
        }
        p.val <- sum(abs(ts.random) >= abs(ts.obs))/num.perm
        c.val <- quantile(ts.random, probs = 0.95)
    }
    
    # 1st value of output is p value, 2nd is 95% critical value, 3rd is the actual test statistic
    if (diag == TRUE) 
        return(list(`p-value` = p.val, `95% crit val` = c.val, `Obs. TS` = ts.obs, ts.dist = ts.random)) else {
        return(list(`p-value` = p.val, `95% crit val` = c.val, `Obs. TS` = ts.obs))
    }
}
test.ks.obc <- function(x, y) {
  # Right Now, this function returns the values of KS and OBC and a combined value, which is just the
  # addition of the two numbers
    pv.ks <- perm.test(x, y, ks.res.simp)[[1]]  #Based on Distribution
    pv.obc <- perm.test(x, y, myts)[[1]]
    pv.com <- pv.ks + pv.obc
    z <- list(KS = pv.ks, OBC = pv.obc, Comb = pv.com)
    return(z)
}
new.perm.test <- function(x, y, ..., f, num.perm = 2001, diag = FALSE, exact = TRUE) {
  require(gtools)
  lenx <- length(x)
    # Calculate TS for the ACTUAL data
  ts.obs <- f(x, y,...)
  ts.random <- vector(mode="numeric",length=num.perm)
  # Two sample
  if(is.numeric(y)){
    z <- c(x, y)
    lenz <- length(z)
    if (lenz < 11 & exact == TRUE) {
      all.perm <- permutations(n = lenz, r = lenz, v = z, repeats.allowed = FALSE, set = FALSE)
      all.permx <- all.perm[, 1:lenx]
      all.permy <- all.perm[, (lenx + 1):lenz]
      exact.perm <- dim(all.perm)[1]
      for (i in 1:exact.perm) {
        ts.random[i] <- f(all.permx[i, ], all.permy[i, ])
      }
      p.val <- sum(abs(ts.random) >= abs(ts.obs))/exact.perm
      c.val <- quantile(ts.random, probs = 0.95)
    } else {
      for (i in 1:num.perm) {
        z1 <- sample(z, size = lenz, replace = FALSE)
        a <- z1[1:lenx]
        b <- z1[(lenx + 1):lenz]
        ts.random[i] <- f(a, b)
      }
      p.val <- sum(abs(ts.random) >= abs(ts.obs))/num.perm
      c.val <- quantile(ts.random, probs = 0.95)
    }
    if (diag == TRUE) {
      return(list(`p-value` = p.val, `95% crit val` = c.val, `Obs. TS` = ts.obs, ts.dist = ts.random))
    } else {
      return(list(`p-value` = p.val, `95% crit val` = c.val, `Obs. TS` = ts.obs))
    }
  }  
  
  # One Sample
  #if (is.character(y)){
  #  y <- get(y, mode = "function", envir = parent.frame())
  #}
  #if (is.function(y)){
  if(is.character(y)){
    ry <- dist.conv(funname=y, type="r")
    y <- get(y, mode = "function", envir = parent.frame())
    for (i in 1:num.perm){
      z <- ry(lenx,...)
      ts.random[i] <- f(z,y,...)
    }
    p.val <- sum(abs(ts.random) >= abs(ts.obs))/num.perm
    c.val <- quantile(ts.random, probs = 0.95)
    if (diag == TRUE) {
      return(list(`p-value` = p.val, `95% crit val` = c.val, `Obs. TS` = ts.obs, ts.dist = ts.random))
    } else {
      return(list(`p-value` = p.val, `95% crit val` = c.val, `Obs. TS` = ts.obs))
    }
  }
  
}
# This code snippet allows us to take in a function argument such as qgamma

# 1st value of output is p value, 2nd is 95% critical value, 3rd is the actual test statistic

# Outlier Functions

perm.test.out <- function(x, y, ..., f, num.perm = 2001, diag = FALSE, exact = FALSE) {
    # Takes as input f, a function that returns a statistic
    require(gtools)
    lenx <- length(x)
    # This code snippet allows us to take in a function argument such as qgamma
    if (is.character(y)) 
        y <- get(y, mode = "function", envir = parent.frame())
    if (is.function(y)) 
        y <- y(seq(1/(lenx + 1), lenx/(lenx + 1), length.out = lenx), ...)  #Note: quantiles up for debate
    leny <- length(y)
    # Combine both datasets into one
    z <- c(x, y)
    lenz <- length(z)
    # Calculate TS for the ACTUAL data
    ts.obs <- f(x, y)[[1]]
    ts.random <- c(NULL)
    ### 
    if (lenz < 11 & exact == TRUE) {
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
        quan.mat <- matrix(data = NA, nrow = num.perm, ncol = 2)
        for (i in 1:num.perm) {
            z1 <- sample(z, size = lenz, replace = FALSE)
            a <- z1[1:lenx]
            b <- z1[(lenx + 1):lenz]
            a <- sort(a)
            b <- sort(b)
            res <- f(a, b)
            ts.random[i] <- res[[1]]
            ord <- res[[2]]
            quan.mat[i, ] <- c(a[ord], b[ord])
        }
        # Tabulating the values of quan.mat
        tab.quan.mat <- table(quan.mat)
        # Getting the number most repeated
        rep.val <- as(names(which.max(tab.quan.mat)), mode(quan.mat))
        # Getting the number of times it is repeated
        temp <- which.max(tab.quan.mat)[[1]]
        numrept <- tab.quan.mat[temp][[1]]
        # Finds the percent of the time the value repeats
        perc.rep <- numrept/num.perm
        if (perc.rep > 0.9) {
            # value included most
            if (length(which(rep.val - 1e-08 < y & y < rep.val + 1e-08)) > 0) {
                rem <- which(rep.val - 1e-08 < y & y < rep.val + 1e-08)
                y <- y[-rem]
            } else {
                rem <- which(rep.val - 1e-08 < x & x < rep.val + 1e-08)
                x <- x[-rem]
            }
            if (is.integer(rem) == FALSE) {
                stop("Tolerance Failure")
            }
            
            res.rem <- perm.test(x, y, f = myts)
            res.rem[4] <- perc.rep
            res.rem[5] <- c("Outlier Removed")
            res.rem[6] <- rem
            res.rem[7] <- rep.val
            names(res.rem)[4:7] <- c("Percent Reps", "Out Check", "Index", "Repeated Value")
            return(res.rem)
        }
        
        p.val <- sum(abs(ts.random) >= abs(ts.obs))/num.perm
        c.val <- quantile(ts.random, probs = 0.95)
        
    }
    
    # 1st value of output is p value, 2nd is 95% critical value, 3rd is the actual test
    # statistic
    if (diag == TRUE) {
        return(list(`p-value` = p.val, `95% crit val` = c.val, `Obs. TS` = ts.obs, ts.dist = ts.random))
    } else {
        return(list(`p-value` = p.val, `95% crit val` = c.val, `Obs. TS` = ts.obs, `Value Matrix` = quan.mat, 
            `Percent Repeat` = perc.rep, `Remove?` = c("No Outliers Removed")))
    }
}

myts.out <- function(x, y, ..., interp = 4, do.plot = FALSE) {
    # Computes maximum difference of quantiles Args: x: A vector of observations for a
    # R.V. (must be numeric) y: Either (1) Another vector of observations (two sample)
    # (2) A quantile function such as qnorm (one sample) interp: method of interpolation
    # used. For more details, see ?quantile do.plot: Creates a plot illustrating the
    # statistic Returns: The value of the statistic AND associated index of the vector
    
    x <- sort(x)
    # finding lengths
    lenx <- length(x)
    leny <- length(y)
    
    ############ Two Sample
    if (is.numeric(y)) {
        y <- sort(y)
        # finding range values Note: This is controversial. We should consider changing
        # these range values in the future
        x1 <- seq(1/lenx, 1, 1/lenx)
        y1 <- seq(1/leny, 1, 1/leny)
        if (lenx == leny) {
            z <- max(abs(y - x))
            ind <- which.max((abs(y - x)))
            if (do.plot == TRUE) 
                plot.ts.2sam(x, y, x1, y1, lenx, leny)
        } else if (lenx > leny) {
            q1 <- quantile(x, probs = x1, type = interp)
            q_inter <- approxfun(x1, q1, yleft = min(q1), yright = max(q1))
            z <- max(abs(q_inter(y1) - y))
            ind <- which.max((abs(q_inter(y1) - y)))
            if (do.plot == TRUE) 
                plot.ts.2sam(x, y, x1, y1, lenx, leny)
        } else {
            # So length y>x
            q2 <- quantile(y, probs = y1, type = interp)
            q_inter <- approxfun(y1, q2, yleft = min(q2), yright = max(q2))
            z <- max(abs(q_inter(x1) - x))
            ind <- which.max((abs(q_inter(x1) - x)))
            if (do.plot == TRUE) 
                plot.ts.2sam(x, y, x1, y1, lenx, leny)
        }
        return(list(z, ind))
    }
    ############# One Sample
    if (is.function(y)) 
        funname <- as.character(substitute(y))
    if (is.character(y)) {
        funname <- y
    }
    # if (is.character(funname))
    y <- get(funname, mode = "function", envir = parent.frame())
    if (!is.function(y)) 
        stop("'y' must be numeric or a function or a string naming a valid function")
    z <- max(x - y(seq(1/(lenx + 1), lenx/(lenx + 1), length.out = lenx), ...))  #Note: quantiles up for debate
    ind <- which.max(x - y(seq(1/(lenx + 1), lenx/(lenx + 1), length.out = lenx), ...))
    if (do.plot == TRUE) {
        plot.ts.1sam(x, y, ..., funname = funname, lenx = lenx)
    }
    return(list(z, ind))
} 

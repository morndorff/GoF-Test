# Just rewriting the entire outlier code from scratch
# Still using myts.out
new.perm.test.out <- function(x,y, ..., f, fops=NULL, 
                              distops=NULL, num.perm = 2001, diag = FALSE, exact = FALSE) {
  lenx <- length(x)
    # Handing possibly necessary variable transformations
  if (is.function(y)) 
    y <- as.character(substitute(y))

  if (is.character(y)) {
    y <- chartoli(y)
    # Calculating observed test statistic
    if (length(fops) == 0) 
      fops <- NULL
    if (length(distops) == 0) 
      distops <- NULL
    ts.obs <- do.call(f, c(list(x), list(names(y)), distops, fops))[[1]]
  }
  if (is.numeric(y)){
    # Calculating observed test statistic
    ts.obs <- do.call(f, list(x,y, fops))[[1]]
  }
  
  # Calculating distribution of test statistic
  ts.random <- vector(mode = "numeric", length = num.perm)
  
  z <- c(x,y)
  lenz <- length(z)
  
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
  
  #
  # At this point in the function, we have the followigng information:
  # 1) a bunch of test statistics (in ts.random)
  # 2) perc.rep, a variable describing whether to repeat the 
  # test with less data
  
  # If everything is kosher, return the value of new.perm.test
  # using x and y given to the function
  if(perc.rep < .9){
    zout <- do.call(new.perm.test,c(list('x'=x),list('y'=y),list('f'=myts)))#,fops,distops))
    
    return(zout)
  }
  num.loops <- num.loops +1
  # If there is a number of repeated values, call
  # this function again, but with a reduced x and y
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
  num.loops <- num.loops +1
  do.call(new.perm.test.out,c(list(x),list(y),list('f'=f),fops,distops))
}

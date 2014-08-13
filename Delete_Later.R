# Experimentation
perm.test.match.out <- function(x, y, ..., f, fops=NULL, distops= NULL, num.perm = 2001, diag = FALSE, exact = FALSE) {
  # Takes as input f, a function that returns a statistic
  require(gtools)
  lenx <- length(x)
  # Error handling and calculating test statistic -----------------------
  
  # Handling function inputs for y
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
  ts.random <- vector(mode = "numeric", length = num.perm)
  
  # Begin controversial code here -----------------------
    
  z <- c(x,y)
  z_seq <- seq_along(z)
  lenz <- length(z)
  
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
      #z1 <- sample(z, size = lenz, replace = FALSE)
      seq_sam <- sample(z_seq, size = lenz, replace= FALSE)
      sam_x <- z[seq_sam[1:lenx]]
      sam_y <- z[seq_sam[(lenx+1):lenz]]
      sam_x <- sort(sam_x)
      sam_y <- sort(sam_y)
      res <- f(sam_x, sam_y)
      ts.random[i] <- res[[1]]
      ord <- res[[2]]
      x_out <- sam_x[ord]
      y_out <- sam_y[ord]
      x_out_ind <- match(x_out,z)
      y_out_ind <- match(y_out,z)
      quan.mat[i, ] <- c(x_out_ind,y_out_ind)
      #a <- z1[1:lenx]
      #b <- z1[(lenx + 1):lenz]
      #a <- sort(a)
      #b <- sort(b)
      
      #ts.random[i] <- res[[1]]
      #ord <- res[[2]]
      #quan.mat[i, ] <- c(a[ord], b[ord])
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
#       # value included most
#       if (length(which(rep.val - 1e-08 < y & y < rep.val + 1e-08)) > 0) {
#         rem <- which(rep.val - 1e-08 < y & y < rep.val + 1e-08)
#         y <- y[-rem]
#       } else {
#         rem <- which(rep.val - 1e-08 < x & x < rep.val + 1e-08)
#         x <- x[-rem]
#       }
#       if (is.integer(rem) == FALSE) {
#         stop("Tolerance Failure")
#       }
      if(rep.val>lenx){
        y <- y[-(rep.val-lenx)]
      } else{
        x <- x[-rep.val]
      }
      res.rem <- perm.test(x, y, f = myts)
      res.rem[4] <- perc.rep
      res.rem[5] <- c("Outlier Removed")
      res.rem[6] <- rep.val
      res.rem[7] <- z[rep.val]
      #names(res.rem)[4:7] <- c("Percent Reps", "Out Check", "Index", "Repeated Value")
      names(res.rem)[4:7] <- c("Percent Reps", "Out Check", "Index", "Removed Value")
      
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
# Returning to the outlier functions
# attempting to edit so that the behavior is sensical for
# unequal sample sizes. It currently sucks.

perm.test.match.out2 <- function(x, y, ..., f, 
                                 fops=NULL, distops= NULL, 
                                 num.perm = 2001, diag = FALSE, 
                                 exact = FALSE, count=0) {
  require(gtools)
  lenx <- length(x)
  leny <- length(y)
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
  # Equal Sample Sizes
  if(lenx==leny){
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
        count <- count +1 
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
        res.rem[[8]] <- x
        res.rem[[9]] <- y
        res.rem[10] <- count
        names(res.rem)[4:10] <- c("Percent Reps", "Out Check", "Index", 
                                 "Removed Value", "x", "y", "count")
        
        return(res.rem)
      }
    }
  }
  # Unequal sample Sizes
  if(lenx!=leny){
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
      quan.mat <- vector(mode="numeric", length=num.perm)
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
        if(lenx > leny){
          y_out <- sam_y[ord]
          y_out_ind <- match(y_out,z)
          quan.mat[i] <- c(y_out_ind)
        } else {
          x_out <- sam_x[ord]
          x_out_ind <- match(x_out,z)
          quan.mat[i] <- c(x_out_ind)
        }
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
      if (perc.rep > 0.45) {
        count <- count +1 
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
        res.rem[[8]] <- x
        res.rem[[9]] <- y
        res.rem[10] <- count
        names(res.rem)[4:10] <- c("Percent Reps", "Out Check", "Index", 
                                  "Removed Value", "x", "y", "count")        
        return(res.rem)
      }
    }
  }
  
  # If we get here, we haven't returned yet, and we just need
  # to analyze our p-values and output our answers (easy?)
  p.val <- sum(abs(ts.random) >= abs(ts.obs))/num.perm
  c.val <- quantile(ts.random, probs = 0.95)
  
  if (diag == TRUE) {
    return(list(`p-value` = p.val, 
                `95% crit val` = c.val, 
                `Obs. TS` = ts.obs, 
                `ts.dist` = ts.random, 
                `Value Matrix` = quan.mat,
                `Out Check` = c("No Outliers Removed"),
                `x`=x, 
                `y`=y,  
                `count`=count 
    ))
  } else {
    
    return(list(`p-value` = p.val, 
                `95% crit val` = c.val, 
                `Obs. TS` = ts.obs, 
                `Percent Repeat` = perc.rep, 
                `Out Check` = c("No Outliers Removed"), 
                `x`=x, 
                `y`=y,  
                `count`=count ))
  }
  
  
}

perm.test.match.out2.iter <- function(x, y, ..., f, 
                                      fops=NULL, distops= NULL, 
                                      num.perm = 2001, diag = FALSE, 
                                      exact = FALSE, count=0) {
  fcall <- match.call()
  #return(fcall)
  fcall[[1]] <- perm.test.match.out2
  test_results <- eval(fcall, parent.frame())
  if(test_results["Out Check"]=="Outlier Removed"){
    fcall[[1]] <- perm.test.match.out2.iter
    fcall$x <- test_results[["x"]]
    fcall$y <- test_results[["y"]]
    fcall$count <- test_results[["count"]]
    eval(fcall, parent.frame())
  } else{ 
    ind <- which(names(test_results)=="count")
    names(test_results)[ind] <- "Outliers Removed"
    return(test_results)
  } 
}





  
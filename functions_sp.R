# SP Plot Functions
SP_Threshold_Selection <- function(x, y, ...,
                                grid_fine=15, dyadic_after=TRUE, 
                                peek_scaling=TRUE, plotEstFun=F, doplot=F,
                                debug=F, diag=F, plotThresh=F, threshscale=T, includescale=T){
  # includescale: threshold lowest level scaling coefficent?
  # threshscale: Do we want to threshold the detail functions?
  # grid_fine: How fine the cross validation grid is
  # dyadic_after: Making the lambda value smaller or larger 
  # peek_scaling: Effects our scaling. Do we look at all the data when we scale,
  # or do we just look at the in sample data
  # If true, we do, if false we don't
  if(!includescale) stop("have to include scaling right now")
  
  # Only works with dyadic sample size  
  require(wavethresh)
  if(is.numeric(y)) stop("Haven't Done Two Sample Yet")
  
  # One Sample -----
  if (is.list(y)) 
    y <- names(y)
  if (is.function(y)) 
    funname <- as.character(substitute(y))
  if (is.character(y)) 
    funname <- y
  y <- get(funname, mode = "function", envir = parent.frame())
  if (!is.function(y)) 
    stop("'y' must be numeric or a function or a string naming a valid function")
  
  # Calculating Plot Differences -----
  
  x <- sort(x)
  lenx <- length(x)
  i <- 1:lenx
  sp <- (2/pi) * asin(sqrt((i-.5)/lenx))
  si <- (2/pi) * asin(sqrt(y(x, ...))) # Important part
  
  dif <- si-sp
  
  # Split the Samples -------
  if(lenx %% 2==0){
    #print("its even!")
    
    ind <- 1:lenx
    ind_even <- seq(2, lenx, by=2)
    ind_odd <- seq(1, lenx-1, by=2)
    si_even <- dif[ind_even]
    len_even <- lenx / 2
    si_odd <- dif[ind_odd]
    len_odd <- len_even
  } else{
    #print("its odd!")
    ind <- 1:lenx
    ind_even <- seq(2, lenx, by=2)
    ind_odd <- seq(1, lenx-1, by=2)
    
    si_even <- dif[ind_even]
    len_even <- (lenx - 1) / 2
    
    si_odd <- dif[ind_odd]
    len_odd <- len_even + 1
  } 
  # Getting Odd and Even Lambda -----
  Odd_OOS <- Find_Lambda_SP(in_sample_dif = si_even, out_of_sample_dif = si_odd,
                         grid_fine=grid_fine, plotEstFun=plotEstFun, doplot=doplot,
                         insamplefirst=F, plotThresh=plotThresh, threshscale=threshscale,
                         includescale=includescale)
  if(plotThresh) return(Odd_OOS)
  Odd_MSE <- Odd_OOS[["MSE Vector"]]
  Odd_Lambda <- Odd_OOS[["Optimal Lambda"]]
  
  Even_OOS <- Find_Lambda_SP(in_sample_dif = si_odd, out_of_sample_dif =si_even,
                          grid_fine=grid_fine, plotEstFun=plotEstFun, doplot=doplot,
                          insamplefirst=T, plotThresh=plotThresh, threshscale=threshscale,
                          includescale=includescale)
  Even_MSE <- Even_OOS[["MSE Vector"]]
  Even_Lambda <- Even_OOS[["Optimal Lambda"]]

  # Choosing Lambda----------
  
  Chosen_Lambda <- mean(c(Even_Lambda, Odd_Lambda))
  # This Chosen Lambda value needs to be adjusted to 
  # account for the fact that the even and odd samples
  # are of sample size n/2
  Average_Relative_Position <- mean(c(Even_OOS[["Position"]], Odd_OOS[["Position"]]))
  
  
  if(dyadic_after){
    # If the sample is going to be adjusted down or up to a
    # dyadic power of 2 afterwards, we should account for that
    Adjusted_Length <- 2^floor(log2(lenx))
    Adjusted_Lambda <- (1/ sqrt((1 - (log(2)/log(Adjusted_Length))))) * Chosen_Lambda
  } else{
    Adjusted_Lambda <- (1/ sqrt((1 - (log(2)/log(lenx))))) * Chosen_Lambda
  }
  if(diag){
    return(list("Even MSE"=Even_MSE, "Odd MSE"=Odd_MSE, 
                "Even Lambda"=Even_Lambda, "Odd Lambda"= Odd_Lambda,
                "Chosen Threshold (Lambda)"= Chosen_Lambda,
                "Sample Size Adjusted Threshold"=Adjusted_Lambda))
  }
  return(list("Even MSE"=Even_MSE, "Odd MSE"=Odd_MSE, 
              "Even Lambda"=Even_Lambda, "Odd Lambda"= Odd_Lambda,
              "Chosen Threshold (Lambda)"= Chosen_Lambda,
              "Sample Size Adjusted Threshold"=Adjusted_Lambda,
              "Relative Position Chosen"=Average_Relative_Position))
  
}

Find_Lambda_SP <- function(in_sample_dif, out_of_sample_dif, 
                        grid_fine, plotEstFun, doplot, insamplefirst, plotThresh, threshscale, includescale){
  # internal function for cross validation. Not to be used 
  
  # Input: A 'in sample' x and y and an out of sample x and y
  # Output: A list containing the MSE's used as well as the optimal lambda value found

  
  # Decompose input function
  model_wd <- wd(in_sample_dif, filter.number=8, family="DaubLeAsymm") #la8

  if(threshscale){
    Largest_Level <- nlevelsWT(model_wd)
    testfun <- Vectorize(accessD.wd, "level")
    details <- unlist(testfun(model_wd, level=0:(Largest_Level-1)))
    if(includescale) scaling <-  accessC(model_wd, level=0)
    Largest_Coefficient <- max(abs(c(details, scaling))) 
    Smallest_Coefficient <- min(abs(c(details, scaling)))
    lambda_grid <- seq(Smallest_Coefficient, 
                       Largest_Coefficient, length.out=grid_fine)

  } else{
  # Make lambda grid
  Largest_Level <- nlevelsWT(model_wd)
  testfun <- Vectorize(accessD.wd, "level")
  details <- unlist(testfun(model_wd, level=1:(Largest_Level-1)))
  if(includescale) scaling <- accessC(model_wd, level=0)
  #scaling <- accessC(model_wd, level=1)
  Largest_Coefficient <- max(abs(c(details, scaling)))
  Smallest_Coefficient <- min(abs(c(details, scaling)))
  lambda_grid <- seq(Smallest_Coefficient, 
                     Largest_Coefficient, length.out=grid_fine)
  }
  # Making sure last threshold thresholds everything
  # lambda_grid <- append(lambda_grid, Largest_Coefficient + .1*sd(lambda_grid))
  
  # For values in the lambda grid,
  # Perform thresholding, reconstruction, and MSE comparison
  MSE <- vector(length=length(lambda_grid))
  Thresh_list <- vector(mode="list", length=length(lambda_grid))
  for(i in seq_along(lambda_grid)){
    # Thresholding
    if(threshscale){
      thresh_test <- threshold(model_wd, type="hard", policy="manual", 
                               value=lambda_grid[i], levels=0:(nlevelsWT(model_wd)-1) )#, verbose=TRUE)
      if(includescale){
        Lowest_C <- abs(accessC(thresh_test, level=0))
        if(Lowest_C <= lambda_grid[i]){
          thresh_test <- putC(thresh_test, level=0, v=0)
        }
      }
    }
    else{
      thresh_test <- threshold(model_wd, type="hard", policy="manual", 
                               value=lambda_grid[i], levels=1:(nlevelsWT(model_wd)-1)) #, verbose=TRUE)
      if(includescale){
        Lowest_C <- abs(accessC(thresh_test, level=0))
        if(Lowest_C < lambda_grid[i]){
          thresh_test <- putC(thresh_test, level=0, v=0)
        }
      }
    }
    # Reconstruction
    reconstructed_fun <- wr(thresh_test)
    # Here, we estimate the odd from the even, and vice versa
    
    # First, we need to do some endpoint correction
    # Then, we take our estimate as the midpoint between entries
    ma <- function(x,n=2){filter(x,rep(1/n,n), sides=2)}
    if(insamplefirst){
      Dif_First <- reconstructed_fun[1] / 2 # Take midpoint  
      estimated_fun <- ma(reconstructed_fun)
      estimated_fun <- c(Dif_First, head(estimated_fun, -1))
    }else{
      #Interpolate maximum
      Dif_End <- tail(reconstructed_fun , 1) / 2
      estimated_fun <- ma(reconstructed_fun)
      estimated_fun <- c(head(estimated_fun, -1), Dif_End)
    }
    
    # calculating MSE
    errors <- estimated_fun - out_of_sample_dif
    MSE[i] <- sum(errors^2)
    if(plotThresh) Thresh_list[[i]] <- (list("Lambda Value"=paste(i, "of", length(lambda_grid)), 
                               "Thresholded Coefficients"=thresh_test, 
                               "Reconstructed Function"=estimated_fun, 
                               "Unthresholded Wavelet"=model_wd))
    
    if(plotEstFun){
      #layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
      par(mfrow=c(2,1))
      plot(in_sample_dif, type="l",
           main=paste("Lambda Position=", round(i/length(lambda_grid), digits=2), "Estimated In Red"),
           ylab="In Sample Differences")
      points(estimated_fun, type="l", col="red")
      plot(out_of_sample_dif, type="l", ylab="Out of Sample Differences",
           main=paste("MSE is ", round(MSE[i],3),  "Odd in Red"))
      lines(estimated_fun,type="l", col="red")
    }
  }
  
  if(doplot){
    par(mfrow=c(1,1))
    plot(lambda_grid, MSE)
  }
  if(plotThresh) return(Thresh_list)
  Position <- which.min(MSE)
  opt_lambda <- lambda_grid[Position]
  Relative_Position <- Position/length(lambda_grd)
  return(list("MSE Vector"= MSE, "Optimal Lambda"= opt_lambda, "Position"=Relative_Position))
}

SP_Test<- function(x, y, ..., grid=10, Chosen_Threshold, spacing=F, doplot=F, threshscale=T, includescale=T){
  # Wrapper function.
  # uses wd for convenience. Change this later to use other wavelet basis
  
  # Carry out the thresholding, using a supplied value of 
  # Chosen_Threshold
  if(!threshscale) stop("Needs Update!")
  if(!includescale) stop("Needs Update!")
  
  
  if (is.list(y)) 
    y <- names(y)
  if (is.function(y)) 
    funname <- as.character(substitute(y))
  if (is.character(y)) 
    funname <- y
  y <- get(funname, mode = "function", envir = parent.frame())
  if (!is.function(y)) 
    stop("'y' must be numeric or a function or a string naming a valid function")
  
  # Use lambda, threshold appropriately
  x <- sort(x)
  n <- length(x)
  i <- 1:n
  sp <- (2/pi) * asin(sqrt((i-.5)/n))
  si <- (2/pi) * asin(sqrt(y(x, ...)))
  # print(si)
  if(doplot){
    par(mfrow=c(2,1))
    plot(sp,si, xlim=c(0,1), ylim=c(0,1), main="SP-Plot")
    abline(a=0, b=1)
    plot(sp-si)
  }
  z <- sp-si
  #   library(waveslim)
  #   F.dwt <- dwt(F.x(z) - y(z))
  #   coefs <- unlist(F.dwt)

  library(wavethresh)

  model_wd <- wd(z, filter.number=8, family="DaubLeAsymm")

  Largest_Level <- nlevelsWT(model_wd)
  testfun <- Vectorize(accessD.wd, "level")
  details <- unlist(testfun(model_wd, level=0:(Largest_Level-1)))
  scaling <- accessC(model_wd, level=0)
  coefs <- c(details, scaling)
  
  threshold <- sapply(coefs, function(x){
    if(abs(x) < Chosen_Threshold) x <- 0
    x
  })
  STAT <- sum(threshold^2)
  STAT
}
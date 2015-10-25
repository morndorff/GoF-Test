# Cross Validation
rm(list=ls())


Even_Odd_Threshold <- function(x, y, ..., lambda=FALSE, doplot=T, 
                       grid_fine=15, opt="universal"){ 
  if(is.numeric(y)) stop("Haven't Done Two Sample Yet")
  
  # x <- rnorm(40)
  lenx <- length(x)
  x <- sort(x)
  # is even or odd
  if(lenx %% 2==0){
    print("its even!")
    ind <- 1:lenx
    ind_even <- seq(2, lenx, by=2)
    ind_odd <- seq(1, lenx-1, by=2)
    x_even <- x[ind_even]
    len_even <- lenx / 2
    x_odd <- x[ind_odd]
    len_odd <- len_even
  } else{
    print("its odd!")
    ind <- 1:lenx
    ind_even <- seq(2, lenx, by=2)
    ind_odd <- seq(1, lenx-1, by=2)
    x_even <- x[ind_even]
    len_even <- (lenx - 1) / 2
    x_odd <- x[ind_odd]
    len_odd <- len_even + 1
  } 
  # Getting lengths (will need later)

  
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
  if(doplot){
    par(mfrow=c(1,1))
    plot(ecdf(x))
    temp <- seq(min(x), max(x), length.out=200)
    lines(temp, y(temp))
  }
  

  # Create even and odd data sets
  # even
  F.x <- ecdf(x_even)
  F.x <- F.x(x_even)
  F.x <- F.x - (.5 /len_even)
  Dif_Even <- F.x - y(x_even, ...)
  #odd
  F.x <- ecdf(x_odd)
  F.x <- F.x(x_odd)
  F.x <- F.x - (.5 /len_odd)
  Dif_Odd <- F.x - y(x_odd, ...)
  
  even_len <- length(Dif_Even)
  odd_len <- length(Dif_Odd)
  
  
  
  # Transforming to grid (temporary solution)
  x_even_scale <- (x_even - min(x_even)) / (max(x_even) - min(x_even))
  x_odd_scale <- (x_odd - min(x_odd)) / (max(x_odd) - min(x_odd))
  
  even_grid <- makegrid(t=x_even_scale, y=Dif_Even)
  even_grid <- even_grid$gridy
  
  odd_grid <- makegrid(t=x_odd_scale, y=Dif_Odd)
  odd_grid <- odd_grid$gridy
  
  Find_Lambda <- function(in_sample_data, out_of_sample_data){
    # in_sample_data = even_grid
    #out_of_sample_data= odd_grid
    
          # Do wavelet decomp on even entries
    model_wd <- wd(in_sample_data)#, family="DaubExPhase", filter.number=1) # Add more options here later
    
    # model_wd = even_wd
    # Find Universal Threshold
    wavelet_coefs <- NULL
    for(i in 3:(nlevelsWT(model_wd)-1))  
      wavelet_coefs <- c(wavelet_coefs, accessD(model_wd, level=i))
    noise_level <- mad(wavelet_coefs) #noise.level
    universal_threshold <- noise_level * sqrt(2*log(even_len))
    
    # Create lambda grid
    min_lambda <- universal_threshold - 15 * sd(in_sample_data)
    max_lambda <- universal_threshold + 3 * sd(in_sample_data)
    if(min_lambda < 0) min_lambda=0
    lambda_grid <- seq(min_lambda, max_lambda, length.out=grid_fine) # 10 is default
    
    # Thresholding
    MSE <- vector(length=length(lambda_grid))
    for(i in seq_along(lambda_grid)){
      #Threshold
      thresh_test <- threshold(model_wd, type="hard", policy="manual", 
                               value=lambda_grid[i], levels=0:(nlevelsWT(model_wd)-1),
                               verbose=TRUE)

      # Reconstruct
      reconstructed_fun <- wr(thresh_test)
      #par(mfrow=c(2,2))
      errors <- reconstructed_fun - out_of_sample_data
      
      MSE[i] <- sum(errors^2)
      if(doplot){
      layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
      plot(even_grid, type="l")
      plot(reconstructed_fun, type="l")
      plot(odd_grid, type="l", main=paste("MSE is ", round(MSE[i],2),  "Odd in Red"))
      lines(reconstructed_fun,type="l", col="red")
    }
    }
    opt_lambda <- lambda_grid[which.min(MSE)]
    return(list("MSE Vector"= MSE, "Optimal Lambda"= opt_lambda))
  }
  
  Odd_OOS <- Find_Lambda(in_sample_data = even_grid, out_of_sample_data = odd_grid)
  Odd_MSE <- Odd_OOS[["MSE Vector"]]
  Odd_Lambda <- Odd_OOS[["Optimal Lambda"]]
  Even_OOS <- Find_Lambda(in_sample_data = odd_grid, out_of_sample_data = even_grid)
  Even_MSE <- Even_OOS[["MSE Vector"]]
  Even_Lambda <- Even_OOS[["Optimal Lambda"]]
  Chosen_Lambda <- mean(c(Even_Lambda, Odd_Lambda))
  
  if(opt=="universal"){
    lenx <- length(x)
    F.x <- ecdf(x)
    F.x <- F.x(x)
    F.x <- F.x - (.5 /lenx)
    Dif_X <- F.x - y(x, ...)
    # Make the grid
    x_scale <- (x - min(x)) / (max(x) - min(x))
    x_grid <- makegrid(t=x_scale, y=Dif_X)
    x_grid <- x_grid$gridy
    x_wd <- wd(x_grid)
    u_threshold <- threshold(x_wd, return.threshold=TRUE, 
                             type="hard", policy="universal")
    
    return(list("Even MSE"=Even_MSE, "Odd MSE"=Odd_MSE, 
         "Even Lambda"=Even_Lambda, "Odd Lambda"= Odd_Lambda,
         "Chosen Threshold (Lambda)"= Chosen_Lambda, 
         "Universal Threshold" = u_threshold))
  }
  
  
  return(list("Even MSE"=Even_MSE, "Odd MSE"=Odd_MSE, 
              "Even Lambda"=Even_Lambda, "Odd Lambda"= Odd_Lambda,
              "Chosen Threshold (Lambda)"= Chosen_Lambda))
  
  
#   
#   # Do wavelet decomp on even entries
#   even_wd <- wd(even_grid)#, family="DaubExPhase", filter.number=1) # Add more options here later
#   
#   # Find Universal Threshold
#   wavelet_coefs <- NULL
#   for(i in 3:(nlevelsWT(even_wd)-1))  
#     wavelet_coefs <- c(wavelet_coefs, accessD(even_wd, level=i))
#   noise_level <- mad(wavelet_coefs) #noise.level
#   universal_threshold <- noise_level * sqrt(2*log(even_len))
#   
#     # Create lambda grid
#   min_lambda <- universal_threshold - 15 * sd(even_grid)
#   max_lambda <- universal_threshold + 3 * sd(even_grid)
#   print(min_lambda)
#   print(max_lambda)
#   if(min_lambda < 0) min_lambda=0
#   lambda_grid <- seq(min_lambda, max_lambda, length.out=grid_fine) # 10 is default
#   
#   MSE <- vector(length=length(lambda_grid))
#   for(i in seq_along(lambda_grid)){
#     #Threshold
#     thresh_test <- threshold(even_wd, type="hard", policy="manual", 
#                              value=lambda_grid[i], levels=0:(nlevelsWT(even_wd)-1),
#                              verbose=TRUE)
#     #print(lambda_grid[i])
#     #print(thresh_test)
#     # Reconstruct
#     reconstructed_fun <- wr(thresh_test)
#     #par(mfrow=c(2,2))
#     errors <- reconstructed_fun - odd_grid
#     
#     MSE[i] <- sum(errors^2)
#     layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
#     plot(even_grid, type="l")
#     plot(reconstructed_fun, type="l")
#     plot(odd_grid, type="l", main=paste("MSE is ", round(MSE[i],2),  "Odd in Red"))
#     lines(reconstructed_fun,type="l", col="red")
#   }
#   
#   if(doplot){
#     par(mfrow=c(1,2))
#     plot(x_odd, Dif_Odd, type="l")
#     plot(x_even, Dif_Even, type="l")
#   }
#   
#   return(MSE)
  
  # return(list("Even X"=x_even, "Odd X"=x_odd, "Even Dif"=Dif_Even, "Odd Dif"=Dif_Odd))
}

set.seed(5)
x <- rnorm(100)
data <- Even_Odd_Threshold(x,pnorm)
wd_data <- wd(x)
wd_thresh <- threshold(wd_data, type="hard", policy="manual",
                       value=data[["Chosen Threshold (Lambda)"]], 
                       evels=0:(nlevelsWT(wd_data)-1))

# what is the universal threshold for this data?
threshold


# plot(data)




even_dif <- data[["Even Dif"]]
odd_dif <- data[["Odd Dif"]]
even_len <- length(even_dif)
odd_len <- length(odd_dif)
x_even <- data[["Even X"]]
x_odd <- data[["Odd X"]]
# Gridding things



######### Verifying Universal Thesholding
thresh_d <- threshold(even_wd, type="hard", policy="universal")
threshold(even_wd, type="hard", policy="universal"
            return.threshold=TRUE, verbose=TRUE)
thresh_confirm <- NULL; 
# Default thresholding leaves first 2 levels unthresholded
for(i in 3:(nlevelsWT(even_wd)-1))  
  thresh_confirm <- c(testd3, accessD(even_wd, level=i))
MAD_test <- mad(thresh_confirm) #noise.level
thresh_test_value <- MAD_test * sqrt(2*log(even_len))
thresh_test_value
thresh_test <- threshold(even_wd, type="hard", policy="manual", 
                          value=thresh_test_value)
############# 

# Calculate IMSE 
even_threshed <- wr(thresh_d)
errors <- odd_grid - even_threshed
# par(mfrow=c(2,2))
# plot(odd_grid)
# plot(even_threshed)
# plot(errors)
MSE <- sum(errors^2)
















# #universal (don't supply lambda yet)
# denoised <- wr(thresh_d)
# plot(denoised)
# plot(odd_dif)
# MSE <- sum((denoised - odd_grid)^2)
# 
# #testing
# high_level_coef <- accessD(even_wd, level=even_wd$nlevels-1)
# MAD_decomp <- mad(high_level_coef)
# thresh_value <- MAD_decomp*sqrt(2*log(even_len))
# thresh_value2 <- MAD_decomp^2 * sqrt(2*log(even_len))
# thresh_test <- threshold(even_wd, type="hard", policy="manual", 
#                          value=thresh_value)
# testd <- accessD(even_wd, level=4)
# dev(d)
# 
# 
# thresh_test2 <- threshold(even_wd, type="hard", policy="manual", 
#                           value=thresh_value2)
# testd3 <- NULL; 
# for(i in 3:nlevelsWT(even_wd)) 
#   testd3 <- c(testd3, accessD(even_wd, level=(i-1)))
# MAD_decomp3 <- mad(testd3) #noise.level
# thresh_value3 <- MAD_decomp3 * sqrt(2*log(even_len))
# thresh_test3 <- threshold(even_wd, type="hard", policy="manual", 
#                           value=thresh_value3)
# 
# 
# 
# denoised 
# test_denoised <- wr(thresh_test)
# test_denoised
# test_denoised2 <- wr(thresh_test2)
# test_denoised2
# test_denoised3 <- wr(thresh_test3)
# test_denoised3


set.seed(5)
x <- rnorm(40)
x <- sort(x)
lenx <- length(x)
# is even or odd
if(lenx %% 2==0){
  print("its even!")
  ind <- 1:lenx
  ind_even <- seq(2, lenx, by=2)
  ind_odd <- seq(1, lenx-1, by=2)
} else{
  print("its odd!")
  ind <- 1:lenx
  ind_even <- seq(2, lenx, by=2)
  ind_odd <- seq(1, lenx-1, by=2)
} 
y <- pnorm
# We need to get \hat(F)_X(x) - F_X(x)
F.x <- ecdf(x)
F.x <- F.x(x)
F.x <- F.x - .5/lenx
X_Dif <- F.x - pnorm(x)

#OK, we now have the two difference functions
# Let's look at how the threshold does with lambda
wd(X_Dif)








# Tangent 
# Variance of \hat{F}_X(x)- F_X(x)
set.seed(7)
rows <- 100
cols <- 50
x <- matrix(rnorm(cols*rows), nrow=rows, ncol=cols)
x[1, ]
x <- t(apply(x, 1, sort))
origx <- x
x[1, ]
x <- t(apply(x, 1, pnorm))
x[1, ]
y <- matrix(nrow=rows, ncol=cols)
for(i in 1:rows){
  F.x <- ecdf(x[i,])
  y[i, ] <- F.x(x[i,])
}
ex <- y- x

plot(origx[1, ], ex[1,], type="l", xlim=c(-3.2,3.2), ylim=c(-.2,.2))
for(i in 1:rows){
  lines(origx[i,], ex[i,])
}


apply(x,2,lines)

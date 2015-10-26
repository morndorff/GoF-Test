
Even_Odd_Threshold2 <- function(x, y, ..., lambda=FALSE, doplot=T, 
                               grid_fine=10, opt=NULL){ 
  # This version of the function applies the gridding with all the data
  # The previous version gridded even seperately from odd
  
  if(is.numeric(y)) stop("Haven't Done Two Sample Yet")
  # x <- rnorm(40)
  
  lenx <- length(x)
  x <- sort(x)
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

  # Getting the difference between the ecdf(X) and its CDF
  lenx <- length(x)
  x <- sort(x)
  F.x <- ecdf(x)
  F.x <- F.x(x)
  F.x <- F.x - (.5 /lenx)
  Dif_X <- F.x - y(x, ...)
  length(Dif_X)
  
  # Scaling and putting things on a grid
  x_scale <- (x - min(x)) / (max(x)- min(x))
  x_grid <- makegrid(t=x_scale, y=Dif_X)
  dif_grid <- x_grid$gridy
  data_grid <- x_grid$gridt

  plot(data_grid, dif_grid, type="l")
  
  # Create even and odd data sets
  lengrid <- length(dif_grid)
  len_even <- lengrid / 2
  len_odd <- len_even
  ind <- 1:lengrid
  
  ind_even <- seq(2, lengrid, by=2)
  ind_odd <- seq(1, lengrid-1, by=2)
  
  Dif_Even <- dif_grid[ind_even]
  Dif_Odd <- dif_grid[ind_odd]

  Data_Even <- data_grid[ind_even]
  Data_Odd <- data_grid[ind_odd]
  
#   # even
#   F.x <- ecdf(x_even)
#   F.x <- F.x(x_even)
#   F.x <- F.x - (.5 /len_even)
#   Dif_Even <- F.x - y(x_even, ...)
#   #odd
#   F.x <- ecdf(x_odd)
#   F.x <- F.x(x_odd)
#   F.x <- F.x - (.5 /len_odd)
#   Dif_Odd <- F.x - y(x_odd, ...)
  
  # print(Dif_Even)
  # print(Dif_Odd)
  par(mfrow=c(2,1))
  plot(Data_Even, Dif_Even)
  lines(Data_Odd, Dif_Odd, col="red")
  
  Find_Lambda <- function(in_sample_dif, out_of_sample_dif, in_sample_data, out_of_sample_data){
    # in_sample_dif = even_grid
    #out_of_sample_dif= odd_grid
    len_in_samp <- length(in_sample_dif)
    
    # Do wavelet decomp on even entries
    model_wd <- wd(in_sample_dif)#, family="DaubExPhase", filter.number=1) # Add more options here later
    
    # model_wd = even_wd
    # Find Universal Threshold
    wavelet_coefs <- NULL
    for(i in 3:(nlevelsWT(model_wd)-1))  
      wavelet_coefs <- c(wavelet_coefs, accessD(model_wd, level=i))
    noise_level <- mad(wavelet_coefs) #noise.level
    universal_threshold <- noise_level * sqrt(2*log(len_in_samp))
    
    # Create lambda grid from univesal threshold
    # Very heuristic right now
    min_lambda <- universal_threshold - 15 * sd(in_sample_dif)
    max_lambda <- universal_threshold + 3 * sd(in_sample_dif)
    if(min_lambda < 0) min_lambda=0
    lambda_grid <- seq(min_lambda, max_lambda, length.out=grid_fine) # 10 is default
    
    # Thresholding
    MSE <- vector(length=length(lambda_grid))
    for(i in seq_along(lambda_grid)){
      #Threshold
      thresh_test <- threshold(model_wd, type="hard", policy="manual", 
                               value=lambda_grid[i], levels=0:(nlevelsWT(model_wd)-1))#, verbose=TRUE)
      
      # Reconstruct
      reconstructed_fun <- wr(thresh_test)
      #par(mfrow=c(2,2))
      # HERE WE DO THE INTERPOLATION
      if(min(in_sample_data) > min(out_of_sample_data)){
        library(truncdist)
        #Interpolate minimum
        begin_point <- extrunc("norm", b=min(in_sample_data), ...)
        print(paste("Begin Point is", begin_point))
        print(min(in_sample_data))
        Dif_End <- ((1 - .05)/ len_even) - y(begin_point)
        #print(Dif_End)
      }else{
        #Interpolate maximum
        
        #Estimated odd end point
        end_point <- extrunc("norm", a=max(in_sample_data), ...)
        print(paste("End Point is", end_point))
        Dif_End <- ((len_even-.5)/ len_even) - y(end_point)
      }
      
      
      errors <- reconstructed_fun - out_of_sample_dif
      
      MSE[i] <- sum(errors^2)
      if(doplot){
        layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
        plot(in_sample_dif, type="l")
        plot(reconstructed_fun, type="l")
        plot(out_of_sample_dif, type="l", main=paste("MSE is ", round(MSE[i],2),  "Odd in Red"))
        lines(reconstructed_fun,type="l", col="red")
      }
      
      
    }
    plot(lambda_grid, MSE)
    opt_lambda <- lambda_grid[which.min(MSE)]
    return(list("MSE Vector"= MSE, "Optimal Lambda"= opt_lambda))
  }
  print(length(Dif_Even))
  Odd_OOS <- Find_Lambda(in_sample_dif = Dif_Even, out_of_sample_dif = Dif_Odd,
                         in_sample_data = Data_Even, out_of_sample_data = Data_Odd)
  Odd_MSE <- Odd_OOS[["MSE Vector"]]
  Odd_Lambda <- Odd_OOS[["Optimal Lambda"]]
  Even_OOS <- Find_Lambda(in_sample_dif = Dif_Odd, out_of_sample_dif = Dif_Even,
                          in_sample_data = Data_Odd, out_of_sample_data = Data_Even)
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



### OLD CODE






# # Transforming to grid (temporary solution)
# x_even_scale <- (x_even - min(x_even)) / (max(x_even) - min(x_even))
# x_odd_scale <- (x_odd - min(x_odd)) / (max(x_odd) - min(x_odd))
# 
# even_grid <- makegrid(t=x_even_scale, y=Dif_Even)
# even_grid <- even_grid$gridy
# 
# odd_grid <- makegrid(t=x_odd_scale, y=Dif_Odd)
# odd_grid <- odd_grid$gridy
# 
# Find_Lambda <- function(in_sample_data, out_of_sample_dif){
#   # in_sample_data = even_grid
#   #out_of_sample_dif= odd_grid
#   
#   # Do wavelet decomp on even entries
#   model_wd <- wd(in_sample_data)#, family="DaubExPhase", filter.number=1) # Add more options here later
#   
#   # model_wd = even_wd
#   # Find Universal Threshold
#   wavelet_coefs <- NULL
#   for(i in 3:(nlevelsWT(model_wd)-1))  
#     wavelet_coefs <- c(wavelet_coefs, accessD(model_wd, level=i))
#   noise_level <- mad(wavelet_coefs) #noise.level
#   universal_threshold <- noise_level * sqrt(2*log(even_len))
#   
#   # Create lambda grid
#   min_lambda <- universal_threshold - 15 * sd(in_sample_data)
#   max_lambda <- universal_threshold + 3 * sd(in_sample_data)
#   if(min_lambda < 0) min_lambda=0
#   lambda_grid <- seq(min_lambda, max_lambda, length.out=grid_fine) # 10 is default
#   
#   # Thresholding
#   MSE <- vector(length=length(lambda_grid))
#   for(i in seq_along(lambda_grid)){
#     #Threshold
#     thresh_test <- threshold(model_wd, type="hard", policy="manual", 
#                              value=lambda_grid[i], levels=0:(nlevelsWT(model_wd)-1),
#                              verbose=TRUE)
#     
#     # Reconstruct
#     reconstructed_fun <- wr(thresh_test)
#     #par(mfrow=c(2,2))
#     errors <- reconstructed_fun - out_of_sample_dif
#     
#     MSE[i] <- sum(errors^2)
#     if(doplot){
#       layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
#       plot(even_grid, type="l")
#       plot(reconstructed_fun, type="l")
#       plot(odd_grid, type="l", main=paste("MSE is ", round(MSE[i],2),  "Odd in Red"))
#       lines(reconstructed_fun,type="l", col="red")
#     }
#   }
#   opt_lambda <- lambda_grid[which.min(MSE)]
#   return(list("MSE Vector"= MSE, "Optimal Lambda"= opt_lambda))
# }
# 
# Odd_OOS <- Find_Lambda(in_sample_data = even_grid, out_of_sample_dif = odd_grid)
# Odd_MSE <- Odd_OOS[["MSE Vector"]]
# Odd_Lambda <- Odd_OOS[["Optimal Lambda"]]
# Even_OOS <- Find_Lambda(in_sample_data = odd_grid, out_of_sample_dif = even_grid)
# Even_MSE <- Even_OOS[["MSE Vector"]]
# Even_Lambda <- Even_OOS[["Optimal Lambda"]]
# Chosen_Lambda <- mean(c(Even_Lambda, Odd_Lambda))

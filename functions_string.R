# String Manipulation Functions
dist.conv <- function(funname, type){
  # Getting a function name and converting it to a distribution function or 
  # random variable genderation
  # Input
  # funname: a function name, such as qnorm
  # type: either p (returns pnorm) or r (returns rnorm)
  # Returns:
  # 
  if(is.function(funname)){
    funname <- as.character(substitute(funname))
  }
  if(is.character(funname)){
  dist <- strsplit(funname, "q")[[1]][2]
  tdist <- paste(type, dist, sep = "")
  tdist <- get(tdist, mod = "function", envir = parent.frame())
  }
  return(tdist)
}


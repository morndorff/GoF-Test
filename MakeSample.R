# Creating the MakeSample Function
make_sample <- function(n, dist, params) {
    # Args: dist: 'norm', 'unif', etc.  n: sample size params: list(min=0,max=2)
    dist <- eval.parent(dist[[1]], n = 1)
    dist <- as.character(substitute(dist))
    rdist <- paste("r", dist, sep = "")
    rdist <- get(rdist, mode = "function", envir = parent.frame())
    sample <- do.call(rdist, c(list(n = n), params))
    return(sample)
}
sam <- make_sample(n = 100, dist = unif, params = list(min = 0, max = 2))
sam <- make_sample(n = 1000, dist = norm, params = list(mean = 0, sd = 2)) 

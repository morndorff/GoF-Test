# Function archive

perm.test <- function(x, y, ..., f, num.perm = 2001, diag = FALSE, exact = TRUE) {
  # Runs a permutation test for a given function (Note: Input function should JUST
  # output p-value) Takes as input f, a function that returns a statistic Computes the
  # results of a permutation test for a test statistic Args: x: A vector of
  # observations for a R.V. (must be numeric) y: Either (1) Another vector of
  # observations (two sample) (2) A quantile function such as qnorm (one sample) f: a
  # function which outputs a test statistic. MUST output only a numeric test statistic
  # num.perm: Number of permutations. To avoid messiness, is best if not a multiple of
  # 5 diag: if TRUE, also outputs the resulting values of the test statistics from
  # each permutation exact: if TRUE, then for small samples <11 the function
  # calculates every possible permutation Returns: A list containing: [1]: The p-value
  # [2]: The critical value [3]: The test statistic [4](if diag=TRUE): The values of
  # the test statistic in each permutation
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
    for (i in 1:num.perm) {
      z1 <- sample(z, size = lenz, replace = FALSE)
      a <- z1[1:lenx]
      b <- z1[(lenx + 1):lenz]
      ts.random[i] <- f(a, b)
    }
    p.val <- sum(abs(ts.random) >= abs(ts.obs))/num.perm
    c.val <- quantile(ts.random, probs = 0.95)
  }
  
  # 1st value of output is p value, 2nd is 95% critical value, 3rd is the actual test
  # statistic
  if (diag == TRUE) {
    return(list(`p-value` = p.val, `95% crit val` = c.val, `Obs. TS` = ts.obs, ts.dist = ts.random))
  } else {
    return(list(`p-value` = p.val, `95% crit val` = c.val, `Obs. TS` = ts.obs))
  }
}

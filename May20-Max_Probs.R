# May 19 
# This script and the associated R Environment 'max_prob' 
# illustrate the problems of using the
# maximum in several distributions The following objects are imporatnt 
# x - A 80x200 samples from  N(0,3) 
# y - A 80x200 samples from t(3) 
# test3 -results of the max(quantile difs) perm test


#Note: Problem with permutation. Look at the following:
test3 <- perm.test(x[2, ], y[2, ], myts)  #This is from the 2nd sample, 200 draws from N(0,3), t(3)
# Note that this returns a p-value of 1, despite the obvious differences in the distribution
myts(x[2, ], y[2, ], do.plot = TRUE)
# Why is this happening?  For one thing, the biggest differences often occur at the min and the max of
# the distribtion. If we mix using permutation, these maxes are likely to land in different samples,
# so the result is exactly the same as before, and the observed TS is the same So we get a bunch of
# ties
hist(test3[[4]])
table(test3[[4]])
## More on the problem:
probx <- x[2, ]
proby <- y[2, ]
myts(probx, proby)
prob <- perm.test(probx, proby, myts)
prob[1:3]
table(prob[[4]])
# Why does this happen?
sortx <- sort(probx)
sorty <- sort(proby)
difxy <- sortx - sorty
com.xy <- sort(c(probx, proby))
# plot(com.xy) A plot to illustrate some of the problem
plot(ecdf(com.xy))
x1 <- seq(1/length(probx), 1, length.out = length(probx))
points(sortx, x1, col = "red")
y1 <- seq(1/length(proby), 1, length.out = length(proby))
points(sorty, y1, col = "blue")
plot(ecdf(com.xy))
segments(sortx, x1, x1 = sorty, y1 = y1, col = "red")
# want biggest in blue (This code snippet will be useful later)
a <- which.max(abs(difxy))
segments(sortx[a], x1[a], x1 = sorty[a], y1 = y1[a], col = "blue")
# Two notes: (1)Divergence at tails! (2) Even if we do some permutation, we are still screwed, since
# the max distance will be very similar Some graphs of t(3) and N(0,3) , overlayed
grapht.x <- seq(-15, 15, length = 300)
grapht.y <- dt(grapht.x, 3)
graphn.x <- seq(-15, 15, length = 300)
graphn.y <- dnorm(graphn.x, 0, 3)
plot(grapht.x, grapht.y, type = "l", lwd = 2)
lines(graphn.x, graphn.y, col = "red", lwd = 2)
# Zoom in!
grapht.x <- seq(-15, -5, length = 300)
grapht.y <- dt(grapht.x, 3)
graphn.x <- seq(-15, -5, length = 300)
graphn.y <- dnorm(graphn.x, 0, 3)
plot(grapht.x, grapht.y, type = "l", lwd = 2)
lines(graphn.x, graphn.y, col = "red", lwd = 2)
############ To do tomorrow: permute things look at the differences, and see how its the same?  May 19
sm.x <- tail(sortx)
sm.y <- tail(sorty)
plot(ecdf(sm.x))
plot(ecdf(sm.y))
sm.x1 <- seq(1/6, 1, length.out = 6)
sm.y1 <- seq(1/6, 1, length.out = 6)
points(sm.x, sm.x1, col = "red")
sm_test <- perm.test(sm.x, sm.y, myts)
sm_test_ks <- perm.test(sm.x, sm.y, ks.res.simp)
sm_test_ks_trad <- ks.test(sm.x, sm.y)
###################### Toy Example to illustrate problem
toy.x <- c(1, 2, 4, 3, 500)
toy.y <- c(1, 2, 2, 4, 5)
toy.x1 <- seq(1/5, 1, length.out = 5)
toy.y1 <- seq(1/5, 1, length.out = 5)
toy.com <- c(toy.x, toy.y)
sample(toy.com, replace = FALSE)
toy_test <- perm.test(toy.x, toy.y, myts)
toy_test
# We can manipulate this a bit to produce p-value not equal to 1
toy.x <- c(1, 2, 5, 3, 500)
toy.y <- c(1, 2, 2, 4, 4)
toy_test <- perm.test(toy.x, toy.y, myts)
toy_test
perm.test(toy.x, toy.y, ks.res.simp)
ks.test(toy.x, toy.y)
# To be fair, everything else fails in this particular example. But the problem is that our problem
# persists even if other differences appear Max of the quantiles isn't a 'stable' enough statistic to
# be approrpriate? This is shown in the power studies The p-value becomes '1', often when the largest
# value is far away and the 2nd largest value resides in the other sample e.g. in proby, largest value
# is 13.29, 2nd largest is 7.68 in prox, the largest value is 7.94, 2nd largest is 7.56 this becomes
# all that matters, all other information is discarded However, for KS, that is not the case, even
# though it takes a 'maximum' statistic as well The problem doesn't only happen for this kind of thing
# though... 

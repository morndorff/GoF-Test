# Variance of various test statistics under the null hypothesis
source("functions.R")
set.seed(5)


  testmat <- matrix(rnorm(5000*20),nrow=5000,ncol=20)
  ksvar <- apply(testmat, 1, ks.res.simp, y=pnorm, sd=2)
  energyvar <- apply(testmat, 1, wave.energy, y=pnorm, sd=2)
m_ks <- mean(ksvar)
m_en <- mean(energyvar)

sd_ks <- sd(ksvar)
sd_en <- sd(energyvar)

s_ks <- (ksvar-m_ks)/sd_ks
s_en <- (energyvar-m_en)/sd_en

hist(s_ks)
hist(s_en)

hist(ksvar, breaks=25, main="ooc ks")
hist(energyvar, breaks=25, main="ooc energy")


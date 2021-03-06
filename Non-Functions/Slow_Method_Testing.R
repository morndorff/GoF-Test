# "Fast Method" Testing

# November 18, 2014

# rm(list=ls())
source("functions.R")
library(waveslim)
set.seed(5)
Plots <- function (test, weight="Biased") {
  plot(test[[1]], main=paste("Plot of Test Statistic h(t), Weight=",weight), type="l") # Weight = Biased Towards Extreme Obs?
  hist(test[[7]], main="Tau Chosen (Raw)") # Which Tau is Chosen?
  hist(test[[8]], main="Tau Chosen (Chosen/Time)", breaks=20) # Tau Chosen Divided by Time (Scaled 0-1)
  d <- density(test[[8]], adjust=.5)
  plot(d, , "gaussian", xlim=c(0,1), main="Tau/T Chosen at Each Step") #Kernel Density Estimate
  quantile(test[[8]], probs=seq(0,1,.1)) # Quantiles of Tau/T
}
#############################################################
######### CURVE STUFF
# Finding ARLs for energy.curve method
UCL <- seq(.01,.06, length.out=10)
ARL_Curve_SA <- ARL_Proc(UCL=UCL, time=10, 
                         method=Find_IC_RL_Slow, 
                         tstat=wave.energy, 
                         dist="norm",
                         params=list(mean=0, sd=1))
ARL_Curve_SA_Mat <- ARL_Curve_SA[[1]] # Sample Average ARL
ARL <- as.numeric(ARL_Curve_SA_Mat[1,])
f <- approxfun(ARL,UCL)
UCL_SA_200_Curve_Slow <- f(200)


####### Plotting some Sample Paths
set.seed(5)
test <- Find_IC_RL_Slow(num.samp=32, dist="norm", 
                        params=list(mean=0, sd=1), 
                        tstat=wave.energy, 
                        UCL=.07)
###########CHANGE UCL AFTER SIMULATION!!!!


jpeg('ENCURVEIC_Slow2.jpg')
plot(test[[1]], main=paste("Plot of Test Statistic h(t), Slow Method IC N01"), , type="l")
dev.off()

test <- Find_CP_RL_Slow(num.samp=64, dist_one="norm", param_one=list(mean=0, sd=1),
                        dist_two = "norm", param_two = list(mean=0, sd=2),
                        tstat=wave.energy, UCL=.065, cp=1)
jpeg('ENCURVEOOC_Slow.jpg')
plot(test[[1]], main=paste("Plot of Test Statistic h(t) Slow Method IC N01, OOCN02 CP1=1"), type="l")
dev.off()
#########

# AVERAGE RUN LENGTH FOR OOC CONTROL FROM TIME 1

ARL_OOC_CURVE <- Find_ARL_OOC(num.samp=32, dist_one="norm", param_one=list(mean=0, sd=1),
                              dist_two="norm", param_two=list(mean=0, sd=2), cp=1,
                              tstat=wave.energy, UCL=UCL_SA_200_Curve, time=120, method=Find_CP_RL_Fast)
ARL_OOC_CURVE[[1]] # Yikes


###############################################################
######### KS STUFF
set.seed(5)
UCL <- seq(.02,.18, length.out=10)
ARL_Curve_SA <- ARL_Proc(UCL=UCL, time=120, 
                         method=Find_IC_RL_Fast, 
                         tstat=ks.res.simp, 
                         dist="norm",
                         params=list(mean=0, sd=1))
ARL_Curve_SA_Mat <- ARL_Curve_SA[[1]] # Sample Average ARL
ARL <- as.numeric(ARL_Curve_SA_Mat[1,])
f <- approxfun(ARL,UCL)
UCL_SA_200_KS <- f(200)
# Plotting some Sample Paths

set.seed(5)
test <- Find_IC_RL_Fast(num.samp=32, dist="norm", 
                        params=list(mean=0, sd=2), 
                        tstat=ks.res.simp, 
                        UCL=UCL_SA_200_KS,
                        detail=TRUE, weight=FALSE)
jpeg('KSIC.jpg')
plot(test[[1]], main=paste("Plot of Test Statistic h(t) KS IC N01"), , type="l")
dev.off()

test <- Find_CP_RL_Fast(num.samp=32, dist_one="norm", param_one=list(mean=0, sd=1),
                        dist_two = "norm", param_two = list(mean=0, sd=2),
                        tstat=ks.res.simp, 
                        UCL=UCL_SA_200_KS, cp=1,
                        detail=TRUE, weight=FALSE)
jpeg('KSOOC.jpg')
plot(test[[1]], main=paste("Plot of Test Statistic h(t) KS IC N01, OOCN02 "), type="l")
dev.off()

# AVERAGE RUN LENGTH FOR OOC CONTROL FROM TIME 1

ARL_OOC_KS <- Find_ARL_OOC(num.samp=32, dist_one="norm", param_one=list(mean=0, sd=1),
                           dist_two="norm", param_two=list(mean=0, sd=2), cp=1,
                           tstat=ks.res.simp, UCL=UCL_SA_200_KS, time=120, method=Find_CP_RL_Fast)
ARL_OOC_KS[[1]] # Yikes



###### Comparison
ARL_OOC_CURVE[[1]]
ARL_OOC_KS[[1]]
# Conclusions: Algorithm just seems bad. New Task: Test for slow algorithm.

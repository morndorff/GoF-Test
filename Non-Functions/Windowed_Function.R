# January 2, 2015

# Creating Windowed Version

source("functions.R")
test <- Find_IC_RL_Windowed(num.samp=64, dist="norm", 
                    params=list(mean=0, sd=1), 
                    tstat=wave.energy, UCL=.055, WSize=10, throw=TRUE, doplot=TRUE)
plot(test)
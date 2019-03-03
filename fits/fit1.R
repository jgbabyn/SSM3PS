##This is the first trial run...
##SPAM branch 92bbbf4

library(SSM3PS)
source("../data-raw/mdata3Ps.R")

##Careful with the KEYS!
keys <- keyGen(12,4)
keys$keyF <- 0:11
keys$keyQ[1,] <- 0:11
keys$keyQ[2,] <- 12:23
keys$keySD[1,] <- 0:11
keys$keySD[2,] <- 12:23
keys$keySD[3,] <- c(NA,24:34)


setd <- spamSetup(surveys,catch,landings,keys,stock_wt,midy_wt,mat,M=0.2,years=1983:2015,plusGroup=12)

##Watch as it abort traps?

spamFit(setd)

library(SSM3PS)
data(finMod)
modelA2 = finMod$fit

estX <- summary(modelA2$sdr,"random")
C <- solve(modelA2$obj$env$spHess(modelA2$obj$env$last.par.best,random=TRUE))
Xr <- MASS::mvrnorm(1,estX[,1],C)


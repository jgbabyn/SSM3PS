library(SSM3PS)
data(rDat)

dat <- datSetup(rDat$surveys,rDat$landings,rDat$fstock_wt,rDat$fcomm_wt,mat=rDat$fmat,
                age=2:16,years=1959:2016,plusGroup=14)

indi = dat$indices
crls = dat$crls
param <- dat$param
data <- dat$data

##Try with q at 6 plus all equal
qnames <- unique(na.omit(data.frame(iq=data$iq,qname=indi$qname,sur=indi$survey,age=indi$Age)))
qnames$age[qnames$age > 6] = 6
qmap = paste0(qnames$sur,qnames$age)
qmap = factor(qmap)
param$log_qparm = rep(0,length(unique(na.omit(data$iq))))

data$keyF = c(0,1,2,3,4,rep(5,13-5))
#data$keyF = c(0:4,5,5,5,5,5,5,5,5)
data$corflag = 4
data$corflagCRL = 4
param$log_std_logF = rep(0,length(unique(data$keyF)))
param$log_F = matrix(0,nrow=nrow(data$mat),ncol=length(param$log_std_logF))
param$rec_parm = numeric(0);
param$tRhoF = c(0,rep(0.1,5));
##param$param$tRhoF = c(0,0.0000000001,0.000000001,0.00000001)
param$tRhoCRL = c(0,rep(0.1,11))


mapA = list(
  tRhoF = factor(c(NA,0,1,2,3,4)),
  log_qparm = qmap,
  tRhoCRL = factor(c(NA,0,1,2,3,rep(4,12-5))),
  log_std_CRL = factor(c("m","m","m","m","m","m","m","m","m","m","m","m")),
  log_std_landings = factor(NA)

)

modelA1 <- spamFit(data,param,indi,crls,random=c("log_N","log_F"),sdrep=TRUE,map=mapA,control=list(iter.max=500,eval.max=500))
rmA1 <- residuals(modelA1,"oneStepGeneric")

qPA1 <- qPattern(modelA1$sdr,modelA1$indices)

cA1 <- oneStopPlot(rmA1,"catch")
RVOA1 <- oneStopPlot(rmA1,"RVOff")
RVIA1 <- oneStopPlot(rmA1,"RVInOff")
gA1 <- oneStopPlot(rmA1,"geac")
fA1 <- oneStopPlot(rmA1,"french")
lA1 <- oneStopPlot(rmA1,"landings")

cA1r <- oneStopPlot(rmA1,"catch",residuals=TRUE)
RVOA1r <- oneStopPlot(rmA1,"RVOff",residuals=TRUE)
RVIA1r <- oneStopPlot(rmA1,"RVInOff",residuals=TRUE)
gA1r <- oneStopPlot(rmA1,"geac",residuals=TRUE)
fA1r <- oneStopPlot(rmA1,"french",residuals=TRUE)
lA1r <- oneStopPlot(rmA1,"landings",residuals=TRUE)

cA1b <- oneStopPlot(rmA1,"catch",bubble=TRUE)
RVOA1b <- oneStopPlot(rmA1,"RVOff",bubble=TRUE)
RVIA1b <- oneStopPlot(rmA1,"RVInOff",bubble=TRUE)
gA1b <- oneStopPlot(rmA1,"geac",bubble=TRUE)
fA1b <- oneStopPlot(rmA1,"french",bubble=TRUE)
lA1b <- oneStopPlot(rmA1,"landings",bubble=TRUE)

r4pA <- resid4plot(rmA1)
r4pA$yplot
r4pA$cplot
r4pA$aplot
r4pA$eplot

ssb <- tsCIplot("log_ssb",modelA1$sdr,1959:2016,exp=TRUE)

topMod <- list(fit=modelA1,resid=rmA1)

usethis::use_data(topMod)

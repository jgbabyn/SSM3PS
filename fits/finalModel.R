library(SSM3PS)
data(rDat)

dat <- datSetup(rDat$surveys,rDat$landings,rDat$fstock_wt,rDat$fcomm_wt,mat=rDat$fmat,
                age=2:16,years=1959:2016,plusGroup=14)

indi = dat$indices
crls = dat$crls
param <- dat$param
data <- dat$data

##Try with q at 9 plus all equal
qnames <- unique(na.omit(data.frame(iq=data$iq,qname=indi$qname,sur=indi$survey,age=indi$Age)))
qnames$age[qnames$age > 9] = 9
qmap = paste0(qnames$sur,qnames$age)
qmap = factor(qmap)
param$log_qparm = rep(0,length(unique(na.omit(data$iq))))

data$keyF = c(0,1,2,3,4,5,6,7,8,8,8,rep(8,13-11))
data$corflag = 4
data$corflagCRL = 4
data$corflag
param$log_std_logF = rep(0,length(unique(data$keyF)))
param$log_F = matrix(0,nrow=nrow(data$mat),ncol=length(param$log_std_logF))
param$rec_parm = numeric(0);
param$tRhoF = c(0,rep(0.1,max(data$keyF)));
param$tRhoCRL = c(0,rep(0.1,11))


mapA = list(
  tRhoF = factor(c(NA,0,1,2,3,4,5,5,5)),
  log_qparm = qmap,
  tRhoCRL = factor(c(NA,0,1,2,3,4,5,5,5,rep(5,12-9))),
  log_std_CRL = factor(c("xs","s","s","s","m","m","m","l","l","l","l","l")),
  log_std_landings = factor(NA)

)

modelA2 <- spamFit(data,param,indi,crls,random=c("log_N","log_F"),sdrep=TRUE,map=mapA,control=list(iter.max=500,eval.max=500))
rmA2 <- residuals(modelA2,"oneStepGeneric")

ssb2 <- tsCIplot("log_ssb",modelA2$sdr,1959:2016,exp=TRUE)

qPA2 <- qPattern(modelA2$sdr,modelA2$indices)

cA2 <- oneStopPlot(rmA2,"catch")
RVOA2 <- oneStopPlot(rmA2,"RVOff")
RVIA2 <- oneStopPlot(rmA2,"RVInOff")
gA2 <- oneStopPlot(rmA2,"geac")
fA2 <- oneStopPlot(rmA2,"french")
lA2 <- oneStopPlot(rmA2,"landings")

cA2r <- oneStopPlot(rmA2,"catch",residuals=TRUE)
RVOA2r <- oneStopPlot(rmA2,"RVOff",residuals=TRUE)
RVIA2r <- oneStopPlot(rmA2,"RVInOff",residuals=TRUE)
gA2r <- oneStopPlot(rmA2,"geac",residuals=TRUE)
fA2r <- oneStopPlot(rmA2,"french",residuals=TRUE)
lA2r <- oneStopPlot(rmA2,"landings",residuals=TRUE)

cA2b <- oneStopPlot(rmA2,"catch",bubble=TRUE)
RVOA2b <- oneStopPlot(rmA2,"RVOff",bubble=TRUE)
RVIA2b <- oneStopPlot(rmA2,"RVInOff",bubble=TRUE)
gA2b <- oneStopPlot(rmA2,"geac",bubble=TRUE)
fA2b <- oneStopPlot(rmA2,"french",bubble=TRUE)
lA2b <- oneStopPlot(rmA2,"landings",bubble=TRUE)

r4pA2 <- resid4plot(rmA2)
r4pA2$yplot
r4pA2$cplot
r4pA2$aplot
r4pA2$eplot

finMod <- list(fit=modelA2,resid=rmA2)


usethis::use_data(finMod)

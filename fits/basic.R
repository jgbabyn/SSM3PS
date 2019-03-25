library(SSM3PS)
##Load in data read into R package
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
data$corflag = 3
data$corflagCRL = 0
param$log_std_logF = rep(0,length(unique(data$keyF)))
param$log_F = matrix(0,nrow=nrow(data$mat),ncol=length(param$log_std_logF))
param$rec_parm = numeric(0);
param$tRhoF = c(0.1)
##param$param$tRhoF = c(0,0.0000000001,0.000000001,0.00000001)
param$tRhoCRL = numeric(0)


mapL = list(
  #tRhoF = factor(c(NA,1,2,3,4,5)),  
    log_qparm = qmap,
    #tRhoCRL = factor(c(NA,rep(1,11))),
  log_std_CRL = factor(c("m","m","m","m","m","m","m","m","m","m","m","m")),
  log_std_landings = factor(NA)

)

fitt <- spamFit(data,param,indi,crls,random=c("log_N","log_F"),map=mapL,control=list(iter.max=500,eval.max=500))

qPat <- qPattern(fitt$sdr,fitt$indices)

landingsP <- tsCIplot("log_landings_pred",fitt$sdr,1959:2016,exp=FALSE)

resi <- residuals(fitt)
library(ggplot2)

lland <- rDat$landings
lland$landings = lland$landings/1000
lland[,-1] = log(lland[,-1])
landB <- oneStopPlot(resi,"landings")



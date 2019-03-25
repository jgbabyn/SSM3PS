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
data$corflag = 0
data$corflagCRL = 0
param$log_std_logF = rep(0,length(unique(data$keyF)))
param$log_F = matrix(0,nrow=nrow(data$mat),ncol=length(param$log_std_logF))
param$rec_parm = numeric(0);
param$tRhoF = numeric(0);
##param$param$tRhoF = c(0,0.0000000001,0.000000001,0.00000001)
param$tRhoCRL = numeric(0)


mapL = list(
  #tRhoF = factor(c(NA,1,2,3,4,5)),  
    log_qparm = qmap,
    #tRhoCRL = factor(c(NA,rep(1,11))),
  log_std_CRL = factor(c("m","m","m","m","m","m","m","m","m","m","m","m")),
  log_std_landings = factor(NA)

)

qqq <- system.time(fitt <- spamFit(data,param,indi,crls,random=c("log_N","log_F"),map=mapL,control=list(iter.max=500,eval.max=500)))

qPat <- qPattern(fitt$sdr,fitt$indices)

landingsP <- tsCIplot("log_landings_pred",fitt$sdr,1959:2016,exp=FALSE)

resi <- residuals(fitt)
library(ggplot2)

lland <- rDat$landings
lland$landings = lland$landings/1000
lland[,-1] = log(lland[,-1])
landB <- oneStopPlot(resi,"landings")


modelOpts <- expand.grid(cfF=0:4,cfCRL=0:4,rf=0:3)
modelOpts$time <- NA
modelOpts$AIC <- NA
modelOpts$BIC <- NA
modelOpts$logLik <- NA
modelOpts$message <- NA

for(i in 75:nrow(modelOpts)){
    tempD <- data
    tempP <- param
    tempM <- mapL

    cfF = modelOpts$cfF[i]
    cfCRL = modelOpts$cfCRL[i]
    rf = modelOpts$rf[i]
    
    tempD$corflag = cfF
    tempD$corflagCRL = cfCRL
    tempD$recflag = rf

    ##Parmeters for F correlation
    if(cfF == 0 | cfF == 1){
        tempP$tRhoF = numeric(0)
    }
    if(cfF == 2 | cfF == 3){
        tempP$tRhoF = 0.1
    }
    if(cfF == 4){
        tempP$tRhoF = c(0,rep(0.1,5))
        tempM$tRhoF = factor(c(NA,0,1,2,3,4))
    }

    ##CRL corr params
    if(cfCRL == 0 | cfCRL == 1){
        tempP$tRhoCRL = numeric(0)
    }
    if(cfCRL == 2 | cfCRL == 3){
        tempP$tRhoCRL = 0.1
    }
    if(cfCRL == 4){
        tempP$tRhoCRL = c(0,rep(0.1,11))
        tempM$tRhoCRL = factor(c(NA,0,1,2,3,rep(4,12-5)))
    }

    ##rec param
    if(rf != 0){
        tempP$rec_parm = c(0,0)
    }
    
    time <- system.time(tempF <- spamFit(tempD,tempP,indi,crls,random=c("log_N","log_F"),sdrep=FALSE,map=tempM,control=list(iter.max=500,eval.max=500)))

    modelOpts$time[i] = time[3]
    modelOpts$AIC[i] = AIC(tempF)
    modelOpts$BIC[i] = BIC(tempF)
    modelOpts$logLik[i] = logLik(tempF)
    modelOpts$message[i] = tempF$opt$message
    print(i)
}

RCmo <- modelOpts[modelOpts$message == "relative convergence (4)",]
RCmo <- dplyr::arrange(RCmo,AIC,BIC)
FCmo <- modelOpts[modelOpts$message == "false convergence (8)",]
SCmo <- modelOpts[modelOpts$message == "singular convergence (7)",]

usethis::use_data(modelOpts)

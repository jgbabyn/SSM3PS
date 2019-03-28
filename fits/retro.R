library(SSM3PS)
data(rDat)

years <- 1959:2016

##perform the retro
retro <- list()
for(i in 2007:2016){
    j = i - 2006

    dat <- datSetup(rDat$surveys,rDat$landings,rDat$fstock_wt,rDat$fcomm_wt,mat=rDat$fmat,
                    age=2:16,years=1959:i,plusGroup=14)

    
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

    mod <- spamFit(data,param,indi,crls,random=c("log_N","log_F"),sdrep=TRUE,map=mapA,control=list(iter.max=500,eval.max=500))
    retro[[j]] <- mod
}

##Get recruitment from sdr, sdr does column wise matrix unpacking
sum10 <- summary(retro[[10]]$sdr)
ind = which(rownames(sum10) == "log_N")
sum10 = sum10[ind[1:58],]
trec = list()
trec$value = sum10[,1]
trec$sd = sum10[,2]

##Make plots from the retro
ssbRet <- tsCIplot("log_ssb",retro[[10]]$sdr,1959:2016,exp=TRUE)
bioRet <- tsCIplot("log_biomass",retro[[10]]$sdr,1959:2016,exp=TRUE)
aveF46 <- tsCIplot("log_aveFXY",retro[[10]]$sdr,1959:2016,exp=TRUE)
recRet <- tsCIplot("log_N",trec,1959:2016,exp=TRUE)
for(i in 2007:2015){
    j = i - 2006
    colr = sample(colors(),1)

    summ <- summary(retro[[j]]$sdr)
    ind = which(rownames(summ) == "log_N")
    summ = summ[ind[1:(48+j)],]
    ttrec = list()
    ttrec$value = summ[,1]
    ttrec$sd = summ[,2]

    tSSB <- tsCIplot("log_ssb",retro[[j]]$sdr,1959:i,exp=TRUE,add=TRUE,color=colr,fcolor=colr)
    tBio <- tsCIplot("log_biomass",retro[[j]]$sdr,1959:i,exp=TRUE,add=TRUE,color=colr,fcolor=colr)
    tF46 <- tsCIplot("log_aveFXY",retro[[j]]$sdr,1959:i,exp=TRUE,add=TRUE,color=colr,fcolor=colr)
    tRec <- tsCIplot("log_N",ttrec,1959:i,exp=TRUE,add=TRUE,color=colr,fcolor=colr)
    
    ssbRet = ssbRet + tSSB$line
    bioRet = bioRet + tBio$line
    aveF46 = aveF46 + tF46$line
    recRet = recRet + tRec$line
}

retL <- list(retro=retro,ssbRet=ssbRet,bioRet=bioRet,aveRet=aveF46,recRet=recRet)

usethis::use_data(retL)

##Leave one out fleet runs

##Leave out GEAC
surveys = rDat$surveys
surveys$geac = NULL
dat <- datSetup(surveys,rDat$landings,rDat$fstock_wt,rDat$fcomm_wt,mat=rDat$fmat,
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

   ngeac <- spamFit(data,param,indi,crls,random=c("log_N","log_F"),sdrep=TRUE,map=mapA,control=list(iter.max=500,eval.max=500))

##Leave out french

surveys = rDat$surveys
surveys$french = NULL
dat <- datSetup(surveys,rDat$landings,rDat$fstock_wt,rDat$fcomm_wt,mat=rDat$fmat,
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

   nfrench <- spamFit(data,param,indi,crls,random=c("log_N","log_F"),sdrep=TRUE,map=mapA,control=list(iter.max=500,eval.max=500))

##Leave out RVOff
surveys = rDat$surveys
surveys$RVOff = NULL
dat <- datSetup(surveys,rDat$landings,rDat$fstock_wt,rDat$fcomm_wt,mat=rDat$fmat,
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

   nrvoff <- spamFit(data,param,indi,crls,random=c("log_N","log_F"),sdrep=TRUE,map=mapA,control=list(iter.max=500,eval.max=500))

##Leave out RVInOff
surveys = rDat$surveys
surveys$RVInOff = NULL
dat <- datSetup(surveys,rDat$landings,rDat$fstock_wt,rDat$fcomm_wt,mat=rDat$fmat,
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

   nrvinoff <- spamFit(data,param,indi,crls,random=c("log_N","log_F"),sdrep=TRUE,map=mapA,control=list(iter.max=500,eval.max=500))

cva <- list(ngeac,nfrench,nrvoff,nrvinoff)

ssblfo <- tsCIplot("log_ssb",retro[[10]]$sdr,1959:2016,exp=TRUE)
biolfo <- tsCIplot("log_biomass",retro[[10]]$sdr,1959:2016,exp=TRUE)
avelfo <- tsCIplot("log_aveFXY",retro[[10]]$sdr,1959:2016,exp=TRUE)
reclfo <- tsCIplot("log_N",trec,1959:2016,exp=TRUE)
colls <- c("green","red","blue","orange")

for(i in 1:4){
    colr = colls[i]

    summ <- summary(cva[[i]]$sdr)
    ind = which(rownames(summ) == "log_N")
    summ = summ[ind[1:58],]
    ttrec = list()
    ttrec$value = summ[,1]
    ttrec$sd = summ[,2]

    tSSB <- tsCIplot("log_ssb",cva[[i]]$sdr,1959:2016,exp=TRUE,add=TRUE,color=colr,fcolor=colr)
    tBio <- tsCIplot("log_biomass",cva[[i]]$sdr,1959:2016,exp=TRUE,add=TRUE,color=colr,fcolor=colr)
    tF46 <- tsCIplot("log_aveFXY",cva[[i]]$sdr,1959:2016,exp=TRUE,add=TRUE,color=colr,fcolor=colr)
    tRec <- tsCIplot("log_N",ttrec,1959:2016,exp=TRUE,add=TRUE,color=colr,fcolor=colr)
    
    ssblfo = ssblfo + tSSB$line
    biolfo = biolfo + tBio$line
    avelfo = avelfo + tF46$line
    reclfo = reclfo + tRec$line

}

lfo = list(ssblfo=ssblfo,biolfo=biolfo,avelfo=avelfo,reclfo=reclfo)
usethis::use_data(lfo,overwrite=TRUE)

##Bounds testing
landMed <- as.data.frame(rDat$landings)
landMed[,2] <- landMed[,2]/1000
landMed$low <- landMed[,2]*landMed[,4]
landMed$up <- landMed[,2]*landMed[,3]
landMed$pred <- retro[[10]]$rep$landings_pred
library(ggplot2)

landBoundP <- ggplot(landMed,aes(Year,landings)) + geom_line(aes(Year,landings)) + geom_ribbon(aes(ymin=low,ymax=up),fill="green",alpha=0.3)
landBoundP <- landBoundP + geom_line(aes(Year,pred),data=landMed,color="green")

landTi <- as.data.frame(rDat$landings)
landTi$UB[landTi$Year %in% 1959:1993] = 1.25
landTi$UB[landTi$Year %in% 1994:1996] = 1.05
landTi$UB[landTi$Year %in% 1997:1999] = 1.50
landTi$UB[landTi$Year %in% 2000:2016] = 1.15

landTi$LB[landTi$Year %in% 1959:1993] = 0.75
landTi$LB[landTi$Year %in% 1994:1996] = 0.95
landTi$LB[landTi$Year %in% 1997:1999] = 0.5
landTi$LB[landTi$Year %in% 2000:2016] = 0.85


landWi <- as.data.frame(rDat$landings)
landWi$UB[landWi$Year %in% 1959:1993] = 1.75
landWi$UB[landWi$Year %in% 1994:1996] = 1.20
landWi$UB[landWi$Year %in% 1997:1999] = 2.0
landWi$UB[landWi$Year %in% 2000:2016] = 1.5

landWi$LB[landWi$Year %in% 1959:1993] = 0.35
landWi$LB[landWi$Year %in% 1994:1996] = 0.80
landWi$LB[landWi$Year %in% 1997:1999] = 0.25
landWi$LB[landWi$Year %in% 2000:2016] = 0.5


dat <- datSetup(rDat$surveys,landTi,rDat$fstock_wt,rDat$fcomm_wt,mat=rDat$fmat,
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

tibound <- spamFit(data,param,indi,crls,random=c("log_N","log_F"),sdrep=TRUE,map=mapA,control=list(iter.max=500,eval.max=500))


dat <- datSetup(rDat$surveys,landWi,rDat$fstock_wt,rDat$fcomm_wt,mat=rDat$fmat,
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

wibound <- spamFit(data,param,indi,crls,random=c("log_N","log_F"),sdrep=TRUE,map=mapA,control=list(iter.max=500,eval.max=500))

landTi[,2] <- landTi[,2]/1000
landTi$low <- landTi[,2]*landTi[,4]
landTi$up <- landTi[,2]*landTi[,3]
landTi$pred <- tibound$rep$landings_pred

landWi[,2] <- landWi[,2]/1000
landWi$low <- landWi[,2]*landWi[,4]
landWi$up <- landWi[,2]*landWi[,3]
landWi$pred <- wibound$rep$landings_pred



landBoundP <- ggplot(landTi,aes(Year,landings)) + geom_line(aes(Year,landings)) + geom_ribbon(aes(ymin=low,ymax=up),fill="blue",alpha=0.3)
landBoundP <- landBoundP + geom_line(aes(Year,pred),data=landTi,color="blue")


landBoundP <- landBoundP + geom_line(aes(Year,landings),landMed) + geom_ribbon(aes(ymin=low,ymax=up),landMed,fill="green",alpha=0.3)
landBoundP <- landBoundP + geom_line(aes(Year,pred),data=landMed,color="green")

landBoundP <- landBoundP + geom_line(aes(Year,landings),landWi) + geom_ribbon(aes(ymin=low,ymax=up),landWi,fill="red",alpha=0.2)
landBoundP <- landBoundP + geom_line(aes(Year,pred),data=landWi,color="red") + ylab("Landings (Kt)")

usethis::use_data(landBoundP,overwrite=TRUE)

mboundSSB <- tsCIplot("log_ssb",retro[[10]]$sdr,1959:2016,exp=TRUE,color="green",fcolor="green")
boundSSB <- tsCIplot("log_ssb",tibound$sdr,1959:2016,exp=TRUE,color="blue",fcolor="blue",add=TRUE)
wboundSSB <- tsCIplot("log_ssb",wibound$sdr,1959:2016,exp=TRUE,color="red",fcolor="red",add=TRUE)
mboundSSB = mboundSSB + boundSSB$line + wboundSSB$line + ylab("SSB (Kt)")

usethis::use_data(mboundSSB)

##Weights test
stkwt = matrix(NA,nrow=58,ncol=16)
rownames(stkwt) = 1959:2016
ms = apply(rDat$stock_wt[rDat$stock_wt$Year %in% 1983:1987,-1],2,median)
stkwt[as.character(1959:1982),] = rep(ms,each=length(1959:1982))
stkwt[as.character(1983:2016),] = as.matrix(rDat$stock_wt[,-1])
stkwt = cbind(Year=1959:2016,as.data.frame(stkwt))
names(stkwt) = c("Year",paste0("Age",1:16))

midwt = matrix(NA,nrow=58,ncol=16)
rownames(midwt) = 1959:2016
ms = apply(rDat$midy_wt[rDat$midy_wt$Year %in% 1983:1987,-1],2,median)
midwt[as.character(1959:1982),] = rep(ms,each=length(1959:1982))
midwt[as.character(1983:2016),] = as.matrix(rDat$midy_wt[,-1])
midwt = cbind(Year=1959:2016,as.data.frame(midwt))
names(midwt) = c("Year",paste0("Age",1:16))


dat <- datSetup(rDat$surveys,rDat$landings,stkwt,midwt,mat=rDat$mat,
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

fwts <- spamFit(data,param,indi,crls,random=c("log_N","log_F"),sdrep=TRUE,map=mapA,control=list(iter.max=500,eval.max=500))


dat <- datSetup(rDat$surveys,rDat$landings,stkwt,midwt,mat=rDat$fmat,
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

fwts2 <- spamFit(data,param,indi,crls,random=c("log_N","log_F"),sdrep=TRUE,map=mapA,control=list(iter.max=500,eval.max=500))


##Eyeballed points from 2017 3Ps relSSB resdoc
eyepts <- c(1.5,1.55,1.7,1.65,1.6,1.7,1.95,1.8,1.45,1.25,1.05,1,1.4,1.45,1.25,1.3,1.5,1.6,1.65,1.95,2.25,2.26,1.95,1.4,1.25,1.05,0.95,1.1,1.6,1.75,1.4,1.1,1.2)
yearE = c(1983:2005,2007:2016)
relpts = data.frame(eyepts,Year=yearE)
plot(yearE,eyepts,ylim=c(0,3.5),type="l")
abline(h=1,lty=2)
abline(h=2,lty=2)

relSSB <- tsCIplot("log_ssb",retro[[10]]$sdr,1959:2016,relYear=1994)
relSSBo <- tsCIplot("log_ssb",fwts$sdr,1959:2016,relYear=1994,color="red",fcolor="red",add=TRUE)
relSSBo2 <- tsCIplot("log_ssb",fwts2$sdr,1959:2016,relYear=1994,color="blue",fcolor="blue",add=TRUE)
relSSB <- relSSB + relSSBo$line + relSSBo2$line
relSSB <- relSSB + geom_line(aes(Year,eyepts,color="orange"),relpts) + geom_point(aes(Year,eyepts,color="orange"),relpts) +geom_hline(yintercept=1,linetype="dashed")
relSSB <- relSSB + geom_hline(yintercept=2,linetype="dashed") + ylab("Relative SSB to 1994") + guides(fill=FALSE)
relSSB <- relSSB + scale_color_manual(label=c("FPCA wts. + FPCA mat.","Provided Wts. + FPCA mat.","DFO SURBA","Provided Wts. + Mat."),
                                      values=c("red","blue","black","orange"))

usethis::use_data(relSSB,overwrite=TRUE)

absSSB <- tsCIplot("log_ssb",retro[[10]]$sdr,1959:2016,exp=TRUE)
absSSBo <- tsCIplot("log_ssb",fwts$sdr,1959:2016,exp=TRUE,color="red",fcolor="red",add=TRUE)
absSSBo2 <- tsCIplot("log_ssb",fwts2$sdr,1959:2016,exp=TRUE,color="blue",fcolor="blue",add=TRUE)
absSSB <- absSSB + absSSBo$line + absSSBo2$line + ylab("SSB (Kt)") + guides(fill=FALSE)
absSSB <- absSSB + scale_color_manual(label=c("FPCA wts. + FPCA mat.","Provided Wts. + FPCA mat.","Provided Wts. + Mat."),
                                      values=c("red","blue","orange"))


usethis::use_data(absSSB)

fwgtss <- list(f1=fwts,f2=fwts2)
usethis::use_data(fwgtss)

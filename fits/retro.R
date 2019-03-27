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

##Jittering the parameters

final <- retro[[10]]$opt$par
obj <- retro[[10]]$obj
jittRun <- list()
for(i in 1:100){
    obj$par <- runif(63,-1,1)
    names(obj$par) = names(final)
    run <- nlminb(obj$par,obj$fn,obj$gr)
    jittRun[[i]] = all.equal(final,run$par)
}

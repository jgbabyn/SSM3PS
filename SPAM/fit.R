
library(TMB)

setwd("~/Dropbox/workingOk")
#compile("fit.cpp","-O0 -g")
compile("fit.cpp")
##dyn.unload(dynlib("fit.cpp"))
dyn.load(dynlib("fit"))

library(SSM3PS)
source("~/Dropbox/workingOk/SSM3PS/data-raw/mdata3Ps.R")

#surveys$catch = catch
surveys = list(RVOff = surveys$RVOff,catch = surveys$catch)
#surveys$catch[,-1] = surveys$catch[,-1]/1000
dat <- datSetup(surveys,landings,fstock_wt,fcomm_wt,mat=fmat,
                age=3:16,years=1983:2016,plusGroup=14)

indd <- dat$indices
##param <- paramSetup(dat)
param <- dat$param
## param$param$log_std_CRL = rep(0,ncol(dat$data$mat)-1)
## param$param$log_std_landings = log(0.02)
## param$param$log_sdS = 0
## param$param$log_N = matrix(0,nrow=nrow(dat$data$mat),ncol=ncol(dat$data$mat))

dat <- dat$data

#dat$lowerMult = rep(0.75,nrow(indd))
#dat$upperMult = rep(1.5,nrow(indd))
#dat$use_cb = NULL
dat$keyF = c(0,1,rep(2,13-2))
#dat$recflag = 0
#dat$recflag = 0
#dat$corflag = 0
#dat$corflagCRL = 0
#dat$gammaSq = 0.1
## dat$crlVec = dat$log_index[1:374]
## dat$log_index = dat$log_index[-c(1:374)]
## dat$i_zero = dat$i_zero[-c(1:374)]
## dat$iyear = dat$iyear[-c(1:374)]
## dat$iage = dat$iage[-c(1:374)]
## dat$isurvey = dat$isurvey[-c(1:374)]
## dat$iq = dat$iq[-c(1:374)]
## dat$fs = dat$fs[-c(1:374)]
## dat$ft = dat$ft[-c(1:374)]
#dat$aveFages = c(1,3)
qparm_name = levels(as.factor(dat$iq))
param$param$log_std_logF = rep(0,3)
param$param$log_F = matrix(0,nrow=nrow(dat$mat),ncol=3)
param$param$rec_parm = numeric(0);
param$param$tRhoF = numeric(0)
param$param$tRhoCRL = numeric(0)

                                 
##check sorting?

##Map for fitting landings right now
mapL = list( 
  log_qparm = factor(c(0,1,2,rep(3,10))),
  log_std_CRL = factor(c("4+","4+","4+","4+","4+","4+","4+","4+","4+","4+","4+")),
  log_std_cye=factor(NA),  
  logit_ar_cye_year = factor(NA),  
  cye=factor(matrix(NA,nrow=nrow(dat$mat),ncol=ncol(dat$mat))),
  log_std_landings = factor(NA)

)

dat$index_censor = 2
obj <- MakeADFun(dat,param$param,random=c("log_N","log_F"),DLL="fit",map=mapL)



ff <- obj$gr(obj$par)
names(ff) <- names(obj$par)


opt <- nlminb(obj$par,obj$fn,obj$gr,control=list(iter.max=500,eval.max=500))

tein = which(dat$i_zero == 1)
osp <- oneStepPredict(obj,"log_index","keep",discrete=FALSE,"oneStepGeneric")

ospCRL <- oneStepPredict(obj,"crlVec","keepCRL",discrete=FALSE)
rep2 <- obj$report()
sd.rep2 <- sdreport(obj)

##Trying to solve this problem
ind2 <- indd
ind2$log_index = log(ind2$index)
ind2 = cbind(ei=rep2$Elog_index,ind2)

cind = ind2[ind2$survey == "catch",c("ei","log_index")]
cind$eei = exp(cind$ei)
cind$index = exp(cind$log_index)/1000
cind$diff = cind$eei -cind$index

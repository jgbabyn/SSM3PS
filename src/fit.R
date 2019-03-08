library(TMB)

setwd("~/Documents/GitHub/SSM3PS/src")
compile("fit.cpp","-O0 -g")
#dyn.unload(dynlib("fit.cpp"))
dyn.load(dynlib("fit"))

library(SSM3PS)
source("../data-raw/mdata3PS.R")

#dat <- datSetup(surveys,catch$catch,landings,stock_wt,
#                midy_wt,mat,ages=3:16,years=1983:2015,plusGroup=12,match3NOdims=TRUE)

dat <- datSetup(surveys[-2],catch$catch,landings,stock_wt,midy_wt,mat,
                age=2:16,years=1983:2015,plusGroup=14,naz.rm=TRUE)


param <- paramSetup(dat)

qparm_name = levels(as.factor(dat$iq))

mapq = qparm_name
agep = str_pad(6:14 , 2, pad = "0")

ind = mapq %in% c(paste("OFF:",agep,sep=''),paste("IO:",agep,sep=''))

mapq[ind] = paste("OFF:",rep('6+',9),sep='')


mapN = list( 
  log_qparm = factor(mapq),
  log_std_logF = factor(c("2","3","4+","4+","4+","4+","4+","4+","4+","4+","4+","4+","4+")),
  log_std_log_C = factor(c("2","3","4+","4+","4+","4+","4+","4+","4+","4+","4+","4+","4+")),          
  logit_ar_logF_age = factor(NA),              
  logit_ar_logF_year = factor(NA), 
  log_std_pe=factor(NA),  
  logit_ar_pe_year = factor(NA),     
  logit_ar_pe_age = factor(NA),
  pe=factor(matrix(NA,nrow=nrow(dat$mat),ncol=ncol(dat$mat))),
  log_std_cye=factor(NA),  
  logit_ar_cye_year = factor(NA),  
  cye=factor(matrix(NA,nrow=nrow(dat$C),ncol=ncol(dat$C))) 
)

dat$index_censor = 2
dat$catch_censor = 2
obj <- MakeADFun(dat,param$param,random=c("log_Rec_dev","log_F","pe"),DLL="fit",map=mapN)

obj$gr(obj$par)

opt <- nlminb(obj$par,obj$fn,obj$gr,control=list(iter.max=500,eval.max=500))

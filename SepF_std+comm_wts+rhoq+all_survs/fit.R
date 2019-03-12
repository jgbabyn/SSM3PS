
setwd("C:\\Users\\mroberts\\OneDrive - Memorial University of Newfoundland\\Marine Institute\\Classes\\FISH_6005\\Project 1\\Project 1\\Working model\\SepF_std+comm_wts+rhoq+all_survs")

library(TMB)

rm(list = ls())
dyn.unload("fit");

load("tmb.RData")

str(tmb.data)
#trying to match dimensions based on 3NO SSM
#tmb.data$C<-tmb.data$C[1:57,]
#tmb.data$log_C<-tmb.data$log_C[1:57,]
#tmb.data$C_zero<-tmb.data$C_zero[1:57,]
#tmb.data$landings<-tmb.data$landings[1:57]
#tmb.data$log_landings<-tmb.data$log_landings[1:57]

setparm = function(opt,rep){
  
  parm = parameters
  uname = names(parameters)
  names.par = names(opt$par)
  for(i in 1:length(uname)){
    ind = names.par==uname[i]
    parm[[i]] = opt$par[ind]
  }         
  parm$log_Rec_dev = rep$log_Rec_dev 
  parm$log_F = rep$log_F
  parm$pe = rep$pe     
  parm$cye = rep$cye               
  return(parm)
}

library(TMB) 
library(stringr)

compile("fit.cpp")  

dyn.load("fit") 
#dyn.unload("fit");

qparm_name = levels(as.factor(indices$qname))

parameters <- list( 
  log_No=rep(log(15),tmb.data$A-1),
  log_Rec_mean=5,      
  log_std_log_R=-1, 
  log_qparm=rep(0,length(unique(tmb.data$iq))),
  log_rho=rep(0,length(unique(tmb.data$ia))),
  log_std_index=rep(log(0.3),length(unique(tmb.data$isurvey))), 
  log_std_rho=rep(log(0.3),length(unique(tmb.data$ia))), 
  #log_std_ratio=rep(log(0.3),length(unique(tmb.data$isurvey))), 
  log_std_logF = rep(log(0.1),tmb.data$A),
  log_std_logF_new = rep(log(0.1),tmb.data$A),
  log_std_pe = log(0.1),      
  log_std_cye = log(0.5),                 
  logit_ar_logF_year = log(0.99/0.01), 
  logit_ar_logF_age = -10,            
  logit_ar_pe_year = -10, 
  logit_ar_pe_age = -10, 
  logit_ar_cye_year = 0, 
  log_std_log_C = rep(log(0.3),tmb.data$A), 
  log_Rec_dev=rep(5,tmb.data$Y),     
  log_F=matrix(log(0.3),nrow=tmb.data$Y,ncol=tmb.data$A),
  pe=matrix(0,nrow=tmb.data$Y,ncol=tmb.data$A),
  cye=matrix(0,nrow=tmb.data$Y,ncol=tmb.data$A)
)

parameters.L <- list( 
  log_No=rep(-Inf,tmb.data$A-1),
  log_Rec_mean=0,    
  log_std_log_R=-5,   
  log_qparm=rep(-Inf,length(unique(tmb.data$iq))),
  log_rho=rep(-Inf,length(unique(tmb.data$ia))),
  log_std_index=rep(log(0.01),length(unique(tmb.data$isurvey))),
  log_std_rho=rep(log(0.01),length(unique(tmb.data$ia))),
  #log_std_ratio=rep(log(0.01),length(unique(tmb.data$isurvey))),
  log_std_logF = rep(-Inf,tmb.data$A), 
  log_std_logF_new = rep(-Inf,tmb.data$A), 
  log_std_pe = -Inf,      
  log_std_cye = -Inf,                 
  logit_ar_logF_year = -Inf, 
  logit_ar_logF_age = -Inf,            
  logit_ar_pe_year = -Inf, 
  logit_ar_pe_age = -Inf, 
  logit_ar_cye_year = -Inf, 
  log_std_log_C = rep(-Inf,tmb.data$A)          
)            

parameters.U <- list(
  log_No=rep(log(4000),tmb.data$A-1),
  log_Rec_mean=50,    
  log_std_log_R=5,    
  log_qparm=rep(Inf,length(unique(tmb.data$iq))),
  log_rho=rep(Inf,length(unique(tmb.data$ia))),
  log_std_index=rep(log(2),length(unique(tmb.data$isurvey))), 
  log_std_rho=rep(log(2),length(unique(tmb.data$ia))),
  #log_std_ratio=rep(log(2),length(unique(tmb.data$isurvey))), 
  log_std_logF = rep(Inf,tmb.data$A),  
  log_std_logF_new = rep(-Inf,tmb.data$A), 
  log_std_pe = log(0.5),      
  log_std_cye = Inf,                 
  logit_ar_logF_year = log(0.99/0.01), 
  logit_ar_logF_age = log(0.99/0.01),            
  logit_ar_pe_year = log(0.99/0.01), 
  logit_ar_pe_age = log(0.99/0.01), 
  logit_ar_cye_year = log(0.99/0.01), 
  log_std_log_C = rep(Inf,tmb.data$A)      
)

lower = unlist(parameters.L);
upper = unlist(parameters.U);

mapq = qparm_name

age1 = str_pad(4:14 , 2, pad = "0")
age2 = str_pad(4:14 , 2, pad = "0")
age3 = str_pad(4:14 , 2, pad = "0")

ind = mapq %in% c(paste("BOTH:",age1,sep=''),paste("FRANCE:",age2,sep=''),paste("GEAC:",age3,sep=''))

mapq[ind] = c(paste("BOTH:",rep('6+',11),sep=''), paste("FRANCE:",rep('6+',11),sep=''), paste("GEAC:",rep('6+',11),sep=''))


map = list( 
  log_qparm = factor(mapq),
  #log_std_rho=factor(unique(tmb.data$ia)),
  log_std_logF = factor(c("1", "2","3","4+","4+","4+","4+","4+","4+","4+","4+","4+","4+")),
  log_std_logF_new = factor(c("1","2","3","4+","4+","4+","4+","4+","4+","4+","4+","4+","4+")),
  log_std_log_C = factor(c("1","2","3","4+","4+","4+","4+","4+","4+","4+","4+","4+","4+")),  
  logit_ar_logF_age = factor(NA),              
  logit_ar_logF_year = factor(NA), 
  log_std_pe=factor(NA),  
  logit_ar_pe_year = factor(NA),     
  logit_ar_pe_age = factor(NA),
  pe=factor(matrix(NA,nrow=tmb.data$Y,ncol=tmb.data$A)),
  log_std_cye=factor(NA),  
  logit_ar_cye_year = factor(NA),  
  cye=factor(matrix(NA,nrow=tmb.data$Y,ncol=tmb.data$A)) 
)

tp=parameters.L;  
tp$log_std_logF=rep(-Inf,4);
tp$log_std_logF_new=rep(-Inf,4);
tp$log_std_log_C=rep(-Inf,4); 
tp$log_qparm=rep(-Inf,length(unique(mapq))); 
#tp$log_std_rho=rep(-Inf, length(unique(tmb.data$ia)));
tp$log_std_pe=NULL;   
tp$logit_ar_logF_year = NULL 
tp$logit_ar_logF_age = NULL  
tp$logit_ar_pe_year = NULL 
tp$logit_ar_pe_age = NULL  
tp$log_std_cye=NULL;
tp$logit_ar_cye_year = NULL
lower = unlist(tp); 

tp=parameters.U;      
tp$log_std_logF=rep(Inf,4); 
tp$log_std_logF_new=rep(Inf,4); 
tp$log_std_log_C=rep(Inf,4); 
tp$log_qparm=rep(Inf,length(unique(mapq)));
#tp$log_std_rho=rep(Inf, length(unique(tmb.data$ia)));
tp$log_std_pe=NULL;     
tp$logit_ar_logF_year = NULL 
tp$logit_ar_logF_age = NULL  
tp$logit_ar_pe_year = NULL 
tp$logit_ar_pe_age = NULL   
tp$log_std_cye=NULL;
tp$logit_ar_cye_year = NULL
upper = unlist(tp);

## first produce starting values by fitting to upper bounds;

tmb.data$index_censor = 2; 
tmb.data$catch_censor = 2; 
tmb.data$use_pe = 0;       
tmb.data$use_cye = 0; 


obj <- MakeADFun(tmb.data,parameters,map=map,
                 random=c("log_Rec_dev","log_F"), 
                 DLL="fit",inner.control=list(maxit=500,trace=F))  



sp=obj$par

length(lower)
length(upper)
length(sp)

obj$gr(sp)  

opt<-nlminb(sp,obj$fn,obj$gr,upper=upper,lower=lower,
                  control = list(trace=10,eval.max=2000,iter.max=1000))

obj$gr(opt$par) 

rep = obj$report()

#rep.start=rep

sd.rep<-sdreport(obj)

save.image(file='fit.RData')

#names(obj$env)
#sp=obj$env$last.par[1:length(obj$par)]



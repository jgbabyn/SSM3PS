## explore saveral model for L@A 3PS cod
### kincrement - Noel Cadigan
##try to improve the fit

rm(list=ls())
library(tidyverse);library(plyr);library(reshape2); library(dplyr);library(ggplot2)
library(lattice);library(stats4);library(MASS);library(car);library(TMB);library("xtable") 
library(latticeExtra)


source("C:/Users/ChampagnatJ/Dropbox/IFREMER/3Ps_cod/Cours/Cours_Noel_Cadigan/F6005_winter_2019/3Pscod/growth_new/growth/resid.txt")
source("C:/Users/ChampagnatJ/Dropbox/IFREMER/3Ps_cod/Cours/Cours_Noel_Cadigan/F6005_winter_2019/3Pscod/growth_new/growth/pm.txt")

setwd("C:/Users/ChampagnatJ/Dropbox/IFREMER/3Ps_cod/Cours/Cours_Noel_Cadigan/F6005_winter_2019/3Pscod/growth_new/growth/k_increments")

len <- read.delim("..//data//lengths.txt",header=T)

survey.dates <- read.delim("..//data//dates.txt",header=T)
start = as.Date(as.character(survey.dates$Start.Date), "%d-%b-%y")
end = as.Date(as.character(survey.dates$End.Date), "%d-%b-%y")

start.doy=as.numeric(strftime(start, format = "%j"))
end.doy=as.numeric(strftime(end, format = "%j"))

sf = ((start.doy+end.doy)/2)/365



year=len$year
age=1:12
n.age=length(age)
n.year=length(year)
agem = matrix(age,nrow=n.year,ncol=n.age,byrow=T)
yearm = matrix(year,nrow=n.year,ncol=n.age,byrow=F) 
sfm = matrix(sf,nrow=n.year,ncol=n.age,byrow=F)

lts = data.frame(len=as.numeric(unlist(len[,2:(n.age+1)])),
                 age = as.numeric(unlist(agem)),year = as.numeric(unlist(yearm)),sf=as.numeric(unlist(sfm)))
lts$cohort = lts$year - lts$age  
#lts$age = lts$age + lts$sf
ltsp = subset(lts,!is.na(len))

ltsp$ia = as.numeric(as.factor(ltsp$age))-1 
ltsp$iy = as.numeric(as.factor(ltsp$year))-1  
ltsp$ic = as.numeric(as.factor(ltsp$cohort))-1

#compile("vonb1.cpp")
dyn.load(dynlib("vonb")) 
#dyn.unload("vonb")


tmb.data = list(
  nobs = nrow(ltsp),
  nc = length(unique(ltsp$cohort)),
  na = length(unique(ltsp$age)),
  y = ltsp$len,    
  sf = ltsp$sf,
  cohort = ltsp$cohort,
  age = ltsp$age,
  ic = ltsp$ic, 
  iy = ltsp$iy,
  ia = ltsp$ia,
  REk = 0, #0 for no, 1 for yes;       
  RELinf = 0, #0 for no, 1 for yes;      
  REpo = 0 #0 for no, 1 for yes;
)

parameters <- list( 
  log_Linf = log(120),
  log_k_main=log(0.04),  
  logit_po=-5,  
  log_std_log_Linf = log(0.1),
  log_std_log_k = log(0.1),     
  log_std_lpo = log(0.1), 
  log_std_me = log(0.2),
  rlog_Linf = rep(0,tmb.data$nc),
  rlog_k_dev = matrix(0,tmb.data$na,tmb.data$nc),
  rlogit_po = rep(0,tmb.data$nc)
)

upper.L =  list( 
  log_Linf = log(220),
  log_k_main=log(0.4), 
  logit_po=1, 
  log_std_log_Linf = log(1),
  log_std_log_k = log(1),     
  log_std_lpo = log(0.25),  
  log_std_me = log(0.8)
) 
lower.L =  list(
  log_Linf = log(80),
  log_k_main=log(1e-7),  
  logit_po=-Inf,
  log_std_log_Linf = log(1e-6), 
  log_std_log_k = log(1e-6),      
  log_std_lpo = log(1e-6), 
  log_std_me = log(0.01)
) 

upper = unlist(upper.L)
lower = unlist(lower.L)

#####
###1 no cohort effects for k or Linf ################;
#####
map1 = list(
  log_std_log_Linf=factor(NA),
  log_std_log_k=factor(NA),    
  log_std_lpo=factor(NA),
  rlog_Linf=factor(rep(NA,length(parameters$rlog_Linf))),
  rlog_k_dev=factor(matrix(NA,tmb.data$na,tmb.data$nc)),
  rlogit_po=factor(rep(NA,length(parameters$rlogit_po)))
)

obj <- MakeADFun(tmb.data, parameters,map=map1, DLL="vonb") 

low=lower[match(names(obj$par),names(lower))]
upp=upper[match(names(obj$par),names(upper))]

length(obj$par)
length(low)
length(upp)

obj$gr(obj$par)

opt1<-nlminb(obj$par,obj$fn,obj$gr,lower=low,upper=upp)

rep1 = obj$report()   
sd.rep1<-sdreport(obj)
rep1$nparm=length(obj$par)
stats1 = list(Dev = 2*opt1$objective,nll = opt1$objective,nparm = rep1$nparm)
stats1$AIC = stats1$Dev + 2*rep1$nparm
stats1$BIC= stats1$Dev + log(tmb.data$nobs)*rep1$nparm


#####
### 2 cohort effects for k and Linf ################;
#####

parameters$log_std_log_Linf = log(0.2)
parameters$log_std_log_k = log(0.2)     
parameters$log_std_lpo = log(0.2)  
tmb.data$REk=1     
tmb.data$RELinf=1  
tmb.data$REpo=0

map2 = list(
  #log_std_log_Linf=factor(NA),
  #log_std_log_k=factor(NA),    
  log_std_lpo=factor(NA),
  #rlog_Linf=factor(rep(NA,length(parameters$rlog_Linf))),
  #rlog_k_dev=factor(matrix(NA,tmb.data$na,tmb.data$nc)),
  rlogit_po=factor(rep(NA,length(parameters$rlogit_po)))
)


obj <- MakeADFun(tmb.data, parameters,map=map2, random=c("rlog_Linf","rlog_k_dev"), DLL="vonb1") 

low=lower[match(names(obj$par),names(lower))]
upp=upper[match(names(obj$par),names(upper))]

length(obj$par)
length(low)
length(upp) 


obj$gr(obj$par)

opt2<-nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper)

rep2 = obj$report()   
sd.rep2<-sdreport(obj)
rep2$nparm=length(obj$par)
stats2 = list(Dev = 2*opt2$objective,nll = opt2$objective,nparm = rep2$nparm)
stats2$AIC =  2*rep2$nparm + stats2$Dev
stats2$BIC= stats2$Dev + log(tmb.data$nobs)*rep2$nparm

#####
###3 cohort effects for Linf ################;
#####

tmb.data$REk=0     
tmb.data$RELinf=1  
tmb.data$REpo=0

map3 = list(
  #log_std_log_Linf=factor(NA),
  log_std_log_k=factor(NA),    
  log_std_lpo=factor(NA),
  #rlog_Linf=factor(rep(NA,length(parameters$rlog_Linf))),
  rlog_k_dev=factor(matrix(NA,tmb.data$na,tmb.data$nc)),
  rlogit_po=factor(rep(NA,length(parameters$rlogit_po)))
)

obj <- MakeADFun(tmb.data, parameters,map=map3, random=c("rlog_Linf"), DLL="vonb1")

low=lower[match(names(obj$par),names(lower))]
upp=upper[match(names(obj$par),names(upper))]

length(obj$par)
length(low)
length(upp) 

obj$gr(obj$par)

opt3<-nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper)

rep3 = obj$report()   
sd.rep3<-sdreport(obj)
rep3$nparm=length(obj$par)
stats3 = list(Dev = 2*opt3$objective,nll = opt3$objective,nparm = rep3$nparm)
stats3$AIC = stats3$Dev + 2*rep3$nparm
stats3$BIC= stats3$Dev + log(tmb.data$nobs)*rep3$nparm


#####
###4 cohort effects for k ################;
#####

tmb.data$REk=1     
tmb.data$RELinf=0  
tmb.data$REpo=0

map4 = list(
  log_std_log_Linf=factor(NA),
  #log_std_log_k=factor(NA),    
  log_std_lpo=factor(NA),
  rlog_Linf=factor(rep(NA,length(parameters$rlog_Linf))),
  #rlog_k_dev=factor(matrix(NA,tmb.data$na,tmb.data$nc)),
  rlogit_po=factor(rep(NA,length(parameters$rlogit_po)))
)

obj <- MakeADFun(tmb.data, parameters,map=map4, random=c("rlog_k_dev"), DLL="vonb")

low=lower[match(names(obj$par),names(lower))]
upp=upper[match(names(obj$par),names(upper))]

length(obj$par)
length(low)
length(upp) 

obj$gr(obj$par)

opt4<-nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper)

rep4 = obj$report()   
sd.rep4<-sdreport(obj)
rep4$nparm=length(obj$par)
stats4 = list(Dev = 2*opt4$objective,nll = opt4$objective,nparm = rep4$nparm)
stats4$AIC = stats4$Dev + 2*rep4$nparm
stats4$BIC= stats4$Dev + log(tmb.data$nobs)*rep4$nparm 


#######
### 5 cohort effects for k + stddev k by age block
#####

#compile("vonb2.cpp")
dyn.load(dynlib("vonb2")) 
#dyn.unload("vonb")
 
tmb.data$REk=1     
tmb.data$RELinf=0  
tmb.data$REpo=0

#map r_log_kdev to separate estimates by ages
parameters$log_std_log_k <- rep(log(0.1),tmb.data$na)
upper.L$log_std_log_k <- rep(log(1),tmb.data$na)
lower.L$log_std_log_k <- rep(log(1e-6),tmb.data$na)


upper = unlist(upper.L)
lower = unlist(lower.L)

map5 = list(
  log_std_log_Linf=factor(NA),
  log_std_log_k=factor(c(1,2,3,4,5,6,7,8,9,10,11,12)),    ##change la std en 3 groupes d'age 1-2 , 3-9 ,10-12
  log_std_lpo=factor(NA),
  rlog_Linf=factor(rep(NA,length(parameters$rlog_Linf))),
  #rlog_k_dev=factor(NA),
  rlogit_po=factor(rep(NA,length(parameters$rlogit_po)))
)


obj <- MakeADFun(tmb.data, parameters,map=map5, random=c("rlog_k_dev"), DLL="vonb2")

low=lower[match(names(obj$par),names(lower))]
upp=upper[match(names(obj$par),names(upper))]

length(obj$par)
length(low)
length(upp) 

obj$gr(obj$par)

opt5<-nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper,control=list(iter.max=10000,eval.max=10000))

rep5 = obj$report()   
sd.rep5<-sdreport(obj)
rep5$nparm=length(obj$par)
stats5 = list(Dev = 2*opt5$objective,nll = opt5$objective,nparm = rep5$nparm)
stats5$AIC = stats5$Dev + 2*rep5$nparm
stats5$BIC= stats5$Dev + log(tmb.data$nobs)*rep5$nparm 



################################################
##################Results plotting##############
################################################

#set up the run you wanna plot
rep <- rep5
sd.rep <- sd.rep5
stats <- stats5
opt <- opt5


##output text saving

out.tab = matrix(NA,(5+stats$nparm),1)
out.tab[1,1] = stats$Dev 
out.tab[2,1] = stats$nparm
out.tab[3,1] = stats$AIC
out.tab[4,1] = stats$BIC 
out.tab[5,1] = mean(rep$resid**2) 
out.tab[6:nrow(out.tab),1] = opt$par  
out.tab[c(6:7,9:nrow(out.tab)),]=exp(out.tab[c(6:7,9:nrow(out.tab)),])  

rownames(out.tab)=c("Dev","Nparm","AIC","BIC","MSE",gsub("log_","",names(opt$par)))
colnames(out.tab) <- "Model name"

ctext = 'put some caption text here'
print(xtable(out.tab,digits=c(0,3),caption=ctext,
             align=c('l','r')),type='html',
      file='model_output.doc',caption.placement="top")



## plotting ####;

ind = names(sd.rep$value)=='log_Linf_c'
Linf = exp(sd.rep$value[ind])
Linf.se = sd.rep$sd[ind]  
Linf.L = exp(sd.rep$value[ind] - 1.96*Linf.se)
Linf.U = exp(sd.rep$value[ind] + 1.96*Linf.se)

ucohort <- sort(unique(tmb.data$cohort)) 

gname <- 'Linf_est.jpeg'
jpeg(file=gname,width=3,height=3,units='in',res=300) 

par(mar=c(3.5,3.5,0.2,0.2))
xt = range(ucohort) 
yt=c(min(Linf.L),max(Linf.U))
plot(ucohort,Linf,xlab="",ylab="",type='l',lty=1,lwd=3,ylim=yt,las=1)
segments(ucohort,Linf.L,ucohort,Linf.U,col='grey',lty=1)
lines(ucohort,Linf,lty=1,lwd=3)
abline(h=mean(Linf),lty=2)
mtext(side=2,line=2.2,outer=F,substitute(paste(L[infinity]," (cm)")),las=0)  
dev.off()  

##################Standard diags###########################;

fitdata = data.frame(year=ltsp$year,age=floor(ltsp$age),cohort = ltsp$cohort,
                     predict=rep$pred,std.resid = rep$std_resid)

gname <- 'resid5p.jpeg'
jpeg(file=gname,width=4,height=5,units='in',res=300) 

resid4p(fitdata)  
dev.off()


######## mat plot  ###############;

gname <- 'matplot_year.jpeg'
jpeg(file=gname,width=4,height=4,units='in',res=300) 

par(mar=c(3,3,0.5,1.5))
pmplot(fitdata$age,fitdata$year,fitdata$std.resid,scale=5)

mtext(line=2,outer=F,"Age",side=2)
mtext(line=2,outer=F,"Year",side=1)

dev.off()     



#############################  fits by age #################################################;


my.padding <- list(layout.heights = list( 
  top.padding = 0, 
  main.key.padding = 0, 
  key.axis.padding = 0, 
  axis.xlab.padding = 0, 
  xlab.key.padding = 0, 
  key.sub.padding = 0), 
  layout.widths = list( 
    left.padding = 0, 
    key.ylab.padding = 0, 
    ylab.axis.padding = 0, 
    axis.key.padding = 0, 
    right.padding = 0) 
) 


fitdata = data.frame(year=ltsp$year,age=floor(ltsp$age),cohort = ltsp$cohort,
                     predict=rep$pred,std.resid = rep$std_resid)

uage=1:12
fitdata$len=tmb.data$y  
xt <- range(fitdata$year)
yt <- range(fitdata$len)

gname <- 'len_by_age.jpeg'
jpeg(file=gname,width=6,height=6,units='in',res=300)

ttl <- list(c("Length at age by year"), cex=1)
xttl <- list(c("Year"), cex=1)
yttl <- list(c("Length (cm)"), cex=1)
stripttl <- list(cex=0.2,font=1)
ax <- list(cex=0.5,relation="free")

p1=xyplot(predict~year|factor(age), data=fitdata, type=c("l"), main=ttl, xlab=xttl,lwd=2,
          ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(2,6,1),las=1,
          xlim=c(1983,2016),
          col=2,par.strip.text = list(cex=0.6),
          par.settings = list(layout.heights = list(strip = .8)))

p2=xyplot(len~year|factor(age), data=fitdata, type=c("p"), main=ttl, xlab=xttl,lwd=2,
          ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(2,6,1),las=1,
          xlim=c(1983,2016),
          col=1,par.strip.text = list(cex=0.6),
          par.settings = list(layout.heights = list(strip = .8),
                              plot.symbol = list(col = "black",fill = "black",cex = 0.5,pch = 21)))

print(p2 + as.layer(p1),  packet.panel = function(layout, row, column, ...) { 
  layout <- layout[c(2, 1, 3)]
  packet.panel.default(layout = layout, row = column, column = row, ...) 
})


dev.off()

gname <- 'len_by_cohort.jpeg'
jpeg(file=gname,width=5,height=7,units='in',res=300)

ttl <- list(c("Length at age by cohort"), cex=1)
xttl <- list(c("Year"), cex=1)
yttl <- list(c("Length (cm)"), cex=1)
stripttl <- list(cex=0.2,font=1)
ax <- list(cex=0.5,relation="free")

p1=xyplot(predict~age|factor(cohort), data=fitdata, type=c("l"), main=ttl, xlab=xttl,lwd=2,
          ylab=yttl, strip.text=stripttl, as.table=TRUE, layout=c(5,9,1),las=1,
          xlim=c(1,13),ylim=c(10,100),
          col=2,par.strip.text = list(cex=0.6),
          par.settings = list(layout.heights = list(strip = .8)))

p2=xyplot(len~age|factor(cohort), data=fitdata, type=c("p"), main=ttl, xlab=xttl,lwd=2,
          ylab=yttl, strip.text=stripttl, as.table=TRUE, layout=c(5,9,1),las=1,
          xlim=c(1,13),ylim=c(10,100),
          col=1,par.strip.text = list(cex=0.6),
          par.settings = list(layout.heights = list(strip = .8),
                              plot.symbol = list(col = "black",fill = "black",cex = 0.5,pch = 21)))

print(p1 + as.layer(p2))

dev.off()


datk = expand.grid(cohort=sort(unique(tmb.data$cohort)), age=1:12)
datk$dk = exp(as.vector(t(rep$log_k)))
datk$year = datk$cohort+datk$age

datk = subset(datk,(year>=1983)&(year<=2016))


gname <- 'k_3D.jpeg'
jpeg(file=gname,width=5,height=5,units='in',res=300)

wireframe(dk ~ year*age, data = datk,zoom=0.95,
          xlab = "Year", ylab = "Age",
          main = "k",
          drape = TRUE,
          colorkey = TRUE,
          screen = list(z = -60, x = -60)
)

dev.off()


gname <- 'k_contour.jpeg'
jpeg(file=gname,width=5,height=5,units='in',res=300)

contourplot(dk ~ year*age, data=datk,region=T, 
            at=seq(0.03,0.063,by=0.002),
            colorkey=list(at=seq(0.03,0.063,by=0.002)), 
            col.regions=terrain.colors(17),
            panel = latticeExtra::panel.2dsmoother) 

dev.off()


####
#windows();
ggplot(fitdata)+geom_line(aes(year,predict,col=factor(age)))+ theme(legend.position = "bottom")
min <- NULL
max<-NULL
for (i in 1:length(unique(fitdata$age))){  
  agei <- unique(fitdata$age)[i]
  min[i]<-min(fitdata[fitdata$age==agei,]$predict)
  max[i]<-max(fitdata[fitdata$age==agei,]$predict)
}
cv <- max-min


########################
######MODEL APPLICATION = L@A & W@A computation
#########################

#save final parameter
rep <- rep4
sd.rep <- sd.rep4 
final <- list(rho_c = rep$po_c ,Linf_c = rep$Linf_c,k_c = rep$ck)

############  Jan 1 lts  ############################; NOEL CODE, DIDN'T UPDATE IT YET

uage=1:20 
ucohort <- c(sort(unique(fitdata$cohort)),2016:2018)
n.age = length(uage)
n.cohort = length(ucohort)

pred.data=data.frame(age=rep(uage,times=n.cohort),cohort = rep(ucohort,each=n.age))
pred.data$year = pred.data$age + pred.data$cohort
pred.data$cfac=factor(pred.data$cohort)


pred.data$length <-



np = length(rep$Linf_c)
Linf = c(rep$Linf_c,rep(rep$Linf_c[np],3)) 
kv = c(rep$k_c,rep(rep$k_c[np],3))  
pred.data$predict = as.numeric(vonb(Linf,parameters$b,kv,exp(opt4$par[3]),pred.data))

temp = subset(pred.data,(pred.data$year>=1983)&(pred.data$year<=2020))
temp$cfac=NULL
temp = temp[order(temp$year,temp$age),]

uyear = sort(unique(temp$year))
nyear=length(uyear)
uage = sort(unique(temp$age))
nage=length(uage)
lmat = matrix(NA,nrow=length(uyear),ncol=length(uage))

for(y in uyear){
  ind = temp$year==y
  ly = temp$predict[ind]
  ai = temp$age[ind]
  lmat[y-1982,ai]=ly
}

for(a in uage){
  ind = is.na(lmat[,a])
  if(any(ind)){
    pmax = max((1:nyear)[!ind]) 
    pmin = min((1:nyear)[!ind])
    
    py = mean(uyear[ind])
    if(py>mean(uyear)){ind1 = (pmax-2):pmax}
    if(py<mean(uyear)){ind1 = pmin:(pmin+2)}
    lmat[ind,a]=mean(lmat[ind1,a])
  }}

colnames(lmat) = as.character(uage) 
rownames(lmat) = as.character(uyear)

write.table(lmat,file='pred.lts.dat')

xdigits = c(1,rep(2,nage)); #a digit for rowname, and all rows;
xalign = c('c',rep('r',nage)); # alignment 'l','c', or 'r';


ctext = 'put some caption text here'
print(xtable(lmat,digits=xdigits,caption=ctext,align=xalign),type='html',
      file='predicted_len.html',caption.placement="top")


############  midyear lts  ############################;

uage=1:20 + 0.5
ucohort <- c(sort(unique(fitdata$cohort)),2016:2018)
n.age = length(uage)
n.cohort = length(ucohort)

pred.data=data.frame(age=rep(uage,times=n.cohort),cohort = rep(ucohort,each=n.age))
pred.data$year = pred.data$age + pred.data$cohort
pred.data$cfac=factor(pred.data$cohort)

np = length(rep$Linf_c)
Linf = c(rep$Linf_c,rep(rep$Linf_c[np],3)) 
kv = c(rep$k_c,rep(rep$k_c[np],3))  
pred.data$predict = as.numeric(vonb(Linf,parameters$b,kv,exp(opt3$par[3]),pred.data))

temp = subset(pred.data,(pred.data$year>=1983)&(pred.data$year<=2020))
temp$cfac=NULL
temp = temp[order(temp$year,temp$age),]

uyear = sort(unique(temp$year))
nyear=length(uyear)
uage = sort(unique(temp$age))
nage=length(uage)
lmat = matrix(NA,nrow=length(uyear),ncol=length(uage))

for(y in uyear){
  ind = temp$year==y
  ly = temp$predict[ind]
  ai = temp$age[ind]
  lmat[y-1982,ai]=ly
}

for(a in uage){
  ind = is.na(lmat[,a])
  if(any(ind)){
    pmax = max((1:nyear)[!ind]) 
    pmin = min((1:nyear)[!ind])
    
    py = mean(uyear[ind])
    if(py>mean(uyear)){ind1 = (pmax-2):pmax}
    if(py<mean(uyear)){ind1 = pmin:(pmin+2)}
    lmat[ind,a]=mean(lmat[ind1,a])
  }}

colnames(lmat) = as.character(uage) 
rownames(lmat) = as.character(uyear)

write.table(lmat,file='pred_midy_lts.dat')

xdigits = c(1,rep(2,nage)); #a digit for rowname, and all rows;
xalign = c('c',rep('r',nage)); # alignment 'l','c', or 'r';


ctext = 'put some caption text here'
print(xtable(lmat,digits=xdigits,caption=ctext,align=xalign),type='html',
      file='predicted__midy_len.html',caption.placement="top")




require(reshape2)
library(lattice)
library(latticeExtra)
library('ggplot2')    
library(lme4)          
library(cAIC4)  

#setwd("C:\\home\\CADIGAN\\GradProgram\\2019\\F6005\\project1\\data_EDA\\weights")  

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


vec_func = function(mdat){
  vdat = melt(mdat,id=c("Year"))
  vdat$Age = as.numeric(substr(vdat$variable,4,5))
  vdat$variable=NULL
  vdat$index=vdat$value
  vdat$value=NULL
  return(vdat)
}

## Estimate RV length weight relationship;


mlen = read.table(file='../RV_mlen@age.txt',header=TRUE)
mwt = read.table(file='../RV_mwt@age.txt',header=TRUE)
year=mlen[,1]

dat = data.frame(log_length = unlist(log(mlen[,2:13])),
 log_weight = unlist(log(mwt[,2:13])), 
 year = as.vector(matrix(year,nrow=length(year),ncol=12)))
 
dat = subset(dat,!is.na(dat$log_length))
dat$Length = exp(dat$log_length)
dat$Weight = exp(dat$log_weight)
dat$fyear = as.factor(dat$year)

#gname <- 'length_weight.jpeg'
#jpeg(file=gname,width=3,height=3,units='in',res=300)

par(mar=c(3,3,0.2,0.2),mgp=c(2,1,0))
plot(dat$Length,dat$Weight,pch=19,xlab='Length (cm)',ylab='Weight (kg)') 

#dev.off()

mixfit <- lmer(log_weight ~ log_length + (log_length|fyear),REML = FALSE, data=dat)

summary(mixfit)

logLik(mixfit) 
AIC(mixfit)
cAIC(mixfit)

fit <- glm(log_weight ~ log_length, data=dat)

summary(fit)

logLik(fit) 
AIC(fit)

## AIC for the linear model with random slopes and intercepts was substantially lower;
## but for simplicity I have not pursued using an annual length-weight relationship.

# Mixfit annual length-weight relationship coefficients
parms = fixef(mixfit)
rparms = ranef(mixfit)

mdat = data.frame(year = as.numeric(rownames(rparms$fyear)),
  a = exp(rparms$fyear[,1]+parms[1]),
  b = rparms$fyear[,2]+parms[2])

a = exp(coef(fit)[1])
b = coef(fit)[2]    

par(mfrow=c(2,1),mgp=c(2,1,0),mar=c(3,3,0.2,0.2))  
plot(a~year,data=mdat,typ='l',ylab='a',xlab='year')
abline(h=a,lty=2)
plot(b~year,data=mdat,typ='l',ylab='b',xlab='year')  
abline(h=b,lty=2)

#glmmTMB(log_weight ~ log_length + ar1(log_length|fyear), data=dat,family=gaussian)


## make weights from VonB model lengths and do some comparison plotting with DFO

cnames = c('Year',paste('Age',1:16,sep=""))
len.dat = read.table(file='pred.lts.dat',header=FALSE,col.names=cnames,skip=1)
midy.len.dat = read.table(file='pred_midy_lts.dat',header=FALSE,col.names=cnames,skip=1) 

wtm = a*(as.matrix(len.dat[,2:17])**b) 
midy_wtm = a*(as.matrix(midy.len.dat[,2:17])**b)

RV.wt = read.table(file='..\\RV_mwt@age.txt',header=TRUE)
Comm.wt = read.table(file='..\\Comm_wt.txt',header=TRUE)  
Comm.wt = Comm.wt[!is.na(match(Comm.wt$Year,year)),]

temp = cbind(year,wtm)
colnames(temp) = c('Year',paste('Age',1:16,sep=''))
write.table(temp,file='..//stock_wt.dat')

wt.dat = vec_func(as.data.frame(temp)) 
colnames(wt.dat) = c('Year','Age','vb_wt')
 
temp = cbind(year,midy_wtm)
colnames(temp) = c('Year',paste('Age',1:16,sep=''))  
write.table(temp,file='..//midy_wt.dat')

wt.dat$vb_midy_wt = vec_func(as.data.frame(temp))[,3]

temp = vec_func(as.data.frame(RV.wt))
colnames(temp) = c('Year','Age','RV_wt')

wt.dat = merge(wt.dat,temp,by.x=c('Year','Age'),all.x=TRUE)

temp = vec_func(as.data.frame(Comm.wt))
colnames(temp) = c('Year','Age','CM_wt')  
  
wt.dat = merge(wt.dat,temp,by.x=c('Year','Age'),all.x=TRUE)

# names(wt.dat)
#"Year"       "Age"        "vb_wt"      "vb_midy_wt" "RV_wt"      "CM_wt"

  xttl <- list(c("Year"), cex=1)
  yttl <- list(c("Weight (kg)"), cex=1)
  stripttl <- list(cex=0.2,font=1)
  ax <- list(cex=0.5,y=list(relation=c("free")))
  my.key <- list(space = "top",columns=3,
                border = FALSE,
                size = 5,
                lines = list(lty = 1,col=1:3,lwd=2),
                text = list(c("Jan1","Mid-year","RV Mean")))

                
png(file='compare_wt_rv.png',width=5,height=5,units='in',res=300)

  print(xyplot(vb_wt+vb_midy_wt+RV_wt~Year|factor(Age), data=wt.dat,
  type="l", xlab=xttl,lwd=2,lty = 1,
  ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(4,4,1),las=1,
  key=my.key,col=1:3,par.strip.text = list(cex=0.6),
  par.settings = list(layout.heights = list(strip = .8))))
  
dev.off()

  my.key <- list(space = "top",columns=3,
                border = FALSE,
                size = 5,
                lines = list(lty = 1,col=1:3,lwd=2),
                text = list(c("Jan1","Mid-year","Comm")))

                
png(file='compare_wt_cm.png',width=5,height=5,units='in',res=300)

  print(xyplot(vb_wt+vb_midy_wt+CM_wt~Year|factor(Age), data=wt.dat,
  type="l", xlab=xttl,lwd=2,lty = 1,
  ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(4,4,1),las=1,
  key=my.key,col=1:3,par.strip.text = list(cex=0.6),
  par.settings = list(layout.heights = list(strip = .8))))
  
dev.off()



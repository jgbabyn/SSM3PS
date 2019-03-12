
setwd("C:\\Users\\mroberts\\OneDrive - Memorial University of Newfoundland\\Marine Institute\\Classes\\FISH_6005\\Project 1\\Project 1\\Working model\\SepF_std+comm_wts+rhoq+all_survs")

my.padding <- list(layout.heights = list( 
                        top.padding = 0, 
                        main.key.padding = 1, 
                        key.axis.padding = 0, 
                        axis.xlab.padding = 0, 
                        xlab.key.padding = 0, 
                        key.sub.padding = 0), 
                layout.widths = list( 
                        left.padding = 1,
                        key.ylab.padding = 0, 
                        ylab.axis.padding = 1, 
                        axis.key.padding = 0, 
                        right.padding = 1) 
                ) 

require(lattice)
library(latticeExtra)
require(reshape2)  
library(lme4)    
library(ggplot2)  
library(psych)    
library(stargazer)  
library(numDeriv)   
library(splines)
require(corrplot)  
require(xtable)
 
load('tmb.RData')
load('fit.RData')

setwd(paste(getwd(),"/figs",sep=""))

# get log indices and wt=0 for index=0
indices$YC = indices$Year-indices$Age 


setwd("C:\\Users\\mroberts\\OneDrive - Memorial University of Newfoundland\\Marine Institute\\Classes\\FISH_6005\\Project 1\\Project 1\\SAM")

sam.F =  read.table("..\\SAM\\fbar_updated.dat",header=T,sep = " ")
#year<-seq(1959, 2015, by=1)

#assess.ssb = read.table("..\\data_EDA\\assmt_SSB.txt",
#header=F,col.names=c('year','ssb'),sep = "\t")
#assess.ssb$ssb = assess.ssb$ssb/1000

sam.ssb = read.table("..\\SAM\\ssb_updated.dat",header=T,sep = " ")

setwd("C:\\Users\\mroberts\\OneDrive - Memorial University of Newfoundland\\Marine Institute\\Classes\\FISH_6005\\Project 1\\Project 1\\Working model\\SepF_std+comm_wts+rhoq+all_survs\\figs")


gname <- 'exp_pred_rho.jpeg'
jpeg(file=gname,width=4,height=4,units='in',res=300)
expected_ratio<-c(1.1778062, 1.0216776, 0.9474586, 0.9517923, 0.9413553, 0.9336698, 0.9330784, 0.8977035, 0.8858351, 0.8584615, 0.8738739, 0.8620690, 0.9230769)
pred_ratio<-exp(rep$log_rho)
plot(expected_ratio-pred_ratio~seq(from=2, to=14, by=1), type="l", lwd=2, xlab="Ages", ylab="Expected-predicted rho")
abline(h=0, lty=2, col="red")
dev.off()

year=seq(from=1959, to=2016, by=1)

value.names = names(sd.rep$value)

ind = value.names == "log_aveF_46";
est_46 = exp(sd.rep$value[ind]);
stdi =  sd.rep$sd[ind]
low_46 = exp(sd.rep$value[ind] - qnorm(0.975)*stdi); 
high_46 = exp(sd.rep$value[ind] + qnorm(0.975)*stdi); 

ind = value.names == "log_aveF_69";
est_69 = exp(sd.rep$value[ind]);
stdi =  sd.rep$sd[ind]
low_69 = exp(sd.rep$value[ind] - qnorm(0.975)*stdi); 
high_69 = exp(sd.rep$value[ind] + qnorm(0.975)*stdi);

ylim=range(low_46,low_69,high_69,high_46)
#ylim[2]=2.5
 
gname <- 'aveF.jpeg'
jpeg(file=gname,width=4,height=4,units='in',res=300)
  par(mfcol=c(1,1),oma=c(0,1,0.5,0.5),mar=c(3,2.5,0,0),las=1)
  plot(year,est_46,type='n',lwd=2,xlab='',ylab='',xlim=c(1959,2016),ylim=ylim)
  
  sx = year
  polygon(c(sx,rev(sx)),c(low_46,rev(high_46)),col='grey',border = NA)       
  box(lty=1)   
  lines(year,est_46,type='l',lty=1,lwd=2)
  #lines(assess.F$year,assess.F$F4_6,type='l',lty=3,lwd=2)  
  lines(assess.F$year,sam.F$Estimate,type='l',lty=1,lwd=2,col='red')    
  legend("topright",lwd=c(2,2),lty=c(1,1),bty='n',
    col=c('black','red'),legend=c('SSM',"SAM"),cex=0.75) 
  
  mtext(side=2,outer=F,line=2.5,'Ave F (Ages 4-6)',las=0)        
  #plot(year,est_69,type='n',lwd=2,xlab='',ylab='',xlim=c(1959,2016),ylim=ylim)
  #polygon(c(sx,rev(sx)),c(low_69,rev(high_69)),col='grey',border = NA)       
  #box(lty=1)   
  #lines(year,est_69,type='l',lty=1,lwd=2)   
  #lines(assess.F$year,assess.F$F6_9,type='l',lty=3,lwd=2)
  #mtext(side=2,outer=F,line=2.5,'Ave F (Ages 6-9)',las=0)
    
  mtext(side=1,outer=F,line=2,'Year') 
dev.off()

year=seq(from=1959, to=2016, by=1)

ind = value.names == "log_ssb";
est = exp(sd.rep$value[ind]);
stdi =  sd.rep$sd[ind]
low = exp(sd.rep$value[ind] - qnorm(0.975)*stdi); 
high = exp(sd.rep$value[ind] + qnorm(0.975)*stdi);

ylim=range(est,low,high,sam.ssb$Estimate)
#ylim[2]=250

gname <- 'ssb.jpeg'
jpeg(file=gname,width=4,height=3,units='in',res=300)  
  par(mar=c(3,3.5,0.5,4),las=1)

  plot(year,est,type='n',lwd=2,xlab='',ylab='',ylim=c(0,120))
  
  sx = year
  polygon(c(sx,rev(sx)),c(low,rev(high)),col='grey',border = NA)       
  box(lty=1)   
  lines(year,est,type='l',lty=1,lwd=2)
  #lines(year,assess.ssb$ssb,type='l',lty=3,lwd=2)
  lines(year,sam.ssb$Estimate,type='l',lty=1,lwd=2,col='red') 
  mtext(side=2,outer=F,line=2.5,'SSB (Kt)',las=0)      
  mtext(side=1,outer=F,line=2,'Year')   
    
  par(new=T)
  plot(year,stdi,axes=F,xlab="",ylab="",type='l',lty=1,lwd=2,col=grey(0.6))
  axis(side=4,las=1,hadj=0)
  mtext("Standard Error",side=4,line=3,las=0) 

  legend("topright",lwd=c(2,2,2),lty=c(1,1,1),bty='n',
    col=c('black',grey(0.6),'red'),legend=c('SSM','SSM_SE',"SAM"),cex=0.75)
 
  
dev.off()

df<-data.frame(year, est)
rel_ssb<-est/est[36]

plot(year, rel_ssb, type="l", xlim=c(1984,2016), lwd=2, ylim=c(0.5, 7), ylab="Relative SSB", xlab="Year")
abline(h=1, lty=2, col="lightgrey")
abline(h=3, lty=2, col="lightgrey")

fit = data.frame(resid = rep$resid_index,resid_std=rep$std_resid_index)

fit$Age = indices$Age    
fit$Year = indices$Year
fit$Cohort = fit$Year-fit$Age
fit$survey = indices$survey
fit$pred = rep$Elog_index
fit$index = indices$index 
fit$log_index = log(indices$index)

fit$colr = 'deepskyblue'
fit$colr[fit$resid>0]='firebrick2' 
fit$pch=21       
fit$colr[(indices$i_zero==1)&(fit$resid>0)]='azure2'
fit$colr[(indices$i_zero==1)&(fit$resid<=0)]='dimgray'

my.padding$layout.heights$main.key.padding = 0
#fit$pch[indices$i_zero==1]=25

jDarkGray <- 'grey20'
jPch <- 21
resid.plot1 = xyplot(factor(Age)~Year|survey, data=fit, 
       ylab='Age',main='Residuals',
       cex = 1.5*sqrt(abs(fit$resid_std)/pi), fill.color = fit$colr,
       col = jDarkGray,
       par.settings = my.padding,
       par.strip.text = list(cex=0.5),
       panel = function(x, y, ..., cex, fill.color, pch, subscripts) { 
         panel.abline(v=seq(1983,2016,by=5), col.line=grey(0.9))
         panel.xyplot(x, y, cex = cex[subscripts],
                       pch=fit$pch, fill = fill.color[subscripts], ...)
         })
         
print(resid.plot1)
                  
gname = c("resid_matrix_bubbles.jpeg")
jpeg(file=gname,width=5,height=5,units='in',res=300)
  print(resid.plot1)
dev.off()


par.names = names(sd.rep$par.fixed)
#ind = par.names == "log_qparm"
ind = value.names == "log_qparm";

qest = aggregate(indices$iq,list(age=indices$Age,survey=indices$survey),unique)
qest$iq=qest$x+1  
temp = as.numeric(table(qest$iq))
qest$nx = temp[qest$iq] 

qest$Eest = exp(sd.rep$value[ind]);
qest$sd =  sd.rep$sd[ind]
qest$E.L = exp(sd.rep$value[ind] - qnorm(0.975)*qest$sd); 
qest$E.U = exp(sd.rep$value[ind] + qnorm(0.975)*qest$sd);

#qest$est = sd.rep$par.fixed[ind][qest$iq]
#qest$sd = sqrt(diag(sd.rep$cov.fixed))[ind][qest$iq]
#qest$Eest = exp(qest$est)
#qest$E.L = exp(qest$est + qnorm(0.025)*qest$sd) 
#qest$E.U = exp(qest$est + qnorm(0.975)*qest$sd)

qest$colr='grey'
qest$colr[qest$nx>1]='black'


temp = tapply(qest$Eest,qest$survey,max)
qest$max_Eest = temp[match(qest$survey,names(temp))]

qest$sEest = qest$Eest/qest$max_Eest
qest$sE.L = qest$E.L/qest$max_Eest 
qest$sE.U = qest$E.U/qest$max_Eest

qest$survey1 = paste(qest$survey,' max Q = ',round(qest$max_Eest,digits=3),' (x10000)',sep='')

qest  = qest[order(qest$survey1,qest$age),]

ind1 = qest$survey %in% names(temp)


jpeg(file='q.jpeg',width=4,height=4,units='in',res=300)

  xttl <- list(c("Age Class"), cex=1)
  yttl <- list(c("Survey Q Pattern"), cex=1)
  stripttl <- list(cex=0.1,font=1)
  my.key <- list(space = "top",columns=3,
                border = FALSE,
                size = 3,between=0.1,
                lines = list(lty = c(1,1,1),col=c('black','grey','grey'),lwd=c(2,1,1)),
                text = list(c("Est","L","U")))
  
  ax <- list(cex=0.5,relation="free",tck=0.3)
  p1 = xyplot(sEest+sE.L+sE.U~age|factor(survey1), data=qest,type='l',
  xlab=xttl,lwd=c(2,1,1),lty=c(1,1,1),ylim=c(0,1.8),xlim=range(qest$age),
  ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(1,3,1),las=1,
  key=my.key,col=c('black','grey','grey'),par.strip.text = list(cex=0.4),
  par.settings = my.padding)
  
  p2 = xyplot(sEest~age|factor(survey1), data=qest,type='p',fill.color = qest$colr,
  xlab=xttl,lwd=c(2,1,1),lty=c(1,1,1),ylim=c(0,1.8),xlim=range(qest$age),col='black',
  ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(1,3,1),las=1,
  key=my.key,par.strip.text = list(cex=0.4),
  par.settings = my.padding,
       panel = function(x, y, ..., cex, fill.color, subscripts) {
         panel.xyplot(x, y, pch = 21, fill = fill.color[subscripts], ...)
         })
           
  print(p1 + as.layer(p2) + as.layer(p1))
  print(p1)
dev.off()

dsd = sqrt(diag(sd.rep$cov.fixed))
corr.matrix = diag(1/dsd)%*%sd.rep$cov.fixed%*%diag(1/dsd)
rownames(corr.matrix) = names(dsd)                
colnames(corr.matrix) = names(dsd)

 
jpeg(file='corr_parm.jpeg',width=10,height=10,units='in',res=300)
par(mar=c(0,5,0,1))

corrplot.mixed(corr.matrix, tl.pos='lt',tl.srt=45,tl.cex=0.25,cl.cex=0.5,tl.col='black',number.cex=0.25)

dev.off() 

################# pop plot ############################;

yscale=1

ind = value.names == "log_Rec";
est = sd.rep$value[ind]        
sd.est = sd.rep$sd[ind]

Popest = list(recruit = exp(est)/yscale,
   recruit.L95 = exp(est - qnorm(0.975)*sd.est)/yscale, 
   recruit.U95 = exp(est + qnorm(0.975)*sd.est)/yscale) 

ind = value.names == "log_biomass";   
  Popest$biomass = exp(sd.rep$value[ind])/yscale;
  Popest$biomass.L95 = exp(sd.rep$value[ind] - qnorm(0.975)*sd.rep$sd[ind])/yscale; 
  Popest$biomass.U95 = exp(sd.rep$value[ind] + qnorm(0.975)*sd.rep$sd[ind])/yscale;

ind = value.names == "log_ssb";   
  Popest$ssb = exp(sd.rep$value[ind])/yscale;
  Popest$ssb.L95 = exp(sd.rep$value[ind] - qnorm(0.975)*sd.rep$sd[ind])/yscale; 
  Popest$ssb.U95 = exp(sd.rep$value[ind] + qnorm(0.975)*sd.rep$sd[ind])/yscale;   

ind = value.names == "log_aveF_46";   
  Popest$aveF_46 = c(exp(sd.rep$value[ind]),NA);
  Popest$aveF_46.L95 = c(exp(sd.rep$value[ind] - qnorm(0.975)*sd.rep$sd[ind]),NA); 
  Popest$aveF_46.U95 = c(exp(sd.rep$value[ind] + qnorm(0.975)*sd.rep$sd[ind]),NA);  

ind = value.names == "log_aveF_69";   
  Popest$aveF_69 = c(exp(sd.rep$value[ind]),NA);
  Popest$aveF_69.L95 = c(exp(sd.rep$value[ind] - qnorm(0.975)*sd.rep$sd[ind]),NA); 
  Popest$aveF_69.U95 = c(exp(sd.rep$value[ind] + qnorm(0.975)*sd.rep$sd[ind]),NA);

Popest$year=year 
Popest = data.frame(Popest) 

mean.line=TRUE
rec.age=2
stock.name='3Ps Cod'
ydim1 <- c("(millions)")
ydim2 <- c("(000s tonnes)")
  
jpeg(file='pop.jpeg',width=3,height=4,units='in',res=300)

par(mfcol=c(3,1),oma=c(3,4.5,1,1),mar=c(0.5,0,1,0),las=1,cex.lab=0.8,cex.axis=0.8)

x <- c(min(Popest$year),max(Popest$year))
yc <- Popest$year
y <- c(min(Popest$recruit.L95),max(Popest$recruit.U95))
y[2]=240
plot(x,y,xlab='',ylab='',type='n',xaxt='n',lwd=2)

low=Popest$recruit.L95
high=Popest$recruit.U95
polygon(c(yc,rev(yc)),c(low,rev(high)),col='grey',border = NA)       
box(lty=1)   
lines(yc,Popest$recruit,type='l',lwd=2) 
  
mtext(side=3,line=0.5,outer=F,stock.name,las=0,cex=0.8)
mtext(side=2,line=3.5,outer=F,'Recruitment',las=0,cex=0.8)
ytxt <- paste('age',rec.age,' ',ydim1,sep="")
mtext(side=2,line=2.5,outer=F,ytxt,las=0,cex=0.8)
if(mean.line){abline(h=mean(Popest$recruit),lty=2)}

ylim <- c(0,max(Popest$biomass.U95,Popest$ssb.U95))
ylim[2]=300

plot(Popest$year,Popest$biomass,xlab='',ylab='',type='n',lty=1,lwd=2,ylim=ylim,xaxt='n')
low=Popest$biomass.L95
high=Popest$biomass.U95
polygon(c(yc,rev(yc)),c(low,rev(high)),col = rgb(1,0,0,0.3),border = NA)       
box(lty=1)    
abline(h=mean(Popest$biomass),lty=2) 

low=Popest$ssb.L95
high=Popest$ssb.U95
polygon(c(yc,rev(yc)),c(low,rev(high)),col = rgb(0,0,1,0.3),border = NA)       
box(lty=1)   
lines(yc,Popest$ssb,type='l',lwd=2,col='blue')  
lines(yc,Popest$biomass,type='l',lwd=2,col='red')
abline(h=mean(Popest$ssb),lty=2,col='blue') 
legend("topright",lty=c(1,1),lwd=c(2,2),col=c('red','blue'),
bty='n',c('Total','SSB'),cex=0.6)     

mtext(side=2,line=2.5,outer=F,"Biomass (Mt)",las=0,cex=0.8)  

ylim <- c(0,max(Popest$aveF_46.U95,Popest$aveF_69.U95,na.rm=T))
ylim[2]=2

ind = !is.na(Popest$aveF_46)
x = Popest$year[ind]
y = Popest$aveF_46[ind]
low=Popest$aveF_46.L95[ind]
high=Popest$aveF_46.U95[ind]

plot(x,y,xlab='',ylab='',type='n',lty=1,lwd=2,ylim=ylim)
polygon(c(x,rev(x)),c(low,rev(high)),col = rgb(1,0,0,0.3),border = NA)       
box(lty=1)    
abline(h=mean(y),lty=2) 

y1 = Popest$aveF_69[ind]
low=Popest$aveF_69.L95[ind]
high=Popest$aveF_69.U95[ind]

polygon(c(x,rev(x)),c(low,rev(high)),col = rgb(0,0,1,0.3),border = NA)       
box(lty=1)   
lines(x,y1,type='l',lwd=2,col='blue')  
lines(x,y,type='l',lwd=2,col='red')
abline(h=mean(y1),lty=2,col='blue') 
legend("topright",lty=c(1,1),lwd=c(2,2),col=c('red','blue'),
bty='n',c('4-6','6-9'),cex=0.6)     

mtext(side=2,line=2.5,outer=F,'Average F',las=0,cex=0.8)   
mtext(side=1,line=2.2,outer=F,'Year',las=0,cex=0.8)

dev.off()

ssb.matrix = rep$SSB
uage=sort(unique(fit$Age))

temp = expand.grid(Year=year, Age=uage)
temp$ssb = as.vector(rep$SSB_matrix)
            

gname = "SSB_3D.jpeg"
jpeg(file=gname,width=5,height=5,units='in',res=300)

  par(mar=c(0,2,1,0),mgp=c(2,0.7,0))
  persp(year,uage,ssb.matrix, 
      main="SSB-at-age Surface",
      xlab='Year',ylab='Age',zlab='Weight (Kg/tow)',
      col="lightgreen", 
      theta=130, phi=40,
      r=50,
      d=0.1,
      expand=0.5,
      ltheta=90, lphi=180,
      shade=0.25,
      ticktype="detailed",
      nticks=5) # produces the 3-D plot
      
dev.off() 

gname = "N_3D.jpeg"
jpeg(file=gname,width=5,height=5,units='in',res=300)

  par(mar=c(0,2,1,0),mgp=c(2,0.7,0))
persp(year,uage,rep$N, 
      main="Numbers-at-age Surface",
      xlab='Year',ylab='Age',zlab='Number',
      col="lightgreen", 
      theta=130, phi=40,
      r=50,
      d=0.1,
      expand=0.5,
      ltheta=90, lphi=180,
      shade=0.25,
      ticktype="detailed",
      nticks=5) # produces the 3-D plot
      
dev.off()

yearm1=year[1:(length(year)-1)]   

gname = "F_3D.jpeg"
jpeg(file=gname,width=5,height=5,units='in',res=300)

  par(mar=c(0,2,1,0),mgp=c(2,0.7,0))
persp(yearm1,uage,rep$F, 
      main="F-at-age Surface",
      xlab='Year',ylab='Age',zlab='F',
      col="lightgreen", 
      theta=130, phi=40,
      r=50,
      d=0.1,
      expand=0.5,
      ltheta=90, lphi=180,
      shade=0.25,
      ticktype="detailed",
      nticks=5) # produces the 3-D plot
      
dev.off()

age=uage
    
jpeg(file='prop_age.jpeg',width=6,height=3,units='in',res=300)

par(mfrow=c(1,2),mar=c(3,1,1,1),las=1,mgp=c(2,1,0),oma=c(0,3,0,0))

Bs = t(rep$N)
Bs = Bs%*%diag(1/apply(Bs,2,sum))
colnames(Bs) = as.character(year)
rownames(Bs) = as.character(age)
 
col1 <- colorRampPalette(rev(c("red","green","yellow","cyan", "blue")))  
barplot(Bs, col=col1(nrow(Bs)) , border="white", space=0.04, cex.axis=1,
xlab="Year",ylab='',main='')
mtext(side=2,line=2,"Proportion at age 2-12",las=0,outer=T)  
mtext(side=3,line=0,"Numbers at age",las=0)

Bs = t(rep$B)
Bs = Bs%*%diag(1/apply(Bs,2,sum))
colnames(Bs) = as.character(year)
rownames(Bs) = as.character(age)
 
col1 <- colorRampPalette(rev(c("red","green","yellow","cyan", "blue")))  
barplot(Bs, col=col1(nrow(Bs)) , border="white", space=0.04, cex.axis=1,
xlab="Year",ylab='',main='')
legend("topright", rownames(Bs), fill = col1(nrow(Bs)), bty = "n",
xjust=1,inset=c(-0.09,0),xpd=NA,cex=0.7,x.intersp=0.1)
mtext(side=3,line=0,"Biomass at age",las=0)

dev.off()

## predicted and observed total catches;

ind = !is.na(Popest$aveF_46)
catch.tot =  data.frame(year=year[ind],obs = tmb.data$landings/yscale,
  pred = rep$landings_pred/yscale)
  #,std.res = rep$std_landings_resid)
catch.tot$pdiff = 100*(catch.tot$obs/catch.tot$pred-1)

ylim=range(catch.tot[,2:3])

jpeg(file='obs_pred_landings.jpeg',width=4,height=3,units='in',res=300)

par(mar=c(3,4,0.5,3),las=1,mgp=c(2,1,0))

plot(catch.tot$year,catch.tot$obs,type='l',las=1,xlab='Year',ylab='',lwd=2,ylim=ylim)
lines(catch.tot$year,catch.tot$pred,lty=1,type='l',lwd=2,col='red')
 
yname = c("Landings (000 tonnes)")
mtext(side=2,line=2.5,yname,las=0,cex=1.2)

par(new=T)
plot(catch.tot$year,catch.tot$pdiff,axes=F,xlab="",ylab="",type='l',lty=2,lwd=2,col=grey(0.6))
abline(h=0,col=grey(0.6))
axis(side=4,las=1,hadj=0)
mtext("Difference (%)",side=4,line=2,las=0) 

legend('topright',bty='n',
c('Reported','Predicted','Difference'),col=c('black','red',grey(0.6)),
lwd=c(2,2,2),lty=c(1,1,2),cex=0.7)

dev.off()

#jpeg(file='landings_resid_ts.jpeg',width=4,height=3,units='in',res=300)

#par(mar=c(3,4,0.5,3),las=1,mgp=c(2,1,0))

#plot(catch.tot$year,catch.tot$std.res,type='l',las=1,xlab='Year',ylab='',lwd=2)
#abline(h=0,lty=2,lwd=2)
 
#yname = c("Std Landings Residuals")
#mtext(side=2,line=3,yname,las=0,cex=1.2)   

#dev.off()


######################catch age comparison #########################################;

vec_func = function(mdat){
  vdat = melt(mdat,id=c("Year"))
  vdat$Age = as.numeric(substr(vdat$variable,4,5))
  vdat$variable=NULL
  vdat$index=vdat$value
  vdat$value=NULL
  return(vdat)
}

X=data.frame(tmb.data$log_C)
X$Year = year[ind]
X.vec = vec_func(X)
mdat = log(rep$EC)
colnames(mdat)=colnames(tmb.data$C) 
X=data.frame(mdat)
X$Year = year[ind]
X.vec$pred = vec_func(X)$index   
mdat = tmb.data$C_zero
colnames(mdat)=colnames(tmb.data$C) 
X=data.frame(mdat)
X$Year = year[ind] 
X.vec$C_zero = vec_func(X)$index   
mdat = rep$std_C_resid
colnames(mdat)=colnames(tmb.data$C) 
X=data.frame(mdat)
X$Year = year[ind] 
X.vec$std_C_resid = vec_func(X)$index


jpeg(file='obs_pred_catch_age.jpeg',width=5,height=5,units='in',res=300)

ttl <- list(c("Catch Composition, Ages 2-12"), cex=1)
  xttl <- list(c("Year"), cex=1)
  yttl <- list(c("Log Number (10^6)"), cex=1)
  stripttl <- list(cex=0.2,font=1)
  ax <- list(cex=0.5,relation="free",rot=0)    
ax <- list(cex=0.5,tck=0.3,y=list(rot=0,relation=c("free")))
  my.key <- list(space = "top",columns=2,
                border = FALSE,
                size = 5,
                lines = list(lty = 1,col=c('red','black'),lwd=2),
                text = list(c("Model","Sampled")))
  
  p1=xyplot(index+pred~Year|factor(Age), data=X.vec, 
  main=ttl, xlab=xttl,lwd=c(2,2),type=c('l'),
  ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(3,4,1),las=1,
  key=my.key,col=c('black','red'),par.strip.text = list(cex=0.6), fill.color = X.vec$color,   
  par.settings = my.padding)
  
  p2=xyplot(index~Year|factor(Age), data=X.vec,groups = factor(C_zero), 
  main=ttl, xlab=xttl,lwd=c(2,2),type=c('p'),pch=19,
  ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(3,4,1),las=1,
  key=my.key,col=c('white','grey'),par.strip.text = list(cex=0.6),   
  par.settings = my.padding)
       
  print(p1 + as.layer(p2) + as.layer(p1))
  
dev.off()

jpeg(file='std_resid_catch_age.jpeg',width=5,height=5,units='in',res=300)

ttl <- list(c("Catch Composition, Ages 2-12"), cex=1)
  xttl <- list(c("Year"), cex=1)
  yttl <- list(c("Standardized Residuals"), cex=1)
  stripttl <- list(cex=0.2,font=1)
  
  p1=xyplot(std_C_resid~Year|factor(Age), data=X.vec,cex=0.5, 
  main=ttl, xlab=xttl,lwd=c(2,2),type=c('p'),pch=19,group=C_zero,
  ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(3,4,1),las=1,
  col=c('black','grey'),par.strip.text = list(cex=0.6),   
  par.settings = my.padding,
    panel = function(...) {
        panel.abline(h=c(0), lty = 2, col = "black",lwd=1)
        panel.xyplot(...)})
       
  print(p1)
  
dev.off()

jpeg(file='catch_age_err_hist.jpeg',width=3,height=3,units='in',res=300)

ind = X.vec$C_zero==0
perc.err = 100*(exp(X.vec$index[ind])-exp(X.vec$pred[ind]))/exp(X.vec$pred[ind])

par(mar=c(3,3,0.5,0.5),las=1,mgp=c(2,1,0))
hist(perc.err,xlab='Relative Difference (%)',breaks=20,main='')
legend("topright",paste('n =',length(perc.err)),bty='n')

dev.off() 

jpeg(file='catch_age_resid_qqplot.jpeg',width=3,height=3,units='in',res=300)

ind = X.vec$C_zero==0

par(mar=c(3,3,2,0.5),las=1,mgp=c(2,1,0))
qqplot(qqnorm(X.vec$std_C_resid[ind]))
abline(a=0,b=1,col='red')

dev.off()

jpeg(file='catch_age_obs_pred.jpeg',width=3,height=3,units='in',res=300)

par(mar=c(3,3,0.5,0.5),las=1,mgp=c(2,0.5,0))

plot(exp(X.vec$pred[ind]),exp(X.vec$index[ind]),
xlab='Predicted (10^6)',ylab='Observed (10^6)')
abline(a=0,b=1,col='grey',lwd=2)

dev.off()

cfit = X.vec
cfit$resid = cfit$std_C_resid
cfit$colr = 'deepskyblue'
cfit$colr[cfit$resid>0]='firebrick2'  
cfit$colr[(cfit$C_zero==1)&(cfit$resid>0)]='azure2'
cfit$colr[(cfit$C_zero==1)&(cfit$resid<=0)]='dimgray'

my.padding$layout.heights$main.key.padding = 0

jDarkGray <- 'grey20'
jPch <- 21
resid.plot1 = xyplot(factor(Age)~Year , data=cfit, 
       ylab='Age',
       cex = 2*sqrt(abs(cfit$resid)/pi), fill.color = cfit$colr,
       col = jDarkGray,
       par.settings = my.padding,
       par.strip.text = list(cex=0.5),
       panel = function(x, y, ..., cex, fill.color, subscripts) { 
         panel.abline(v=seq(1960,2015,by=5), col.line=grey(0.9))
         panel.xyplot(x, y, cex = cex[subscripts],
                      pch = jPch, fill = fill.color[subscripts], ...)
         })
         
print(resid.plot1)
                  
gname = c("C@A_resid_matrix_bubbles.jpeg")
jpeg(file=gname,width=4,height=3,units='in',res=300)
  print(resid.plot1)
dev.off()
catch = subset(cfit,C_zero==0)

jpeg(file='catch_resid_4P.jpeg',width=4,height=5,units='in',res=300)

par(mfcol=c(4,1),oma=c(1,3,1,1),mar=c(3,2,0,2),cex.axis=1.2,cex.lab=1.2,las=1)

x <- catch$Year
y <- catch$std_C_resid

catch$cohort = catch$Year-catch$Age
uage = sort(unique(catch$Age))

plot(x,y,xlab="",ylab="",type='n')
text(x,y,catch$Age,cex=0.5)
mres <- tapply(y,x,"mean")
lines(sort(unique(x)),mres,lty=1,lwd=3,col='red') 
abline(h=0,lty=1)
mtext(side=4,line=0.5,outer=F,'Year',las=0)

x=catch$cohort
plot(x,y,xlab="",ylab="",pch=3)
abline(h=0,lty=1)     
mres <- tapply(y,x,"mean")
lines(sort(unique(x)),mres,lty=1,lwd=3,col='red')
mtext(side=4,line=0.5,outer=F,'Cohort',las=0)

plot(catch$Age,y,xlab="",ylab="",pch=3)
abline(h=0,lty=2)
mtext(side=4,line=0.5,outer=F,'Age',las=0)
lines(uage,tapply(y,catch$Age,mean),lty=1)

plot(catch$pred,y,xlab="",ylab="",pch=3)
abline(h=0,lty=1)
mtext(side=4,line=0.5,outer=F,'Expected',las=0)

ytext <- c("Standardized residuals")
mtext(side=2,ytext,line=1,outer=T,las=0,cex=1.2)

dev.off()

mdat = rep$log_F
colnames(mdat)=colnames(tmb.data$C) 
X=data.frame(mdat)      
X$Year = year[1:nrow(tmb.data$C)] 
age.res = vec_func(X)
names(age.res) = c('year','age','F')

ylim = c(min(age.res$F),max(age.res$F))
ylim[1]=-4

jpeg(file='log_Fa.jpeg',width=4,height=4,units='in',res=300)

  xttl <- list(c("Year"), cex=1)
  yttl <- list(c("Log Fishing Mortality rates"), cex=1)
  stripttl <- list(cex=0.2,font=1)
  
  ax <- list(cex=0.5,relation="free",tck=0.3,rot=0)
  print(xyplot(F~year|factor(age), data=age.res, type="l", xlab=xttl,lwd=2,lty=1,
  ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(3,4,1),
  par.strip.text = list(cex=0.6),ylim=ylim,col=1,
  par.settings = my.padding,
    panel = function(...) {
        panel.abline(h=log(0.5), lty = 1, col = "grey",lwd=1)
        panel.xyplot(...)
}))
  
dev.off() 

mdat = rep$F
colnames(mdat)=colnames(tmb.data$C) 
X=data.frame(mdat)      
X$Year = year[1:nrow(tmb.data$C)] 
age.res = vec_func(X)
names(age.res) = c('year','age','F')

ylim = c(min(age.res$F),max(age.res$F))

jpeg(file='Fa.jpeg',width=4,height=4,units='in',res=300)

  xttl <- list(c("Year"), cex=1)
  yttl <- list(c("Fishing Mortality rates"), cex=1)
  stripttl <- list(cex=0.2,font=1)
  
  ax <- list(cex=0.5,relation="free",tck=0.3,rot=0)
  print(xyplot(F~year|factor(age), data=age.res, type="l", xlab=xttl,lwd=2,lty=1,
  ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(3,4,1),
  par.strip.text = list(cex=0.6),ylim=ylim,col=1,
  par.settings = my.padding,
    panel = function(...) {
        panel.abline(h=c(0.25,0.5,0.75), lty = 1, col = "grey",lwd=1)
        panel.xyplot(...)
}))
  
dev.off() 

bp<-function(x,y,v, scale=3, ...){
  plot(x,y,cex=sqrt(abs(v))*scale, col=ifelse(v<0,'black','grey'), pch=ifelse(v<0,16,16), ...)
  points(x[v>0],y[v>0],cex=sqrt(v[v>0])*scale, col='grey', pch=1, ...)
  points(x[v<0],y[v<0],cex=sqrt(abs(v[v<0]))*scale, col='black', pch=1, ...)
}

n.age <- length(unique(age.res$age))
n.year <- length(unique(age.res$year))

Fm <- rep$F
pr <- age.res$F/rep(tapply(age.res$F,age.res$year,max),times=n.age)

Prm <- matrix(pr,ncol=n.age,byrow=FALSE) 
ave.pr <- tapply(pr,age.res$age,mean)

ave.Prm <- matrix(ave.pr,nrow=n.year,ncol=n.age,byrow=TRUE)
Prdiff <- Prm - ave.Prm;
                                
yval <- unique(age.res$age) 
xval <- unique(age.res$year) 

jpeg(file="pr_F.jpeg",width=5,height=6,units='in',res=300)

par(mfrow=c(2,1),oma=c(1,2,0,0),mar=c(3,1,1,0))
nf <- layout(matrix(c(1,2),2,1,byrow=TRUE), c(7), c(3,6), TRUE)

plot(yval,ave.pr/max(ave.pr),xlab='',ylab='',type='l',lwd=2,las=1,cex=0.8)
mtext(side=3,outer=F,line=0.1,'Average Selectivity Pattern',cex=0.8)
mtext(side=1,outer=F,line=1.8,'Age',cex=1)  

vPrdiff <- as.vector(Prdiff); # this is in age within year order;

scale <- 3
 bp(age.res$year,age.res$age,vPrdiff, ylim=range(age.res$age), xlim=range(age.res$year), 
       las=1, xlab='', ylab='', main="", scale=scale)
    points(age.res$year,age.res$age,cex=sqrt(abs(vPrdiff))*scale, col='black', pch=1)
mtext(side=2,outer=F,line=2,'Age',cex=1) 
mtext(side=1,outer=F,line=2.5,'Year',cex=1) 
mtext(side=3,outer=F,line=0.1,'Change in Selectivity Pattern',cex=0.8)

dev.off()

age.res$pr=pr

jpeg(file='pr_F_levelplot.jpeg',width=4,height=2,units='in',res=300)

par(mar=c(3,3,0.5,0.5),cex.axis=1.2,las=1)

plot1=levelplot(pr~year+age,data=age.res,col.regions = gray.colors(1000),par.settings = my.padding,main=list(c("Fishery Selectivity"), cex=1, font=6))
print(plot1)

dev.off()

############# Survey Fits ##############

fit = indices
sname = unique(fit$survey)
plot_type='jpeg'

fit$year=fit$Year
fit$age=fit$Age     
fit$Elog_index = rep$Elog_index 
fit$Eindex = exp(rep$Elog_index) 
fit$std.resid = rep$std_resid_index
fit$cohort=fit$year-fit$age
fit$est.wt=1
fit$log_index = log(fit$index)   

fit$sa = factor(paste(fit$survey,"_",fit$age,sep=''))
ox = rep(c(4:11,1:3),each=3) + rep(c(0,22,11),times=11)
fit$sa = factor(fit$sa,levels(fit$sa)[ox])

fit.all=fit

fit = subset(fit,i_zero==0)

tran.bias = 1
gname <- 'index.jpeg'
source("C:\\home\\CADIGAN\\Rfunc\\admb_plots\\obs_pred_plot.txt",local=FALSE)
dev.off()


uage = sort(unique(fit$age))
pname='';fig.folder=getwd()
source("C:\\home\\CADIGAN\\Rfunc\\admb_plots\\sres_plot.txt",local=FALSE)


my.key <- list(space = "top",columns=3,
          border = FALSE,
          size = 3,between=0.1,
          points = list(col=c("darkgoldenrod4","red","cyan2"),pch=19),
          text = list(usurvey))
          
my.padding$layout.heights$main.key.padding = 1

jpeg(file='index_fit.jpeg',width=4,height=4,units='in',res=300)

  xttl <- list(c("Predicted Index"), cex=1)
  yttl <- list(c("Observed Index"), cex=1)
  stripttl <- list(cex=0.2,font=1)
  
  ax <- list(cex=0.5,relation="free",tck=0.3,rot=0)
  print(xyplot(log_index~Elog_index|factor(age),group=survey,
  data=fit, type="p", xlab=xttl,pch=19,col=c("darkgoldenrod4","red","cyan2"),
  ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(3,4,1),
  par.strip.text = list(cex=0.6), key=my.key,cex=0.5,
  par.settings = my.padding,
    panel = function(...) {
        panel.abline(a=0,b=1, lty = 1, col = "grey",lwd=1)
        panel.xyplot(...)
}))
  
dev.off() 

my.panel.loess <- function(x, y, span = 0.9, degree = 1, ...) {
    loess.fit <- loess.smooth(x, y, span = span, dgree = degree )
    panel.lines(loess.fit$x, loess.fit$y,col=grey(0.3), ...)
}

jpeg(file='index_resid_ts.jpeg',width=4,height=6,units='in',res=300)

  xttl <- list(c("Year"), cex=1)
  yttl <- list(c("Standardized Residual"), cex=1)
  stripttl <- list(cex=0.2,font=1)
  xlim = range(year)
  ylim=range(fit$std.resid)
  
 # ax <- list(cex=0.5,relation="free",tck=0.3,rot=0)
  print(xyplot(std.resid~year|sa,
  data=fit, type="p", xlab=xttl,pch=19,xlim=xlim,ylim=ylim,
  ylab=yttl, strip.text=stripttl, as.table=TRUE, layout=c(3,11,1),
  par.strip.text = list(cex=0.6), cex=0.5,
  par.settings = my.padding,
    panel = function(...) {
        panel.abline(h=0, lty = 1, col = "grey",lwd=1)   
        my.panel.loess(...) 
        panel.xyplot(...)     
}))
  
dev.off() 

my.key <- list(space = "top",columns=2,
          border = FALSE,
          size = 3,between=0.1,
          lines = list(col=c("black","red"),lty=c(1,1),lwd=c(2,2)),
          text = list(c('Observed','Predicted')))
          
my.padding$layout.heights$main.key.padding = 1


jpeg(file='obs_pred_index.jpeg',width=4,height=6,units='in',res=300)

  xttl <- list(c("Year"), cex=1)
  yttl <- list(c("Observed/Predicted"), cex=1)
  stripttl <- list(cex=0.2,font=1)
  xlim = range(year)
  
 # ax <- list(cex=0.5,relation="free",tck=0.3,rot=0)
  print(xyplot(log_index+Elog_index~year|sa,                          
  data=fit.all, type=c('l','l'), xlab=xttl,lwd=2,xlim=xlim,col=c('black','red'),
  ylab=yttl, strip.text=stripttl, as.table=TRUE, layout=c(3,11,1),
  par.strip.text = list(cex=0.6), cex=0.5,key=my.key,
  par.settings = my.padding))
  
dev.off() 

  xttl <- list(c("Year Class"), cex=1)
  fit1 = subset(fit,(YC>=1965)&(YC<=2013))
  xlim=range(uage)
  ax <- list(cex=0.5,relation=c("same"))
  
jpeg(file='obs_pred_byYC_FRV.jpeg',width=5,height=5,units='in',res=300)
 
  print(xyplot(log_index+Elog_index~age|factor(cohort), data=subset(fit1,survey=='Fall'),xlim=xlim,
  type="l", xlab=xttl,lwd=2,
  ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(6,6,1),las=1,
  key=my.key,col=2:1,par.strip.text = list(cex=0.6),
  par.settings = list(layout.heights = list(strip = .8))))

dev.off()          
  
jpeg(file='obs_pred_byYC_SRV.jpeg',width=5,height=5,units='in',res=300)
 
  print(xyplot(log_index+Elog_index~age|factor(cohort), data=subset(fit1,survey=='Spring'),xlim=xlim,
  type="l", xlab=xttl,lwd=2,
  ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(7,6,1),las=1,
  key=my.key,col=2:1,par.strip.text = list(cex=0.6),
  par.settings = list(layout.heights = list(strip = .8))))

dev.off()

par.names = names(opt$par)
ind = substr(par.names,1,7) == 'log_std'
x = exp(opt$par[ind])   
names(x) = sub("log_std_","",names(x))
xn = names(x)
xn[xn=='index']=usurvey
names(x)=xn

jpeg(file='sd_barplot.jpeg',width=4,height=4,units='in',res=300)

par(mar=c(3,6,1,1),mgp=c(2,0.5,0)) # increase y-axis margin.

barplot(x,names.arg=names(x),horiz=TRUE,xlab='STD',las=1,main='STD')
abline(v=c(0.25,0.5,0.75,1),col='grey',lty=2)

dev.off()

if(tmb.data$use_pe==1){
mdat = rep$pe
colnames(mdat)=colnames(tmb.data$C) 
X=data.frame(mdat)      
X$Year = year 
age.res = vec_func(X)
names(age.res) = c('year','age','pe')

age.res$colr = 'deepskyblue'
age.res$colr[age.res$pe>0]='firebrick2' 
age.res$pch=21     

my.padding$layout.heights$main.key.padding = 0
#fit$pch[indices$i_zero==1]=25

resid.plot1 = xyplot(factor(age)~year, data=age.res, 
       ylab='Age',main='Process Errors',
     #  scales = list(y = 0:5,x=c(2000,2005,2010,2015)),
       cex = 5*sqrt(abs(age.res$pe)/pi), fill.color = age.res$colr,
       col = jDarkGray,
       par.settings = my.padding,
       par.strip.text = list(cex=0.5),
       panel = function(x, y, ..., cex, fill.color, pch, subscripts) { 
         panel.abline(v=seq(1960,2015,by=5), col.line=grey(0.9))
         panel.xyplot(x, y, cex = cex[subscripts],
                       pch=age.res$pch, fill = fill.color[subscripts], ...)
         })
         
#print(resid.plot1)
                  
gname = c("pe_matrix_bubbles.jpeg")
jpeg(file=gname,width=5,height=5,units='in',res=300)
  print(resid.plot1)
dev.off()

ylim = c(min(age.res$pe),max(age.res$pe))

jpeg(file='pe_ts.jpeg',width=4,height=4,units='in',res=300)

  xttl <- list(c("Year"), cex=1)
  yttl <- list(c("Process Errors"), cex=1)
  stripttl <- list(cex=0.2,font=1)
  
  ax <- list(cex=0.5,relation="free",tck=0.3,rot=0)
  print(xyplot(pe~year|factor(age), data=age.res, type="l", xlab=xttl,lwd=2,lty=1,
  ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(3,4,1),
  par.strip.text = list(cex=0.6),ylim=ylim,col=1,
  par.settings = my.padding,
    panel = function(...) {
        panel.abline(h=0, lty = 1, col = "grey",lwd=1)
        panel.xyplot(...)
}))
  
dev.off()  
}

if(tmb.data$use_cye==1){

ye = cbind(catch.tot$year,rep$cye)
colnames(ye) = c('Year',paste('Age',uage,sep=''))

dat.ye = vec_func(as.data.frame(ye))
dat.ye$colr = 'deepskyblue'
dat.ye$colr[dat.ye$index>0]='firebrick2'

resid.plot1 = xyplot(factor(Age)~Year , data=dat.ye, 
       ylab='Age',
       scales = list(y = 0:5,x=c(2000,2005,2010,2015)),
       cex = 2*sqrt(abs(dat.ye$index)/pi), fill.color = dat.ye$colr,
       col = jDarkGray,
       par.settings = my.padding,
       par.strip.text = list(cex=0.5),
       panel = function(x, y, ..., cex, fill.color, subscripts) { 
         panel.abline(h=1:6, v=c(2000,2005,2010,2015), col.line=grey(0.9))
         panel.xyplot(x, y, cex = cex[subscripts],
                      pch = jPch, fill = fill.color[subscripts], ...)
         })
         
#print(resid.plot1)
                  
gname = c("cye_matrix_bubbles.jpeg")
jpeg(file=gname,width=4,height=3,units='in',res=300)
  print(resid.plot1)
dev.off()

mean.ye = apply(rep$cye,1,mean)

gname = c("cye.jpeg")   
jpeg(file=gname,width=3,height=3,units='in',res=300)
 par(mar=c(3,3.5,0.5,0.5),las=1)
  plot(catch.tot$year,mean.ye,type='l',lwd=2,xlab='',ylab='')
  abline(h=0,lty=2,lwd=2) 
  mtext(side=2,line=2.5,'Year Effect Mean',las=0)       
  mtext(side=1,line=2,'Year')
         
dev.off()
}

####################  Table #################################;

vtnames = c('std_log_C','std_logF','std_cye',
'std_log_R','std_index','ar_logF_age','ar_logF_year','ar_pe_year','ar_pe_age','ar_cye_year')

ind = names(sd.rep$value) %in% vtnames

np=length(sd.rep$value[ind])

rnames = names(sd.rep$value[ind])

survey_name=unique(indices$survey)[unique(indices$isurvey)+1]

for(i in 1:length(vtnames)){
  ind1 = rnames==vtnames[i] 
  npi = length(rnames[ind1])
  iname=1:npi
  if(vtnames[i]=='std_index'){iname=survey_name}
  if(npi>1){rnames[ind1] = paste(vtnames[i],iname,sep='')}
}

out.tab = matrix(NA,np,2)
out.tab[,1] = sd.rep$value[ind] 
out.tab[,2] = 100*sd.rep$sd[ind]/sd.rep$value[ind]

colnames(out.tab) = c("Est","CV(%)")
rownames(out.tab)=rnames

t1 = 2*opt$objective+ 2*length(opt$par)
t2 = 2*opt$objective + log(length(tmb.data$iq))*length(opt$par)
tname = paste("Model results, nll = ",round(opt$objective,digits=2),
', AIC = ',round(t1,2),", BIC = ",round(t2,2),sep="") 

print(xtable(out.tab,digits=c(0,3,3),caption=tname,
align=c('l','r','r')),type='html',
  file='model_output.doc',caption.placement="top")





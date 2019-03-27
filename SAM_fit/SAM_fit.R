# SAM AT https://github.com/fishfollower/SAM 
#install.packages("Rcpp")
#devtools::install_github("fishfollower/SAM/stockassessment")

#install.packages("corplot")
#library(corplot)
rm(list=ls())
library(TMB);library(ggplot2);library(gridExtra);library(plyr);library(tidyverse);library(reshape2); library(dplyr);library(ggplot2)
library(lattice);library(stats4);library(MASS);library(car);library("xtable") 
library(stockassessment)

#load data
load("C:/Users/jchampag/Documents/GitHub/SSM3PS/SAM_fit/SAMdata_IO.Rdata")
#load("C:/Users/ChampagnatJ/Documents/GitHub/SSM3PS/SAM_fit/SAMdata_OFF.RData")

dat <- setup.sam.data(
  surveys=surveys,
  residual.fleet=cn, 
  prop.mature=mo, 
  stock.mean.weight=sw, 
  catch.mean.weight=cw, 
  dis.mean.weight=dw, 
  land.mean.weight=lw,
  prop.f=pf, 
  prop.m=pm, 
  natural.mortality=nm, 
  land.frac=lf)

#class(cn);class(mo);class(sw);class(cw);class(dw);
#class(lf);class(pf);class(pm);class(lw);class(nm);
#class(surveys);
#str(surveys)


##########
# Set up configuration
###########
conf <- defcon(dat)

#should be 2
conf$minAge

#should be 14
conf$maxAge

#last age is a plus group so should 1 (no plus group=0)
conf$maxAgePlusGroup

#coupling fishing mortality states, try varying this, if make them all indep then appears to overfit the catch
conf$keyLogFsta
#na <- conf$maxAge-conf$minAge+1 
#conf$keyLogFsta[1,]=rep(0,na) #same for all age

# Correlation of fishing mortality across ages (0 independent, 1 compound symmetry, or 2 AR(1)
conf$corFlag 

# Coupling of the survey catchability parameters (nomally first row is not used, as that is covered by fishing mortality).  
conf$keyLogFpar 

#density dep catchability parms
conf$keyQpow

# Coupling of process variance parameters for log(F)-process (nomally only first row is used)
conf$keyVarF

# Coupling of process variance parameters for log(N)-process
conf$keyVarLogN

#coupling of var parms for the obs
conf$keyVarObs

# Covariance structure for each fleet ("ID" independent, "AR" AR(1), or "US" for unstructured). | Possible values are: "ID" "AR" "US"
conf$obsCorStruct

conf$fbarRange <- c(4,6)

#######
##Fit
#######
par <- defpar(dat,conf)
fit <- sam.fit(dat,conf,par)

#######
##ouptput et plot
#######
setwd('C:/Users/jchampag/Documents/GitHub/SSM3PS/SAM_fit')
saveConf (conf , file ="model.cfg") 
save(fit,file = 'SAM_fit.RData')


fit1 <-fit

#global plot: SSB,R Fbar
#windows();plot(fit1)
jpeg(file='global_fit.jpeg',width=5,height=8,units='in',res=300) ;plot(fit1);dev.off()

##correlation between observation
#windows();corplot(fit1) #idem obscorrplot(fit1)
jpeg(file='correlation_fleet.jpeg',width=8,height=8,units='in',res=300) ;corplot(fit1);dev.off()

#one step ahead residual
#res1 <-residuals(fit1)
#windows();plot(res1)
jpeg(file='residual.step.ahead.jpeg',width=8,height=8,units='in',res=300) ;plot(res1);dev.off()

#join sample residual
#resp <- procres(fit)
#windows();plot(resp)
jpeg(file='residual.resample.jpeg',width=8,height=8,units='in',res=300) ;plot(res1);dev.off()

#catch plot
#windows();catchplot(fit1)
jpeg(file='catch_plot.jpeg',width=8,height=8,units='in',res=300) ;catchplot(fit1);dev.off()

#fleet fit to observation
#windows();fitplot(fit1)
jpeg(file='fit_fleet.jpeg',width=8,height=5,units='in',res=300) ;fitplot(fit1);dev.off()
attributes(surveys)
jpeg(file='fit_fleet_catch.jpeg',width=8,height=5,units='in',res=300) ;fitplot(fit1,fleets=1);dev.off()
jpeg(file='fit_fleet_GEAC.jpeg',width=6,height=5,units='in',res=300) ;fitplot(fit1,fleets = 2);dev.off()
jpeg(file='fit_fleet_ERHAPS.jpeg',width=6,height=5,units='in',res=300) ;fitplot(fit1,fleets = 3);dev.off()
jpeg(file='fit_fleet_IO.jpeg',width=6,height=5,units='in',res=300) ;fitplot(fit1,fleets = 4);dev.off()


##Fishing mortality
#windows();fbarplot(fit1)
jpeg(file='Fbar.jpeg',width=6,height=5,units='in',res=300) ;fbarplot(fit1);dev.off()
F1 <-as.data.frame(faytable(fit1))
F1$years <- as.numeric(rownames(F1))
F1_l <- melt(F1,id.vars='years',value.name = "F",variable.name = "age")
#windows();wireframe(F ~ years*age, data = F1_l,zoom=0.95,xlab = "Year", ylab = "Age",main = "F", drape = TRUE,colorkey = TRUE,screen = list(z = -60, x = -60))
jpeg(file='3D_F.jpeg',width=6,height=5,units='in',res=300) ;wireframe(F ~ years*age, data = F1_l,zoom=0.95,xlab = "Year", ylab = "Age",
                                    main = "F",drape = TRUE,colorkey = TRUE,screen = list(z = -60, x = -60));dev.off()
#windows();ggplot(F1_l)+geom_line(aes(years,F,col=age))
jpeg(file='F.jpeg',width=6,height=5,units='in',res=300) ;ggplot(F1_l)+geom_line(aes(years,F,col=age));dev.off()



####retrospective plot
ny.retro <- 5
#retro1 <- retro(fit1,ny.retro)
#windows();plot(retro1)
jpeg(file='retro_plot.jpeg',width=6,height=8,units='in',res=300) ;plot(retro1);dev.off()

##leave out fleet plot
#lo <- leaveout(fit1)
#windows();plot(lo)
jpeg(file='leaveout_plot.jpeg',width=6,height=8,units='in',res=300) ;plot(lo);dev.off()

#Quatchabilities
#windows();qtableplot(qtable(fit1))
jpeg(file='Q.jpeg',width=6,height=8,units='in',res=300) ;qtableplot(qtable(fit1));dev.off()

#parameter plot => bof
#windows();parplot(fit1)

#didn't those function to work
#window();ssbplot(fit1)
#window();recplot(fit1)
#coef.sam(fit1)

##save stuff
#coefficient
sd.rep <-sdreport(fit$obj)

write.table(summary(fit,digits=5),file='output.dat') 
write.table( modeltable(fit1),file='modelSummary.dat')
write.table( cbind(sd.rep$value,sd.rep$sd),file='raw_output.txt')



#######
#Model comparison
#######
#load different fit
fits <- c(old=fit1 , new = fit2 )
fits <- c(fit1,fit2,fit3 )
modeltable ( fits )
windows();plot(fits , addCI = TRUE )


# SAM AT https://github.com/fishfollower/SAM 
#install.packages("Rcpp")
#devtools::install_github("fishfollower/SAM/stockassessment")

#install.packages("corplot")
#library(corplot)

library(stockassessment)
## read NAFO assessment results

setwd("C:\\Users\\mroberts\\OneDrive - Memorial University of Newfoundland\\Marine Institute\\Classes\\FISH_6005\\Project 1\\Project 1\\SAM\\data_2")

#don't know if these are the best comparisons, be careful when examining plots with these
 assess.ssb = read.table("SSB.txt", header=F, sep = "\t", col.names=c("year", "ssb"))
 assess.ssb$ssb= assess.ssb$ssb/1000
 
 assess.F = read.table("F_bar.dat",header=F, sep = "\t", col.names=c("year",'F6_9','F4_6'))
 #assess.F<-assess.F[-c(57,58),]


#mdata doesn't show all data modifications because I wasn't sure how to code into the ICES format
 
#now includes zeros for age 1 fish, age 14 is a plus group (sum of 14,15,16)
cn <- read.ices("cn.txt")/1000
#age 1 has wgt 0.01, age 2 has 0.1 every year, 2015 and 2016 filled in with median from prior 3 years for every age 
cw <- read.ices("cw.txt")
dw <- read.ices("dw.txt")
lf <- read.ices("lf.txt")
lw <- read.ices("lw.txt")
#removed years outside of the range
mo <- read.ices("mo.txt")
nm <- read.ices("nm.dat")
pf <- read.ices("pf.txt")
pm <- read.ices("pm.txt")
#used DFO stock wts, added 0's for age 1
sw <- read.ices("sw.txt")
#survey.txt is both surveys, survey_OFF.txt is just OFF, survey_IO.txt is just IO
surveys <- read.ices("survey_OFF.txt")

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

conf <- defcon(dat)



#configuration

#should be 1
conf$minAge

#should be 14
conf$maxAge

#last age is a plus group so should be 14
conf$maxAgePlusGroup

#coupling fishing mortality states, try varying this, if make them all indep then appears to overfit the catch
conf$keyLogFsta[1,]=c(-1,-1,0,1,2,3,4,5,5,5,5,5,6,6)

# Correlation of fishing mortality across ages (0 independent, 1 compound symmetry, or 2 AR(1)
conf$corFlag <- 1

# Coupling of the survey catchability parameters (nomally first row is not used, as that is covered by fishing mortality).  
conf$keyLogFpar[2,]=c(-1,-1,0,1,1,1,1,1,1,1,1,1,1,1) 

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

par <- defpar(dat,conf)

#par$logFpar = log(rep(0.3,15))

fit <- sam.fit(dat,conf,par)


F<-faytable(fit)

plot(seq(from=1959, to=2016, by=1), F[,1], ylim=c(0,2), col="white")
for(i in 1:12){
  lines(seq(from=1959, to=2016, by=1), F[,i])
}

library(TMB)
sdreport(fit$obj)

tble = summary(fit,digits=3)

ass.years<-seq(from=1959, to=2016, by=1)

plot(ass.years,tble[,4]/tble[36,4], type="l", ylim=c(0,4), lwd=2)
abline(h=1, lty=2)
abline(h=2, lty=2)

write.table(ssbtable(fit),file='..//ssb_updated.dat') 
write.table(fbartable(fit),file='..//fbar_updated.dat') 


#par(mar=c(3,3,0.5,0.454))
ssbplot(fit)
#lines(assess.ssb$year,assess.ssb$ssb,lwd=2,lty=1,col='red')



cv=fit$sdrep$sd[names(fit$sdrep$value)=="logssb"]
plot(assess.F$year,cv)

fbarplot(fit)   
lines(assess.F$year,assess.F$F4_6,lwd=2,lty=1,col='red')

dev.off()
recplot(fit)

catchplot(fit)

#par(mar=c(2,2,0,0),mgp=c(2,1,0))
fitplot(fit)


res <- residuals(fit)
plot(res)

retro <- retro(fit,year=5)
plot(retro)


require(lattice)
library(latticeExtra)

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


qest = data.frame(survey=rep(c("IO","OFF"),each=14),age=rep(1:14,2))

temp = partable(fit)[1:17,]
t2 = conf$keyLogFpar
temp1 = rbind(temp[t2[2,]+1,],temp[t2[3,]+1,],temp[t2[4,]+1,])
qest = cbind(qest,temp1)[,c(1,2,5,6,7)]
colnames(qest) = c('survey','age','est','Low','High')

#jpeg(file='q.jpeg',width=4,height=4,units='in',res=300)

ylim = range(qest[,3:4])
#ylim[2]=4.5

xttl <- list(c("Age Class"), cex=1)
yttl <- list(c("Survey Q Pattern"), cex=1)
stripttl <- list(cex=0.1,font=1)
my.key <- list(space = "top",columns=3,
               border = FALSE,
               size = 3,between=0.1,
               lines = list(lty = c(1,1,1),col=c('black','grey','grey'),lwd=c(2,1,1)),
               text = list(c("Est","L","U")))

ax <- list(cex=0.5,relation="free",tck=0.3)
print(xyplot(est+Low+High~age|factor(survey), data=qest, type="l",
             xlab=xttl,lwd=c(2,1,1),lty=c(1,1,1),ylim=ylim,xlim=range(qest$age),
             ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(1,3,1),las=1,
             key=my.key,col=c('black','grey','grey'),par.strip.text = list(cex=0.4),
             par.settings = my.padding,
             panel = function(...) {
               panel.abline(h=1, lty = 2, col = "black",lwd=1)
               panel.xyplot(...)}))
#  par.settings = list(layout.heights = list(strip = .8))))

#dev.off()




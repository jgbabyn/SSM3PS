setwd("C:/Users/fhate/documents/Stock Assessment/Fish 6005/lectures/Project 1/data_EDA")
###
library(lattice)

#new.midy=1959:1982
midy.wt.all = read.table(file='midy_wt.dat') 
wtPredNew<- rep(0,24)
for(a in 2:17)
 {
 midy.wt=midy.wt.all[,a]
 midy.wt2 <- midy.wt[length(midy.wt):1]
# fit100=arima(midy.wt, order = c(1, 0, 0),  xreg = NULL, include.mean = TRUE, transform.pars = TRUE, fixed = NULL, init = NULL, method = "ML")
fit101=arima(midy.wt, order = c(1, 0, 1),  xreg = NULL, include.mean = TRUE, transform.pars = TRUE, fixed = NULL, init = NULL, method = "ML")
# fit111=arima(midy.wt, order = c(1, 1, 1),  xreg = NULL, include.mean = TRUE, transform.pars = TRUE, fixed = NULL, init = NULL, method = "ML")
# fit110=arima(midy.wt, order = c(1, 1, 1),  xreg = NULL, include.mean = TRUE, transform.pars = TRUE, fixed = NULL, init = NULL, method = "ML")
# fit001=arima(midy.wt, order = c(1, 1, 1),  xreg = NULL, include.mean = TRUE, transform.pars = TRUE, fixed = NULL, init = NULL, method = "ML")
# fit011=arima(midy.wt, order = c(1, 1, 1),  xreg = NULL, include.mean = TRUE, transform.pars = TRUE, fixed = NULL, init = NULL, method = "ML")
# fit010=arima(midy.wt, order = c(1, 1, 1),  xreg = NULL, include.mean = TRUE, transform.pars = TRUE, fixed = NULL, init = NULL, method = "ML")
# fit100$sigma2
# fit101$sigma2
# fit111$sigma2
# fit110$sigma2
# fit001$sigma2
# fit011$sigma2
# fit010$sigma2
# fit100$aic
# fit101$aic
# fit111$aic
# fit110$aic
# fit001$aic
# fit011$aic
# fit010$aic
 wt.age2=predict(fit101,n.ahead=24)
 wtPred <- wt.age2$pred

#  acf(fit101$residuals)
#  pacf(fit101$residuals)
 wtPredNew <- cbind(wtPredNew,wtPred[24:1])
}
wtPredNew<- wtPredNew[,-1]
wtPredNew<-cbind(1959:1982,wtPredNew)
colnames(wtPredNew)<-colnames(midy.wt.all)
midy_wt.all<-rbind(wtPredNew, midy.wt.all)
write.table(midy_wt.all,file='mid.year59-2016.txt')
#mid = read.table(file='mid.year59-82.txt')
midy.fit <-as.vector(as.matrix(midy_wt.all[,2:17]))

ax <- list(cex=0.5,relation="free",tck=0.3,rot=0)
xyplot(midy.fit ~rep(seq(from=1959, to=2016), times=16)|rep(seq(from=1, to=16), each=58), type="l", scales=ax)

###$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
midy.wt.all = read.table(file='midy_wt.dat',header=TRUE) 
new.wts<-data.frame(Year=seq(from=1959, to=1982, by=1))
new_wts<-cbind(new.wts, t(as.data.frame(apply(midy.wt.all[,c(2:17)],2,median))))
midywt.all<-rbind(new_wts, midy.wt.all)

midy.fit <-as.vector(as.matrix(midywt.all[,2:17]))

ax <- list(cex=0.5,relation="free",tck=0.3,rot=0)
xyplot(midy.fit ~rep(seq(from=1959, to=2016), times=16)|rep(seq(from=1, to=16), each=58), type="l", scales=ax)

#plot(midywt.all$Year,midywt.all$Age3)
####################################################33333
cnames = c('Year',paste('Age',2:16,sep=""))
catch = read.table(file='catch.dat',header=FALSE,col.names=cnames)
catch = catch[,3:14]/1000
#midywt.all<-midywt.all[,3:16]
midy.pred<-midy_wt.all[,4:15]
plandings.pred=apply(catch*midy.pred,1,sum)
#plandings.median=apply(catch*midywt.all,1,sum)
cnames = c('Year','landings')
landings = read.table(file='landings.dat',header=FALSE,col.names=cnames)
#landings$landings[59]=landings$landings[59]+landings$landings[58]   #just added two numbers of landings in 2016
landings=landings[c(-59,-60),]
olandings = landings$landings/1000
y1 = plandings.pred
#y3=plandings.median
y2 = olandings
year=1959:2016
ylim=range(y1,y2)
#ylim=range(y3,y2)

gname = "SoP22.jpeg"
jpeg(file=gname,width=4,height=4,units='in',res=300)    

par(mar=c(3,3.5,0.5,3),las=1)
plot(year,y1,type='l',lwd=2,xlab='',ylab='',ylim=ylim,xlim=range(year)) 
lines(year,y2,type='l',lty=1,lwd=2,col='red')
abline(v=1994.5,lty=1,lwd=3)
lines(year,y2,type='l',lty=1,lwd=2,col='red')
lines(year,y1,type='l',lty=1,lwd=2)
legend("topright",bty='n',legend=c('50x'),cex=0.75) 
legend(1960,260,bty='n',legend=c('Rep','SUM Xprod'),col=c('red','black'),lty=1,lwd=2,cex=0.75)

mtext(side=2,outer=F,line=2.5,'Landings (Kt)',las=0)
mtext(side=1,outer=F,line=2,'Year') 

par(new=T)
pdiff = 100*(olandings/plandings.pred-1)
plot(year,pdiff,axes=F,xlab="",ylab="",type='l',lty=1,lwd=2,col=grey(0.6))
abline(h=mean(pdiff),col=grey(0.6))
axis(side=4,las=1,hadj=0)
mtext("Difference (%)",side=4,line=2,las=0) 

#plot(plandings,olandings,xlab='',ylab='')
#abline(a=0,b=1)   
#mtext(side=2,outer=F,line=2.5,'Landings (Kt)',las=0) 
#mtext(side=1,outer=F,line=2,'Sum C@A*W@A')

dev.off()


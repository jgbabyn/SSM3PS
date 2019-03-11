

setwd("C:/Users/fhate/Documents/Stock Assessment/Fish 6005/lectures/Project 1/data_EDA")

IO<-read.table("RV_IO.dat",header=F, sep = "\t", col.names=c("year",paste(rep("age", times=16), seq(from=1, to=16, by=1), sep="")))
OFF<-read.table("RV_OFF.dat",header=F, sep = "\t", col.names=c("year",paste(rep("age", times=16), seq(from=1, to=16, by=1), sep="")))

library(lattice)
OFF.sub<-subset(OFF, year>=1997.3)


IO.ages<-IO[,-1]
IO.mean<-apply(IO.ages,2, mean, na.rm=TRUE)
IO.old<-IO.ages[-c(21,22),]
IO.mean.old<-apply(IO.old,2, mean, na.rm=TRUE)

OFF.ages<-OFF.sub[,-1]
OFF.mean<-apply(OFF.ages, 2, mean, na.rm=TRUE)
OFF.old<-OFF.ages[-c(21,22),]
OFF.mean.old<-apply(OFF.old,2, mean, na.rm=TRUE)

mixed<-IO.ages/OFF.ages

ages<-rep(colnames(mixed), each=22)
ages<-factor(ages, levels = paste(rep("age", times=16), seq(from=1, to=16, by=1), sep=""))

#update of Figure 22 from CSAS 2017/063
setwd("C:/Users/fhate/Documents/Stock Assessment/Fish 6005/lectures/Project 1/my files")

#jpeg("IO_OFF_mean.jpeg", width=480, height=480)
plot(seq(from=1, to=16, by=1), IO.mean.old/OFF.mean.old, type="l", xlab="Age", ylab="IO/OFF", lwd=2, col="red")
abline(h=1, lty=2)
lines(seq(from=1, to=16, by=1), IO.mean/OFF.mean, type="l",lwd=2, col="blue")
legend(x=10,y=1.25, legend=c("2018", "2016"), lty=1, lwd=2, col=c("blue","red"))
#dev.off()

#same thing but by ages
xyplot(as.vector(as.matrix(mixed))~rep(IO$year, times=16)|ages, type="l", lwd=2,xlab="Year", ylab="IO/OFF")




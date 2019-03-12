#install.packages("reshape2")
##install.packages('lubridate')
require(reshape2)
library(lattice)
library(stringr)
library(lubridate)

assess.year = 1959:2016

##Maybe use relative paths if at all...? Also can use / on Windows too 
##setwd("./data_EDA") 
setwd("C:\\Users\\mroberts\\OneDrive - Memorial University of Newfoundland\\Marine Institute\\Classes\\FISH_6005\\Project 1\\Project 1\\data_EDA")

cnames = c('Year','landings')
landings = read.table(file='landings.dat',header=FALSE,col.names=cnames)
#Apprarently lines 59 and 60 are copy paste mistakes, so we shouldn't use them
#landings[58,2]<-4751+666
landings<-landings[-c(59,60),]

cnames = c('Year',paste('Age',1:16,sep=""))
RV.IO.matrix = read.table(file='RV_IO.dat',header=FALSE,col.names=cnames)

RV.OFF.matrix = read.table(file='RV_OFF.dat',header=FALSE,col.names=cnames) 

#Let's not use age14 weights alone as plus group weight
off.oldages<-RV.OFF.matrix[,15:17] #take off plus group ages
off.oldages$sum<-off.oldages[,1]+off.oldages[,2]+off.oldages[,3]
off.oldages$prop14<-off.oldages[,1]/off.oldages$sum #find the proportion that each age contributes to the plus group
off.oldages$prop15<-off.oldages[,2]/off.oldages$sum 
off.oldages$prop16<-off.oldages[,3]/off.oldages$sum 
#off.oldages<-off.oldages[-c(13),] 
off.oldages<-as.matrix(off.oldages) #is.na functions don't like its current format...
off.oldages[is.nan(off.oldages)]<-0.33 #for NA years, assuming equal proportions for each age..maybe not the best assumption
off.oldages[is.na(off.oldages)]<-0.33
avg.prop_14<-median(off.oldages[,5], na.rm=TRUE) #turn those proportions into an average value so these calculations can be applied to the whole time series
avg.prop_15<-median(off.oldages[,6], na.rm=TRUE) #using median because not normally distributed
avg.prop_16<-median(off.oldages[,7], na.rm=TRUE)

prop_14<-c(rep(avg.prop_14, times=21),off.oldages[,5]) #create vector with averages for years without data and real values for years with data
prop_15<-c(rep(avg.prop_15, times=21),off.oldages[,6])
prop_16<-c(rep(avg.prop_16, times=21),off.oldages[,7])


RV.FRANCE.matrix = read.table(file='erhaps.txt',header=FALSE, sep=";") 
RV.FRANCE.matrix<-RV.FRANCE.matrix[-1,-1]
RV.FRANCE.matrix<-as.data.frame(t(RV.FRANCE.matrix))
RV.FRANCE.matrix$Year<-paste(seq(from=1978, to=1992, by=1), 3, sep=".")
colnames(RV.FRANCE.matrix)<- c(paste('Age',2:15,sep=""), "Year")

RV.GEAC.matrix = read.csv(file='geac.csv',header=TRUE) 
RV.GEAC.matrix$Year<-paste(RV.GEAC.matrix$Year, 12, sep=".")
RV.GEAC.matrix<-RV.GEAC.matrix[,c(1,3:18)]
colnames(RV.GEAC.matrix)<- c('Year',paste('Age',1:16,sep=""))

IO.ages<-RV.IO.matrix[,-1]
IO.mean<-apply(IO.ages,2, mean, na.rm=TRUE)
IO.old<-IO.ages[-c(21,22),]
IO.mean.old<-apply(IO.old,2, mean, na.rm=TRUE)

OFF.sub<-subset(RV.OFF.matrix, Year>=1997.3)
OFF.ages<-OFF.sub[,-1]
OFF.mean<-apply(OFF.ages, 2, mean, na.rm=TRUE)
OFF.old<-OFF.ages[-c(21,22),]
OFF.mean.old<-apply(OFF.old,2, mean, na.rm=TRUE)

cnames = c('Year','MedDate',paste0('Age',1:16))
geac <- read.csv(file="geac.csv",header=TRUE,col.names=cnames)
geac$fs <- lubridate::month(as.Date(geac$MedDate))/10 ##not actually the month, frac of year, so divide by 10 rather than 12

#IO.old[IO.old==0]<-0.00001
#OFF.old[OFF.old==0]<-0.00001
diff_matrix<-as.matrix(IO.old/OFF.old)
diff_matrix[is.nan(diff_matrix)] <- 1 #turning into 1 because both values are zero
diff_matrix[diff_matrix==Inf] <- 1 #turning the Inf into 1 because it was a difference of 0.01 to zero and most other years had a value of 1
diff_matrix<-diff_matrix[,2:16]

ratio<-IO.mean.old/OFF.mean.old

both_RV<-read.table("combined_RV.txt")
colnames(both_RV)[1]<-"Year"

cnames = c('Year',paste('Age',2:16,sep=""))
catch = read.table(file='catch.dat',header=FALSE,col.names=cnames) 
#lets try this with a plus group for age 14
catch[is.na(catch)] <- 0
catch[,14]<-catch[,14]+catch[,15]+catch[,16]
catch<-catch[,c(1:14)]


cnames = c('Year',paste('Age',1:16,sep=""))
mat.all = read.table(file='mat.txt',header=FALSE,col.names=cnames)


#doesnt go all the way back to 1959, for now will use median weights as historical extrapolation
wt.all = read.table(file='stock_wt.dat',header=TRUE)
new.wts<-data.frame(Year=seq(from=1959, to=1982, by=1))
new_wts<-cbind(new.wts, t(as.data.frame(apply(wt.all[,c(2:17)],2,median))))
wt.all<-rbind(new_wts, wt.all)

midy_wt.all = read.table(file='midy_wt.dat',header=TRUE) 
new.wts<-data.frame(Year=seq(from=1959, to=1982, by=1))
new_wts<-cbind(new.wts, t(as.data.frame(apply(midy_wt.all[,c(2:17)],2,median))))
midy_wt.all<-rbind(new_wts, midy_wt.all)


cnames = c('Year',paste('Age',3:14,sep=""))
comm_wt = read.csv(file='comm_wt_mod.csv',header=F) 
colnames(comm_wt)<-cnames

comm_wt2<-cbind(wt.all[,3], comm_wt[,2:13]) #adding on the stock weights for age 2's
comm_wt2<-cbind(comm_wt$Year, comm_wt2)
cnames = c('Year',paste('Age',2:14,sep=""))
colnames(comm_wt2)<-cnames


new_comm_wt<-matrix(ncol=14, nrow=58) #lets create smoothed comm_wts
new_comm_wt[,1]<-comm_wt$Year

for(i in 2:14){
  loess.fit<-loess(comm_wt2[,i]~assess.year, span=0.3) #span affects the alpha parameter (smoothing) could add degree to change the degree of polynomial (1,2)
  new_comm_wt[,i]<-predict(loess.fit) 
}
new_comm_wt<-as.data.frame(new_comm_wt)
colnames(new_comm_wt)<-cnames

plot(new_comm_wt$Age3~assess.year, type="l", col="blue", lwd=2)
lines(comm_wt$Age3~assess.year, lty=2)

midy_wts<-as.vector(as.matrix(midy_wt.all[,2:17]))

ax <- list(cex=0.5,relation="free",tck=0.3,rot=0)
xyplot(midy_wts~rep(seq(from=1959, to=2016), times=16)|rep(seq(from=1, to=16), each=58), type="l", scales=ax)

age<-2:14
A=length(age)

## subset for assessment
#14 is a plus group for RV surveys
RV.IO.matrix2<-RV.IO.matrix
RV.IO.matrix2[is.na(RV.IO.matrix2)] <- 0
RV.IO.matrix2[,15]<-RV.IO.matrix2[,15]+RV.IO.matrix2[,16]+RV.IO.matrix2[,17]
RV.IO.matrix[,15]<-RV.IO.matrix2[,15]
RV.IO.matrix[10,16]<-NA
RV.IO.matrix<-RV.IO.matrix[,c(1,3:15)]


RV.OFF.matrix2<-RV.OFF.matrix
RV.OFF.matrix2[is.na(RV.OFF.matrix2)] <- 0
RV.OFF.matrix2[,15]<-RV.OFF.matrix2[,15]+RV.OFF.matrix2[,16]+RV.OFF.matrix2[,17]
RV.OFF.matrix[,15]<-RV.OFF.matrix2[,15]
RV.OFF.matrix[25,16]<-NA
RV.OFF.matrix<-RV.OFF.matrix[,c(1,3:15)]

RV.FRANCE.matrix2<-RV.FRANCE.matrix
#RV.FRANCE.matrix2[is.na(RV.FRANCE.matrix2)] <- 0
RV.FRANCE.matrix2[,13]<-RV.FRANCE.matrix2[,13]+RV.FRANCE.matrix2[,14]
RV.FRANCE.matrix[,13]<-RV.FRANCE.matrix2[,13]
#RV.FRANCE.matrix[RV.FRANCE.matrix==0]<-NA
RV.FRANCE.matrix<-RV.FRANCE.matrix[,c(1:13,15)]

RV.GEAC.matrix2<-RV.GEAC.matrix
RV.GEAC.matrix2[is.na(RV.GEAC.matrix2)] <- 0
RV.GEAC.matrix2[,15]<-RV.GEAC.matrix2[,15]+RV.GEAC.matrix2[,16]+RV.GEAC.matrix2[,17]
RV.GEAC.matrix[,15]<-RV.GEAC.matrix2[,15]
RV.GEAC.matrix[RV.GEAC.matrix==0]<-NA
RV.GEAC.matrix<-RV.GEAC.matrix[,c(1,3:15)]



#remove duplicates for now from RV.OFF
#RV.OFF.matrix<-RV.OFF.matrix[-c(12),]

#remove extra years
RV.OFF.matrix<-RV.OFF.matrix[1:35,]


landings = subset(landings,Year %in% assess.year)

mat = subset(mat.all,mat.all[,1] %in% seq(from=1959, to=2016, by=1))
mat = mat[,c(3:15)]

wt = subset(wt.all,wt.all[,1] %in% seq(from=1959, to=2016, by=1))
wt$Age14plus<-prop_14*wt$Age14+prop_15*wt$Age15+prop_16*wt$Age16
wt = wt[,c(3:14,18)]              
midy_wt = subset(midy_wt.all,midy_wt.all[,1] %in% seq(from=1959, to=2016, by=1))
midy_wt$Age14plus<-prop_14*midy_wt$Age14+prop_15*midy_wt$Age15+prop_16*midy_wt$Age16
midy_wt = midy_wt[,c(3:14,18)]   

colnames(wt)=colnames(mat) 
colnames(midy_wt)=colnames(mat)

year = catch[,1]
catch = catch[,2:14]/1000

# make vector data frames

vec_func = function(mdat){
  vdat = melt(mdat,id=c("Year"))
  vdat$Age = as.numeric(substr(vdat$variable,4,5))
  vdat$variable=NULL
  vdat$index=vdat$value
  vdat$value=NULL
  return(vdat)
}

IO.vec = vec_func(RV.IO.matrix)
IO.vec$survey = 'IO'     

#change the fractional yr column into yr and month columns
frac_yr<-strsplit(as.character(IO.vec$Year), "[.]")
IO.vec$Year<-as.numeric(sapply(frac_yr, "[[", 1))
IO.vec$month<-as.numeric(sapply(frac_yr, "[[", 2))

OFF.vec = vec_func(RV.OFF.matrix)
OFF.vec$survey = 'OFF'    

#change the fractional yr column into yr and month columns
frac_yr<-strsplit(as.character(OFF.vec$Year), "[.]")
OFF.vec$Year<-as.numeric(sapply(frac_yr, "[[", 1))
OFF.vec$month<-as.numeric(sapply(frac_yr, "[[", 2))

both.vec = vec_func(both_RV)
both.vec$survey = 'BOTH'    

#change the fractional yr column into yr and month columns
frac_yr<-strsplit(as.character(both.vec$Year), "[.]")
both.vec$Year<-as.numeric(sapply(frac_yr, "[[", 1))
both.vec$month<-as.numeric(sapply(frac_yr, "[[", 2))


france.vec = vec_func(RV.FRANCE.matrix)
france.vec$survey = 'FRANCE'    

#change the fractional yr column into yr and month columns
frac_yr<-strsplit(as.character(france.vec$Year), "[.]")
france.vec$Year<-as.numeric(sapply(frac_yr, "[[", 1))
france.vec$month<-as.numeric(sapply(frac_yr, "[[", 2))

xyplot(france.vec$index~france.vec$Year|france.vec$Age, type="l", scales=ax)

geac.vec = vec_func(RV.GEAC.matrix)
geac.vec$survey = 'GEAC'    

#change the fractional yr column into yr and month columns
frac_yr<-strsplit(as.character(geac.vec$Year), "[.]")
geac.vec$Year<-as.numeric(sapply(frac_yr, "[[", 1))
geac.vec$month<-as.numeric(sapply(frac_yr, "[[", 2))

xyplot(geac.vec$index~geac.vec$Year|geac.vec$Age, type="l", scales=ax)

#indices = rbind(IO.vec, OFF.vec)
#can switch this to include as much or as little data as desired
indices<-rbind(both.vec, france.vec, geac.vec)
indices$log_OFF<-c(log(OFF.vec$index), rep(NA, times=length(geac.vec$Year)+length(france.vec$Year)))


#########################################################

temp = cbind(year,catch)
colnames(temp) = c('Year',colnames(catch))
catch.vec = vec_func(temp) 

temp = cbind(year,wt)
colnames(temp) = c('Year',colnames(wt))
wt.vec = vec_func(temp)   

temp = cbind(year,midy_wt)
colnames(temp) = c('Year',colnames(midy_wt))
midy_wt.vec = vec_func(temp)  

temp = cbind(year,new_comm_wt)
colnames(temp) = c('Year',colnames(new_comm_wt))
comm_wt.vec = vec_func(temp)  



temp = cbind(c(year),mat)
colnames(temp) = c('Year',colnames(mat))
mat.vec = vec_func(temp)

#something seems off about the landings data
#landings<-landings[-c(59,60),]

#not really sure what people are using this data for. Add whatever you need
save(assess.year,catch,landings,mat,wt,midy_wt,comm_wt,indices,catch.vec,mat.vec,wt.vec,
     midy_wt.vec,comm_wt.vec,RV.OFF.matrix,RV.IO.matrix,file='3Ps.RData')

indices$fs=NA


indices$fs<-indices$month/10 #not actually a month just the fraction of year, so divide by 10 rather than 12

off.years<-rep(seq(from=1983, to=2016, by=1), times=13)
iyear = as.numeric(factor(mat.vec$Year))-1          
iage = as.numeric(factor(mat.vec$Age))-1

ind_surv<-indices$survey=="BOTH"
indices$irho<-0
indices$irho[ind_surv]<-1


ratio_vec<-rep(log(ratio[2:14]), each=15)
ratio_mat<-matrix(ratio_vec, ncol=13, nrow=15)
all_diffs<-rbind(ratio_mat, log(diff_matrix[,1:13]))

indices$diffs<-NA
indices$diffs[ind_surv]<-as.vector(all_diffs)


indices$iyear = iyear[match(indices$Year,mat.vec$Year)] 
indices$iage = iage[match(indices$Age,mat.vec$Age)]
indices$i_zero=0
ind = indices$index<0.005
indices$i_zero[ind]=1   
indices$index[ind]=0.005
indices = subset(indices,!is.na(index))

indices$isurvey = as.numeric(as.factor(indices$survey))-1
indices$qname = paste(indices$survey,":",str_pad(indices$Age , 2, pad = "0"),sep='')
#ind = indices$Age>=6   
#indices$qname[ind] = paste(indices$survey[ind],":6+",sep='')

#ind = (indices$Age<=5)&(indices$survey=="Fall")&(indices$Year<1995)
#indices$qname[ind] = paste(indices$survey[ind],"E:",indices$Age[ind],sep='')     
#ind = (indices$Age<=5)&(indices$survey=="Spring")&(indices$Year<1996)
#indices$qname[ind] = paste(indices$survey[ind],"E:",indices$Age[ind],sep='')

indices$iq = as.numeric(as.factor(indices$qname))-1 

A = ncol(mat) 
Y = nrow(mat)

M_matrix = matrix(0.2,nrow=Y,ncol=A,byrow=T)


plandings=apply(catch*midy_wt,1,sum, na.rm=TRUE)
olandings = landings$landings/1000
y1 = plandings
y2 = olandings
#ind=year>=1995
#y1[ind]=50*y1[ind] 
#y2[ind]=50*y2[ind]

ylim=range(y1,y2)


gname = "SoP.jpeg"
jpeg(file=gname,width=4,height=4,units='in',res=300)    

par(mar=c(3,3.5,0.5,3),las=1)
plot(year,y1,type='l',lwd=2,xlab='',ylab='',ylim=ylim,xlim=range(year)) 
lines(year,y2,type='l',lty=1,lwd=2,col='red')
abline(v=1994.5,lty=1,lwd=3)
lines(year[!ind],y2[!ind],type='l',lty=1,lwd=2,col='red')
lines(year[!ind],y1[!ind],type='l',lty=1,lwd=2)
legend("topright",bty='n',legend=c('50x'),cex=0.75) 
legend(1960,260,bty='n',legend=c('Rep','SUM Xprod'),col=c('red','black'),lty=1,lwd=2,cex=0.75)

mtext(side=2,outer=F,line=2.5,'Landings (Kt)',las=0)
mtext(side=1,outer=F,line=2,'Year') 

par(new=T)
pdiff = 100*(olandings/plandings-1)
plot(year,pdiff,axes=F,xlab="",ylab="",type='l',lty=1,lwd=2,col=grey(0.6))
abline(h=mean(pdiff),col=grey(0.6))
axis(side=4,las=1,hadj=0)
mtext("Difference (%)",side=4,line=2,las=0) 

#plot(plandings,olandings,xlab='',ylab='')
#abline(a=0,b=1)   
#mtext(side=2,outer=F,line=2.5,'Landings (Kt)',las=0) 
#mtext(side=1,outer=F,line=2,'Sum C@A*W@A')

dev.off()

#changed the dimensions here, because was converting to vector
C_zero = matrix(0,nrow=Y,ncol=A,byrow=T)
C_zero[catch<0.0005]=1
#C_zero[is.na(C_zero)] <- 0
catchp = catch
catchp[catch<0.0005]=0.0005



tmb.data = list(
  M = M_matrix[1:(Y),],
  weight = as.matrix(wt),  
  mat = as.matrix(mat), 
  midy_weight = as.matrix(midy_wt[1:(Y),]),
  comm_wt=as.matrix(new_comm_wt),
  C = as.matrix(catchp),                  
  log_C = log(as.matrix(catchp)),     
  landings = plandings,
  log_landings = log(plandings),
  index = indices$index,
  log_OFF = indices$log_OFF,
  i_zero = indices$i_zero, 
  C_zero = C_zero,   
  iyear = indices$iyear,           
  iage = indices$iage,    
  isurvey = indices$isurvey,  
  iq = indices$iq,
  irho= indices$irho,
  diffs = indices$diffs,
  fs = indices$fs,
  A = A,
  Y = Y
)

names(tmb.data$iq) = indices$qname            
names(tmb.data$isurvey) = indices$survey 

tmb.data$log_index = log(tmb.data$index) 

tmb.data$index_censor = 0;  
tmb.data$catch_censor = 0;  

setwd("C:\\Users\\mroberts\\OneDrive - Memorial University of Newfoundland\\Marine Institute\\Classes\\FISH_6005\\Project 1\\Project 1\\Working model\\SepF_std+comm_wts+rhoq+all_survs")

save(tmb.data,indices,file='tmb.RData')

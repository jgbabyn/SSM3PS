setwd("C:\\home\\CADIGAN\\GradProgram\\2019\\F6005\\project1\\data_EDA")

require(lattice)  
library(latticeExtra)
library('ggplot2') 
require(reshape2) 
library(stringr)

## First Some R code to collect Data

assess.year = 1959:2017
                  
cnames = c('Year','landings')
landings = read.table(file='landings.dat',header=FALSE,col.names=cnames)

cnames = c('Year',paste('Age',1:16,sep=""))       
RV_IO.matrix = read.table(file='RV_IO.dat',header=FALSE,col.names=cnames)
cnames = c('Year',paste('Age',1:16,sep=""))       
RV_OFF.matrix = read.table(file='RV_OFF.dat',header=FALSE,col.names=cnames)

cnames = c('Year',paste('Age',2:16,sep=""))
catch = read.table(file='catch.dat',header=FALSE,col.names=cnames) 

cnames = c('Year',paste('Age',1:16,sep=""))
mat.all = read.table(file='mat.txt',header=FALSE,col.names=cnames)

wt.all = read.table(file='stock_wt.dat',header=TRUE) 
midy_wt.all = read.table(file='midy_wt.dat',header=TRUE) 

mat = subset(mat.all,mat.all[,1] %in% assess.year)
mat = mat[,2:17]

#wt = subset(wt.all,wt.all[,1] %in% assess.year)
wt = wt.all[,2:17] 
## add weight for 2017;
wt = rbind(wt,wt[nrow(wt),])             
#midy_wt = subset(midy_wt.all,midy_wt.all[,1] %in% assess.year)
midy_wt = midy_wt.all[,2:17]

#colnames(wt)=colnames(mat) 
#colnames(midy_wt)=colnames(mat)

year = catch[,1]
catch = catch[,2:16]/1000

# make vector data frames

vec_func = function(mdat){
  vdat = melt(mdat,id=c("Year"))
  vdat$Age = as.numeric(substr(vdat$variable,4,5))
  vdat$variable=NULL
  vdat$index=vdat$value
  vdat$value=NULL
  return(vdat)
}

RV_IO.vec = vec_func(RV_IO.matrix)
RV_IO.vec$survey = 'RV Inshore+Offshore' 
RV_IO.vec$fs = RV_IO.vec$Year - floor(RV_IO.vec$Year)
RV_IO.vec$Year = floor(RV_IO.vec$Year)  
   
RV_OFF.vec = vec_func(RV_OFF.matrix)
RV_OFF.vec$survey = 'RV Offshore'     
RV_OFF.vec$fs = RV_OFF.vec$Year - floor(RV_OFF.vec$Year)
RV_OFF.vec$Year = floor(RV_OFF.vec$Year)  

indices = rbind(RV_IO.vec,RV_OFF.vec)

temp = cbind(year,catch)
colnames(temp) = c('Year',colnames(catch))
catch.vec = vec_func(temp) 

rv.year = subset(wt.all[,1],wt.all[,1] %in% assess.year)
 
temp = cbind(c(rv.year,max(rv.year)+1),wt)
colnames(temp) = c('Year',colnames(wt))
wt.vec = vec_func(temp)   

temp = cbind(rv.year,midy_wt)
colnames(temp) = c('Year',colnames(midy_wt))
midy_wt.vec = vec_func(temp)    

temp = cbind(c(year,max(year)+1),mat)
colnames(temp) = c('Year',colnames(mat))
mat.vec = vec_func(temp)

save(assess.year,rv.year,catch,mat,wt,midy_wt,indices,catch.vec,mat.vec,wt.vec,
midy_wt.vec,RV_IO.matrix,RV_OFF.matrix,file='3Ps.RData')

############ Plot the inputs  ####################;

my.padding <- list(layout.heights = list( 
                        top.padding = 0, 
                        main.key.padding = 0, 
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

year = sort(unique(catch.vec$Year))
age = sort(unique(catch.vec$Age)) 
Y=length(year)
A=length(age)

## make a spay plot

spay_dat = function(X){
 
 Y=nrow(X)
 A=ncol(X)-1
 Ap=A+1
 year=X[,1]
 Xp=X[,2:Ap]
 age = as.numeric(sub('Age','',colnames(X)[2:Ap]))

 Ph = Xp/matrix(apply(Xp,1,sum,na.rm=T),nrow=Y,ncol=A,byrow=F)
 Ph.dev = Ph - matrix(apply(Ph,2,mean,na.rm=T),nrow=Y,ncol=A,byrow=T)
 Ph.dev = Ph.dev/matrix(sqrt(apply(Ph,2,var,na.rm=T)),nrow=Y,ncol=A,byrow=T)

  dat = data.frame(age = as.vector(matrix(age,nrow=Y,ncol=A,byrow=T)),
  year = as.vector(matrix(year,nrow=Y,ncol=A,byrow=F)),
  N = unlist(Xp),spay=unlist(Ph.dev))
  dat$clr = sign(dat$spay)
  dat$cohort = dat$year-dat$age
  return(dat)
}

dat = spay_dat(cbind(year,catch))

png(file='catch_spay.png',width=5,height=5,units='in',res=300)

ggplot(dat, aes(x=year, y=age, color=factor(clr))) + 
    geom_point(aes(size = 1.5*abs(spay)),alpha=0.6, show.legend = FALSE) +
    scale_size_continuous(range = c(0,10)) + 
    theme(plot.margin = margin(2, 2, 2, 2)) +
    ggtitle("Catch SPAY") +
    scale_color_manual(values=c('blue','red'))+ 
    scale_y_discrete(name ="Age",limits=age)+ 
    scale_x_discrete(name ="year",limits=seq(1960,2020,by=5))
    
dev.off()

dat$symbol=19
dat$symbol[dat$N==0]=8
dat$size = 1.5*abs(dat$N) 
dat$size[dat$N==0]=2 

png(file='catch_bubbles.png',width=5,height=5,units='in',res=300)

ggplot(dat, aes(x=year, y=age)) + 
    geom_point(aes(size = size, shape=factor(symbol)),alpha=0.6, show.legend = FALSE) +
    scale_shape_manual(values=c(8,19))+
    scale_size_continuous(range = c(0,10)) +  
    scale_y_continuous(breaks=1:15)             + 
    scale_x_continuous(breaks=seq(1960,2015,by=5))
    
dev.off()

##plot inputs

xttl <- list(c("Year"), cex=1)
yttl <- list(c("Catch (millions)"), cex=1)
stripttl <- list(cex=0.5,font=1)
#ax <- list(cex=1,relation="free",tck=0.3,y=list(rot=0))
ax <- list(cex=0.5,tck=0.3,y=list(rot=0,relation=c("free")))
dat = catch.vec 
nr = ceiling(length(2:16)/3)

P1 = xyplot(index~Year|as.factor(Age),
data=dat, type=c("l"),lwd=2,xlab=xttl, scales=ax,
ylab=yttl, strip.text=stripttl, las=1,col=c('black'),
as.table=TRUE, layout=c(3,nr,1),
par.settings = my.padding,
par.strip.text = list(cex=0.6),
)
gname = "catch_ts.jpeg"
jpeg(file=gname,width=4,height=4,units='in',res=300)       
  print(P1)
dev.off()


## plot RV_IO

xttl <- list(c("Year"), cex=1)
yttl <- list(c("Inshore+Offshore RV Index (mnpt)"), cex=1)
stripttl <- list(cex=0.5,font=1)
ax <- list(cex=0.5,tck=0.3,y=list(rot=0,relation=c("free")))
dat = RV_IO.vec 
nr = ceiling(length(2:12)/3)

P1 = xyplot(index~Year|as.factor(Age),
data=dat, type=c("l"),lwd=2,xlab=xttl, scales=ax,
ylab=yttl, strip.text=stripttl, las=1,col=c('black'),
as.table=TRUE, layout=c(3,nr,1),
par.settings = my.padding,
par.strip.text = list(cex=0.5),
)
gname = "RV_IO_ts.jpeg"
jpeg(file=gname,width=4,height=4,units='in',res=300)       
  print(P1)
dev.off()

dat = spay_dat(RV_IO.matrix)

png(file='RV_IO_spay.png',width=5,height=5,units='in',res=300)

ggplot(dat, aes(x=year, y=age, color=factor(clr))) + 
    geom_point(aes(size = 1.5*abs(spay)),alpha=0.6, show.legend = FALSE) +
    scale_size_continuous(range = c(0,10)) + 
    theme(plot.margin = margin(2, 2, 2, 2)) + 
    ggtitle("Inshore+Offshore RV SPAY") +
    scale_color_manual(values=c('blue','red'))+ 
    scale_y_discrete(name ="Age",limits=age)+ 
    scale_x_discrete(name ="year",limits=seq(1983,2016,by=5))
    
dev.off()

dat$symbol=19
dat$symbol[dat$N==0]=8
dat$size = 1.5*abs(dat$N) 
dat$size[dat$N==0]=2 

png(file='RV_IO_bubbles.png',width=5,height=5,units='in',res=300)

ggplot(dat, aes(x=year, y=age)) + 
    geom_point(aes(size = size, shape=factor(symbol)),alpha=0.6, show.legend = FALSE) +
    scale_shape_manual(values=c(8,19))+
    scale_size_continuous(range = c(0,10)) +  
    scale_y_continuous(breaks=1:15)             + 
    scale_x_continuous(breaks=seq(1983,2016,by=5))
    
dev.off()

## plot Offshore Survey

xttl <- list(c("Year"), cex=1)
yttl <- list(c("Offshore RV Index (mnpt)"), cex=1)
stripttl <- list(cex=0.5,font=1)  
ax <- list(cex=0.5,tck=0.3,y=list(rot=0,relation=c("free")))
dat = RV_OFF.vec 
nr = ceiling(length(2:12)/3)

P1 = xyplot(index~Year|as.factor(Age),
data=dat, type=c("l"),lwd=2,xlab=xttl, scales=ax,
ylab=yttl, strip.text=stripttl, las=1,col=c('black'),
as.table=TRUE, layout=c(3,nr,1),
par.settings = my.padding,
par.strip.text = list(cex=0.5),
)
gname = "RV_OFF_ts.jpeg"
jpeg(file=gname,width=4,height=4,units='in',res=300)       
  print(P1)
dev.off()

dat = spay_dat(RV_OFF.matrix)

png(file='RV_OFF_spay.png',width=5,height=5,units='in',res=300)

ggplot(dat, aes(x=year, y=age, color=factor(clr))) + 
    geom_point(aes(size = 1.5*abs(spay)),alpha=0.6, show.legend = FALSE) +
    scale_size_continuous(range = c(0,10)) + 
    theme(plot.margin = margin(2, 2, 2, 2)) +  
    ggtitle("Offshore RV SPAY") +
    scale_color_manual(values=c('blue','red'))+ 
    scale_y_discrete(name ="Age",limits=age)+ 
    scale_x_discrete(name ="year",limits=seq(1983,2016,by=5))
    
dev.off() 

dat$symbol=19
dat$symbol[dat$N==0]=8
dat$size = 1.5*abs(dat$N) 
dat$size[dat$N==0]=2 

png(file='RV_OFF_bubbles.png',width=5,height=5,units='in',res=300)

ggplot(dat, aes(x=year, y=age)) + 
    geom_point(aes(size = size, shape=factor(symbol)),alpha=0.6, show.legend = FALSE) +
    scale_shape_manual(values=c(8,19))+
    scale_size_continuous(range = c(0,10)) +  
    scale_y_continuous(breaks=1:15)             + 
    scale_x_continuous(breaks=seq(1983,2016,by=5))
    
dev.off()


#Plot mats

xttl <- list(c("Year"), cex=1)
yttl <- list(c("% Mature"), cex=1)
stripttl <- list(cex=0.5,font=1) 
ax <- list(cex=0.5,tck=0.3,y=list(rot=0,relation=c("free")))
dat = mat.vec
dat$index = round(100*dat$index,digits=1) 
nr = ceiling(length(1:15)/3)

P1 = xyplot(index~Year|as.factor(Age),
data=dat, type=c("l"),lwd=2,xlab=xttl, scales=ax,
ylab=yttl, strip.text=stripttl, las=1,col=c('black'),
as.table=TRUE, layout=c(3,nr,1),
par.settings = my.padding,
par.strip.text = list(cex=0.5),
)
gname = "mat_ts.jpeg"
jpeg(file=gname,width=4,height=4,units='in',res=300)       
  print(P1)
dev.off()  

#Plot weights

xttl <- list(c("Year"), cex=1)
yttl <- list(c("Jan1 Body Weight (kg)"), cex=1)
stripttl <- list(cex=0.5,font=1)
ax <- list(cex=0.5,tck=0.3,y=list(rot=0,relation=c("free")))
dat = wt.vec
#dat$index = round(100*dat$index,digits=1) 
nr = ceiling(length(1:16)/4)

P1 = xyplot(index~Year|as.factor(Age),
data=dat, type=c("l"),lwd=2,xlab=xttl, scales=ax,
ylab=yttl, strip.text=stripttl, las=1,col=c('black'),
as.table=TRUE, layout=c(4,nr,1),
par.settings = my.padding,
par.strip.text = list(cex=0.6),
)
gname = "wt_ts.jpeg"
jpeg(file=gname,width=4,height=4,units='in',res=300)       
  print(P1)
dev.off()


xttl <- list(c("Year"), cex=1)
yttl <- list(c("Mid-year Body Weight (kg)"), cex=1)
stripttl <- list(cex=0.5,font=1)
ax <- list(cex=0.5,tck=0.3,y=list(rot=0,relation=c("free")))
dat = midy_wt.vec
#dat$index = round(100*dat$index,digits=1) 
nr = ceiling(length(1:16)/4)

P1 = xyplot(index~Year|as.factor(Age),
data=dat, type=c("l"),lwd=2,xlab=xttl, scales=ax,
ylab=yttl, strip.text=stripttl, las=1,col=c('black'),
as.table=TRUE, layout=c(4,nr,1),
par.settings = my.padding,
par.strip.text = list(cex=0.6),
)
gname = "midy_wt_ts.jpeg"
jpeg(file=gname,width=4,height=4,units='in',res=300)       
  print(P1)
dev.off()


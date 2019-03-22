# SAM AT https://github.com/fishfollower/SAM 
#install.packages("Rcpp")
#devtools::install_github("fishfollower/SAM/stockassessment")

rm(list=ls())
library(TMB);library(plyr);library(tidyverse);library(reshape2); library(dplyr)
library(stockassessment)

#load data => newest dataset from brach:SSM/SepF_std+comm_wts+rhoq+all_survs
load("C:/Users/ChampagnatJ/Documents/GitHub/SSM3PS/data_EDA/3PS.Rdata")
#load("C:/Users/jchampag/Documents/GitHub/SSM3PS/data_EDA/3PS.Rdata")

##plot of weight
windows();ggplot()+geom_line(data=wt.vec,aes(Year,index),col=1)+geom_line(data=midy_wt.vec,aes(Year,index),col=2)+
  geom_line(data=comm_wt.vec,aes(Year,index),col=3)+facet_wrap(~Age,scale='free')
windows();ggplot()+geom_line(data=catch.vec,aes(Year,index),col=1)+facet_wrap(~Age,scale='free')+
  geom_line(data=indices,aes(Year,index,col=survey))
windows();ggplot()+geom_line(data=catch.vec,aes(Year,index,col=factor(Age)))#+facet_wrap(~Age,scale='free')

###DATA SET UP IN SAM FORMAT

#stock weight
names(wt) <- substr(names(wt),4,nchar(names(wt)))
rownames(wt) <- assess.year
sw <- as.matrix(wt)

#catch number
names(catch) <- substr(names(catch),4,nchar(names(catch)))
rownames(catch) <- assess.year
cn <- as.matrix(catch)

#catch weight /!\ age 2 is missing
names(comm_wt) <- substr(names(comm_wt),4,nchar(names(comm_wt)))
rownames(comm_wt) <- assess.year
comm_wt$'2' <-midy_wt[,1] #using midy_sw for age2 for now, Jon extrapolation migth be usefull here
cw <- as.matrix(comm_wt[,c(14,2:13)])

#landing and discard table
dw <- cw #same than catch weight
lf <- cn;lf[] <-1 #landing fraction = 100%
lw <- cw #same than catch weight

#maturity
names(mat) <- substr(names(mat),4,nchar(names(mat)))
rownames(mat) <- assess.year
mo <- as.matrix(mat)

#natural mortality
nm <- cn;nm[] <-0.2

#proportion of mature and dead for ssb calculation = 1st Jan so 0
pf <-  cn;pf[] <-0
pm <-  cn;pm[] <-0



#surveys => annoying process, made it on excel and import it througth :read.ices
# ?/!\ two files can be used surveys-FR-GEAC-IO.txt OR surveys-FR-GEAC-OFF.txt depending what you want

#unique(indices$survey)
#french <- as.data.frame(dcast(indices[indices$survey == "FRANCE",],Year~Age,value.var = 'index'))
#rownames(french) <- french$Year
#geac <- as.data.frame(dcast(indices[indices$survey == "GEAC",],Year~Age,value.var = 'index'))
#rownames(geac) <- geac$Year
#both <- as.data.frame(dcast(indices[indices$survey == "BOTH",],Year~Age,value.var = 'index',fun=mean)) ##not working, probably more than 1 survey
#RV_IO <- as.data.frame(RV.IO.matrix)
#rownames(RV_IO) <- substr(as.character(RV_IO$Year),1,4)
#names(RV_IO) <- substr(names(RV_IO),4,nchar(names(RV_IO)))
#RV_OFF <- as.data.frame(RV.OFF.matrix)
#RV_OFF <- RV_OFF[RV_OFF$Year != '1993.3',] #two survey in 1993, choose to keep 1993.1 or 1993.3
#rownames(RV_OFF) <-  substr(as.character(RV_OFF$Year),1,4)
#names(RV_OFF) <- substr(names(RV_OFF),4,nchar(names(RV_OFF)))
#surveys1 <- list(french=french[,-1],geac=geac[,-1],off=RV_OFF[,-1])#or io=RV_IO[,-1]  )
#write.csv(both,'C:/Users/jchampag/Dropbox/IFREMER/3Ps_cod/Cours/Cours_Noel_Cadigan/F6005_winter_2019/Project1/both.csv')

#choose the survey indices you wanna use
surveys <- read.ices("C:/Users/ChampagnatJ/Documents/GitHub/SSM3PS/data_EDA/surveys-FR-GEAC-IO.txt");survey <- 'IO'
surveys <- read.ices("C:/Users/ChampagnatJ/Documents/GitHub/SSM3PS/data_EDA/surveys-FR-GEAC-OFF.txt");survey <- 'OFF'
#surveys2 <- read.ices("C:/Users/ChampagnatJ/Dropbox/IFREMER/3Ps_cod/Cours/Cours_Noel_Cadigan/F6005_winter_2019/Project1/surveys-FR-GEAC-OFF.txt")

##check consistency between all tables
dim(cn);dim(cw);dim(dw);dim(lf);dim(lw);dim(mo);dim(nm);dim(pf);dim(pm);dim(sw);dim(surveys)
dim(surveys$`GEAC													`);dim(surveys$`IO													`);dim(surveys$`ERHAPS													`);dim(surveys$`OFF													`)

rownames(cn);rownames(cw);rownames(dw);rownames(lf);rownames(lw);
rownames(mo);rownames(nm);rownames(pf);rownames(pm);rownames(sw)

#check if data are accecpted by SAM data set up
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


#saving in a Rdata file
save(surveys,cn,mo,sw,cw,dw,lw,pf,pm,nm,lf,file=paste("C:/Users/ChampagnatJ/Documents/GitHub/SSM3PS/SAM_fit/SAMdata_",survey,".Rdata",sep=""))

##Shouldn't be any real reason to run this...
##This is also just the kludgest cheapest things I could do to get it to run
##I only use 1982 onwards and ages 

##Reading this all in again so I have the bloody years...

setwd("~/Documents/GitHub/SSM3PS/data-raw/data_EDA")

possible.years = 1954:2019
possible.years = factor(possible.years)

cnames = c('Year','landings')
landings = read.table(file='landings.dat',header=FALSE,col.names=cnames)
landings$Year = factor(landings$Year,levels=levels(possible.years))

cnames = c('Year',paste0('Age',1:16))
RV_IO.matrix = read.table(file="RV_IO.dat",header=FALSE,col.names=cnames)
RV_IO.matrix$Year = as.character(RV_IO.matrix$Year)
yfs = do.call(rbind,strsplit(RV_IO.matrix$Year,"\\."))
RV_IO.matrix$Year = yfs[,1]
RV_IO.matrix$fs = as.numeric(yfs[,2])/10
RV_IO.matrix$Year = factor(RV_IO.matrix$Year,levels=levels(possible.years))

cnames = c('Year',paste0('Age',1:16))
RV_OFF.matrix = read.table(file="RV_OFF.dat",header=FALSE,col.names=cnames)
RV_OFF.matrix$Year = as.character(RV_OFF.matrix$Year)
yfs = do.call(rbind,strsplit(RV_OFF.matrix$Year,"\\."))
RV_OFF.matrix$Year = yfs[,1]
RV_OFF.matrix$fs = as.numeric(yfs[,2])/10
RV_OFF.matrix$Year = factor(RV_OFF.matrix$Year,levels=levels(possible.years))

cnames = c('Year',paste0('Age',1:16))
mat = read.table(file="mat.txt",header=FALSE,col.names=cnames)
mat$Year = factor(mat$Year,levels=levels(possible.years))

midy_wt = read.table(file='midy_wt.dat',header=TRUE)
midy_wt$Year = factor(midy_wt$Year,levels=levels(possible.years))

stock_wt = read.table(file='stock_wt.dat',header=TRUE)
stock_wt$Year = factor(stock_wt$Year,levels=levels(possible.years))

cnames = c('Year',paste0('Age',2:16))
catch = read.table(file='catch.dat',header=FALSE,col.names=cnames)
catch$Year = factor(catch$Year,levels=levels(possible.years))

used.years = midy_wt$Year

RV1 <- tidyr::gather(RV_OFF.matrix,key="Age",value="index",Age1:Age16)
RV1$Age = as.numeric(gsub("Age","",RV1$Age))
RV1$survey = "RV Offshore"
RV2 <- tidyr::gather(RV_IO.matrix,key="Age",value="index",Age1:Age16)
RV2$Age = as.numeric(gsub("Age","",RV2$Age))
RV2$survey = "RV Inshore+Offshore"
indices = rbind(RV1,RV2)
indices$survey = factor(indices$survey)

mat.vec <- tidyr::gather(mat,key="Age",value="maturity",Age1:Age16)
mat.vec$Age = as.numeric(gsub("Age","",mat.vec$Age))

catch.vec <- tidyr::gather(catch,key="Age",value="index",Age2:Age16)
catch.vec$Age = as.numeric(gsub("Age","",catch.vec$Age))

##I just can't

##Setting up to try with the fit.cpp...
## used.years = midy_wt$Year

## indices$iyear = as.numeric(indices$Year)-min(as.numeric(used.years))
## indices$iage = indices$Age-1
## indices$i_zero = 0
## ind = indices$index<0.005
## indices$i_zero[ind] =1
## indices$index[ind]=0.005
## indices = subset(indices,!is.na(index))
## indices$isurvey = as.numeric(as.factor(indices$survey))-1
## indices$qname = paste(indices$survey,":",stringr::str_pad(indices$Age , 2, pad = "0"),sep='')
## indices$iq = as.numeric(as.factor(indices$qname))-1

## A = length(2:14)
## Y = length(used.years)

## M_matrix = matrix(0.2,nrow=Y,ncol=A,byrow=TRUE)

## plandings = apply(catch[catch$Year %in% used.years,2:14]*midy_wt[midy_wt$Year %in% used.years,3:15],1,sum)

## C_zero = matrix(0,nrow=Y,ncol=A)
## C_zero[catch[catch$Year %in% used.years,2:14] < 0.0005] =1
## catchp = catch[catch$Year %in% used.years,2:14]
## catchp[catch[catch$Year %in% used.years,2:14]< 0.0005] = 0.0005

## tmb.data = list(
##     M = M_matrix[1:(Y-1),],
##     weight = as.matrix(stock_wt[stock_wt$Year %in% used.years,3:15]),
##     mat = as.matrix(mat[mat$Year %in% used.years,3:15]),
##     midy_weight = as.matrix(midy_wt[midy_wt$Year %in% used.years,3:15]),
##     C = as.matrix(catchp),
##     log_C = log(as.matrix(catchp)),
##     landings = plandings,
##     log_landings = log(plandings),
##     index = indices$index,
##     i_zero = indices$i_zero,
##     C_zero = C_zero,
##     iyear = indices$iyear,
##     iage = indices$iage,
##     isurvey = indices$isurvey,
##     iq = indices$iq,
##     fs = indices$fs,
##     A = A,
##     Y = Y
## )

## names(tmb.data$iq) = indices$qname            
## names(tmb.data$isurvey) = indices$survey 

## tmb.data$log_index = log(tmb.data$index) 


## tmb.data$iyear_index_ye = as.numeric(factor(paste(tmb.data$isurvey,":",tmb.data$iyear,sep='')))-1

## temp = as.numeric(factor(substr(sort(unique(paste(tmb.data$isurvey,":",tmb.data$iyear,sep=''))),1,1)))
## tmb.data$isurvey_index_ye = temp-1

## tmb.data$index_censor = 0;  
## tmb.data$catch_censor = 0;  

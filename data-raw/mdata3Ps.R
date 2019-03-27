library(tidyverse)
library(lubridate)

setwd("../data-raw/data_EDA")

##landings
cnames = c('Year','landings')
landings = read.table(file='landings.dat',header=FALSE,col.names=cnames) %>%
    group_by(Year) %>%
    summarize(landings=sum(landings))

##Asked Danny Ings about bounds on landings broke it down into periods
##"Medium bounds"
##1959-1993, somewhat wide
 landings$UB[landings$Year %in% 1959:1993] = 1.5
 landings$LB[landings$Year %in% 1959:1993] = 0.5
## ##1994-96, bycatch only, narrow
 landings$UB[landings$Year %in% 1994:1996] = 1.10
 landings$LB[landings$Year %in% 1994:1996] = 0.90
## ##1997-98, maybe 99, very wide, wider than 1959-93
 landings$UB[landings$Year %in% 1997:1999] = 1.75
 landings$LB[landings$Year %in% 1997:1999] = 0.35
## ##1999(or 2000) to present, narrower than 1959-1993 but wider than 94-96
 landings$UB[landings$Year %in% 2000:2016] = 1.25
 landings$LB[landings$Year %in% 2000:2016] = 0.75

##Inshore/offshore survey
cnames = c('Year',paste0('Age',1:16))
RV_IO.matrix = read.table(file="RV_IO.dat",header=FALSE,col.names=cnames)
RV_IO.matrix$Year = as.character(RV_IO.matrix$Year)
yfs = do.call(rbind,strsplit(RV_IO.matrix$Year,"\\."))
RV_IO.matrix$Year = yfs[,1]
RV_IO.matrix$fs = as.numeric(yfs[,2])/10
RV_IO.matrix$Year = as.numeric(RV_IO.matrix$Year)
##Remove 2006
RV_IO.matrix =RV_IO.matrix[RV_IO.matrix$Year != "2006",]

##Offshore survey
cnames = c('Year',paste0('Age',1:16))
RV_OFF.matrix = read.table(file="RV_OFF.dat",header=FALSE,col.names=cnames)
RV_OFF.matrix$Year = as.character(RV_OFF.matrix$Year)
yfs = do.call(rbind,strsplit(RV_OFF.matrix$Year,"\\."))
RV_OFF.matrix$Year = yfs[,1]
RV_OFF.matrix$fs = as.numeric(yfs[,2])/10
RV_OFF.matrix$Year = as.numeric(RV_OFF.matrix$Year)
##Remove RV_OFF after 1996 since RV_IO is extension
RV_OFF.matrixR = RV_OFF.matrix[RV_OFF.matrix$Year < "1996",]


##Read in the GEAC survey data
cnames = c('Year','MedDate',paste0('Age',1:16))
geac = read.csv("geac.csv",header=TRUE,col.names=cnames)
geac$fs <- lubridate::month(as.Date(geac$MedDate))/12 ##Yeah I know...
geac = geac[,-2]
##Remove 1997 and 2007
geacR = geac[geac$Year != "1997",]
geacR = geacR[geacR$Year != "2007",]

##Read in the French survey data
frenchS = read.csv("ERHAPS_indices_CAN.csv",sep=";")
frenchSur = data.frame(Year = 1978:1992)
frenchSur = cbind(frenchSur,t(frenchS)[-1,])
names(frenchSur) = c('Year',paste0('Age',1:14))
rownames(frenchSur) = NULL

frenchDates = read.csv("ERHAPS_res.csv")
frenchDates$Start.date = as.Date(frenchDates$Start.date)
frenchDates$End.date = as.Date(frenchDates$End.date)
frenchDates$midDate = frenchDates$Start.date + ((frenchDates$End.date - frenchDates$Start.date)/2)
frenchDates$fs = lubridate::month(frenchDates$midDate)/12
frenchDates = frenchDates[,c("Year","fs")]
frenchSur = dplyr::inner_join(frenchSur,frenchDates)

cnames = c('Year',paste0('Age',1:16))
mat = read.table(file="mat.txt",header=FALSE,col.names=cnames)

midy_wt = read.table(file='midy_wt.dat',header=TRUE)

comm_wt = read.table(file='comm_wt.txt',header=TRUE)

cnames = c('Year',paste0('Age',2:14))
DFO_stock_wt = read.table(file="DFO_stock_wt.txt",header=TRUE,col.names=cnames)

stock_wt = read.table(file='stock_wt.dat',header=TRUE)

cnames = c('Year',paste0('Age',2:16))
catch = read.table(file='catch.dat',header=FALSE,col.names=cnames)


catch$fs = NA
surveys = list(RVOff = RV_OFF.matrixR,RVInOff = RV_IO.matrix,geac=geacR,catch=catch,french=frenchSur)

##Read in FPCA files
fstock_wt = read.csv("fStock_wt.csv",header=TRUE)
fstock_wt = fstock_wt[,-1]
fcomm_wt = read.csv("fComm_wt.csv",header=TRUE)
fcomm_wt[,-1] = fcomm_wt[,-1]
fmat = read.csv("fMat.csv",header=TRUE)
fmat[,-1] = fmat[,-1]

rDat = list(surveys=surveys,fstock_wt=fstock_wt,fcomm_wt=fcomm_wt,fmat=fmat,
            catch=catch,DFO_stock_wt=DFO_stock_wt,stock_wt=stock_wt,midy_wt=midy_wt,mat=mat,
            RVOff=RV_OFF.matrix,RVInOff=RV_IO.matrix,geac=geac,french=frenchSur,landings=landings,comm_wt=comm_wt)
usethis::use_data(rDat,overwrite=TRUE)


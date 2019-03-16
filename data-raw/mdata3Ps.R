library(tidyverse)
library(lubridate)
##Reading this all in again so I have the bloody years...

setwd("~/Documents/GitHub/SSM3PS/data-raw/data_EDA")

possible.years = 1954:2019
possible.years = factor(possible.years)

##landings
cnames = c('Year','landings')
landings = read.table(file='landings.dat',header=FALSE,col.names=cnames) %>%
    group_by(Year) %>%
    summarize(landings=sum(landings))
##Upper and lower multipliers for censored bounds of landings
landings$UB = rep(1.5,nrow(landings))
landings$LB = rep(0.75,nrow(landings))

##Inshore/offshore survey
cnames = c('Year',paste0('Age',1:16))
RV_IO.matrix = read.table(file="RV_IO.dat",header=FALSE,col.names=cnames)
RV_IO.matrix$Year = as.character(RV_IO.matrix$Year)
yfs = do.call(rbind,strsplit(RV_IO.matrix$Year,"\\."))
RV_IO.matrix$Year = yfs[,1]
RV_IO.matrix$fs = as.numeric(yfs[,2])/10
RV_IO.matrix$Year = as.numeric(RV_IO.matrix$Year)

##Offshore survey
cnames = c('Year',paste0('Age',1:16))
RV_OFF.matrix = read.table(file="RV_OFF.dat",header=FALSE,col.names=cnames)
RV_OFF.matrix$Year = as.character(RV_OFF.matrix$Year)
yfs = do.call(rbind,strsplit(RV_OFF.matrix$Year,"\\."))
RV_OFF.matrix$Year = yfs[,1]
RV_OFF.matrix$fs = as.numeric(yfs[,2])/10
RV_OFF.matrix$Year = as.numeric(RV_OFF.matrix$Year)

##Read in the GEAC survey data
cnames = c('Year','MedDate',paste0('Age',1:16))
geac = read.csv("geac.csv",header=TRUE,col.names=cnames)
geac$fs <- lubridate::month(as.Date(geac$MedDate))/12 ##Yeah I know...
geac = geac[,-2]

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

used.years = midy_wt$Year

catch$fs = NA
surveys = list(RVOff = RV_OFF.matrix,RVInOff = RV_IO.matrix,geac=geac,catch=catch,french=frenchSur)


##Weird old FDA tangent because
library(fdapace)
tcomm = t(comm_wt[comm_wt$Year %in% 1978:2016,-1])


ids = rep(1:ncol(tcomm),each=nrow(tcomm))
s <- 3:14
tv <- rep(s,ncol(tcomm))

fdaIn <- MakeFPCAInputs(IDs=ids,tVec=tv,tcomm)
cWClust <- FClust(fdaIn$Ly,fdaIn$Lt,k=3,optnsFPCA=list(FVEthreshold=0.9))
names(cWClust$cluster) = 1978:2014

tcomm = t(comm_wt[,-1])


ids = rep(1:ncol(tcomm),each=nrow(tcomm))
s <- 3:14
tv <- rep(s,ncol(tcomm))

fdaIn <- MakeFPCAInputs(IDs=ids,tVec=tv,tcomm)
FPCAWeights <- FPCA(fdaIn$Ly,fdaIn$Lt,optns=list(dataType="Sparse"))

##Instead of ages 3 to 14, I need 1 to 16
newGrid <- 1:16
nMu <- spline(FPCAWeights$workGrid,FPCAWeights$mu,xout=newGrid)$y
nPhi <- apply(FPCAWeights$phi,2,function(y) spline(FPCAWeights$workGrid,y,xout=newGrid)$y)
nFPCA <- list(cumFVE = FPCAWeights$cumFVE,xiEst=FPCAWeights$xiEst,mu=nMu,phi=nPhi)

##Try to impute ages 1,2,15,16 imputing at bounds is sort of gross, but the mean function is so simple this might actually
##be OKAY. Worth a shot anyways.

##Remove noise by extrapolating on the fitted values
ncomm = ConvertSupport(FPCAWeights$workGrid,3:14,phi=t(fitted(FPCAWeights,K=1)))
##Spots for missing ages
aM = rep(NA,ncol(ncomm))
ncomm = rbind(aM,ncomm)
ncomm = rbind(aM,ncomm)
ncomm = rbind(ncomm,aM)
ncomm = rbind(ncomm,aM)

##I don't think 2 year old fish weigh nothing
tstock = t(DFO_stock_wt[DFO_stock_wt$Year %in% 1978:2019,-(1:2)])

ids = rep(1:ncol(tstock),each=nrow(tstock))
s <- 3:14
tv <- rep(s,ncol(tstock))

fdaIn <- MakeFPCAInputs(IDs=ids,tVec=tv,tstock)
DSClust <- FClust(fdaIn$Ly,fdaIn$Lt,optnsFPCA=list(FVEthreshold=0.9))
names(DSClust$cluster) = 1978:2017

tstock = t(DFO_stock_wt[,-(1:2)])

ids = rep(1:ncol(tstock),each=nrow(tstock))
s <- 3:14
tv <- rep(s,ncol(tstock))

fdaIn <- MakeFPCAInputs(IDs=ids,tVec=tv,tstock)
FPCASWeights <- FPCA(fdaIn$Ly,fdaIn$Lt,optns=list(dataType="Sparse"))

##Repeat steps for imputing ages
nstock = ConvertSupport(FPCASWeights$workGrid,3:14,phi=t(fitted(FPCASWeights,K=1)))
##Spots for missing ages
aM = rep(NA,ncol(nstock))
nstock = rbind(aM,nstock)
nstock = rbind(aM,nstock)
nstock = rbind(nstock,aM)
nstock = rbind(nstock,aM)

nMu <- spline(FPCASWeights$workGrid,FPCASWeights$mu,xout=newGrid)$y
nPhi <- apply(FPCASWeights$phi,2,function(y) spline(FPCASWeights$workGrid,y,xout=newGrid)$y)
nSFPCA <- list(cumFVE = FPCASWeights$cumFVE,xiEst=FPCASWeights$xiEst,mu=nMu,phi=nPhi)


##LITERALLY DID NOT THINK ID BE USING THIS AGAIN
FPCAImpute <- function(FPCAobj,data,tau){
    L <- min(which(FPCAobj$cumFVE >= tau*100))
    impute <- matrix(0,nrow=nrow((FPCAobj$xiEst)),ncol=nrow(FPCAobj$phi))
    for(i in 1:length(FPCAobj$xiEst[,1])){
     impute[i,] =  FPCAobj$xiEst[i,1:L]%*%t(FPCAobj$phi[,1:L]) + FPCAobj$mu
    }
    for(i in 1:length(data))
        {
            if(is.na(data[i])){
                data[i] <- impute[i]
                }
        }
    data
}

fpcaStock_wt <- FPCAImpute(nSFPCA,t(nstock),0.9)
fpcaComm_wt <- FPCAImpute(nFPCA,t(ncomm),0.9)

fStock_wt <- data.frame(Year=DFO_stock_wt$Year)
fStock_wt <- cbind(fStock_wt,fpcaStock_wt)
names(fStock_wt) = c("Year",paste0("Age",1:16))

fComm_wt <- data.frame(Year=comm_wt$Year)
fComm_wt <- cbind(fComm_wt,fpcaComm_wt)
names(fComm_wt) = c("Year",paste0("Age",1:16))

##Now for Maturities!
##Squeeze to apply qlogis, so I can plogis on the return and ensure 0-1 bounds
mats = apply(mat[,-1],1,function(x) mice::squeeze(x,bounds=c(0.00000000001,0.99999999999)))
mats = qlogis(mats)

ids = rep(1:nrow(mats),each=ncol(mats))
s <- mat$Year
tv <- rep(s,nrow(mats))

fdaIn <- MakeFPCAInputs(IDs=ids,tVec=tv,t(mats))
FPCAMat <- FPCA(fdaIn$Ly,fdaIn$Lt,optns=list(methodMuCovEst="smooth",userBwMu=14))
fMaty <- fitted(FPCAMat,K=1)

fpMat <- plogis(fMaty)
fMat <- mat$Year
fMat <- cbind(fMat,t(fpMat))
fMat <- as.data.frame(fMat)

names(fMat) <- c("Year",paste0("Age",1:16))

matplot(t(fMat[fMat$Year %in% 1954:1978,-1]),type="l",col="blue")
matplot(t(fMat[fMat$Year %in% 1978:2019,-1]),type="l",col="red",add=TRUE)

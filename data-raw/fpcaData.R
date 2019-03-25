setwd("~/Documents/GitHub/SSM3PS/data-raw/data_EDA")
cnames = c('Year',paste0('Age',2:14))
DFO_stock_wt = read.table(file="DFO_stock_wt.txt",header=TRUE,col.names=cnames)
comm_wt = read.table(file='comm_wt.txt',header=TRUE)
cnames = c('Year',paste0('Age',1:16))
mat = read.table(file="mat.txt",header=FALSE,col.names=cnames)

##Weird old FDA tangent because
library(fdapace)

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

tcomm = t(comm_wt[comm_wt$Year %in% 1978:2016,-1])


ids = rep(1:ncol(tcomm),each=nrow(tcomm))
##Mid-year weights so everybody aged .5
s <- seq(3.5,14.5,by=1)
tv <- rep(s,ncol(tcomm))

fdaIn <- MakeFPCAInputs(IDs=ids,tVec=tv,tcomm)
cWClust <- FClust(fdaIn$Ly,fdaIn$Lt,k=3,optnsFPCA=list(FVEthreshold=0.9))
names(cWClust$cluster) = 1978:2014

tcomm = t(comm_wt[,-1])


ids = rep(1:ncol(tcomm),each=nrow(tcomm))
s <- seq(3.5,14.5,by=1)
tv <- rep(s,ncol(tcomm))

fdaIn <- MakeFPCAInputs(IDs=ids,tVec=tv,tcomm)
FPCAWeights <- FPCA(fdaIn$Ly,fdaIn$Lt,optns=list(dataType="Sparse"))

##Instead of ages 3.5 to 14.5, I need 1.5 to 16.5 for commercial weights and 1.0 to 16.0 for stock
newGrid <- seq(0.5,16.5,by=1)
nMu <- spline(FPCAWeights$workGrid,FPCAWeights$mu,xout=newGrid)$y
nPhi <- apply(FPCAWeights$phi,2,function(y) spline(FPCAWeights$workGrid,y,xout=newGrid)$y)
nFPCA <- list(cumFVE = FPCAWeights$cumFVE,xiEst=FPCAWeights$xiEst,mu=nMu,phi=nPhi)


##Remove noise by extrapolating on the fitted values
ncomm = ConvertSupport(FPCAWeights$workGrid,s,phi=t(fitted(FPCAWeights,K=1)))
##Spots for missing ages
aM = rep(NA,ncol(ncomm))
ncomm = rbind(aM,ncomm)
ncomm = rbind(aM,ncomm)
ncomm = rbind(aM,ncomm)
ncomm = rbind(ncomm,aM)
ncomm = rbind(ncomm,aM)

fpcaComm_wt <- FPCAImpute(nFPCA,t(ncomm),0.9)

fComm_wt <- data.frame(Year=comm_wt$Year)
fComm_wt <- cbind(fComm_wt,fpcaComm_wt)
names(fComm_wt) = c("Year",paste0("Age",newGrid))
scm = tidyr::gather(fComm_wt,key="Age",value="Weight",-Year)
scm$Age = as.numeric(gsub("Age","",scm$Age))

                     
fc2016 = apply(fComm_wt[fComm_wt$Year %in% 2012:2014,-1],2,median)
ffComm_wt = rbind(fComm_wt,c(2015,fc2016))
ffComm_wt = rbind(ffComm_wt,c(2016,fc2016))
ffComm_wt = ffComm_wt[,-2]
names(ffComm_wt) = c("Year",paste0("Age",1:16))
write.csv(ffComm_wt,"fComm_wt.csv")

##FakeSW are basically just a shift down the curve half a year
fakeSW = t(ConvertSupport(newGrid,1:16,phi=t(fpcaComm_wt)))
fakeSW = cbind(fComm_wt$Year,fakeSW)
fakeSW = as.data.frame(fakeSW)
names(fakeSW) = c("Year",paste0("Age",1:16))
fs2016 = apply(fakeSW[fakeSW$Year %in% 2012:2014,-1],2,median)
ffakeSW = rbind(fakeSW,c(2015,fs2016))
ffakeSW = rbind(ffakeSW,c(2016,fs2016))
write.csv(ffakeSW,"fStock_wt.csv")

ssw = tidyr::gather(fakeSW,"Age","Weight",-Year)
ssw$Age = as.numeric(gsub("Age","",ssw$Age))

#gplot <- ggplot(scm,aes(Age,Weight)) + geom_point(aes(Age,Weight,group=Year),col="red") + geom_point(aes(Age,Weight,group=Year),data=ssw,col="blue")

##Stock weights are supposed to be January weights, extrapolate back to Jan.
##Using FPCA on DFO stock weights odd, so try back-calculating this way 

##Try to impute ages 1,2,15,16 imputing at bounds is sort of gross, but the mean function is so simple this might actually
##be OKAY. Worth a shot anyways.


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




fpcaStock_wt <- FPCAImpute(nSFPCA,t(nstock),0.9)

fStock_wt <- data.frame(Year=DFO_stock_wt$Year)
fStock_wt <- cbind(fStock_wt,fpcaStock_wt)
names(fStock_wt) = c("Year",paste0("Age",1:16))




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
write.csv(fMat,file="fMat.csv")

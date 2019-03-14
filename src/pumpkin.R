library(TMB)

setwd("~/Documents/GitHub/SSM3PS/src")
compile("fit.cpp","-O0 -g")
#dyn.unload(dynlib("fit.cpp"))
dyn.load(dynlib("fit"))

library(SSM3PS)
source("../data-raw/mdata3PS.R")

dat <- datSetup(surveys,catch$catch,landings,stock_wt,
                midy_wt,mat,ages=3:16,years=1983:2015,plusGroup=12,match3NOdims=FALSE)

param <- paramSetup(dat)


obj <- MakeADFun(dat,param$param,random=c("log_Rec_dev","log_F","pe"),DLL="fit")

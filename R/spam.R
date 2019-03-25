#'Fit the SPAM model
#'
#' Basic wrapper so I get a class to make things like AIC eaiser.
#'
#' @param dat list of data and parameters from spamSetup
#'
#'
#' @useDynLib 'fit'
#' 
#' @export
#'
spamFit <- function(data,parameters,indices,cVec,random=NULL,map=list(),control=list(),silent=TRUE,...){

    obj <- TMB::MakeADFun(data,parameters,map=map,random=random,silent=silent,...)

    opt <- nlminb(obj$par,obj$fn,obj$gr,control=control)

    fit <- list()

    fit$dat <- data
    fit$parameters <- parameters
    fit$map = map
    fit$obj <- obj
    fit$opt <- opt
    fit$rep <- obj$report()
    fit$sdr <- TMB::sdreport(fit$obj)
    fit$indices = indices
    fit$cVec = cVec
    

    class(fit) <- "spam"
    print(opt$message)
    fit
}


#'Returns the loglikelihood of a spam object
#'
#' @param object the spam object to get the loglikelihood of
#' @export
logLik.spam <- function(object,...)
{
    val = -object$opt$objective
    attr(val,"df") = length(object$opt$par)
    attr(val,"nobs") = length(c(object$dat$log_index,object$dat$crlVec))
    class(val) = "logLik"
    val
}

#'Returns the one step ahead residuals of a spam object and some other stuff for plotting
#'
#' @param object the spam object
#' @param method the one step ahead residuals to use
#' @param ... options to give oneStepPredict
#'
#' @export
residuals.spam <- function(object,method="oneStepGaussianOffMode",...){
    ospIL = TMB::oneStepPredict(object$obj,"log_index","keep",discrete=FALSE,method=method,...)
    ospCRL = TMB::oneStepPredict(object$obj,"crlVec","keepCRL",discrete=FALSE,method=method,...)

    ##Add crud for plotting
    ospIL = cbind(ospIL,Age=object$indices$Age,survey=object$indices$survey,Year=object$indices$Year,EV=object$rep$Elog_index,
                  i_zero=object$indices$i_zero)
    ospCRL = cbind(ospCRL,Age=object$cVec$Age,survey=object$cVec$survey,Year=object$cVec$Year,EV=as.vector(t(object$rep$ECRL)),
                   i_zero=object$cVec$i_zero)
    ret = rbind(ospIL,ospCRL)
    ret
}



#'Helper function for generating plusGroup information
pGroup <- function(x,plusGroup,sOrM="survey",ages,years){
        x <- x[x$Age %in% ages,]
        x <- x[x$Year %in% years,]
        if(!is.null(plusGroup)){
            x <- dplyr::mutate(x,Age = ifelse(Age > plusGroup,plusGroup,Age))
        }       
        if(sOrM == "survey"){
            x <- dplyr::group_by(x,Year,Age,survey,fs)
            if(!is.null(plusGroup))
               x <- dplyr::summarize(x,index = sum(index,na.rm=TRUE))
            x <- dplyr::ungroup(x)
            x <- dplyr::arrange(x,survey,fs,Year,Age)
            x <- as.data.frame(x)
        }
        else if(sOrM == "catch"){
            x <- dplyr::group_by(x,Year,Age)
            if(!is.null(plusGroup))
                x <- dplyr::summarize(x,index = sum(index,na.rm=TRUE))
            x <- dplyr::ungroup(x)
            x <- dplyr::arrange(x,Year,Age)
            x <- tidyr::spread(x,key=Age,value=index)
            x <- dplyr::select(x,-Year)
            x <- as.matrix(x)
        }
        else{
            x <- dplyr::group_by(x,Year,Age)
            if(!is.null(plusGroup))
               x <- dplyr::summarize(x,index = mean(index,na.rm=TRUE))
            x <- dplyr::ungroup(x)
            x <- dplyr::arrange(x,Year,Age)
            x <- tidyr::spread(x,key=Age,value=index)
            x <- dplyr::select(x,-Year)
            x <- as.matrix(x)
        }
        
        x
}


#'Setup data for SPAM model
#'
#'Constructs data and parameters to be usable by SPAM.
#' This accepts the data setup by mdata3Ps.R in the SPAM data-raw folder.
#' 
#'
#' @param surveys named list of surveys/catch to include, catch must be named catch in list
#' @param landings data.frame of landings to include with column for upper and lower censored bound multipliers
#' @param stock_wt data.frame of stock weights
#' @param midy_wt data.frame of midyear weights
#' @param mat data.frame of maturity
#' @param recType recruitment option
#' @param ages optional vector of ages to include, if using plus group, keep ages above plusGroup in vector
#' @param years optional vector of years to include
#' @param plusGroup optional age to specify as a plus group, indices are summed by that age and weights taken to be the mean of ages in plus group
#' @param idetect The survey detection limit
#'
#' @importFrom Rcpp sourceCpp
#'
#' @useDynLib SSM3PS
#' 
#' @export
#'
datSetup <- function(surveys,landings,stock_wt,midy_wt,mat,M=0.2,ages=NULL,years=NULL,plusGroup=NULL,idetect=0.0005,aveFages=c(1,3)){
    if(is.null(ages)){
        ageL <- grep("Age",names(mat))
        ages <- as.numeric(gsub("Age","",names(mat)[ageL]))
    }

    if(is.null(years)){
        years = min(mat$Year):max(mat$Year)
    }

    if(length(M) == 1){
        if(!is.null(plusGroup)){
            maxAge = plusGroup
        }else{
            maxAge = max(ages)
        }
        M = matrix(M,nrow=length(years),ncol=length(min(ages):maxAge))
    }else{
        stopifnot(nrow(M) != length(years) & ncol(M) != length(ages),"M is not right size matrix")
    }
    
    
    vMake <- function(x){
        ageL <- grep("Age",names(x))
        v <- tidyr::gather(x,key="Age",value="index",ageL)
        v$Age <- as.numeric(gsub("Age","",v$Age))
        v
    }
                
    ##Handle the surveys into the indices
    surVec <- lapply(seq_along(surveys),function(x,nam,i){
        ageL <- grep("Age",names(x[[i]]))
        surV <- tidyr::gather(x[[i]],key="Age",value="index",ageL)
        surV$Age <- as.numeric(gsub("Age","",surV$Age))
        surV$survey = nam[[i]]
        surV
    },x=surveys,nam=names(surveys))
    surVec <- do.call(rbind,surVec)
    surVec <- pGroup(surVec,plusGroup,sOrM="survey",ages,years)
    surVec$iage <- surVec$Age - min(ages)
    surVec$iyear <- surVec$Year -min(years)
    surVec$i_zero <- ifelse(surVec$index < idetect | is.na(surVec$index),1L,0L)
    surVec$qname <- paste0(surVec$survey,":",stringr::str_pad(surVec$Age,2,pad="0",))
    surVec$qname <- ifelse(surVec$survey == "catch" | surVec$survey == "landings",NA,surVec$qname)
    surVec$iq <- as.numeric(factor(surVec$qname,levels=sort(unique(surVec$qname))))-1
    surVec$ft <- dplyr::case_when(surVec$survey == "catch" ~ 1L,
                                  surVec$survey == "landings" ~ 2L,
                                  TRUE ~ 0L) ##Otherwise it's a survey
    surVec$isurvey <- ifelse(surVec$survey == "catch" | surVec$survey == "landings",NA,surVec$survey)
    surVec$isurvey <- as.numeric(factor(surVec$isurvey))-1

    ##Setup for the plus group...
    matVec <- vMake(mat)
    stockVec <- vMake(stock_wt)
    midyVec <- vMake(midy_wt)

    pMat <- pGroup(matVec,plusGroup,sOrM="notS",ages,years)
    pStock <- pGroup(stockVec,plusGroup,sOrM="notS",ages,years)
    pMidy <- pGroup(midyVec,plusGroup,sOrM="notS",ages,years)

    landings$landings <- landings$landings/1000
    landings <- landings[landings$Year %in% years,]
    UB = landings$UB
    LB = landings$LB
    landings <- data.frame(Year=landings$Year,Age=rep(NA,nrow(landings)),
                           survey="landings",fs=rep(NA,nrow(landings)),
                           index = landings$landings,iage=rep(NA,nrow(landings)),
                           iyear=landings$Year-min(years),i_zero=0,
                           qname=rep(NA,nrow(landings)),iq=rep(NA,nrow(landings)),
                           ft=2,isurvey=rep(NA,nrow(landings)))
    surVec <- rbind(surVec,landings)
    
    surVec$index[is.na(surVec$index)] = 0
    surVec$index[surVec$index < idetect] = idetect
    surVec$index[surVec$survey == "catch"] = surVec$index[surVec$survey == "catch"]/1000

    cVec = surVec[surVec$survey == "catch",]
    surVec = surVec[surVec$survey != "catch",]
    catt = cVec[,c("index","Age","Year")]
    catt = tidyr::spread(catt,"Age","index")
    cprop = t(apply(catt[,-1],1,catchProp))
    crls = makeCRLs(cprop)
    crls = as.data.frame(cbind(years,crls))
    a2 = sort(unique(cVec$Age))
    a2 = a2[1:(length(a2)-1)]
    names(crls) = c("Year",paste0("Age",a2))
    crls = tidyr::gather(crls,"Age","index",-Year)
    crls$Age = as.numeric(gsub("Age","",crls$Age))
    cVec$index = NULL
    crls = dplyr::left_join(crls,cVec)
    crls = dplyr::arrange(crls,iyear,iage) ##Sorted for my shitty hack
    ##surVec = rbind(crls,surVec)
    
    tmb.data <- list(
        M = M,
        weight = pStock,
        mat = pMat,
        midy_weight = pMidy,
        log_index = log(surVec$index),
        i_zero = surVec$i_zero,
        iyear = surVec$iyear,
        iage = surVec$iage,
        isurvey = surVec$isurvey,
        iq = surVec$iq,
        fs = surVec$fs,
        ft = surVec$ft,
        lowerMult = LB,
        upperMult = UB,
        keyF = rep(0,length(unique(na.omit(surVec$iage)))),
        recflag = 0,
        corflag = 0,
        corflagCRL = 0,
        gammaSq = 0.1,
        crlVec = log(crls$index),
        aveFages = aveFages
        )
    

    dat <- list()
    dat$data <- tmb.data
    dat$param <- paramSetup(tmb.data)
    dat$indices <- surVec
    dat$crls <- crls
    dat
}

#'Basic Setup of parameters
#'
#' @export
paramSetup <- function(dat){
    A = ncol(dat$mat)
    Y = nrow(dat$mat)

    parm <- list(
        log_std_log_R = 0,
        log_qparm=rep(0,length(unique(na.omit(dat$iq)))),
        log_std_index=rep(log(0.3),length(unique(na.omit(dat$isurvey)))),
        log_std_logF = rep(log(0.1),length(unique(dat$keyF))),
        log_std_CRL = rep(0,A-1),
        log_std_landings = log(0.02),
        log_sdS = 0,
        rec_parm = numeric(0),
        log_F=matrix(log(0.3),nrow=Y,ncol=length(unique(dat$keyF))),
        log_N = matrix(0,nrow=Y,ncol=A)
    )
    

    parm
}

    

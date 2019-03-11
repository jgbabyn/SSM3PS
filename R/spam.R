
#'create blank key structure
#'
#' Generate the blank key matrixces for filling
#'
#' @param numAges number of ages
#' @param numFleets number of fleets
#'
#' @export
keyGen <- function(numAges,numFleets){
    ret <- list()
    ret$keyF <- rep(NA,numAges)
    ret$keyQ <- matrix(NA,nrow=numFleets,ncol=numAges)
    ret$keySD <- matrix(NA,nrow=numFleets,ncol=numAges)
    ret
}


    
#'Fit the SPAM model
#'
#'
#' @param dat list of data and parameters from spamSetup
#'
#'
#' @useDynLib 'SSM3PS'
#' 
#' @export
#'
spamFit <- function(data,parameters,map=list(),control=list(),lower=-Inf,upper=Inf,...){

    obj <- TMB::MakeADFun(data,parameters,map=map,random=dat$random,...)

    opt <- nlminb(obj$par,obj,obj$gr,control=control,lower=lower,upper=upper)

    fit <- list()

    fit$dat <- data
    fit$obj <- obj
    fit$opt <- opt
    fit$rep <- obj$report()
    fit$sdr <- TMB::sdreport()

    class(fit) <- "spam"
    fit
}


#'Returns the loglikelihood of a spam object
#'
#' @param object the spam object to get the loglikelihood of
#' 
logLik.spam <- function(object,...)
{
    val = -object$opt$objective
    attr(val,"df") = length(object$opt$par)
    attr(val,"nobs") = length(object$dat$data$obs)
    class(val) = "logLik"
    val
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
#' @param surveys named list of surveys to include
#' @param catch data.frame of catch to include
#' @param landings data.frame of landings to include
#' @param stock_wt data.frame of stock weights
#' @param midy_wt data.frame of midyear weights
#' @param mat data.frame of maturity
#' @param recType recruitment option
#' @param ages optional vector of ages to include, if using plus group, keep ages above plusGroup in vector
#' @param years optional vector of years to include
#' @param plusGroup optional age to specify as a plus group, indices are summed by that age and weights taken to be the mean of ages in plus group
#' @param idetect The survey detection limit
#' @param cdetect The catch dectection limit
#' @param match3NOdims remove final row from C,C_zero,landings
#' @param naz.rm replaces NAs and zeros with dectection limits
#' 
#'
#' @export
#'
datSetup <- function(surveys,catch,landings,stock_wt,midy_wt,mat,M=0.2,ages=NULL,years=NULL,plusGroup=NULL,idetect=0.0005,
                     cdetect=0.0005,match3NOdims=TRUE,naz.rm=TRUE){
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
    surVec$iq <- as.numeric(factor(surVec$qname,levels=sort(unique(surVec$qname))))-1

    ##Setup for the plus group...
    matVec <- vMake(mat)
    catchVec <- vMake(catch)
    stockVec <- vMake(stock_wt)
    midyVec <- vMake(midy_wt)

    pMat <- pGroup(matVec,plusGroup,sOrM="notS",ages,years)
    pCatch <- pGroup(catchVec,plusGroup,sOrM="catch",ages,years)
    pStock <- pGroup(stockVec,plusGroup,sOrM="notS",ages,years)
    pMidy <- pGroup(midyVec,plusGroup,sOrM="notS",ages,years)

    landings <- landings[landings$Year %in% years,]/1000

    if(match3NOdims){
        pMidy = pStock[-nrow(pMidy),]
        pCatch = pCatch[-nrow(pCatch),]
        landings = landings[-nrow(landings),]
        surVec = surVec[surVec$iyear < max(surVec$iyear),]
    }
    C_zero =  matrix(ifelse(pCatch < cdetect | is.na(pCatch),1,0),nrow=nrow(pCatch),ncol=ncol(pCatch))

    if(naz.rm==TRUE){
        pCatch[is.na(pCatch)] = 0
        surVec$index[is.na(surVec$index)] = 0
        pCatch[pCatch < cdetect] =cdetect
        surVec$index[surVec$index < idetect] = idetect
        }
        
    

    tmb.data <- list(
        M = M,
        weight = pStock,
        mat = pMat,
        midy_weight = pMidy,
        C = pCatch,
        landings = landings$landings,
        index = surVec$index,
        i_zero = surVec$i_zero,
        C_zero = C_zero,
        iyear = surVec$iyear,
        iage = surVec$iage,
        isurvey = as.numeric(factor(surVec$survey))-1,
        iq = surVec$iq,
        fs = surVec$fs,
        index_censor = 0,
        catch_censor = 0,
        use_pe = 0,
        use_cye = 0
    )

    dat <- list()
    dat$data <- tmb.data
    dat$param <- paramSetup(tmb.data)
    dat$indices <- surVec
    dat
}

#'Setup parameters
#'
#' @export
paramSetup <- function(dat){
    A = ncol(dat$mat)
    Y = nrow(dat$mat)

    parm <- list(
        log_No=rep(log(1),A-1),
        log_Rec_mean = 5,
        log_std_log_R = 0,
        log_qparm=rep(0,length(unique(dat$iq))),
        log_std_index=rep(log(0.3),length(unique(dat$isurvey))),
        log_std_logF = rep(log(0.1),A),
        log_std_pe = log(0.1),
        log_std_cye = log(0.5),
        logit_ar_logF_year = log(0.99/0.01),
        logit_ar_logF_age = -10,
        logit_ar_pe_year = -10,
        logit_ar_pe_age = -10,
        logit_ar_cye_year = 0,
        log_std_log_C = rep(log(0.3),A),
        log_Rec_dev=rep(5,Y),
        log_F=matrix(log(0.3),nrow=Y,ncol=A),
        pe=matrix(0,nrow=Y,ncol=A),
        cye=matrix(0,nrow=Y,ncol=A)
    )

    parameters.L <- list( 
        log_No=rep(-Inf,A),
        log_Rec_mean=0,    
        log_std_log_R=-5,   
        log_qparm=rep(-Inf,length(unique(dat$iq))),
        log_std_index=rep(log(0.01),length(unique(dat$isurvey))), 
        log_std_logF = rep(-Inf,A),  
        log_std_pe = -Inf,      
        log_std_cye = -Inf,                 
        logit_ar_logF_year = -Inf, 
        logit_ar_logF_age = -Inf,            
        logit_ar_pe_year = -Inf, 
        logit_ar_pe_age = -Inf, 
        logit_ar_cye_year = -Inf, 
        log_std_log_C = rep(-Inf,A)          
    )            

    parameters.U <- list(
        log_No=rep(log(4000),A-1),
        log_Rec_mean=50,    
        log_std_log_R=5,    
        log_qparm=rep(Inf,length(unique(dat$iq))),
        log_std_index=rep(log(2),length(unique(dat$isurvey))), 
        log_std_logF = rep(Inf,A),  
        log_std_pe = log(0.5),      
        log_std_cye = Inf,                 
        logit_ar_logF_year = log(0.99/0.01), 
        logit_ar_logF_age = log(0.99/0.01),            
        logit_ar_pe_year = log(0.99/0.01), 
        logit_ar_pe_age = log(0.99/0.01), 
        logit_ar_cye_year = log(0.99/0.01), 
        log_std_log_C = rep(Inf,A)      
    )

    ret <- list(
        param = parm,
        Lower = parameters.L,
        Upper = parameters.U
    )
}

    

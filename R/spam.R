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

    obj <- TMB::MakeADFun(data,parameters,map=map,random=data$random,...)

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
#' @export
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
#' @param cdetect The catch dectection limit
#' @param naz.rm replaces NAs and zeros with dectection limits
#' @param fit_landings fit landings and catch proportions instead of catch at age
#' @param use_pe use process error?
#' @param use_cye use catch year effects?
#' @param use_Fcorr use correlation in Fs?
#' 
#' @export
#'
datSetup <- function(surveys,landings,stock_wt,midy_wt,mat,M=0.2,ages=NULL,years=NULL,plusGroup=NULL,idetect=0.0005,cdetect=0.0005,naz.rm=TRUE,fit_landings=FALSE,use_pe=FALSE,use_cye=FALSE,use_Fcorr=FALSE){
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


    if(naz.rm==TRUE){
        surVec$index[is.na(surVec$index)] = 0
        surVec$index[surVec$index < idetect] = idetect
        }
        
    

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
        index_censor = 2,
        use_pe = ifelse(use_pe,1,0),
        use_cye = ifelse(use_cye,1,0),
        fit_land = ifelse(fit_landings,1,0),
        lowerMult = LB,
        upperMult = UB,
        use_cb = 3
    )

    dat <- list()
    dat$data <- tmb.data
    dat$param <- paramSetup(tmb.data,use_Fcorr)
    dat$indices <- surVec
    dat
}

#'Setup parameters
#'
#' @export
paramSetup <- function(dat,use_Fcorr=FALSE){
    A = ncol(dat$mat)
    Y = nrow(dat$mat)

    parm <- list(
        log_No=rep(log(1),A-1),
        log_Rec_mean = 5,
        log_std_log_R = 0,
        log_qparm=rep(0,length(unique(na.omit(dat$iq)))),
        log_std_index=rep(log(0.3),length(unique(na.omit(dat$isurvey)))),
        log_std_logF = rep(log(0.1),A),
        log_std_pe = log(0.1),
        log_std_cye = log(0.5),
        logit_ar_logF_year = log(0.99/0.01),
        logit_ar_logF_age = -10,
        logit_ar_pe_year = -10,
        logit_ar_pe_age = -10,
        logit_ar_cye_year = 0,
        log_std_log_C = rep(log(0.3),A),
        log_std_CRL = numeric(0),
        log_std_landings = numeric(0),
        log_Rec_dev=rep(5,Y),
        log_F=matrix(log(0.3),nrow=Y,ncol=A),
        pe=matrix(0,nrow=Y,ncol=A),
        cye=matrix(0,nrow=Y,ncol=A)
    )
    

    qparm_name = levels(as.factor(na.omit(dat$iq)))

    mapq = qparm_name
    agep = str_pad(6:14 , 2, pad = "0")

    ind = mapq %in% c(paste("OFF:",agep,sep=''),paste("IO:",agep,sep=''))

    mapq[ind] = paste("OFF:",rep('6+',9),sep='')

    parm$log_qparm <- log(c(0.1,0.3,0.7,1,1,1,1,1,1,1,1,1,1))

    map = list( 
        log_qparm = factor(mapq),
        log_std_logF = factor(c(0,1,rep(2,times=A-2))),
        log_std_log_C = factor(c(0,1,rep(2,times=A-2))),          
        logit_ar_logF_age = factor(NA),              
        logit_ar_logF_year = factor(NA), 
        log_std_pe=factor(NA),  
        logit_ar_pe_year = factor(NA),     
        logit_ar_pe_age = factor(NA),
        pe=factor(matrix(NA,nrow=nrow(dat$mat),ncol=ncol(dat$mat))),
        log_std_cye=factor(NA),  
        logit_ar_cye_year = factor(NA),  
        cye=factor(matrix(NA,nrow=nrow(dat$mat),ncol=ncol(dat$mat)))
    )

    
    
    if(dat$fit_land == 1){
        parm$log_std_CRL = rep(0,A-1)
        parm$log_std_landings = log(0.02)

        map = list( 
            log_qparm = factor(mapq),
            log_std_logF = factor(c(0,1,rep(2,times=A-2))),
            log_std_log_C = factor(rep(NA,A)),
            log_std_CRL = factor(c(0,1,rep(2,times=A-3))),
            logit_ar_logF_age = factor(NA),              
            logit_ar_logF_year = factor(NA), 
            log_std_pe=factor(NA),  
            logit_ar_pe_year = factor(NA),     
            logit_ar_pe_age = factor(NA),
            pe=factor(matrix(NA,nrow=nrow(dat$mat),ncol=ncol(dat$mat))),
            log_std_cye=factor(NA),  
            logit_ar_cye_year = factor(NA),  
            cye=factor(matrix(NA,nrow=nrow(dat$mat),ncol=ncol(dat$mat))),
            log_std_landings = factor(NA)
        )

        }

    if(dat$use_cye == 1){
        map$log_std_cye = NULL
        map$logit_ar_cye_year = NULL
        map$cye = NULL
    }

    if(dat$use_pe == 1){
        map$log_std_pe = NULL
        map$logit_ar_pe_year = NULL
        map$logit_ar_pe_age = NULL
        map$pe = NULL
    }

    if(use_Fcorr == TRUE){
        map$logit_ar_logF_age = NULL
        map$logit_ar_logF_year = NULL
    }
    
            
    

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
        Upper = parameters.U,
        map = map
    )
}

    

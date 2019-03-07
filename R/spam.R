#'Setup data for SPAM model
#'
#'Constructs data and parameters to be usable by SPAM.
#' Assumes data is in the format of ../data-raw/mdata3Ps.R with years. 
#'
#' @param surveys named list of surveys to include
#' @param catch named list of catch to include
#' @param landings data.frame of landings to include
#' @param stock_wt data.frame of stock weights
#' @param midy_wt data.frame of midyear weights
#' @param mat data.frame of maturity
#' @param recType recruitment option
#' @param ages optional vector of ages to include, if using plus group, keep ages above plusGroup in vector
#' @param years optional vector of years to include
#' @param plusGroup optional age to specify as a plus group, indices are summed by that age and weights taken to be the weighted mean of ages in plus group by index
#'
#' @export
#'
spamSetup <- function(surveys,catch,landings,keys,stock_wt,midy_wt,mat,recType=0,corrFFlag=0,M,years=NULL,ages=NULL,plusGroup=NULL){
    obs <- obsSetup(surveys,catch,landings,plusGroup,years=years,ages=ages)

    if(missing(years))
        years = sort(unique(obs$aux[,1]))

    if(missing(ages))
        ages = sort(unique(obs$aux[,2]))

    if(!is.null(plusGroup)){
        ages = min(ages):plusGroup
    }
    
    

    Mmat = matrix(M,nrow=length(years),ncol=length(ages))

    cAge = paste0("Age",ages)    
    sw = stock_wt[stock_wt$Year %in% years,]
    my = midy_wt[midy_wt$Year %in% years,]
    matty = mat[mat$Year %in% years,]
    
    obs$keyF <- keys$keyF
    obs$keyQ <- keys$keyQ
    obs$keySD <- keys$keySD

    obs$corrFFlag = corrFFlag

    obs$recType = recType
    obs$M = Mmat
    obs$stock_wt = as.matrix(sw)
    obs$midy_wt = as.matrix(my)
    obs$mat = as.matrix(matty)
    obs$minYear = min(years)
    obs$minAge = min(ages)

    parm <- list()
    parm$logN = matrix(0,nrow=length(years),ncol=length(ages))
    parm$logF = matrix(0,nrow=length(years),ncol=(max(obs$keyF)+1))
    parm$logSDrec = 0
    parm$logSDsur = 0
    parm$tRhoF = ifelse(corrFFlag >= 4,rep(0.1,ages),0.1)
    parm$logsdF <- numeric(max(obs$keyF,na.rm=TRUE)+1)
    parm$logQ <- numeric(max(obs$keyQ,na.rm=TRUE)+1)
    parm$obsSD <- numeric(max(obs$keySD,na.rm=TRUE)+1)

    ret <- list()
    ret$data <- obs
    ret$parm <- parm

    ret
}

    
    
    
    


#'Setup observation vector, aux stuff and idx
#'
#' @param surveys list of survey fleets
#' @param catch list of catch fleets
#' @param landings landings statistics
#' @param plusGroup Form a plus group at some year?
#'
obsSetup <- function(surveys,catch,landings,plusGroup,years,ages){
    ##Converting data to vector form

    surVec <- lapply(seq_along(surveys),function(x,nam,i){
        ageL <- grep("Age",names(x[[i]]))
        surV <- tidyr::gather(x[[i]],key="Age",value="index",ageL)
        surV$Age <- as.numeric(gsub("Age","",surV$Age))
        surV$fleet = i
        surV$survey = nam[[i]]
        surV
    },x=surveys,nam=names(surveys))
    surVec <- do.call(rbind,surVec)

    catchVec <- lapply(seq_along(catch),function(x,nam,i,lsur){
        x[[i]]$fs <- NA ##For compatibility with surVec and rbind
        ageL <- grep("Age",names(x[[i]]))
        catV <- tidyr::gather(x[[i]],key="Age",value="index",ageL)
        catV$Age <- as.numeric(gsub("Age","",catV$Age))
        catV$fleet = i + lsur
        catV$survey = NA ##Again compat, maybe if there were multiple catch types or fleets?
        catV
    },x=catch,nam=names(catch),lsur=length(surveys))
    catchVec <- do.call(rbind,catchVec)

    ##Convert landings to be similar to surVec and catchVec for rbind
   landVec <- data.frame(Year = landings$Year,fs=rep(NA,nrow(landings)),Age=rep(NA,nrow(landings)),index=landings$landings,
                            fleet = 1 + length(surveys)+length(catch),survey = rep(NA,nrow(landings)))

    obs <- do.call(rbind,list(surVec,catchVec,landVec))

    ##Plus group summarization
    ##Is weighted mean doing what I think it is doing? Yes
    if(!is.null(plusGroup)){
        obs$pg = FALSE
        obs <- dplyr::mutate(obs,Age = ifelse(Age > plusGroup,plusGroup,Age))
        obs <- dplyr::group_by(obs,Year,Age,fs,fleet,survey)
        obs <- dplyr::summarize(obs,index_n=sum(index,na.rm=TRUE))
        obs <- dplyr::mutate(obs,pg = Age == plusGroup,
                             pg = ifelse(is.na(pg),FALSE,pg))
        obs <- dplyr::arrange(obs,fleet,survey,fs,Year,Age) ##Get back in order!
        obs <- as.data.frame(obs) ##It's OK to not be a tibble anymore
    }else{
        obs$pg = FALSE
    }

    if(!is.null(years))
        obs <- obs[obs$Year %in% years,]

    if(!is.null(ages))
        obs <- obs[obs$Age %in% ages,]
    
    ret <- list()
    ret$obs <- obs$index
    ret$aux <- as.matrix(data.frame(Year=obs$Year,Age = obs$Age,fleet=obs$fleet))
    ret$fs <- obs$fs
    ret$plusGroup <- ifelse(!is.null(plusGroup),1,0)
    ret$fleetTypes <- c(rep(1,length(surveys)),rep(0,length(catch)),2)

    mmfun <- function(f,y,ff){idx<-which(ret$aux[,1]==y & ret$aux[,2] == f); ifelse(length(idx)==0,NA,
                                                                                ff(idx)-1)}
    ret$idx <- outer(1:length(ret$fleetTypes),sort(unique(ret$aux[,1])), Vectorize(mmfun,
                                                                                c("f","y")),ff=min)

    ret
    
}

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
#' @useDynLib 'SSM3Ps'
#' 
#' @export
#'
spamFit <- function(dat,map=list(),...){

    obj <- TMB::MakeADFun(dat$data,dat$parm,map=map,random=dat$random,...)

    opt <- nlminb(obj$par,obj,obj$gr,obj$he)

    fit <- list()

    fit$dat <- dat
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


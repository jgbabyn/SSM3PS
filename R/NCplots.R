#'Plot line graph with confidence interval from sd.report object
#'
#' Given a name of something in an sd.report object and the number of corresponding
#' years, plot the time series and the corresponding confidence intervals and return a ggplot object
#'
#' @param sName the name of the series to plot
#' @param sd.rep the sd.report object
#' @param years vector of corresponding years
#' @param alpha alpha value for confidence intervals
#' @param exp take the exponent of the value?
#' @param add return a list of gproto objects to add together plots?
#' @param color the color of the line
#' @param fcolor the color of the confidence interval
#'
#' @export
#'
tsCIplot <- function(sName,sd.rep,years,alpha=0.05,exp=FALSE,add=FALSE,color="black",fcolor="grey70"){
    crit = qnorm(alpha/2,lower.tail=FALSE)

    ind = names(sd.rep$value) == sName
    stdi = sd.rep$sd[ind]
    if(exp == TRUE){
        low = exp(sd.rep$value[ind] - crit*stdi)
        val = exp(sd.rep$value[ind])
        high = exp(sd.rep$value[ind] + crit*stdi)
    }else{
        low = sd.rep$value[ind] - crit*stdi
        val = sd.rep$value[ind]
        high = sd.rep$value[ind] + crit*stdi
    }
    dat = data.frame(val=val,year=years,high=high,low=low)

    pRet = ggplot2::ggplot() + ggplot2::geom_line(ggplot2::aes(x=year,y=val),dat,color=color) +
        ggplot2::geom_ribbon(ggplot2::aes(x=year,ymin=low,ymax=high),dat,fill=fcolor,alpha=0.3)

    if(add==TRUE){
        pRet = list(line = ggplot2::geom_line(ggplot2::aes(x=year,y=val),dat,color=color),
                    ribbon = ggplot2::geom_ribbon(ggplot2::aes(x=year,ymin=low,ymax=high),dat,alpha=0.3),fill=fcolor)
    }
    
    
    
    pRet
}


#'bubble plot of residuals
#'
#'Bubble plots of the residuals using either oneStepPredict residuals or not
#'
#' @param rep the report object of the fit
#' @param indices the indices data frame used in the data
#' @param oneStep optional oneStepPredict residuals
#'
#' @export
bubblePlot <- function(rep,indices,oneStep=NULL,pch=21){
    dat = data.frame(Age = indices$Age, Year = indices$Year,
                     Cohort = indices$Year-indices$Age,
                     survey = indices$survey,
                     pred = rep$Elog_index,
                     index = indices$index,
                     log_index = log(indices$index),
                     resid = rep$resid_index,
                     resid_std = rep$std_resid_index)

    if(!is.null(oneStep)){
        dat$resid_std = oneStep$residual
    }
    
    dat$colr = "deepskyblue"
    dat$colr[dat$resid >0]='firebrick2'
    dat$pch = pch
    dat$colr[(indices$i_zero==1)&(dat$resid>0)]='azure2'
    dat$colr[(indices$i_zero==1)&(dat$resid<=0)]='dimgray'

    jDarkGray <- 'grey20'
    pRet = xyplot(factor(Age) ~ Year | survey,data=dat,
                  ylab="Age",main="Residuals",
                  cex = 1.5*sqrt(abs(dat$resid_std)/pi),fill.color=dat$colr,
                  col=jDarkGray,
                  par.strip.text=list(cex=0.5),
                  panel=function(x,y,...,cex,fill.color,pch,subscripts){
                      panel.abline(v=seq(min(dat$Year),max(dat$Year),by=5),
                                   col.line=grey(0.9))
                      panel.xyplot(x,y,cex=cex[subscripts],
                                   pch=dat$pch,fill=fill.color[subscripts],...)
                  })
    pRet
}

#'Plot estimated q patterns for surveys
#'
#' Function to plot the estimated q catchabilty patterns for surveys.
#' Returns a ggplot2 plot. Values are standardized to the max of the q.
#'
#' @param sd.rep the sd.rep object
#' @param indices the indices object
#' @param alpha alpha value for the confidence intervals
#' @param qpname name of q-parameters in report if different
#'
#' @export
qPattern <- function(sd.rep,indices,alpha=0.05,qpname="log_qparm"){
    crit = qnorm(alpha/2,lower.tail=FALSE)

    ind = names(sd.rep$value) == qpname

    indy = indices[!is.na(indices$iq),] ##Needed because catch is now together with surveys 
    qest = aggregate(indy$iq,list(age=indy$Age,survey=indy$survey),unique)
    qest$iq = qest$x+1
    temp = as.numeric(table(qest$iq))
    qest$nx = temp[qest$iq]

    qest$Eest = exp(sd.rep$value[ind]);
    qest$sd =  sd.rep$sd[ind]
    qest$E.L = exp(sd.rep$value[ind] - crit*qest$sd); 
    qest$E.U = exp(sd.rep$value[ind] + crit*qest$sd);

    qest$colr='grey'
    qest$colr[qest$nx>1]='black'


    temp = tapply(qest$Eest,qest$survey,max)
    qest$max_Eest = temp[match(qest$survey,names(temp))]

    qest$sEest = qest$Eest/qest$max_Eest
    qest$sE.L = qest$E.L/qest$max_Eest 
    qest$sE.U = qest$E.U/qest$max_Eest

    qest$survey1 = paste(qest$survey,' max Q = ',round(qest$max_Eest,digits=3),' (x10000)',sep='')

    qest  = qest[order(qest$survey1,qest$age),]

    ind1 = qest$survey %in% names(temp)

    pRet = ggplot2::ggplot(qest,ggplot2::aes(x=age,y=sEest)) + ggplot2::geom_line() + ggplot2::geom_point() +
        ggplot2::geom_ribbon(ggplot2::aes(ymin=sE.L,ymax=sE.U),alpha=0.3) + facet_wrap(~survey1) + ggplot2::xlab("Age Class") +
        ggplot2::ylab("Survey Q Pattern")

    pRet
}

#'Corrplot of elements
#'
#' Corrplot of elements, uses ggcorrplot because I prefer storable plots...
#'
#' @param sd.rep
#'
#' @export
#'
corrParm <- function(sd.rep){
    dsd = sqrt(diag(sd.rep$cov.fixed))
    corr.matrix = diag(1/dsd)%*%sd.rep$cov.fixed%*%diag(1/dsd)
    ##Patch for all the indentical names
    dsdU = make.unique(names(dsd)) 
    rownames(corr.matrix) = dsdU
    colnames(corr.matrix) = dsdU

    p.mat = ggcorrplot::cor_pmat(corr.matrix)

    pRet = ggcorrplot::ggcorrplot(corr.matrix,method="circle",p.mat=p.mat,insig="blank",colors=c("red","white","blue"))
    pRet
}

#'Survey/Catch composition plot
#'
#' Compare the predicted to actual value of the survey or catch index by age or not.
#' Returns a ggplot plot. This works for both catch and surveys because everything is done in one
#' vector now.
#'
#' @param rep the report object containing the predicted indices
#' @param indices the indices list of info
#' @param indName optional the name of the catch,survey,landing etc. to look at
#' @param facetAges facet wrap the plots by age? Takes sum of index and mean of residuals over ages if false 
#' @param residuals plot the residuals instead of the line plots
#' @param oneStep optional oneStepPredict residuals to supply
#' 
#' @export
surCompPlot <- function(rep,indices,indName=NULL,facetAges=TRUE,residuals=FALSE,oneStep=NULL,loess=FALSE){
    indy = cbind(ei=rep$Elog_index,indices)
    indy$logInd = log(indy$index)
    indy$Age = as.factor(indy$Age)
    indy$zero = factor(indy$i_zero,labels=c("nonzero","zero"),levels=c(0,1))

    if(!is.null(oneStep)){
        indy$std_resid_index = oneStep$residual
    }else{
        indy$std_resid_index = rep$std_resid_index
    }
    
    if(!is.null(indName)){
        indy = indy[indy$survey %in% indName,]
    }

    

    if(facetAges == FALSE){
        indy = dplyr::group_by(indy,Year,survey)
        indy = summarize(indy,logInd=sum(logInd),ei=sum(ei),std_resid_index=mean(std_resid_index))
        facForm = as.formula(~ survey)
    }else{
        facForm = as.formula(~ survey + Age)
    }

    if(!("zero" %in% colnames(indy))){
        indy$zero = 0
        indy$zero = factor(indy$zero,labels=c("nonzero","zero"),levels=c(0,1))
    }
    
    
    if(residuals == FALSE){        
        pRet = ggplot2::ggplot(indy,ggplot2::aes(x=Year,y=logInd)) + ggplot2::geom_line(ggplot2::aes(x=Year,y=logInd))+ ggplot2::geom_line(ggplot2::aes(x=Year,y=ei),color="red") + ggplot2::xlab("Year") + ggplot2::ylab("Log Index") + ggplot2::facet_wrap(facForm) + ggplot2::geom_point(ggplot2::aes(x=Year,y=logInd,col=zero)) + scale_color_manual(values=c(zero="orange",nonzero="black")) 
    }else{
        pRet = ggplot2::ggplot(indy,ggplot2::aes(x=Year,y=std_resid_index)) + ggplot2::geom_point(ggplot2::aes(x=Year,y=std_resid_index,col=zero)) + ggplot2::xlab("Year") + ggplot2::ylab("Standardized Residuals") + ggplot2::facet_wrap(facForm) + scale_color_manual(values=c(zero="orange",nonzero="black"))
        if(loess==TRUE){
            pRet = pRet + ggplot2::geom_smooth()
        }
        
    }
    
    pRet
}

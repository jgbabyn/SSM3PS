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

    



noel.padding <- function(){
    my.padding <- list(layout.heights = list( 
                        top.padding = 0, 
                        main.key.padding = 1, 
                        key.axis.padding = 0, 
                        axis.xlab.padding = 0, 
                        xlab.key.padding = 0, 
                        key.sub.padding = 0), 
                layout.widths = list( 
                        left.padding = 1,
                        key.ylab.padding = 0, 
                        ylab.axis.padding = 1, 
                        axis.key.padding = 0, 
                        right.padding = 1) 
                )
    my.padding
}


#'Spay plot of model fit
#'
#' @param rep the report object containing the residuals
#' @param indices the indices data object containing the indices setup data
#' @param type residuals for survey or catch at age
#' 
#' @export
spayPlot <- function(rep,indices,type="survey"){
    resdat = data.frame(resid=rep$resid_index,resid_std=rep$std_resid_index,
                        Age = indices$Age,Year=indices$Year,Cohort=indices$Year-indices$Age,survey=indices$survey,
                        pred=rep$Elog_index,index=indices$index,log_index=log(indices$index))

    resdat$colr = 'deepskyblue'
    resdat$colr[resdat$resid > 0] = 'firebrick2'
    resdat$pch=21
    resdat$colr[(indices$i_zero==1)&(resdat$resid>0)]='azure2'
    resdat$colr[(indices$i_zero==1)&(resdat$resid<=0)] = 'dimgray'

    my.padding <- noel.padding()
    my.padding$layout.heights$main.key.padding = 0

    jDarkGray <- 'grey20'
    jPch <- 21
    residPlot = lattice::xyplot(factor(Age)~Year|survey,data=resdat,
                                ylab='Age',main='Residuals',
                                cex= 1.5*sqrt(abs(resdat$resid_std)))
}

    
    

    

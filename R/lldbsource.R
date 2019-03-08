#'Debug TMB object code through LLDB
#'
#' This can be used to debug C++ from within R, in the same fashion as TMB's gdbsource,
#' except it replaces gdb with the lldb debugger. I wrote this as I could not get gdbsource to work on OS X 10.14.
#'
#' @param file the file to debug
#' @param interactive whether to interactively control the lldb session.
#' @param object.name optional name of the compiled object if different from file name.
#'
#' @export
lldbsource <- function(file,interactive=FALSE,object.name=NULL){
    if(!file.exists(file))
        stop("File '",file, "'not found")

    lldbscript <- tempfile()
    if(interactive){
        lldbcmd <- c(paste("process launch -- --vanilla -f",file),"thread backtrace","bt")
        lldbcmd <- paste(lldbcmd,"\n",collapse="")
        cat(lldbcmd,file=lldbscript)
        cmd <- paste("R -d lldb --debugger-args=\"-s",lldbscript,"\"")
        system(cmd,intern=FALSE,ignore.stdout=FALSE,ignore.stderr=TRUE)
        return(NULL)
    }
    else{
        lldbcmd <- c(paste("process launch -- --vanilla -f",file),"thread backtrace","bt","quit")
        lldbcmd <- paste(lldbcmd,"\n",collapse="")
        cat(lldbcmd,file=lldbscript)
        cmd <- paste("R -d lldb --debugger-args=\"-b -s",lldbscript,"\"")
        txt <- system(cmd,intern=TRUE,ignore.stdout=FALSE,
                      ignore.stderr=TRUE)
        if(is.null(object.name))
            attr(txt,"file") <- file
        else
            attr(txt,"file") <- object.name
        class(txt) <- "backtrace"
        return(txt)
    }
}

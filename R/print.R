## Print function for ssanova objects
print.ssanova <- function(obj)
{
    ## call
    cat("\nCall:\n",deparse(obj$call),"\n\n",sep="")
    ## terms
    cat("Terms:\n")
    print.default(obj$terms$labels)
    cat("\n")
    ## terms overview
    cat("Number of fixed and random effects:\n\n")
    print.default(obj$desc)
    cat("\n")
    if (obj$method=="v") Method <- "GCV.\n"
    if (obj$method=="m") Method <- "GML.\n"
    if (obj$method=="u") Method <- "Mallows CL.\n"
    cat("Smoothing parameters are selected by",Method)
    cat("\n")
    ## the rest are suppressed
    invisible()
}

## Print function for summary.ssanova objects
print.summary.ssanova <- function (obj,digits=6)
{
    ## call
    cat("\nCall:\n",deparse(obj$call),"\n",sep="")
    cat("\nEstimate of error standard deviation:",obj$sigma,"\n")
    ## residuals
    res <- obj$res
    cat("\nResiduals:\n")
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    rq <- structure(quantile(res), names = nam)
    print(rq,digits=digits)
    cat("Residual sum of squares:",obj$rss)
    cat("\nR square:",obj$r.squared)
    ## selected summaries
    cat("\n\nPenalty associated with the fit:",obj$pen)
    cat("\n\n")
    invisible()
}

## Print function for summary.gssanova objects
print.summary.gssanova <- function (obj,digits=6)
{
    ## call
    cat("\nCall:\n",deparse(obj$call),"\n",sep="")
    if (obj$method=="u")
        cat("\n(Dispersion parameter for ",obj$family,
            " family taken to be ",format(obj$dispersion),")\n\n",sep="")
    if (obj$method=="v")
        cat("\n(Dispersion parameter for ",obj$family,
            " family estimated to be ",format(obj$dispersion),")\n\n",sep="")
    ## residuals
    res <- obj$res
    cat("Working residuals:\n")
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    rq <- structure(quantile(res), names = nam)
    print(rq,digits=digits)
    cat("Residual sum of squares:",obj$rss,"\n")
    ## deviance residuals
    res <- obj$dev.res
    cat("\nDeviance residuals:\n")
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    rq <- structure(quantile(res), names = nam)
    print(rq,digits=digits)
    cat("Deviance:",obj$deviance)
    cat("\nNull deviance:",obj$dev.null)
    ## selected summaries
    cat("\n\nPenalty associated with the fit:",obj$pen)
    cat("\n\nNumber of performance-oriented iterations:",obj$iter)
    cat("\n\n")
    invisible()
}

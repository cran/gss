print.ssanova <-
function(obj) {
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
  if (obj$method=="m") Method <- "Type-II ML.\n"
  if (obj$method=="u") Method <- "Mallows CL.\n"
  cat("Smoothing parameters are selected by",Method)
  cat("\n")
  ## the rest are suppressed
  invisible()
}

print.summary.ssanova <-
function (x,digits=6) {
  ## call
  cat("\nCall:\n",deparse(x$call),"\n",sep="")
  cat("\nEstimate of error standard deviation:",x$sigma,"\n")
  ## residuals
  res <- x$res
  cat("\nResiduals:\n")
  nam <- c("Min", "1Q", "Median", "3Q", "Max")
  rq <- structure(quantile(res), names = nam)
  print(rq,digits=digits)
  cat("Residual sum of squares:",x$rss)
  cat("\nR square:",x$r.squared)
  ## selected summaries
  cat("\n\nPenalty associated with the fit:",x$pen)
#  cat("\n\nNumber of Observations:",length(x$res))
  cat("\n\n")
  invisible()
}

print.summary.gssanova <-
function (x,digits=6) {
  ## call
  cat("\nCall:\n",deparse(x$call),"\n",sep="")
  if (x$method=="u")
    cat("\n(Dispersion parameter for ",x$family,
        " family taken to be ",format(x$dispersion),")\n\n",sep="")
  if (x$method=="v")
    cat("\n(Dispersion parameter for ",x$family,
        " family estimated to be ",format(x$dispersion),")\n\n",sep="")
  ## residuals
  res <- x$res
  cat("Working residuals:\n")
  nam <- c("Min", "1Q", "Median", "3Q", "Max")
  rq <- structure(quantile(res), names = nam)
  print(rq,digits=digits)
  cat("Residual sum of squares:",x$rss,"\n")
  ## deviance residuals
  res <- x$dev.res
  cat("\nDeviance residuals:\n")
  nam <- c("Min", "1Q", "Median", "3Q", "Max")
  rq <- structure(quantile(res), names = nam)
  print(rq,digits=digits)
  cat("Deviance:",x$deviance)
  cat("\nNull deviance:",x$dev.null)
  ## selected summaries
  cat("\n\nPenalty associated with the fit:",x$pen)
  cat("\n\nNumber of performance-oriented iterations:",x$iter)
  cat("\n\n")
  invisible()
}

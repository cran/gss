getsms <- ## Obtain var-cov matrix for fixed effects
function(obj) {
  ## Check input
  if (!any(class(obj)=="ssanova")) {
    stop("gss error in getsms: inputs are of wrong types")
  }
  nobs <- length(obj$c)
  nnull <- length(obj$d)
  ## Call RKPACK ulitity DSMS
  z <- .Fortran("dsms",
                as.double(obj$swk), as.integer(nobs),
                as.integer(nobs), as.integer(nnull),
                as.integer(obj$jpvt),
                as.double(obj$qwk), as.integer(nobs),
                as.double(obj$nlambda),
                sms=double(nnull*nnull), as.integer(nnull),
                double(2*nobs), integer(1))["sms"]
  ## Return the nnull-by-nnull matrix
  matrix(z$sms,nnull,nnull)
}

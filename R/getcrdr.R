getcrdr <- ## Obtain c & d for new y's
function(obj,r) {
  ## Check inputs
  if (is.vector(r)) r <- as.matrix(r)
  if (!(any(class(obj)=="ssanova")&is.matrix(r))) {
    stop("gss error in getcrdr: inputs are of wrong types")
  }
  nobs <- length(obj$c)
  nnull <- length(obj$d)
  nr <- dim(r)[2]
  if (!((dim(r)[1]==nobs)&(nr>0))) {
    stop("gss error in getcrdr: inputs have wrong dimensions")
  }
  ## Call RKPACK ulitity DCRDR
  z <- .Fortran("dcrdr",
                as.double(obj$swk), as.integer(nobs),
                as.integer(nobs), as.integer(nnull),
                as.double(obj$qraux), as.integer(obj$jpvt),
                as.double(obj$qwk), as.integer(nobs),
                as.double(obj$nlambda),
                as.double(r), as.integer(nobs), as.integer(nr),
                cr=double(nobs*nr), as.integer(nobs),
                dr=double(nnull*nr), as.integer(nnull),
                double(2*nobs), integer(1))[c("cr","dr")]
  ## Return cr and dr
  z$cr <- matrix(z$cr,nobs,nr)
  z$dr <- matrix(z$dr,nnull,nr)
  z
}

mspreg <- ## Fit Multiple Smoothing Parameter REGression
function(s,q,y,method="v",varht=1,prec=1e-7,maxiter=30) {
  ## Check inputs
  if (is.vector(s)) s <- as.matrix(s)
  if (!(is.matrix(s)&is.array(q)&(length(dim(q))==3)
        &is.vector(y)&is.character(method))) {
    stop("gss error in mspreg: inputs are of wrong types")
  }
  nobs <- length(y)
  nnull <- dim(s)[2]
  nq <- dim(q)[3]
  if (!((dim(s)[1]==nobs)&(dim(q)[1]==nobs)&(dim(q)[2]==nobs)
        &(nobs>=nnull)&(nnull>0)&(nq>1))) {
    stop("gss error in mspreg: inputs have wrong dimensions")
  }
  ## Set method for smoothing parameter selection
  code <- (1:3)[c("v","m","u")==method]
  if (!length(code)) {
    stop("gss error: unsupported method for smoothing parameter selection")
  }
  ## Call RKPACK driver DMUDR
  z <- .Fortran("dmudr0",
                as.integer(code),
                as.double(s),           # s
                as.integer(nobs), as.integer(nobs), as.integer(nnull),
                as.double(q),           # q
                as.integer(nobs), as.integer(nobs), as.integer(nq),
                as.double(y),           # y
                as.double(0), as.integer(0),
                as.double(prec), as.integer(maxiter),
                theta=double(nq), nlambda=double(1),
                score=double(1), varht=as.double(varht),
                c=double(nobs), d=double(nnull),
                double(nobs*nobs*(nq+2)),
                info=integer(1))[c("theta","info")]
  ## Check info for error
  if (info<-z$info) {               
    if (info>0)
      stop("gss error in mspreg: matrix s is rank deficient")
    if (info==-2)
      stop("gss error in mspreg: matrix q is indefinite")
    if (info==-1)
      stop("gss error in mspreg: input data have wrong dimensions")
    if (info==-3)
      stop("gss error in mspreg: unknown method for smoothing parameter selection.")
    if (info==-4)
      stop("gss error in mspreg: iteration fails to converge, try to increase maxiter")
    if (info==-5)
      stop("gss error in mspreg: iteration fails to find a reasonable descent direction")
  }
  qwk <- 10^z$theta[1]*q[,,1]
  for (i in 2:nq) qwk <- qwk + 10^z$theta[i]*q[,,i]
  ## Call RKPACK driver DSIDR
  zz <- .Fortran("dsidr0",
                 as.integer(code),
                 swk=as.double(s), as.integer(nobs),
                 as.integer(nobs), as.integer(nnull),
                 as.double(y),
                 qwk=as.double(qwk), as.integer(nobs),
                 as.double(0), as.integer(0), double(2),
                 nlambda=double(1), score=double(1), varht=as.double(varht),
                 c=double(nobs), d=double(nnull),
                 qraux=double(nnull), jpvt=integer(nnull),
                 double(3*nobs),
                 info=integer(1))
  ## Return the fit
  c(list(method=method,theta=z$theta),
    zz[c("c","d","nlambda","score","varht","swk","qraux","jpvt","qwk")])
}

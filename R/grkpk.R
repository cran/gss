## Fit Single Smoothing Parameter REGression by Performance-Oriented Iteration
sspregpoi <- function(family,s,q,y,wt,offset,method="u",
                      varht=1,prec=1e-7,maxiter=30)
{
    ## Check inputs
    if (is.vector(s)) s <- as.matrix(s)
    if (!(is.matrix(s)&is.matrix(q)&is.character(method))) {
        stop("gss error in sspregpoi: inputs are of wrong types")
    }
    nobs <- dim(s)[1]
    nnull <- dim(s)[2]
    if (!((dim(s)[1]==nobs)&(dim(q)[1]==nobs)&(dim(q)[2]==nobs)
          &(nobs>=nnull)&(nnull>0))) {
        stop("gss error in sspregpoi: inputs have wrong dimensions")
    }
    ## Set method for smoothing parameter selection
    code <- (1:3)[c("v","m","u")==method]
    if (!length(code)) {
        stop("gss error: unsupported method for smoothing parameter selection")
    }
    eta <- rep(0,nobs)
    nla0 <- log10(mean(abs(diag(q))))
    limnla <- nla0+c(-.5,.5)
    iter <- 0
    alpha <- NULL
    repeat {
        iter <- iter+1
        dat <- switch(family,
                      binomial=mkdata.binomial(y,eta,wt,offset),
                      nbinomial=mkdata.nbinomial(y,eta,wt,offset,alpha),
                      poisson=mkdata.poisson(y,eta,wt,offset),
                      inverse.gaussian=mkdata.inverse.gaussian(y,eta,wt,offset),
                      Gamma=mkdata.Gamma(y,eta,wt,offset))
        alpha <- dat$alpha
        w <- as.vector(sqrt(dat$wt))
        ywk <- w*dat$ywk
        swk <- w*s
        qwk <- w*t(w*q)
        ## Call RKPACK driver DSIDR
        z <- .Fortran("dsidr0",
                      as.integer(code),
                      swk=as.double(swk), as.integer(nobs),
                      as.integer(nobs), as.integer(nnull),
                      as.double(ywk),
                      qwk=as.double(qwk), as.integer(nobs),
                      as.double(0), as.integer(-1), as.double(limnla),
                      nlambda=double(1), score=double(1), varht=as.double(varht),
                      c=double(nobs), d=double(nnull),
                      qraux=double(nnull), jpvt=integer(nnull),
                      double(3*nobs),
                      info=integer(1))
        ## Check info for error
        if (info<-z$info) {               
            if (info>0)
                stop("gss error in sspregpoi: matrix s is rank deficient")
            if (info==-2)
                stop("gss error in sspregpoi: matrix q is indefinite")
            if (info==-1)
                stop("gss error in sspregpoi: input data have wrong dimensions")
            if (info==-3)
                stop("gss error in sspregpoi: unknown method for smoothing parameter selection.")
        }
        eta.new <- (ywk-10^z$nlambda*z$c)/w
        if (!is.null(offset)) eta.new <- eta.new + offset
        disc <- sum(dat$wt*((eta-eta.new)/(1+abs(eta)))^2)/sum(dat$wt)
        limnla <- pmax(z$nlambda+c(-.5,.5),nla0-5)
        if (disc<prec) break
        if (iter>=maxiter) {
            warning("gss warning: performance-oriented iteration fails to converge")
            break
        }
        eta <- eta.new
    }
    ## Return the fit
    c(list(method=method,theta=0,w=as.vector(dat$wt),
           eta=as.vector(eta),iter=iter,alpha=alpha),
      z[c("c","d","nlambda","score","varht","swk","qraux","jpvt","qwk")])
}

## Fit Multiple Smoothing Parameter REGression by Performance-Oriented Iteration
mspregpoi <- function(family,s,q,y,wt,offset,method="u",
                      varht=1,prec=1e-7,maxiter=30)
{
    ## Check inputs
    if (is.vector(s)) s <- as.matrix(s)
    if (!(is.matrix(s)&is.array(q)&(length(dim(q))==3)
          &is.character(method))) {
        stop("gss error in mspregpoi: inputs are of wrong types")
    }
    nobs <- dim(s)[1]
    nnull <- dim(s)[2]
    nq <- dim(q)[3]
    if (!((dim(s)[1]==nobs)&(dim(q)[1]==nobs)&(dim(q)[2]==nobs)
          &(nobs>=nnull)&(nnull>0))) {
        stop("gss error in sspregpoi: inputs have wrong dimensions")
    }
    ## Set method for smoothing parameter selection
    code <- (1:3)[c("v","m","u")==method]
    if (!length(code)) {
        stop("gss error: unsupported method for smoothing parameter selection")
    }
    eta <- rep(0,nobs)
    init <- 0
    theta <- rep(0,nq)
    iter <- 0
    alpha <- NULL
    qwk <- array(0,c(nobs,nobs,nq))
    repeat {
        iter <- iter+1
        dat <- switch(family,
                      binomial=mkdata.binomial(y,eta,wt,offset),
                      nbinomial=mkdata.nbinomial(y,eta,wt,offset,alpha),
                      poisson=mkdata.poisson(y,eta,wt,offset),
                      inverse.gaussian=mkdata.inverse.gaussian(y,eta,wt,offset),
                      Gamma=mkdata.Gamma(y,eta,wt,offset))
        alpha <- dat$alpha
        w <- as.vector(sqrt(dat$wt))
        ywk <- w*dat$ywk
        swk <- w*s
        for (i in 1:nq) qwk[,,i] <- w*t(w*q[,,i])
        ## Call RKPACK driver DMUDR
        z <- .Fortran("dmudr0",
                      as.integer(code),
                      as.double(swk),   # s
                      as.integer(nobs), as.integer(nobs), as.integer(nnull),
                      as.double(qwk),   # q
                      as.integer(nobs), as.integer(nobs), as.integer(nq),
                      as.double(ywk),   # y
                      as.double(0), as.integer(init),
                      as.double(prec), as.integer(maxiter),
                      theta=as.double(theta), nlambda=double(1),
                      score=double(1), varht=as.double(varht),
                      c=double(nobs), d=double(nnull),
                      double(nobs*nobs*(nq+2)),
                      info=integer(1))[c("theta","nlambda","c","info")]
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
        eta.new <- (ywk-10^z$nlambda*z$c)/w
        if (!is.null(offset)) eta.new <- eta.new + offset
        disc <- sum(dat$wt*((eta-eta.new)/(1+abs(eta)))^2)/sum(dat$wt)
        if (disc<prec) break
        if (iter>=maxiter) {
            warning("gss warning: performance-oriented iteration fails to converge")
            break
        }
        init <- 1
        theta <- z$theta
        eta <- eta.new
    }
    qqwk <- 10^z$theta[1]*qwk[,,1]
    for (i in 2:nq) qqwk <- qqwk + 10^z$theta[i]*qwk[,,i]
    ## Call RKPACK driver DSIDR
    z <- .Fortran("dsidr0",
                  as.integer(code),
                  swk=as.double(swk), as.integer(nobs),
                  as.integer(nobs), as.integer(nnull),
                  as.double(ywk),
                  qwk=as.double(qqwk), as.integer(nobs),
                  as.double(0), as.integer(0), double(2),
                  nlambda=double(1), score=double(1), varht=as.double(varht),
                  c=double(nobs), d=double(nnull),
                  qraux=double(nnull), jpvt=integer(nnull),
                  double(3*nobs),
                  info=integer(1))
    ## Check info for error
    if (info<-z$info) {               
        if (info>0)
            stop("gss error in sspregpoi: matrix s is rank deficient")
        if (info==-2)
            stop("gss error in sspregpoi: matrix q is indefinite")
        if (info==-1)
            stop("gss error in sspregpoi: input data have wrong dimensions")
        if (info==-3)
            stop("gss error in sspregpoi: unknown method for smoothing parameter selection.")
    }
    ## Return the fit
    c(list(method=method,theta=theta,w=as.vector(dat$wt),
           eta=as.vector(eta),iter=iter,alpha=alpha),
      z[c("c","d","nlambda","score","varht","swk","qraux","jpvt","qwk")])
}

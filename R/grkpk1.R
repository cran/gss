## Fit Single Smoothing Parameter Non-Gaussian REGression
sspngreg1 <- function(family,s,r,q,y,wt,offset,alpha,nu,random)
{
    nobs <- nrow(r)
    nxi <- ncol(r)
    if (!is.null(s)) {
        if (is.vector(s)) nnull <- 1
        else nnull <- ncol(s)
    }
    else nnull <- 0
    if (!is.null(random)) nz <- ncol(as.matrix(random$z))
    else nz <- 0
    nxiz <- nxi + nz
    nn <- nxiz + nnull
    ## cv function
    cv <- function(lambda) {
        if (is.null(random)) q.wk <- 10^(lambda+theta)*q
        else {
            q.wk <- matrix(0,nxiz,nxiz)
            q.wk[1:nxi,1:nxi] <- 10^(lambda[1]+theta)*q
            q.wk[(nxi+1):nxiz,(nxi+1):nxiz] <- random$sigma(lambda[-1])
        }
        alpha.wk <- max(0,log.la0-lambda[1]-5)*(3-alpha) + alpha
        alpha.wk <- min(alpha.wk,3)
        z <- ngreg1(dc,family,cbind(s,10^theta*r),q.wk,y,wt,offset,nu,alpha.wk)
        assign("dc",z$dc,inherit=TRUE)
        assign("fit",z[c(1:3,5:10)],inherit=TRUE)
        z$score
    }
    ## initialization
    tmp <- sum(r^2)
    if (is.null(s)) theta <- 0
    else theta <- log10(sum(s^2)/nnull/tmp*nxi) / 2
    log.la0 <- log10(tmp/sum(diag(q))) + theta
    if (!is.null(random)) r <- cbind(r,10^(-theta)*random$z)
    ## lambda search
    dc <- rep(0,nn)
    fit <- NULL
    if (is.null(random)) la <- log.la0
    else la <- c(log.la0,random$init)
    if (length(la)-1) {
        counter <- 0
        repeat {
            zz <- nlm(cv,la,stepmax=1,ndigit=7)
            if (zz$code<=3) break
            la <- zz$est
            counter <- counter + 1
            if (counter>=5) {
                warning("gss warning in ssanova1: iteration for model selection fails to converge")
                break
            }
        }
    }
    else {
        repeat {
            mn <- la-1
            mx <- la+1
            zz <- nlm0(cv,c(mn,mx))
            if (min(zz$est-mn,mx-zz$est)>=1e-3) break
            else la <- zz$est
        }
    }
    ## return
    jk <- cv(zz$est)
    if (is.null(random)) q.wk <- 10^theta*q
    else {
        q.wk <- matrix(0,nxiz,nxiz)
        q.wk[1:nxi,1:nxi] <- 10^theta*q
        q.wk[(nxi+1):nxiz,(nxi+1):nxiz] <- 10^(-zz$est[1])*random$sigma(zz$est[-1])
    }
    zzz <- La.eigen(q.wk,TRUE)
    rkq <- min(fit$rkv-nnull,sum(zzz$val/zzz$val[1]>sqrt(.Machine$double.eps)))
    val <- zzz$val[1:rkq]
    vec <- zzz$vec[,1:rkq,drop=FALSE]
    qinv <- vec%*%diag(1/val,rkq)%*%t(vec)
    se.aux <- t(fit$w*cbind(s,10^theta*r))%*%(10^theta*r)%*%qinv
    c <- fit$dc[nnull+(1:nxi)]
    if (nnull) d <- fit$dc[1:nnull]
    else d <- NULL
    if (nz) b <- fit$dc[nnull+nxi+(1:nz)]
    else b <- NULL
    c(list(theta=theta,c=c,d=d,b=b,nlambda=zz$est[1],zeta=zz$est[-1]),
      fit[-1],list(qinv=qinv,se.aux=se.aux))
}

## Fit Multiple Smoothing Parameter Non-Gaussian REGression
mspngreg1 <- function(family,s,r,q,y,wt,offset,alpha,nu,random)
{
    nobs <- nrow(r)
    nxi <- ncol(r)
    if (!is.null(s)) {
        if (is.vector(s)) nnull <- 1
        else nnull <- ncol(s)
    }
    else nnull <- 0
    if (!is.null(random)) nz <-ncol(as.matrix(random$z))
    else nz <- 0
    nxiz <- nxi + nz
    nn <- nxiz + nnull
    nq <- dim(q)[3]
    ## cv function
    cv <- function(theta) {
        r.wk <- qq.wk <- 0
        for (i in 1:nq) {
            r.wk <- r.wk + 10^theta[i]*r[,,i]
            qq.wk <- qq.wk + 10^theta[i]*q[,,i]
        }
        if (is.null(random)) q.wk <- 10^nlambda*qq.wk
        else {
            r.wk <- cbind(r.wk,random$z)
            q.wk <- matrix(0,nxiz,nxiz)
            q.wk[1:nxi,1:nxi] <- 10^nlambda*qq.wk
            q.wk[(nxi+1):nxiz,(nxi+1):nxiz] <- random$sigma(theta[-(1:nq)])
        }
        alpha.wk <- max(0,theta[1:nq]-log.th0-5)*(3-alpha) + alpha
        alpha.wk <- min(alpha.wk,3)
        z <- ngreg1(dc,family,cbind(s,r.wk),q.wk,y,wt,offset,nu,alpha.wk)
        assign("dc",z$dc,inherit=TRUE)
        assign("fit",z[c(1:3,5:10)],inherit=TRUE)
        z$score
    }
    ## initialization
    theta <- -log10(apply(q,3,function(x)sum(diag(x))))
    r.wk <- q.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
        q.wk <- q.wk + 10^theta[i]*q[,,i]
    }
    ## theta adjustment
    z <- sspngreg1(family,s,r.wk,q.wk,y,wt,offset,alpha,nu,random)
    theta <- theta + z$theta
    r.wk <- q.wk <- 0
    for (i in 1:nq) {
        theta[i] <- 2*theta[i] + log10(t(z$c)%*%q[,,i]%*%z$c)
        r.wk <- r.wk + 10^theta[i]*r[,,i]
        q.wk <- q.wk + 10^theta[i]*q[,,i]
    }
    log.la0 <- log10(sum(r.wk^2)/sum(diag(q.wk)))
    log.th0 <- theta-log.la0
    ## lambda search
    z <- sspngreg1(family,s,r.wk,q.wk,y,wt,offset,alpha,nu,random)
    nlambda <- z$nlambda
    log.th0 <- log.th0 + z$lambda
    theta <- theta + z$theta
    ## theta search
    dc <- rep(0,nn)
    fit <- NULL
    if (!is.null(random)) theta <- c(theta,z$zeta)
    counter <- 0
    repeat {
        zz <- nlm(cv,theta,stepmax=1,ndigit=7)
        if (zz$code<=3)  break
        theta <- zz$est        
        counter <- counter + 1
        if (counter>=5) {
            warning("gss warning in gssanova1: iteration for model selection fails to converge")
            break
        }
    }
    ## return
    jk <- cv(zz$est)
    r.wk <- qq.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^zz$est[i]*r[,,i]
        qq.wk <- qq.wk + 10^zz$est[i]*q[,,i]
    }
    if (is.null(random)) q.wk <- qq.wk
    else {
        r.wk <- cbind(r.wk,random$z)
        q.wk <- matrix(0,nxiz,nxiz)
        q.wk[1:nxi,1:nxi] <- qq.wk
        q.wk[(nxi+1):nxiz,(nxi+1):nxiz] <- 10^(-nlambda)*random$sigma(zz$est[-(1:nq)])
    }
    zzz <- La.eigen(q.wk,TRUE)
    rkq <- min(fit$rkv-nnull,sum(zzz$val/zzz$val[1]>sqrt(.Machine$double.eps)))
    val <- zzz$val[1:rkq]
    vec <- zzz$vec[,1:rkq,drop=FALSE]
    qinv <- vec%*%diag(1/val,rkq)%*%t(vec)
    se.aux <- t(fit$w*cbind(s,r.wk))%*%r.wk%*%qinv
    c <- fit$dc[nnull+(1:nxi)]
    if (nnull) d <- fit$dc[1:nnull]
    else d <- NULL
    if (nz) b <- fit$dc[nnull+nxi+(1:nz)]
    else b <- NULL
    c(list(theta=zz$est[1:nq],c=c,d=d,b=b,nlambda=nlambda,zeta=zz$est[-(1:nq)]),
      fit[-1],list(qinv=qinv,se.aux=se.aux))
}

## Non-Gaussian regression with fixed smoothing parameters
ngreg1 <- function(dc,family,sr,q,y,wt,offset,nu,alpha)
{
    nobs <- nrow(sr)
    nn <- ncol(sr)
    nxi <- nrow(q)
    nnull <- nn - nxi
    ## initialization
    eta <- sr%*%dc
    if (!is.null(offset)) eta <- eta + offset
    if (family=="nbinomial") nu <- NULL
    else nu <- list(nu,is.null(nu))
    iter <- 0
    flag <- 0
    ## Newton iteration
    repeat {
        iter <- iter+1
        dat <- switch(family,
                      binomial=mkdata.binomial(y,eta,wt,offset),
                      nbinomial=mkdata.nbinomial(y,eta,wt,offset,nu),
                      poisson=mkdata.poisson(y,eta,wt,offset),
                      inverse.gaussian=mkdata.inverse.gaussian(y,eta,wt,offset),
                      Gamma=mkdata.Gamma(y,eta,wt,offset),
                      weibull=mkdata.weibull(y,eta,wt,offset,nu),
                      lognorm=mkdata.lognorm(y,eta,wt,offset,nu),
                      loglogis=mkdata.loglogis(y,eta,wt,offset,nu))
        nu <- dat$alpha
        ## weighted least squares fit
        w <- as.vector(sqrt(dat$wt))
        ywk <- w*dat$ywk
        srwk <- w*sr
        if (!is.finite(sum(w,ywk,srwk))) {
            if (flag) stop("gss error in gssanova1: Newton iteration diverges")
            eta <- rep(0,nobs)
            iter <- 0
            flag <- 1
            next
        }
        z <- .Fortran("reg",
                      as.double(srwk), as.integer(nobs), as.integer(nnull),
                      as.double(q), as.integer(nxi), as.double(ywk),
                      as.integer(4),
                      double(1), double(1), double(1), dc=double(nn),
                      as.double(.Machine$double.eps),
                      double(nn*nn), double(nn),
                      as.integer(c(rep(1,nnull),rep(0,nxi))),
                      double(max(nobs,nn)), integer(1), integer(1),
                      PACKAGE="gss")["dc"]
        eta.new <- sr%*%z$dc
        if (!is.null(offset)) eta.new <- eta.new + offset
        disc <- sum(dat$wt*((eta-eta.new)/(1+abs(eta)))^2)/sum(dat$wt)
        if (!is.finite(disc)) {
            if (flag) stop("gss error in gssanova1: Newton iteration diverges")
            eta <- rep(0,nobs)
            iter <- 0
            flag <- 1
            next
        }
        eta <- eta.new
        if (disc<1e-7) break
        if (iter<=30) next
        if (!flag) {
            eta <- rep(0,nobs)
            iter <- 0
            flag <- 1
        }
        else {
            warning("gss warning in gssanova1: Newton iteration fails to converge")
            break
        }
    }
    ## calculate cv
    dat <- switch(family,
                  binomial=mkdata.binomial(y,eta,wt,offset),
                  nbinomial=mkdata.nbinomial(y,eta,wt,offset,nu),
                  poisson=mkdata.poisson(y,eta,wt,offset),
                  inverse.gaussian=mkdata.inverse.gaussian(y,eta,wt,offset),
                  Gamma=mkdata.Gamma(y,eta,wt,offset),
                  weibull=mkdata.weibull(y,eta,wt,offset,nu),
                  lognorm=mkdata.lognorm(y,eta,wt,offset,nu),
                  loglogis=mkdata.loglogis(y,eta,wt,offset,nu))
    nu <- dat$alpha
    ## weighted least squares fit
    w <- as.vector(sqrt(dat$wt))
    ywk <- w*dat$ywk
    srwk <- w*sr
    z <- .Fortran("reg",
                  as.double(srwk), as.integer(nobs), as.integer(nnull),
                  as.double(q), as.integer(nxi), as.double(ywk),
                  as.integer(5),
                  double(1), double(1), double(1), dc=double(nn),
                  as.double(.Machine$double.eps),
                  chol=double(nn*nn), double(nn),
                  jpvt=as.integer(c(rep(1,nnull),rep(0,nxi))),
                  hat=double(max(nobs+1,nn)), rkv=integer(1), integer(1),
                  PACKAGE="gss")[c("dc","chol","jpvt","hat","rkv")]
    cv <- switch(family,
                 binomial=cv.binomial(y,eta,wt,z$hat[1:nobs],alpha),
                 poisson=cv.poisson(y,eta,wt,z$hat[1:nobs],alpha,sr,q),
                 Gamma=cv.Gamma(y,eta,wt,z$hat[1:nobs],z$hat[nobs+1],alpha))
    c(z,cv,list(eta=eta,nu=nu[[1]]))
}

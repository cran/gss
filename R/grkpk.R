## Fit Single Smoothing Parameter Non-Gaussian REGression
sspngreg <- function(family,s,r,q,y,wt,offset,alpha,nu,random)
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
        if (nu[[2]]) {
            la.wk <- lambda[-2]
            nu.wk <- list(exp(lambda[2]),FALSE)
        }
        else {
            la.wk <- lambda
            nu.wk <- nu
        }
        if (is.null(random)) q.wk <- 10^(la.wk+theta)*q
        else {
            q.wk <- matrix(0,nxiz,nxiz)
            q.wk[1:nxi,1:nxi] <- 10^(la.wk[1]+theta)*q
            q.wk[(nxi+1):nxiz,(nxi+1):nxiz] <-
                10^(2*ran.scal)*random$sigma$fun(la.wk[-1],random$sigma$env)
        }
        alpha.wk <- max(0,log.la0-la.wk[1]-5)*(3-alpha) + alpha
        alpha.wk <- min(alpha.wk,3)
        z <- ngreg(dc,family,cbind(s,10^theta*r),q.wk,y,wt,offset,nu.wk,alpha.wk)
        assign("dc",z$dc,inherit=TRUE)
        assign("fit",z[c(1:3,5:10)],inherit=TRUE)
        z$score
    }
    cv.wk <- function(lambda) cv.scale*cv(lambda)+cv.shift
    ## initialization
    tmp <- sum(r^2)
    if (is.null(s)) theta <- 0
    else theta <- log10(sum(s^2)/nnull/tmp*nxi) / 2
    log.la0 <- log10(tmp/sum(diag(q))) + theta
    if (!is.null(random)) {
        ran.scal <- theta - log10(sum(random$z^2)/nz/tmp*nxi) / 2
        r <- cbind(r,10^(ran.scal-theta)*random$z)
    }
    else ran.scal <- NULL
    if (nu[[2]]&is.null(nu[[1]])) {
        eta <- rep(0,nobs)
        wk <- switch(family,
                      nbinomial=mkdata.nbinomial(y,eta,wt,offset,NULL),
                      weibull=mkdata.weibull(y,eta,wt,offset,nu),
                      lognorm=mkdata.lognorm(y,eta,wt,offset,nu),
                      loglogis=mkdata.loglogis(y,eta,wt,offset,nu))
        nu[[1]] <- wk$nu[[1]]
    }
    ## lambda search
    dc <- rep(0,nn)
    fit <- NULL
    la <- log.la0
    if (nu[[2]]) la <- c(la, log(nu[[1]]))
    if (!is.null(random)) la <- c(la,random$init)
    if (length(la)-1) {
        counter <- 0
        ## scale and shift cv
        tmp <- abs(cv(la))
        cv.scale <- 1
        cv.shift <- 0
        if (tmp<1&tmp>10^(-4)) {
            cv.scale <- 10/tmp
            cv.shift <- 0
        }
        if (tmp<10^(-4)) {
            cv.scale <- 10^2
            cv.shift <- 10
        }
        repeat {
            zz <- nlm(cv.wk,la,stepmax=1,ndigit=7)
            if (zz$code<=3) break
            la <- zz$est
            counter <- counter + 1
            if (counter>=5) {
                warning("gss warning in ssanova: iteration for model selection fails to converge")
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
    if (nu[[2]]) {
        nu.wk <- exp(zz$est[2])
        zz$est <- zz$est[-2]
    }
    else nu.wk <- NULL
    if (is.null(random)) q.wk <- 10^theta*q
    else {
        q.wk <- matrix(0,nxiz,nxiz)
        q.wk[1:nxi,1:nxi] <- 10^theta*q
        q.wk[(nxi+1):nxiz,(nxi+1):nxiz] <-
            10^(2*ran.scal-zz$est[1])*random$sigma$fun(zz$est[-1],random$sigma$env)
    }
    qinv <- eigen(q.wk,TRUE)
    se.aux <- t(fit$w*cbind(s,10^theta*r))%*%(10^theta*r)%*%qinv$vec
    c <- fit$dc[nnull+(1:nxi)]
    if (nnull) d <- fit$dc[1:nnull]
    else d <- NULL
    if (nz) b <- 10^(ran.scal)*fit$dc[nnull+nxi+(1:nz)]
    else b <- NULL
    c(list(theta=theta,ran.scal=ran.scal,c=c,d=d,b=b,nlambda=zz$est[1],
           zeta=zz$est[-1],nu=nu.wk),fit[-1],list(qinv=qinv,se.aux=se.aux))
}

## Fit Multiple Smoothing Parameter Non-Gaussian REGression
mspngreg <- function(family,s,r,q,y,wt,offset,alpha,nu,random)
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
        if (nu[[2]]) {
            the.wk <- theta[-(nq+1)]
            nu.wk <- list(exp(theta[nq+1]),FALSE)
        }
        else {
            the.wk <- theta
            nu.wk <- nu
        }
        r.wk <- qq.wk <- 0
        for (i in 1:nq) {
            r.wk <- r.wk + 10^the.wk[i]*r[,,i]
            qq.wk <- qq.wk + 10^the.wk[i]*q[,,i]
        }
        if (is.null(random)) q.wk <- 10^nlambda*qq.wk
        else {
            r.wk <- cbind(r.wk,10^(ran.scal)*random$z)
            q.wk <- matrix(0,nxiz,nxiz)
            q.wk[1:nxi,1:nxi] <- 10^nlambda*qq.wk
            q.wk[(nxi+1):nxiz,(nxi+1):nxiz] <-
                10^(2*ran.scal)*random$sigma$fun(the.wk[-(1:nq)],random$sigma$env)
        }
        alpha.wk <- max(0,the.wk[1:nq]-log.th0-5)*(3-alpha) + alpha
        alpha.wk <- min(alpha.wk,3)
        z <- ngreg(dc,family,cbind(s,r.wk),q.wk,y,wt,offset,nu.wk,alpha.wk)
        assign("dc",z$dc,inherit=TRUE)
        assign("fit",z[c(1:3,5:10)],inherit=TRUE)
        z$score
    }
    cv.wk <- function(theta) cv.scale*cv(theta)+cv.shift
    ## initialization
    theta <- -log10(apply(q,3,function(x)sum(diag(x))))
    r.wk <- q.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
        q.wk <- q.wk + 10^theta[i]*q[,,i]
    }
    ## theta adjustment
    z <- sspngreg(family,s,r.wk,q.wk,y,wt,offset,alpha,nu,random)
    if (nu[[2]]) nu[[1]] <- z$nu
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
    z <- sspngreg(family,s,r.wk,q.wk,y,wt,offset,alpha,nu,random)
    if (nu[[2]]) nu[[1]] <- z$nu
    nlambda <- z$nlambda
    log.th0 <- log.th0 + z$lambda
    theta <- theta + z$theta
    if (!is.null(random)) ran.scal <- z$ran.scal
    ## theta search
    dc <- rep(0,nn)
    fit <- NULL
    if (nu[[2]]) theta <- c(theta, log(nu[[1]]))
    if (!is.null(random)) theta <- c(theta,z$zeta)
    counter <- 0
    tmp <- abs(cv(theta))
    cv.scale <- 1
    cv.shift <- 0
    if (tmp<1&tmp>10^(-4)) {
        cv.scale <- 10/tmp
        cv.shift <- 0
    }
    if (tmp<10^(-4)) {
        cv.scale <- 10^2
        cv.shift <- 10
    }
    repeat {
        zz <- nlm(cv.wk,theta,stepmax=1,ndigit=7)
        if (zz$code<=3)  break
        theta <- zz$est        
        counter <- counter + 1
        if (counter>=5) {
            warning("gss warning in gssanova: iteration for model selection fails to converge")
            break
        }
    }
    ## return
    jk <- cv(zz$est)
    if (nu[[2]]) {
        nu.wk <- exp(zz$est[nq+1])
        zz$est <- zz$est[-(nq+1)]
    }
    else nu.wk <- NULL
    r.wk <- qq.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^zz$est[i]*r[,,i]
        qq.wk <- qq.wk + 10^zz$est[i]*q[,,i]
    }
    if (is.null(random)) q.wk <- qq.wk
    else {
        r.wk <- cbind(r.wk,10^(ran.scal)*random$z)
        q.wk <- matrix(0,nxiz,nxiz)
        q.wk[1:nxi,1:nxi] <- qq.wk
        q.wk[(nxi+1):nxiz,(nxi+1):nxiz] <-
            10^(2*ran.scal-nlambda)*random$sigma$fun(zz$est[-(1:nq)],random$sigma$env)
    }
    qinv <- eigen(q.wk,TRUE)
    se.aux <- t(fit$w*cbind(s,r.wk))%*%r.wk%*%qinv$vec
    c <- fit$dc[nnull+(1:nxi)]
    if (nnull) d <- fit$dc[1:nnull]
    else d <- NULL
    if (nz) b <- 10^(ran.scal)*fit$dc[nnull+nxi+(1:nz)]
    else b <- NULL
    c(list(theta=zz$est[1:nq],c=c,d=d,b=b,nlambda=nlambda,zeta=zz$est[-(1:nq)],nu=nu.wk),
      fit[-1],list(qinv=qinv,se.aux=se.aux))
}

## Non-Gaussian regression with fixed smoothing parameters
ngreg <- function(dc,family,sr,q,y,wt,offset,nu,alpha)
{
    nobs <- nrow(sr)
    nn <- ncol(sr)
    nxi <- nrow(q)
    nnull <- nn - nxi
    ## initialization
    cc <- dc[nnull+(1:nxi)]
    eta <- sr%*%dc
    if (!is.null(offset)) eta <- eta + offset
    if ((family=="nbinomial")&is.vector(y)) y <- cbind(y,nu[[1]])
    iter <- 0
    flag <- 0
    dev <- switch(family,
                  binomial=dev.resid.binomial(y,eta,wt),
                  nbinomial=dev.resid.nbinomial(y,eta,wt),
                  poisson=dev.resid.poisson(y,eta,wt),
                  Gamma=dev.resid.Gamma(y,eta,wt),
                  weibull=dev.resid.weibull(y,eta,wt,nu[[1]]),
                  lognorm=dev0.resid.lognorm(y,eta,wt,nu[[1]]),
                  loglogis=dev0.resid.loglogis(y,eta,wt,nu[[1]]))
    dev <- sum(dev) + t(cc)%*%q%*%cc
    ## Newton iteration
    repeat {
        iter <- iter+1
        dat <- switch(family,
                      binomial=mkdata.binomial(y,eta,wt,offset),
                      nbinomial=mkdata.nbinomial(y,eta,wt,offset,nu),
                      poisson=mkdata.poisson(y,eta,wt,offset),
                      Gamma=mkdata.Gamma(y,eta,wt,offset),
                      weibull=mkdata.weibull(y,eta,wt,offset,nu),
                      lognorm=mkdata.lognorm(y,eta,wt,offset,nu),
                      loglogis=mkdata.loglogis(y,eta,wt,offset,nu))
        ## weighted least squares fit
        w <- as.vector(sqrt(dat$wt))
        ywk <- w*dat$ywk
        srwk <- w*sr
        if (!is.finite(sum(w,ywk,srwk))) {
            if (flag) stop("gss error in gssanova: Newton iteration diverges")
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
        dc.diff <- z$dc-dc
        adj <- 0
        repeat {
            dc.new <- dc + dc.diff
            cc <- dc.new[nnull+(1:nxi)]
            eta.new <- sr%*%dc.new
            if (!is.null(offset)) eta.new <- eta.new + offset
            dev.new <- switch(family,
                              binomial=dev.resid.binomial(y,eta.new,wt),
                              nbinomial=dev.resid.nbinomial(y,eta.new,wt),
                              poisson=dev.resid.poisson(y,eta.new,wt),
                              Gamma=dev.resid.Gamma(y,eta.new,wt),
                              weibull=dev.resid.weibull(y,eta.new,wt,nu[[1]]),
                              lognorm=dev0.resid.lognorm(y,eta.new,wt,nu[[1]]),
                              loglogis=dev0.resid.loglogis(y,eta.new,wt,nu[[1]]))
            dev.new <- sum(dev.new) + t(cc)%*%q%*%cc
            if (!is.finite(dev.new)) dev.new <- Inf
            if (dev.new-dev<(1+abs(dev))*1e-1) break
            adj <- 1
            dc.diff <- dc.diff/2
        }
        disc <- sum(dat$wt*((eta-eta.new)/(1+abs(eta)))^2)/sum(dat$wt)
        if (!is.finite(disc)) {
            if (flag) stop("gss error in gssanova: Newton iteration diverges")
            eta <- rep(0,nobs)
            iter <- 0
            flag <- 1
            next
        }
        dc <- dc.new
        eta <- eta.new
        dev <- dev.new
        if (adj) next
        if (disc<1e-7) break
        if (iter<=30) next
        if (!flag) {
            eta <- rep(0,nobs)
            iter <- 0
            flag <- 1
        }
        else {
            warning("gss warning in gssanova: Newton iteration fails to converge")
            break
        }
    }
    ## calculate cv
    dat <- switch(family,
                  binomial=mkdata.binomial(y,eta,wt,offset),
                  nbinomial=mkdata.nbinomial(y,eta,wt,offset,nu),
                  poisson=mkdata.poisson(y,eta,wt,offset),
                  Gamma=mkdata.Gamma(y,eta,wt,offset),
                  weibull=mkdata.weibull(y,eta,wt,offset,nu),
                  lognorm=mkdata.lognorm(y,eta,wt,offset,nu),
                  loglogis=mkdata.loglogis(y,eta,wt,offset,nu))
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
                 Gamma=cv.Gamma(y,eta,wt,z$hat[1:nobs],z$hat[nobs+1],alpha),
                 nbinomial=cv.nbinomial(y,eta,wt,z$hat[1:nobs],alpha),
                 weibull=cv.weibull(y,eta,wt,z$hat[1:nobs],nu[[1]],alpha),
                 lognorm=cv.lognorm(y,eta,wt,z$hat[1:nobs],nu[[1]],alpha),
                 loglogis=cv.loglogis(y,eta,wt,z$hat[1:nobs],nu[[1]],alpha))
    c(z,cv,list(eta=eta))
}

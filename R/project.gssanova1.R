## Calculate Kullback-Leibler projection from gssanova1 objects
project.gssanova1 <- function(object,include,...)
{
    nobs <- nrow(object$mf)
    nxi <- length(object$id.basis)
    ## evaluate full model
    family <- object$family
    eta <- object$eta
    y <- model.response(object$mf,"numeric")
    wt <- model.weights(object$mf)
    offset <- model.offset(object$mf)
    nu <- object$nu
    dat <- switch(family,
                  binomial=proj0.binomial(y,eta,wt,offset),
                  poisson=proj0.poisson(y,eta,wt,offset),
                  Gamma=proj0.Gamma(y,eta,wt,offset))
    fit0 <- dat[c("mu","theta","b")]
    y0 <- fit0$mu
    if (family=="binomial") {
        if (is.null(wt)) wt <- rep(1,nobs)
        if (!is.vector(y)) wt <- as.vector(wt*(y[,1]+y[,2]))
    }
    if (family=="nbinomial") {
        if (!is.vector(y)) y0 <- cbind(y0,y[,2])
    }
    kl0 <- switch(family,
                  binomial=dev.null.binomial(y0,wt,offset),
                  nbinomial=dev.null.nbinomial(y0,wt,offset,nu),
                  poisson=dev.null.poisson(y0,wt,offset),
                  inverse.gaussian=dev.null.inverse.gaussian(y0,wt,offset),
                  Gamma=dev.null.Gamma(y0,wt,offset),
                  weibull=dev.null.weibull(y0,wt,offset,nu),
                  lognorm=dev.null.lognorm(y0,wt,offset,nu),
                  loglogis=dev.null.loglogis(y0,wt,offset,nu))/2/nobs
    ## extract terms in subspace
    s <- matrix(1,nobs,1)
    philist <- object$term[["1"]]$iphi
    r <- NULL
    theta <- NULL
    nq.wk <- nq <- 0
    for (label in object$terms$labels) {
        if (label=="1") next
        x <- object$mf[,object$term[[label]]$vlist]
        x.basis <- object$mf[object$id.basis,object$term[[label]]$vlist]
        nphi <- object$term[[label]]$nphi
        nrk <- object$term[[label]]$nrk
        if (nphi) {
            phi <- object$term[[label]]$phi
            for (i in 1:nphi) {
                if (!any(label==include)) next
                philist <- c(philist,object$term[[label]]$iphi+(i-1))
                s <- cbind(s,phi$fun(x,nu=i,env=phi$env))
            }
        }
        if (nrk) {
            rk <- object$term[[label]]$rk
            for (i in 1:nrk) {
                nq.wk <- nq.wk + 1
                if (!any(label==include)) next
                nq <- nq + 1
                theta <- c(theta,object$theta[nq.wk])
                r <- array(c(r,rk$fun(x,x.basis,nu=i,env=rk$env,out=TRUE)),
                           c(nobs,nxi,nq))
            }
        }
    }
    if (any(include=="partial")) {
        nphi <- object$term$partial$nphi
        for (i in 1:nphi)
            philist <- c(philist,object$term$partial$iphi+(i-1))
        s <- cbind(s,object$mf$partial)
    }
    if (!is.null(object$random)) s <- cbind(s,object$random$z)
    ## calculate projection
    my.wls <- function(theta1=NULL) {
        theta.wk <- 1:nq
        theta.wk[fix] <- theta[fix]
        if (nq-1) theta.wk[-fix] <- theta1
        sr <- 0
        for (i in 1:nq) sr <- sr + 10^theta.wk[i]*r[,,i]
        q <- sr[object$id.basis,]
        sr <- cbind(s,sr)
        z <- ngreg.proj(dc,family,sr,q,fit0,wt,offset,nu)
        assign("dc",z$dc,inherit=TRUE)
        assign("fit1",z[c("mu","theta","b")],inherit=TRUE)
        mean(wt*(fit0$mu*(fit0$theta-fit1$theta)+fit1$b-fit0$b))
    }
    cv.wk <- function(theta) cv.scale*my.wls(theta)+cv.shift
    ## initialization
    r.wk <- 0
    for (i in 1:nq) r.wk <- r.wk + 10^theta[i]*r[,,i]
    if (is.null(s)) theta.wk <- 0
    else theta.wk <- log10(sum(s^2)/ncol(s)/sum(r.wk^2)*nxi) / 2
    theta <- theta + theta.wk
    tmp <- NULL
    for (i in 1:nq) tmp <- c(tmp,10^theta[i]*sum(r[cbind(object$id.basis,1:nxi,i)]))
    fix <- rev(order(tmp))[1]
    ## projection
    dc <- c(object$d[philist],object$b,10^(-theta.wk)*object$c)
    fit1 <- NULL
    if (nq-1) {
        ## scale and shift cv
        tmp <- abs(my.wls(theta[-fix]))
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
        zz <- nlm(cv.wk,theta[-fix],stepmax=.5,ndigit=7)
        if (zz$code>3)
            warning("gss warning in project.gssanova1: theta iteration fails to converge")
        kl <- my.wls(zz$est)
    }
    else kl <- my.wls()
    ## check
    kl1 <- switch(family,
                  binomial=dev.null.binomial(fit1$mu,wt,offset),
                  nbinomial=dev.null.nbinomial(fit1$mu,wt,offset,nu),
                  poisson=dev.null.poisson(fit1$mu,wt,offset),
                  inverse.gaussian=dev.null.inverse.gaussian(fit1$mu,wt,offset),
                  Gamma=dev.null.Gamma(fit1$mu,wt,offset),
                  weibull=dev.null.weibull(fit1$mu,wt,offset,nu),
                  lognorm=dev.null.lognorm(fit1$mu,wt,offset,nu),
                  loglogis=dev.null.loglogis(fit1$mu,wt,offset,nu))/2/nobs
    list(ratio=kl/kl0,kl=kl,check=(kl+kl1)/kl0)
}

## KL projection with Non-Gaussian regression
ngreg.proj <- function(dc,family,sr,q,fit0,wt,offset,nu)
{
    ## initialization
    y <- fit0$mu
    q <- 10^(-5)*q
    eta <- sr%*%dc
    if (!is.null(offset)) eta <- eta + offset
    fit1 <- switch(family,
                  binomial=proj0.binomial(y,eta,wt,offset),
                  poisson=proj0.poisson(y,eta,wt,offset),
                  Gamma=proj0.Gamma(y,eta,wt,offset))
    kl <- mean(wt*(fit0$mu*(fit0$theta-fit1$theta)+fit1$b-fit0$b))
    nobs <- length(eta)
    nn <- ncol(as.matrix(sr))
    nxi <- ncol(q)
    nnull <- nn-nxi
    iter <- 0
    flag <- 0
    ## Newton iteration
    repeat {
        iter <- iter+1
        ## weighted least squares fit
        if (!is.finite(sum(fit1$wt,fit1$ywk))) {
            if (flag) stop("gss error in project.gssanova1: Newton iteration diverges")
            eta <- rep(0,nobs)
            fit1 <- switch(family,
                           binomial=proj0.binomial(y,eta,wt,offset),
                           poisson=proj0.poisson(y,eta,wt,offset),
                           Gamma=proj0.Gamma(y,eta,wt,offset))
            kl <- mean(wt*(fit0$mu*(fit0$theta-fit1$theta)+fit1$b-fit0$b))
            iter <- 0
            flag <- 1
            next
        }
        mumax <- max(abs(t(sr)%*%fit1$u+c(rep(0,nnull),q%*%dc[-(1:nnull)])))
        w <- sqrt(as.vector(fit1$wt))
        z <- .Fortran("reg",
                      as.double(w*sr), as.integer(nobs), as.integer(nnull),
                      as.double(q), as.integer(nxi), as.double(w*fit1$ywk),
                      as.integer(4),
                      double(1), double(1), double(1), dc=double(nn),
                      as.double(.Machine$double.eps),
                      double(nn*nn), double(nn), as.integer(rep(0,nn)),
                      double(max(nobs,nn)), integer(1), integer(1),
                      PACKAGE="gss")["dc"]
        dc <- z$dc
        eta.new <- sr%*%dc
        if (!is.null(offset)) eta.new <- eta.new + offset
        fit1 <- switch(family,
                      binomial=proj0.binomial(y,eta.new,wt,offset),
                      poisson=proj0.poisson(y,eta.new,wt,offset),
                      Gamma=proj0.Gamma(y,eta.new,wt,offset))
        klnew <- mean(wt*(fit0$mu*(fit0$theta-fit1$theta)+fit1$b-fit0$b))
        disc0 <- max((mumax/(1+kl))^2,abs(klnew-kl)/(1+kl))
        disc <- sum(fit1$wt*((eta-eta.new)/(1+abs(eta)))^2)/sum(fit1$wt)
        if (is.nan(disc)) {
            if (flag) stop("gss error in project.gssanova1: Newton iteration diverges")
            eta <- rep(0,nobs)
            fit1 <- switch(family,
                           binomial=proj0.binomial(y,eta,wt,offset),
                           poisson=proj0.poisson(y,eta,wt,offset),
                           Gamma=proj0.Gamma(y,eta,wt,offset))
            kl <- mean(wt*(fit0$mu*(fit0$theta-fit1$theta)+fit1$b-fit0$b))
            iter <- 0
            flag <- 1
            next
        }
        eta <- eta.new
        kl <- klnew
        if (disc0<1e-5) break
        if (disc<1e-5) break
        if (iter<=30) next
        if (!flag) {
            eta <- rep(0,nobs)
            fit1 <- switch(family,
                           binomial=proj0.binomial(y,eta,wt,offset),
                           poisson=proj0.poisson(y,eta,wt,offset),
                           Gamma=proj0.Gamma(y,eta,wt,offset))
            kl <- mean(wt*(fit0$mu*(fit0$theta-fit1$theta)+fit1$b-fit0$b))
            iter <- 0
            flag <- 1
        }
        else {
            warning("gss warning in project.gssanova1: Newton iteration fails to converge")
            break
        }
    }
    fit1 <- switch(family,
                  binomial=proj0.binomial(y,eta,wt,offset),
                  poisson=proj0.poisson(y,eta,wt,offset),
                  Gamma=proj0.Gamma(y,eta,wt,offset))
    c(list(dc=dc),fit1[c("mu","b","theta")])
}

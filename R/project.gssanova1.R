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
    if(is.null(wt)) wt <- rep(1,nobs)
    offset <- model.offset(object$mf)
    if (!is.null(object$random)) {
        if (is.null(offset)) offset <- 0
        offset <- offset + object$random$z%*%object$b
    }
    nu <- object$nu
    y0 <- switch(family,
                 binomial=y0.binomial(y,eta,wt),
                 poisson=y0.poisson(eta),
                 Gamma=y0.Gamma(eta),
                 nbinomial=y0.nbinomial(y,eta,nu),
                 weibull=y0.weibull(y,eta,nu),
                 lognorm=y0.lognorm(y,eta,nu),
                 loglogis=y0.loglogis(y,eta,nu))
    # calculate constant fit
    cfit <- switch(family,
                   binomial=cfit.binomial(y,wt,offset),
                   poisson=cfit.poisson(y,wt,offset),
                   Gamma=cfit.Gamma(y,wt,offset),
                   nbinomial=cfit.nbinomial(y,eta,wt,nu),
                   weibull=cfit.weibull(y,wt,offset,nu),
                   lognorm=cfit.lognorm(y,wt,offset,nu),
                   loglogis=cfit.loglogis(y,wt,offset,nu))
    # calculate total entropy
    kl0 <- switch(family,
                  binomial=kl.binomial(eta,cfit,y0$wt),
                  poisson=kl.poisson(eta,cfit,wt),
                  Gamma=kl.Gamma(eta,cfit,wt),
                  nbinomial=kl.nbinomial(eta,cfit,wt,y0$nu),
                  weibull=kl.weibull(eta,cfit,wt,nu,y0$int),
                  lognorm=kl.lognorm(eta,cfit,wt,nu,y0),
                  loglogis=kl.loglogis(eta,cfit,wt,nu,y0))
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
    ## calculate projection
    my.wls <- function(theta1=NULL) {
        if (!nq) {
            q <- matrix(0)
            sr <- cbind(s,0)
            z <- ngreg.proj(dc,family,sr,q,y0,wt,offset,nu)
        }
        else {
            theta.wk <- 1:nq
            theta.wk[fix] <- theta[fix]
            if (nq-1) theta.wk[-fix] <- theta1
            sr <- 0
            for (i in 1:nq) sr <- sr + 10^theta.wk[i]*r[,,i]
            q <- sr[object$id.basis,]
            sr <- cbind(s,sr)
            z <- ngreg.proj(dc,family,sr,q,y0,wt,offset,nu)
        }
        assign("dc",z$dc,inherit=TRUE)
        assign("eta1",z$eta,inherit=TRUE)
        z$kl
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
    if (nq) dc <- c(object$d[philist],10^(-theta.wk)*object$c)
    else dc <- c(object$d[philist],0)
    eta1 <- NULL
    if (nq>1) {
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
        zz <- nlm(cv.wk,theta[-fix],stepmax=1,ndigit=7)
        if (zz$code>3)
            warning("gss warning in project.gssanova1: theta iteration fails to converge")
        kl <- my.wls(zz$est)
    }
    else kl <- my.wls()
    ## check
    kl1 <- switch(family,
                  binomial=kl.binomial(eta1,cfit,y0$wt),
                  poisson=kl.poisson(eta1,cfit,wt),
                  Gamma=kl.Gamma(eta1,cfit,wt),
                  nbinomial=kl.nbinomial(eta1,cfit,wt,y0$nu),
                  weibull=kl.weibull(eta1,cfit,wt,nu,y0$int),
                  lognorm=kl.lognorm(eta1,cfit,wt,nu,y0),
                  loglogis=kl.loglogis(eta1,cfit,wt,nu,y0))
    list(ratio=kl/kl0,kl=kl,check=(kl+kl1)/kl0)
}

## KL projection with Non-Gaussian regression
ngreg.proj <- function(dc,family,sr,q,y0,wt,offset,nu)
{
    ## initialization
    q <- 10^(-5)*q
    eta <- sr%*%dc
    nobs <- length(eta)
    nn <- ncol(as.matrix(sr))
    nxi <- ncol(q)
    nnull <- nn-nxi
    if (!is.null(offset)) eta <- eta + offset
    iter <- 0
    flag <- 0
    adj <- 0
    fit1 <- switch(family,
                   binomial=proj0.binomial(y0,eta,offset),
                   poisson=proj0.poisson(y0,eta,wt,offset),
                   Gamma=proj0.Gamma(y0,eta,wt,offset),
                   nbinomial=proj0.nbinomial(y0,eta,wt,offset),
                   weibull=proj0.weibull(y0,eta,wt,offset,nu),
                   lognorm=proj0.lognorm(y0,eta,wt,offset,nu),
                   loglogis=proj0.loglogis(y0,eta,wt,offset,nu))
    kl <- fit1$kl
    ## Newton iteration
    repeat {
        if (!adj) iter <- iter+1
        ## weighted least squares fit
        if (!is.finite(sum(fit1$wt,fit1$ywk))) {
            if (flag) stop("gss error in project.gssanova1: Newton iteration diverges")
            eta <- rep(0,nobs)
            fit1 <- switch(family,
                           binomial=proj0.binomial(y0,eta,offset),
                           poisson=proj0.poisson(y0,eta,wt,offset),
                           Gamma=proj0.Gamma(y0,eta,wt,offset),
                           nbinomial=proj0.nbinomial(y0,eta,wt,offset),
                           weibull=proj0.weibull(y0,eta,wt,offset,nu),
                           lognorm=proj0.lognorm(y0,eta,wt,offset,nu),
                           loglogis=proj0.loglogis(y0,eta,wt,offset,nu))
            kl <- fit1$kl
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
        dc.diff <- z$dc-dc
        adj <- 0
        repeat {
            dc.new <- dc + dc.diff
            eta.new <- sr%*%dc.new
            if (!is.null(offset)) eta.new <- eta.new + offset
            fit1 <- switch(family,
                           binomial=proj0.binomial(y0,eta.new,offset),
                           poisson=proj0.poisson(y0,eta.new,wt,offset),
                           Gamma=proj0.Gamma(y0,eta.new,wt,offset),
                           nbinomial=proj0.nbinomial(y0,eta.new,wt,offset),
                           weibull=proj0.weibull(y0,eta.new,wt,offset,nu),
                           lognorm=proj0.lognorm(y0,eta.new,wt,offset,nu),
                           loglogis=proj0.loglogis(y0,eta.new,wt,offset,nu))
            kl.new <- fit1$kl
            if (!is.finite(kl.new)) kl.new <- Inf
            if (kl.new-kl<(1e-4+abs(kl))*1e-1) break
            adj <- 1
            dc.diff <- dc.diff/2
        }
        disc0 <- max((mumax/(1+kl))^2,abs(kl.new-kl)/(1+kl))
        disc <- sum(fit1$wt*((eta-eta.new)/(1+abs(eta)))^2)/sum(fit1$wt)
        if (is.nan(disc)) {
            if (flag) stop("gss error in project.gssanova1: Newton iteration diverges")
            eta <- rep(0,nobs)
            fit1 <- switch(family,
                           binomial=proj0.binomial(y0,eta,offset),
                           poisson=proj0.poisson(y0,eta,wt,offset),
                           Gamma=proj0.Gamma(y0,eta,wt,offset),
                           nbinomial=proj0.nbinomial(y0,eta,wt,offset),
                           weibull=proj0.weibull(y0,eta,wt,offset,nu),
                           lognorm=proj0.lognorm(y0,eta,wt,offset,nu),
                           loglogis=proj0.loglogis(y0,eta,wt,offset,nu))
            kl <- fit1$kl
            iter <- 0
            flag <- 1
            next
        }
        dc <- dc.new
        eta <- eta.new
        kl <- kl.new
        if (adj) next
        if (disc0<1e-5) break
        if (disc<1e-5) break
        if (iter<=30) next
        warning("gss warning in project.gssanova1: Newton iteration fails to converge")
        break
    }
    fit1 <- switch(family,
                   binomial=proj0.binomial(y0,eta,offset),
                   poisson=proj0.poisson(y0,eta,wt,offset),
                   Gamma=proj0.Gamma(y0,eta,wt,offset),
                   nbinomial=proj0.nbinomial(y0,eta,wt,offset),
                   weibull=proj0.weibull(y0,eta,wt,offset,nu),
                   lognorm=proj0.lognorm(y0,eta,wt,offset,nu),
                   loglogis=proj0.loglogis(y0,eta,wt,offset,nu))
    kl <- fit1$kl
    list(dc=dc,eta=eta,kl=kl)
}

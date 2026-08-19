## Calculate Kullback-Leibler projection from gssanova objects
project1.gssanova <- function(object,include,nrep=200,type=1,...)
{
    if (class(object)[1]=="gssanova0")
        stop("gss error: Kullback-Leibler projection is not implemented for gssanova0")
    nobs <- nrow(object$mf)
    nxi <- length(object$id.basis)
    labels.p <- object$lab.p
    ## evaluate full model
    family <- object$family
    eta <- object$eta
    if (object$family=="polr") {
        y <- model.response(object$mf)
        if (!is.factor(y))
            stop("gss error in gssanova1: need factor response for polr family")
        lvls <- levels(y)
        if (nlvl <- length(lvls)<3)
            stop("gss error in gssanova1: need at least 3 levels to fit polr family")
        y <- outer(y,lvls,"==")
    }
    else y <- model.response(object$mf,"numeric")
    wt <- model.weights(object$mf)
    if(is.null(wt)) wt <- rep(1,nobs)
    offset <- model.offset(object$mf)
    if (!is.null(object$random)) {
        if (is.null(offset)) offset <- 0
        offset <- offset + object$random$z%*%object$b
    }
    ## extract terms in subspace
    s <- matrix(1,nobs,1)
    philist <- object$term[["1"]]$iphi
    r <- NULL
    idx.th <- NULL
    nq.wk <- nq <- 0
    for (label in object$terms$labels) {
        if (label=="1") next
        if (label%in%labels.p) next
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
                idx.th <- c(idx.th,nq.wk)
                r <- array(c(r,rk$fun(x,x.basis,nu=i,env=rk$env,out=TRUE)),
                           c(nobs,nxi,nq))
            }
        }
    }
    if (!is.null(object$partial)) {
        nu <- length(object$d)-length(object$lab.p)
        matx.p <- model.matrix(object$partial$mt,object$mf)[,-1,drop=FALSE]
        matx.p <- scale(matx.p)
        for (label in labels.p) {
            nu <- nu+1
            if (!any(label==include)) next
            philist <- c(philist,nu)
            s <- cbind(s,matx.p[,label])
        }
    }
    ## calculate projection
    my.proj <- function(y,eta,offset,theta) {
        y0 <- switch(family,
                     binomial=y0.binomial(y,eta,wt),
                     poisson=y0.poisson(eta),
                     Gamma=y0.Gamma(eta),
                     inverse.gaussian=y0.inverse.gaussian(eta),
                     nbinomial=y0.nbinomial(y,eta,nu),
                     polr=y0.polr(list(eta=eta,nu=nu)),
                     weibull=y0.weibull(y,eta,nu),
                     lognorm=y0.lognorm(y,eta,nu),
                     loglogis=y0.loglogis(y,eta,nu))
        ## calculate constant fit
        cfit <- switch(family,
                       binomial=cfit.binomial(y,wt,offset),
                       poisson=cfit.poisson(y,wt,offset),
                       Gamma=cfit.Gamma(y,wt,offset),
                       inverse.gaussian=cfit.inverse.gaussian(y,wt,offset),
                       nbinomial=cfit.nbinomial(y,wt,offset,nu),
                       polr=cfit.polr(y,wt,offset,nu),
                       weibull=cfit.weibull(y,wt,offset,nu),
                       lognorm=cfit.lognorm(y,wt,offset,nu),
                       loglogis=cfit.loglogis(y,wt,offset,nu))
        ## calculate total entropy
        kl0 <- switch(family,
                      binomial=kl.binomial(eta,cfit,y0$wt),
                      poisson=kl.poisson(eta,cfit,wt),
                      Gamma=kl.Gamma(eta,cfit,wt),
                      inverse.gaussian=kl.inverse.gaussian(eta,cfit,wt),
                      nbinomial=kl.nbinomial(eta,cfit,wt,y0$nu),
                      polr=kl.polr(list(eta=eta,nu=nu),cfit,wt),
                      weibull=kl.weibull(eta,cfit,wt,nu,y0$int),
                      lognorm=kl.lognorm(eta,cfit,wt,nu,y0),
                      loglogis=kl.loglogis(eta,cfit,wt,nu,y0))
        ## projection
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
            assign("dc",z$dc,inherits=TRUE)
            assign("eta1",z$eta,inherits=TRUE)
            z$kl
        }
        cv.wk <- function(theta) cv.scale*my.wls(theta)+cv.shift
        ## initialization
        if (nq) {
            r.wk <- 0
            for (i in 1:nq) r.wk <- r.wk + 10^theta[i]*r[,,i]
            if (is.null(s)) theta.wk <- 0
            else theta.wk <- log10(sum(wt*s^2)/ncol(s)/sum(wt*r.wk^2)*nxi) / 2
            theta <- theta + theta.wk
            tmp <- NULL
            for (i in 1:nq) tmp <- c(tmp,10^theta[i]*sum(r[cbind(object$id.basis,1:nxi,i)]))
            fix <- rev(order(tmp))[1]
        }
        ## projection
        if (nq) dc <- c(object$d[philist],10^(-theta.wk)*object$c)
        else dc <- c(object$d[philist],0)
        eta1 <- NULL
        if (nq>1) {
            if (object$skip.iter) kl <- my.wls(theta[-fix])
            else {
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
                kl <- my.wls(zz$est)
            }
        }
        else kl <- my.wls()
        ## check
        kl1 <- switch(family,
                      binomial=kl.binomial(eta1,cfit,y0$wt),
                      poisson=kl.poisson(eta1,cfit,wt),
                      Gamma=kl.Gamma(eta1,cfit,wt),
                      inverse.gaussian=kl.inverse.gaussian(eta1,cfit,wt),
                      nbinomial=kl.nbinomial(eta1,cfit,wt,y0$nu),
                      polr=kl.polr(list(eta=eta1,nu=nu),cfit,wt),
                      weibull=kl.weibull(eta1,cfit,wt,nu,y0$int),
                      lognorm=kl.lognorm(eta1,cfit,wt,nu,y0),
                      loglogis=kl.loglogis(eta1,cfit,wt,nu,y0))
        list(ratio=kl/kl0,kl=kl,check=(kl+kl1)/kl0,eta1=eta1)
    }
    ## initialization
    if (nq) dc <- c(object$d[philist],object$c)
    else dc <- c(object$d[philist],0)
    nu <- object$nu
    proj0 <- my.proj(y,eta,offset,object$theta[idx.th])
    ## cloning
    cln <- clone.gssanova(object,nrep,type,,proj0$eta1)
    ratio.wk <- NULL
    offset <- model.offset(object$mf)
    for (i in 1:nrep) {
        if (is.null(object$random)) offset.wk <- offset
        else {
            if (is.null(offset)) offset.wk <- object$random$z%*%cln$b[,i]
            else offset.wk <- offset + object$random$z%*%cln$b[,i]
        }
        proj.wk <- suppressWarnings(my.proj(cln$y[[i]],cln$eta[[i]],
                                            offset.wk,cln$theta[idx.th,i]))
        ratio.wk <- c(ratio.wk,proj.wk$ratio)
    }
    pvalue <- mean(proj0$ratio<ratio.wk)
    list(ratio=proj0$ratio,kl=proj0$kl,check=proj0$check,pvalue=pvalue)
}

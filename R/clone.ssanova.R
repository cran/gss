## S3 method
clone <- function (object,...) UseMethod("clone")
## Fit ssanova model
clone.ssanova <- function(object,nrep=1000,type=1,seed=NULL,eta=NULL,fix=FALSE,...)
{
    ## Extract relevant infomation from object
    mf <- object$mf
    term <- object$terms
    offset <- model.offset(mf)
    if (is.null(offset)) offset <- 0
    wt <- model.weights(mf)
    if (!is.null(wt)) wt <- sqrt(wt)
    ## Generate s and r
    s <- r <- NULL
    nobs <- dim(mf)[1]
    nxi <- length(object$id.basis)
    nq <- 0
    for (label in term$labels) {
        if (label=="1") {
            s <- cbind(s,rep(1,len=nobs))
            next
        }
        x <- mf[,term[[label]]$vlist]
        x.basis <- mf[object$id.basis,term[[label]]$vlist]
        nphi <- term[[label]]$nphi
        nrk <- term[[label]]$nrk
        if (nphi) {
            phi <- term[[label]]$phi
            for (i in 1:nphi)
                s <- cbind(s,phi$fun(x,nu=i,env=phi$env))
        }
        if (nrk) {
            rk <- term[[label]]$rk
            for (i in 1:nrk) {
                nq <- nq+1
                r <- array(c(r,rk$fun(x,x.basis,nu=i,env=rk$env,out=TRUE)),c(nobs,nxi,nq))
            }
        }
    }
    ## Partial term in s
    if (!is.null(object$partial)) {
        matx.p <- model.matrix(object$partial$mt.p)[,-1,drop=FALSE]
        matx.p <- scale(matx.p)
        s <- cbind(s,matx.p)
    }
    if (!is.null(s)) {
        if (is.vector(s)) nnull <- 1
        else nnull <- ncol(s)
    }
    else nnull <- 0
    ## Aggregate r
    if (is.matrix(r)) r <- array(r,c(nobs,nxi,1))
    rwk <- 0
    for (i in 1:nq) rwk <- rwk+10^object$theta[i]*r[,,i]
    ## fitted fixed effects
    if (is.null(eta)) eta <- as.vector(s%*%object$d+rwk%*%object$c)+offset
    ## Random term
    if (!is.null(random<-object$random)) {
        nz <- ncol(as.matrix(random$z))
        sig.z <- random$sigma$fun(object$zeta,random$sigma$env)
        q.z <- 10^(2*object$ran.scal)*sig.z
        rwk.z <- 10^object$ran.scal*random$z
        sig.z <- t(chol(sig.z))
    }
    else nz <- 0
    ##
    nxiz <- nxi + nz
    nn <- nxiz + nnull
    if (!is.null(wt)) swk <- wt*s
    else swk <- s
    if (fix) {
        q <- rwk[object$id.basis,]
        if (!is.null(wt)) rwkwk <- wt*rwk
        else rwkwk <- rwk
    }
    ## Component roughness and weights
    qq <- r[object$id.basis,,,drop=FALSE]
    rho0 <- pwt <- NULL
    for (i in 1:nq) {
        rho0 <- c(rho0,log10(sum(object$c*as.vector(qq[,,i]%*%object$c))))
        pwt <- c(pwt,var(r[,,i]%*%object$c)*10^(2*object$theta[i]))
    }
    pwt <- pwt/sum(pwt)
    ## fixed theta
    fix1 <- function(nlam) {
        if (is.null(random)) qwk <- 10^(nlam)*q
        else {
            rwk <- cbind(rwk,rwk.z)
            qwk <- matrix(0,nxiz,nxiz)
            qwk[1:nxi,1:nxi] <- 10^(nlam)*q
            qwk[(nxi+1):nxiz,(nxi+1):nxiz] <- q.z
        }
        z <- .Fortran("reg",
                      as.double(cbind(swk,rwkwk)), as.integer(nobs), as.integer(nnull),
                      as.double(qwk), as.integer(nxiz), as.double(y),
                      as.integer(4), as.double(1), varht=as.double(1),
                      score=double(1), dc=double(nn),
                      as.double(.Machine$double.eps),
                      chol=double(nn*nn), double(nn),
                      jpvt=as.integer(c(rep(1,nnull),rep(0,nxiz))),
                      wk=double(3*nobs+nnull+nz), rkv=integer(1), info=integer(1),
                      PACKAGE="gss")["dc"]
        c <- z$dc[nnull+(1:nxi)]
        d <- z$dc[1:nnull]
        eta1 <- as.vector(s%*%d+rwk%*%c)+offset
        assign("dc.wk",z$dc,inherits=TRUE)
        assign("eta.wk",eta1,inherits=TRUE)
        if (type-1) {
            if (is.null(wt)) mean((eta-eta1)^2)
            else sum((wt*(eta-eta1))^2)/sum(wt^2)
        }
        else {
            rho <- NULL
            for (i in 1:nq) rho <- c(rho,log10(sum(c*as.vector(qq[,,i]%*%c))))
            sum(pwt*(rho-rho0)^2)
        }
    }
    ## varying theta
    fix0 <- function(theta) {
        ## Aggregate r
        nq <- length(theta)
        rwk <- 0
        for (i in 1:nq) rwk <- rwk+10^theta[i]*r[,,i]
        q <- rwk[object$id.basis,]
        if (is.null(random)) qwk <- 10^(nlam)*q
        else {
            rwk <- cbind(rwk,rwk.z)
            qwk <- matrix(0,nxiz,nxiz)
            qwk[1:nxi,1:nxi] <- 10^(nlam)*q
            qwk[(nxi+1):nxiz,(nxi+1):nxiz] <- q.z
        }
        if (!is.null(wt)) rwkwk <- wt*rwk
        else rwkwk <- rwk
        z <- .Fortran("reg",
                      as.double(cbind(swk,rwkwk)), as.integer(nobs), as.integer(nnull),
                      as.double(qwk), as.integer(nxiz), as.double(y),
                      as.integer(4), as.double(1), varht=as.double(1),
                      score=double(1), dc=double(nn),
                      as.double(.Machine$double.eps),
                      chol=double(nn*nn), double(nn),
                      jpvt=as.integer(c(rep(1,nnull),rep(0,nxiz))),
                      wk=double(3*nobs+nnull+nz), rkv=integer(1), info=integer(1),
                      PACKAGE="gss")["dc"]
        c <- z$dc[nnull+(1:nxi)]
        d <- z$dc[1:nnull]
        eta1 <- as.vector(s%*%d+rwk%*%c)+offset
        assign("dc.wk",z$dc,inherits=TRUE)
        assign("eta.wk",eta1,inherits=TRUE)
        if (type-1) {
            if (is.null(wt)) mean((eta-eta1)^2)
            else sum((wt*(eta-eta1))^2)/sum(wt^2)
        }
        else {
            rho <- NULL
            for (i in 1:nq) rho <- c(rho,log10(sum(c*as.vector(qq[,,i]%*%c))))
            sum(pwt*(rho-rho0+2*(theta-object$theta))^2)
        }
    }
    ## cloning
    dc <- mse <- theta <- NULL
    sig <- sqrt(object$varht)
    nlam <- object$nlambda
    if (is.null(seed)) {
        seed <- abs(object$c[1])+abs(rev(object$c)[1])
        seed <- round((seed/ceiling(seed))^2,6)*1000000
        seed <- seed%%10000
    }
    set.seed(seed)
    for (i in 1:nrep) {
        dc.wk <- eta.wk <- NULL
        ee <- rnorm(eta)*sig
        if (!is.null(wt)) ee <- ee/wt
        if (!is.null(random)) ee <- ee+random$z%*%(sig.z%*%rnorm(nz))
        y <- eta+ee-offset
        if (!is.null(wt)) y <- wt*y
        if (fix) z <- nlm(fix1,nlam,stepmax=0.5)
        else z <- nlm(fix0,object$theta,stepmax=0.5)
        if (type-1) mse.wk <- as.vector(z[[1]])
        else {
            if (is.null(wt)) mse.wk <- mean((eta-eta.wk)^2)
            else mse.wk <- sum((wt*(eta-eta.wk))^2)/sum(wt^2)
        }
        if (!fix) theta <- cbind(theta,as.vector(z[[2]]))
        mse <- c(mse,mse.wk)
        dc <- cbind(dc,dc.wk)
    }
    c <- dc[nnull+(1:nxi),]
    if (nnull) d <- dc[1:nnull,,drop=FALSE]
    else d <- NULL
    if (nz) b <- 10^(object$ran.scal)*dc[nnull+nxi+(1:nz),,drop=FALSE]
    else b <- NULL
    obj <- list(fit=object,nrep=nrep,c=c,d=d,b=b,kl=mse,type=type,theta=theta)
    class(obj) <- "clone"
    obj
}

## Calculate Kullback-Leibler projection from ssanova1 objects
project.ssanova1 <- function(object,include,...)
{
    nobs <- nrow(object$mf)
    nxi <- length(object$id.basis)
    ## evaluate full model
    mf <- object$mf
    yy <- predict(object,mf)
    wt <- model.weights(object$mf)
    offset <- model.offset(object$mf)
    if (!is.null(offset)) yy <- yy - offset
    ## extract terms in subspace
    s <- matrix(1,nobs,1)
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
    if (any(include=="partial")) s <- cbind(s,object$mf$partial)
    ## calculate projection
    my.ls <- function(theta1=NULL) {
        if (!nq) {
            q <- matrix(0)
            sr <- cbind(s,0)
        }
        else {
            theta.wk <- 1:nq
            theta.wk[fix] <- theta[fix]
            if (nq-1) theta.wk[-fix] <- theta1
            sr <- 0
            for (i in 1:nq) sr <- sr + 10^theta.wk[i]*r[,,i]
            q <- 10^(-5)*sr[object$id.basis,]
            sr <- cbind(s,sr)
        }
        nn <- ncol(as.matrix(sr))
        nnull <- nn-nxi
        if (!is.null(wt)) {
            wt <- sqrt(wt)
            sr <- wt*sr
            yy <- wt*yy
        }
        z <- .Fortran("reg",
                      as.double(sr), as.integer(nobs), as.integer(nnull),
                      as.double(q), as.integer(nxi), as.double(yy),
                      as.integer(4),
                      double(1), double(1), double(1), dc=double(nn),
                      as.double(.Machine$double.eps),
                      double(nn*nn), double(nn), as.integer(rep(0,nn)),
                      double(max(nobs,nn)), integer(1), integer(1),
                      PACKAGE="gss")["dc"]
        assign("yhat",sr%*%z$dc,inherit=TRUE)
        mean((yy-yhat)^2)
    }
    cv.wk <- function(theta) cv.scale*my.ls(theta)+cv.shift
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
    yhat <- NULL
    if (nq>1) {
        ## scale and shift cv
        tmp <- abs(my.ls(theta[-fix]))
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
            warning("gss warning in project.ssanova1: theta iteration fails to converge")
        kl <- my.ls(zz$est)
    }
    else kl <- my.ls()
    kl0 <- mean((yy-mean(yy))^2)
    kl <- mean((yy-yhat)^2)
    kl1 <- mean((mean(yy)-yhat)^2)
    list(ratio=kl/kl0,kl=kl,check=(kl+kl1)/kl0)
}

## Fit Single Smoothing Parameter (Gaussian) REGression
sspreg1 <- function(s,r,q,y,method,alpha,varht,random)
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
        z <- .Fortran("reg",
                      as.double(cbind(s,10^theta*r)), as.integer(nobs), as.integer(nnull),
                      as.double(q.wk), as.integer(nxiz), as.double(y),
                      as.integer(switch(method,"u"=1,"v"=2,"m"=3)),
                      as.double(alpha), varht=as.double(varht),
                      score=double(1), dc=double(nn),
                      as.double(.Machine$double.eps),
                      chol=double(nn*nn), double(nn),
                      jpvt=as.integer(c(rep(1,nnull),rep(0,nxiz))),
                      wk=double(nobs+nnull+nz), rkv=integer(1), info=integer(1),
                      PACKAGE="gss")[c("score","varht","dc","chol","jpvt","wk","rkv","info")]
        if (z$info) stop("gss error in ssanova: evaluation of GML score fails")
        assign("fit",z[c(1:5,7)],inherit=TRUE)
        score <- z$score
        alpha.wk <- max(0,log.la0-lambda[1]-5)*(3-alpha) + alpha
        alpha.wk <- min(alpha.wk,3)
        if (alpha.wk>alpha) {
            if (method=="u") score <- score + (alpha.wk-alpha)*2*varht*z$wk[2]
            if (method=="v") score <- z$wk[1]/(1-alpha.wk*z$wk[2])^2
        }
        score
    }
    ## initialization
    tmp <- sum(r^2)
    if (is.null(s)) theta <- 0
    else theta <- log10(sum(s^2)/nnull/tmp*nxi) / 2
    log.la0 <- log10(tmp/sum(diag(q))) + theta
    if (!is.null(random)) r <- cbind(r,10^(-theta)*random$z)
    ## lambda search
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
    jk1 <- cv(zz$est)
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
    se.aux <- t(cbind(s,10^theta*r))%*%(10^theta*r)%*%qinv
    c <- fit$dc[nnull+(1:nxi)]
    if (nnull) d <- fit$dc[1:nnull]
    else d <- NULL
    if (nz) b <- fit$dc[nnull+nxi+(1:nz)]
    else b <- NULL
    c(list(method=method,theta=theta,c=c,d=d,b=b,nlambda=zz$est[1],zeta=zz$est[-1]),
      fit[-3],list(qinv=qinv,se.aux=se.aux))
}

## Fit Multiple Smoothing Parameter (Gaussian) REGression
mspreg1 <- function(s,r,q,y,method,alpha,varht,random)
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
        z <- .Fortran("reg",
                      as.double(cbind(s,r.wk)), as.integer(nobs), as.integer(nnull),
                      as.double(q.wk), as.integer(nxiz), as.double(y),
                      as.integer(switch(method,"u"=1,"v"=2,"m"=3)),
                      as.double(alpha), varht=as.double(varht),
                      score=double(1), dc=double(nn),
                      as.double(.Machine$double.eps),
                      chol=double(nn*nn), double(nn),
                      jpvt=as.integer(c(rep(1,nnull),rep(0,nxiz))),
                      wk=double(nobs+nnull+nz), rkv=integer(1), info=integer(1),
                      PACKAGE="gss")[c("score","varht","dc","chol","jpvt","wk","rkv","info")]
        if (z$info) stop("gss error in ssanova: evaluation of GML score fails")
        assign("fit",z[c(1:5,7)],inherit=TRUE)
        score <- z$score
        alpha.wk <- max(0,theta[1:nq]-log.th0-5)*(3-alpha) + alpha
        alpha.wk <- min(alpha.wk,3)
        if (alpha.wk>alpha) {
            if (method=="u") score <- score + (alpha.wk-alpha)*2*varht*z$wk[2]
            if (method=="v") score <- z$wk[1]/(1-alpha.wk*z$wk[2])^2
        }
        score
    }
    ## initialization
    theta <- -log10(apply(q,3,function(x)sum(diag(x))))
    r.wk <- q.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
        q.wk <- q.wk + 10^theta[i]*q[,,i]
    }
    ## theta adjustment
    z <- sspreg1(s,r.wk,q.wk,y,method,alpha,varht,random)
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
    z <- sspreg1(s,r.wk,q.wk,y,method,alpha,varht,random)
    nlambda <- z$nlambda
    log.th0 <- log.th0 + z$lambda
    theta <- theta + z$theta
    ## theta search
    fit <- NULL
    if (!is.null(random)) theta <- c(theta,z$zeta)
    counter <- 0
    repeat {
        zz <- nlm(cv,theta,stepmax=1,ndigit=7)
        if (zz$code<=3)  break
        theta <- zz$est        
        counter <- counter + 1
        if (counter>=5) {
            warning("gss warning in ssanova1: iteration for model selection fails to converge")
            break
        }
    }
    ## return
    jk1 <- cv(zz$est)
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
    se.aux <- t(cbind(s,r.wk))%*%r.wk%*%qinv
    c <- fit$dc[nnull+(1:nxi)]
    if (nnull) d <- fit$dc[1:nnull]
    else d <- NULL
    if (nz) b <- fit$dc[nnull+nxi+(1:nz)]
    else b <- NULL
    c(list(method=method,theta=zz$est[1:nq],c=c,d=d,b=b,nlambda=nlambda,zeta=zz$est[-(1:nq)]),
      fit[-3],list(qinv=qinv,se.aux=se.aux))
}

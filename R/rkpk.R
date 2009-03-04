## Fit Single Smoothing Parameter (Gaussian) REGression
sspreg1 <- function(s,r,q,y,method,alpha,varht,random)
{
    qr.trace <- FALSE
    if ((alpha<0)&(method%in%c("u","v"))) qr.trace <- TRUE
    alpha <- abs(alpha)
    ## get dimensions
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
            q.wk[(nxi+1):nxiz,(nxi+1):nxiz] <-
                10^(2*ran.scal)*random$sigma$fun(lambda[-1],random$sigma$env)
        }
        if (qr.trace) {
            qq.wk <- chol(q.wk,pivot=TRUE)
            sr <- cbind(s,10^theta*r[,attr(qq.wk,"pivot")])
            sr <- rbind(sr,cbind(matrix(0,nxiz,nnull),qq.wk))
            sr <- qr(sr,tol=0)
            rss <- mean(qr.resid(sr,c(y,rep(0,nxiz)))[1:nobs]^2)
            trc <- sum(qr.Q(sr)[1:nobs,]^2)/nobs
            if (method=="u") score <- rss + alpha*2*varht*trc
            if (method=="v") score <- rss/(1-alpha*trc)^2
            alpha.wk <- max(0,log.la0-lambda[1]-5)*(3-alpha) + alpha
            alpha.wk <- min(alpha.wk,3)
            if (alpha.wk>alpha) {
                if (method=="u") score <- score + (alpha.wk-alpha)*2*varht*trc
                if (method=="v") score <- rss/(1-alpha.wk*trc)^2
            }
            if (return.fit) {
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
                z$score <- score
                assign("fit",z[c(1:5,7)],inherit=TRUE)
            }
        }
        else {
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
        }
        score
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
    ## lambda search
    return.fit <- FALSE
    fit <- NULL
    if (is.null(random)) la <- log.la0
    else la <- c(log.la0,random$init)
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
    return.fit <- TRUE
    jk1 <- cv(zz$est)
    if (is.null(random)) q.wk <- 10^theta*q
    else {
        q.wk <- matrix(0,nxiz,nxiz)
        q.wk[1:nxi,1:nxi] <- 10^theta*q
        q.wk[(nxi+1):nxiz,(nxi+1):nxiz] <-
            10^(2*ran.scal-zz$est[1])*random$sigma$fun(zz$est[-1],random$sigma$env)
    }
    zzz <- eigen(q.wk,TRUE)
    rkq <- min(fit$rkv-nnull,sum(zzz$val/zzz$val[1]>sqrt(.Machine$double.eps)))
    val <- zzz$val[1:rkq]
    vec <- zzz$vec[,1:rkq,drop=FALSE]
    qinv <- vec%*%diag(1/val,rkq)%*%t(vec)
    if (nnull) {
        qr.s <- qr(s)
        wk1 <- (qr.qty(qr.s,10^theta*r))[-(1:nnull),]
    }
    else wk1 <- 10^theta*r
    wk2 <- wk1%*%qinv%*%t(wk1)
    diag(wk2) <- diag(wk2) + 10^zz$est[1]
    wk2 <- chol(wk2)
    se.aux0 <- backsolve(wk2,wk1%*%qinv,trans=TRUE)
    se.aux <- t(cbind(s,10^theta*r))%*%(10^theta*r)%*%qinv
    c <- fit$dc[nnull+(1:nxi)]
    if (nnull) d <- fit$dc[1:nnull]
    else d <- NULL
    if (nz) b <- 10^(ran.scal)*fit$dc[nnull+nxi+(1:nz)]
    else b <- NULL
    c(list(method=method,theta=theta,ran.scal=ran.scal,c=c,d=d,b=b,
           nlambda=zz$est[1],zeta=zz$est[-1]),
      fit[-3],list(qinv=qinv,se.aux=list(se.aux,se.aux0)))
}

## Fit Multiple Smoothing Parameter (Gaussian) REGression
mspreg1 <- function(s,r,q,y,method,alpha,varht,random)
{
    qr.trace <- FALSE
    if ((alpha<0)&(method%in%c("u","v"))) qr.trace <- TRUE
    alpha <- abs(alpha)
    ## get dimensions
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
            r.wk <- cbind(r.wk,10^(ran.scal)*random$z)
            q.wk <- matrix(0,nxiz,nxiz)
            q.wk[1:nxi,1:nxi] <- 10^nlambda*qq.wk
            q.wk[(nxi+1):nxiz,(nxi+1):nxiz] <-
                10^(2*ran.scal)*random$sigma$fun(theta[-(1:nq)],random$sigma$env)
        }
        if (qr.trace) {
            qq.wk <- chol(q.wk,pivot=TRUE)
            sr <- cbind(s,r.wk[,attr(qq.wk,"pivot")])
            sr <- rbind(sr,cbind(matrix(0,nxiz,nnull),qq.wk))
            sr <- qr(sr,tol=0)
            rss <- mean(qr.resid(sr,c(y,rep(0,nxiz)))[1:nobs]^2)
            trc <- sum(qr.Q(sr)[1:nobs,]^2)/nobs
            if (method=="u") score <- rss + alpha*2*varht*trc
            if (method=="v") score <- rss/(1-alpha*trc)^2
            alpha.wk <- max(0,theta[1:nq]-log.th0-5)*(3-alpha) + alpha
            alpha.wk <- min(alpha.wk,3)
            if (alpha.wk>alpha) {
                if (method=="u") score <- score + (alpha.wk-alpha)*2*varht*trc
                if (method=="v") score <- rss/(1-alpha.wk*trc)^2
            }
            if (return.fit) {
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
                z$score <- score
                assign("fit",z[c(1:5,7)],inherit=TRUE)
            }
        }
        else {
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
        }
        score
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
    return.fit <- FALSE
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
    log.th0 <- log.th0 + z$nlambda
    theta <- theta + z$theta
    if (!is.null(random)) ran.scal <- z$ran.scal
    ## theta search
    fit <- NULL
    if (!is.null(random)) theta <- c(theta,z$zeta)
    counter <- 0
    ## scale and shift cv
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
            warning("gss warning in ssanova: iteration for model selection fails to converge")
            break
        }
    }
    ## return
    return.fit <- TRUE
    jk1 <- cv(zz$est)
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
    zzz <- eigen(q.wk,TRUE)
    rkq <- min(fit$rkv-nnull,sum(zzz$val/zzz$val[1]>sqrt(.Machine$double.eps)))
    val <- zzz$val[1:rkq]
    vec <- zzz$vec[,1:rkq,drop=FALSE]
    qinv <- vec%*%diag(1/val,rkq)%*%t(vec)
    if (nnull) {
        qr.s <- qr(s)
        wk1 <- (qr.qty(qr.s,r.wk))[-(1:nnull),]
    }
    else wk1 <- r.wk
    wk2 <- wk1%*%qinv%*%t(wk1)
    diag(wk2) <- diag(wk2) + 10^nlambda
    wk2 <- chol(wk2)
    se.aux0 <- backsolve(wk2,wk1%*%qinv,trans=TRUE)
    se.aux <- t(cbind(s,r.wk))%*%r.wk%*%qinv
    c <- fit$dc[nnull+(1:nxi)]
    if (nnull) d <- fit$dc[1:nnull]
    else d <- NULL
    if (nz) b <- 10^(ran.scal)*fit$dc[nnull+nxi+(1:nz)]
    else b <- NULL
    c(list(method=method,theta=zz$est[1:nq],c=c,d=d,b=b,nlambda=nlambda,zeta=zz$est[-(1:nq)]),
      fit[-3],list(qinv=qinv,se.aux=list(se.aux,se.aux0)))
}

## Fit single smoothing parameter density
sspdsty <- function(s,r,q,cnt,qd.s,qd.r,qd.wt,prec,maxiter,alpha)
{
    nxi <- dim(r)[1]
    nobs <- dim(r)[2]
    nqd <- length(qd.wt)
    if (!is.null(s)) nnull <- dim(s)[1]
    else nnull <- 0
    nxis <- nxi+nnull
    if (is.null(cnt)) cnt <- 0
    ## cv function
    cv <- function(lambda) {
        fit <- .Fortran("dnewton",
                        cd=as.double(cd), as.integer(nxis),
                        as.double(10^lambda*q), as.integer(nxi),
                        as.double(rbind(r,s)), as.integer(nobs),
                        as.double(sum(cnt)), as.integer(cnt),
                        as.double(t(rbind(qd.r,qd.s))), as.integer(nqd),
                        as.double(qd.wt),
                        as.double(prec), as.integer(maxiter),
                        as.double(.Machine$double.eps),
                        wk=double(2*(nqd+nobs)+nxis*(nxis+4)+max(nxis,3)),
                        info=integer(1),PACKAGE="gss")
        if (fit$info==1) stop("gss error in ssden: Newton iteration diverges")
        if (fit$info==2) warning("gss warning in ssden: Newton iteration fails to converge")
        assign("cd",fit$cd,inherit=TRUE)
        assign("int",fit$wk[3],inherit=TRUE)
        cv <- alpha*fit$wk[2]-fit$wk[1]
        alpha.wk <- max(0,log.la0-lambda-3.5)*(3-alpha) + alpha
        alpha.wk <- min(alpha.wk,3)
        adj <- ifelse (alpha.wk>alpha,(alpha.wk-alpha)*fit$wk[2],0)
        cv+adj
    }
    ## initialization
    mu <- apply(qd.wt*t(qd.r),2,sum)/sum(qd.wt)
    v <- apply(qd.wt*t(qd.r^2),2,sum)/sum(qd.wt)
    log.la0 <- log10(sum(v-mu^2)/sum(diag(q)))
    ## lambda search
    cd <- rep(0,nxi+nnull)
    int <- NULL
    la <- log.la0
    counter <- 0
    repeat {
        fit <- nlm(cv,la,stepmax=1,ndigit=7)
        if (fit$code<=3) break
        la <- fit$est
        counter <- counter + 1
        if (counter>=5) {
            warning("gss warning in ssden: CV iteration fails to converge")
            break
        }
    }
    ## return
    jk1 <- cv(fit$est)
    c <- cd[1:nxi]
    if (nnull) d <- cd[nxi+(1:nnull)]
    else d <- NULL
    list(lambda=fit$est,theta=0,c=c,d=d,int=int,cv=fit$min)
}

## Fit multiple smoothing parameter density
mspdsty <- function(s,r,q,cnt,qd.s,qd.r,qd.wt,prec,maxiter,alpha)
{
    nxi <- dim(r)[1]
    nobs <- dim(r)[2]
    nqd <- length(qd.wt)
    nq <- dim(q)[3]
    if (!is.null(s)) nnull <- dim(s)[1]
    else nnull <- 0
    nxis <- nxi+nnull
    if (is.null(cnt)) cnt <- 0
    ## cv function
    cv <- function(theta) {
        r.wk <- q.wk <- qd.r.wk <- 0
        for (i in 1:nq) {
            r.wk <- r.wk + 10^theta[i]*r[,,i]
            q.wk <- q.wk + 10^theta[i]*q[,,i]
            qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[,,i]
        }
        fit <- .Fortran("dnewton",
                        cd=as.double(cd), as.integer(nxis),
                        as.double(q.wk), as.integer(nxi),
                        as.double(rbind(r.wk,s)), as.integer(nobs),
                        as.double(sum(cnt)), as.integer(cnt),
                        as.double(t(rbind(qd.r.wk,qd.s))), as.integer(nqd),
                        as.double(qd.wt),
                        as.double(prec), as.integer(maxiter),
                        as.double(.Machine$double.eps),
                        wk=double(2*(nqd+nobs)+nxis*(nxis+4)+max(nxis,3)),
                        info=integer(1),PACKAGE="gss")
        if (fit$info==1) stop("gss error in ssden: Newton iteration diverges")
        if (fit$info==2) warning("gss warning in ssden: Newton iteration fails to converge")
        assign("cd",fit$cd,inherit=TRUE)
        assign("int",fit$wk[3],inherit=TRUE)
        cv <- alpha*fit$wk[2]-fit$wk[1]
        alpha.wk <- max(0,theta-log.th0-3.5)*(3-alpha) + alpha
        alpha.wk <- min(alpha.wk,3)
        adj <- ifelse (alpha.wk>alpha,(alpha.wk-alpha)*fit$wk[2],0)
        cv+adj
    }
    ## initialization
    theta <- -log10(apply(q,3,function(x)sum(diag(x))))
    r.wk <- q.wk <- qd.r.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
        q.wk <- q.wk + 10^theta[i]*q[,,i]
        qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[,,i]
    }
    mu <- apply(qd.wt*t(qd.r.wk),2,sum)/sum(qd.wt)
    v <- apply(qd.wt*t(qd.r.wk^2),2,sum)/sum(qd.wt)
    log.la0 <- log10(sum(v-mu^2)/sum(diag(q.wk)))
    log.th0 <- theta-log.la0
    ## lambda search
    z <- sspdsty(s,r.wk,q.wk,cnt,qd.s,qd.r.wk,qd.wt,prec,maxiter,alpha)
    theta <- theta - z$lambda
    cd <- c(z$c*10^z$lambda,z$d)
    int <- z$int
    ## theta search
    counter <- 0
    repeat {
        fit <- nlm(cv,theta,stepmax=1,ndigit=7)
        if (fit$code<=3)  break
        theta <- fit$est        
        counter <- counter + 1
        if (counter>=5) {
            warning("gss warning in ssden: CV iteration fails to converge")
            break
        }
    }
    ## return
    jk1 <- cv(fit$est)
    c <- cd[1:nxi]
    if (nnull) d <- cd[nxi+(1:nnull)]
    else d <- NULL
    list(lambda=0,theta=fit$est,c=c,d=d,int=int,cv=fit$min)
}

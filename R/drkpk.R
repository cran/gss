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
                        as.double(10^(lambda+theta)*q), as.integer(nxi),
                        as.double(rbind(10^theta*r,s)), as.integer(nobs),
                        as.integer(sum(cnt)), as.integer(cnt),
                        as.double(t(rbind(10^theta*qd.r,qd.s))), as.integer(nqd),
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
        alpha.wk <- max(0,log.la0-lambda-5)*(3-alpha) + alpha
        alpha.wk <- min(alpha.wk,3)
        adj <- ifelse (alpha.wk>alpha,(alpha.wk-alpha)*fit$wk[2],0)
        cv+adj
    }
    ## initialization
    mu.r <- apply(qd.wt*t(qd.r),2,sum)/sum(qd.wt)
    v.r <- apply(qd.wt*t(qd.r^2),2,sum)/sum(qd.wt)
    mu.s <- apply(qd.wt*t(qd.s),2,sum)/sum(qd.wt)
    v.s <- apply(qd.wt*t(qd.s^2),2,sum)/sum(qd.wt)
    if (is.null(s)) theta <- 0
    else theta <- log10(sum(v.s-mu.s^2)/nnull/sum(v.r-mu.r^2)*nxi) / 2
    log.la0 <- log10(sum(v.r-mu.r^2)/sum(diag(q))) + theta
    ## lambda search
    cd <- rep(0,nxi+nnull)
    int <- NULL
    la <- log.la0
    repeat {
        mn <- la-1
        mx <- la+1
        zz <- nlm0(cv,c(mn,mx))
        if (min(zz$est-mn,mx-zz$est)>=1e-3) break
        else la <- zz$est
    }
    ## return
    jk1 <- cv(zz$est)
    c <- cd[1:nxi]
    if (nnull) d <- cd[nxi+(1:nnull)]
    else d <- NULL
    list(lambda=zz$est,theta=theta,c=c,d=d,int=int,cv=zz$min)
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
                        as.double(10^lambda*q.wk), as.integer(nxi),
                        as.double(rbind(r.wk,s)), as.integer(nobs),
                        as.integer(sum(cnt)), as.integer(cnt),
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
        alpha.wk <- max(0,theta-log.th0-5)*(3-alpha) + alpha
        alpha.wk <- min(alpha.wk,3)
        adj <- ifelse (alpha.wk>alpha,(alpha.wk-alpha)*fit$wk[2],0)
        cv+adj
    }
    cv.wk <- function(theta) cv.scale*cv(theta)+cv.shift
    ## initialization
    theta <- -log10(apply(q,3,function(x)sum(diag(x))))
    r.wk <- q.wk <- qd.r.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
        q.wk <- q.wk + 10^theta[i]*q[,,i]
        qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[,,i]
    }
    ## theta adjustment
    z <- sspdsty(s,r.wk,q.wk,cnt,qd.s,qd.r.wk,qd.wt,prec,maxiter,alpha)
    theta <- theta + z$theta
    for (i in 1:nq) {
        theta[i] <- 2*theta[i] + log10(t(z$c)%*%q[,,i]%*%z$c)
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
    lambda <- z$lambda
    log.th0 <- log.th0 + z$lambda
    theta <- theta + z$theta
    cd <- c(z$c,z$d)
    int <- z$int
    ## theta search
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
            warning("gss warning in ssden: CV iteration fails to converge")
            break
        }
    }
    ## return
    jk1 <- cv(zz$est)
    c <- cd[1:nxi]
    if (nnull) d <- cd[nxi+(1:nnull)]
    else d <- NULL
    list(lambda=lambda,theta=zz$est,c=c,d=d,int=int,cv=zz$min)
}

## Calculate Kullback-Leibler projection from ssden objects
project.ssden <- function(object,include,mesh=FALSE,...)
{
    qd.pt <- object$quad$pt
    qd.wt <- object$quad$wt
    ## evaluate full model
    mesh0 <- dssden(object,qd.pt)
    ## extract terms in subspace
    nobs <- dim(object$mf)[1]
    nqd <- length(qd.wt)
    nxi <- length(object$id.basis)
    qd.s <- qd.r <- q <- NULL
    theta <- d <- NULL
    n0.wk <- nq.wk <- nq <- 0
    for (label in object$terms$labels) {
        x.basis <- object$mf[object$id.basis,object$term[[label]]$vlist]
        qd.x <- qd.pt[,object$term[[label]]$vlist]
        nphi <- object$term[[label]]$nphi
        nrk <- object$term[[label]]$nrk
        if (nphi) {
            phi <- object$term[[label]]$phi
            for (i in 1:nphi) {
                n0.wk <- n0.wk + 1
                if (!any(label==include)) next
                d <- c(d,object$d[n0.wk])
                qd.s <- cbind(qd.s,phi$fun(qd.x,nu=i,env=phi$env))
            }
        }
        if (nrk) {
            rk <- object$term[[label]]$rk
            for (i in 1:nrk) {
                nq.wk <- nq.wk + 1
                if (!any(label==include)) next
                nq <- nq + 1
                theta <- c(theta,object$theta[nq.wk])
                qd.r <- array(c(qd.r,rk$fun(x.basis,qd.x,nu=i,env=rk$env,out=TRUE)),
                              c(nxi,nqd,nq))
                q <- cbind(q,rk$fun(x.basis,x.basis,nu=i,env=rk$env,out=FALSE))
            }
        }
    }
    if (!is.null(qd.s)) {
        nn <- nxi + ncol(qd.s)
        qd.s <- t(qd.s)
    }
    else nn <- nxi
    ## calculate projection
    rkl <- function(theta1=NULL) {
        theta.wk <- 1:nq
        theta.wk[fix] <- theta[fix]
        if (nq-1) theta.wk[-fix] <- theta1
        qd.rs <- 0
        for (i in 1:nq) qd.rs <- qd.rs + 10^theta.wk[i]*qd.r[,,i]
        qd.rs <- rbind(qd.rs,qd.s)
        z <- .Fortran("drkl",
                      cd=as.double(cd), as.integer(nn),
                      as.double(t(qd.rs)), as.integer(nqd), as.double(qd.wt),
                      mesh=as.double(mesh0), as.double(.Machine$double.eps),
                      double(nqd), double(nqd), double(nn), double(nn*nn),
                      integer(nn), double(nn), double(nn), double(nqd),
                      as.double(1e-6), as.integer(30),
                      info=integer(1), PACKAGE="gss")
        if (z$info==1)
            stop("gss error in project.ssden: Newton iteration diverges")
        if (z$info==2)
            warning("gss warning in project.ssden: Newton iteration fails to converge")
        assign("cd",z$cd,inherit=TRUE)
        assign("mesh1",z$mesh,inherit=TRUE)
        sum(qd.wt*log(mesh0/mesh1)*mesh0)
    }
    ## initialization
    qd.r.wk <- 0
    for (i in 1:nq) qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[,,i]
    mu.r <- apply(qd.wt*t(qd.r.wk),2,sum)/sum(qd.wt)
    v.r <- apply(qd.wt*t(qd.r.wk^2),2,sum)/sum(qd.wt)
    mu.s <- apply(qd.wt*t(qd.s),2,sum)/sum(qd.wt)
    v.s <- apply(qd.wt*t(qd.s^2),2,sum)/sum(qd.wt)
    if (is.null(qd.s)) theta.wk <- 0
    else theta.wk <- log10(sum(v.s-mu.s^2)/(nn-nxi)/sum(v.r-mu.r^2)*nxi) / 2
    theta <- theta + theta.wk
    tmp <- NULL
    for (i in 1:nq) tmp <- c(tmp,10^theta[i]*sum(q[,i]))
    fix <- rev(order(tmp))[1]
    ## projection
    cd <- c(10^(-theta.wk)*object$c,d)
    mesh1 <- NULL
    if (nq-1) {
        if (nq-2) {
            zz <- nlm(rkl,theta[-fix],stepmax=.5,ndigit=7)
            if (zz$code>3)
                warning("gss warning in project.ssden: theta iteration fails to converge")
        }
        else {
            the.wk <- theta[-fix]
            repeat {
                mn <- the.wk-1
                mx <- the.wk+1
                zz <- nlm0(rkl,c(mn,mx))
                if (min(zz$est-mn,mx-zz$est)>=1e-3) break
                else the.wk <- zz$est
            }
        }
        kl <- rkl(zz$est)
    }
    else kl <- rkl()
    kl0 <- sum(qd.wt*log(mesh0)*mesh0) + log(sum(qd.wt))
    obj <- list(ratio=kl/kl0,kl=kl)
    if (mesh) obj$mesh <- mesh1
    obj
}

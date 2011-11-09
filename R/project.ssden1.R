## Calculate square error projection from ssden1 objects
project.ssden1 <- function(object,include,...)
{
    ## calculate log(rho) and integrals
    rho1 <- sum(object$rho.int)
    rho2 <- rho1^2-sum(object$rho.int^2)+sum(object$rho.int2)
    ## calculate cross integrals of rho, phi, and rk
    s <- object$int$s
    r <- object$int$r
    s.rho <- object$int$s.rho - s*rho1
    r.rho <- object$int$r.rho - r*rho1
    ss <- object$int$ss
    sr <- object$int$sr
    rr <- object$int$rr
    ## evaluate full model
    d <- object$d
    c <- object$c
    theta <- object$theta
    nq <- length(theta)
    s.eta <- ss%*%d
    r.eta <- tmp <- NULL
    r.wk <- r.rho.wk <- sr.wk <- rr.wk <- 0
    for (i in 1:nq) {
        tmp <- c(tmp,10^(2*theta[i])*sum(diag(rr[,,i,i])))
        s.eta <- s.eta + 10^theta[i]*sr[,,i]%*%c
        r.eta.wk <- t(sr[,,i])%*%d
        r.wk <- r.wk + 10^theta[i]*r[,i]
        r.rho.wk <- r.rho.wk + 10^theta[i]*r.rho[,i]
        sr.wk <- sr.wk + 10^theta[i]*sr[,,i]
        for (j in 1:nq) {
            r.eta.wk <- r.eta.wk + 10^theta[j]*rr[,,i,j]%*%c
            rr.wk <- rr.wk + 10^(theta[i]+theta[j])*rr[,,i,j]
        }
        r.eta <- cbind(r.eta,r.eta.wk)
    }
    s.eta <- s.eta - s*(sum(s*d)+sum(r.wk*c))
    r.eta <- r.eta - r*(sum(s*d)+sum(r.wk*c))
    ss <- ss - outer(s,s,"*")
    sr.wk <- sr.wk - outer(s,r.wk,"*")
    rr.wk <- rr.wk - outer(r.wk,r.wk,"*")
    rho.eta <- sum(s.rho*d) + sum(r.rho.wk*c)
    eta2 <- sum(c*(rr.wk%*%c)) + sum(d*(ss%*%d)) + 2*sum(d*(sr.wk%*%c))
    mse <- eta2 + rho2-rho1^2 + 2*rho.eta
    ## extract terms in subspace
    include <- union(names(object$mf),include)
    id.s <- id.q <- NULL
    for (label in include) {
        if (!any(label==object$terms$labels)) next
        term <- object$terms[[label]]
        if (term$nphi>0) id.s <- c(id.s,term$iphi+(1:term$nphi)-2)
        if (term$nrk>0) id.q <- c(id.q,term$irk+(1:term$nrk)-1)
    }
    ## calculate projection
    rkl <- function(theta1=NULL) {
        theta.wk <- 1:nq
        theta.wk[fix] <- theta[fix]
        theta.wk[id.q0] <- theta1
        ##
        ss.wk <- ss[id.s,id.s]
        r.eta.wk <- r.wk <- sr.wk <- rr.wk <- 0
        for (i in id.q) {
            r.eta.wk <- r.eta.wk + 10^theta.wk[i]*r.eta[,i]
            r.wk <- r.wk + 10^theta.wk[i]*r[,i]
            sr.wk <- sr.wk + 10^theta.wk[i]*sr[id.s,,i]
            for (j in id.q) {
                rr.wk <- rr.wk + 10^(theta.wk[i]+theta.wk[j])*rr[,,i,j]
            }
        }
        sr.wk <- sr.wk - outer(s[id.s],r.wk,"*")
        rr.wk <- rr.wk - outer(r.wk,r.wk,"*")
        v <- cbind(rbind(ss.wk,t(sr.wk)),rbind(sr.wk,rr.wk))
        mu <- c(s.eta[id.s],r.eta.wk)
        nn <- length(mu)
        z <- .Fortran("dchdc",
                      v=as.double(v), as.integer(nn), as.integer(nn),
                      double(nn),
                      jpvt=as.integer(rep(0,nn)),
                      as.integer(1),
                      rkv=integer(1),
                      PACKAGE="base")
        m.eps <- .Machine$double.eps
        v <- matrix(z$v,nn,nn)
        rkv <- z$rkv
        while (v[rkv,rkv]<2*sqrt(m.eps)*v[1,1]) rkv <- rkv - 1
        if (rkv<nn) v[(1:nn)>rkv,(1:nn)>rkv] <- diag(v[1,1],nn-rkv)
        mu <- backsolve(v,mu[z$jpvt],tran=TRUE)
        eta2 - sum(mu[1:rkv]^2)
    }
    cv.wk <- function(theta) cv.scale*rkl(theta)+cv.shift
    ## initialization
    tmp[-id.q] <- 0
    fix <- rev(order(tmp))[1]
    id.q0 <- id.q[id.q!=fix]
    ## projection
    skip.iter <- TRUE
    if (skip.iter) se <- rkl(theta[id.q0])
    else {
        nq0 <- length(id.q)
        if (nq0-2) {
            ## scale and shift cv
            tmp <- abs(rkl(theta[id.q0]))
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
            zz <- nlm(cv.wk,theta[id.q0],stepmax=.5,ndigit=7)
        }
        else {
            the.wk <- theta[id.q0]
            repeat {
                mn <- the.wk-1
                mx <- the.wk+1
                zz <- nlm0(rkl,c(mn,mx))
                if (min(zz$est-mn,mx-zz$est)>=1e-3) break
                else the.wk <- zz$est
            }
        }
        se <- rkl(zz$est)
    }
    list(ratio=se/mse,se=se)
}

## Calculate Kullback-Leibler projection from sshzd objects
project.sshzd <- function(object,include,mesh=FALSE,...)
{
    if (!(object$tname%in%include))
        stop("gss error in project.sshzd: time main effect missing in included terms")
    quad.pt <- object$quad$pt
    quad.wt <- object$quad$wt
    nx <- dim(object$qd.wt)[2]
    nbasis <- length(object$id.basis)
    mesh0 <- object$mesh0
    ## extract terms in subspace
    nqd <- length(quad.pt)
    nxi <- length(object$id.basis)
    d <- qd.s <- q <- theta <- NULL
    qd.r <- as.list(NULL)
    n0.wk <- nu <- nq.wk <- nq <- 0
    for (label in object$terms$labels) {
        vlist <- object$terms[[label]]$vlist
        x.list <- object$xnames[object$xnames%in%vlist]
        xy.basis <- object$mf[object$id.basis,vlist]
        qd.xy <- data.frame(matrix(0,nqd,length(vlist)))
        names(qd.xy) <- vlist
        if (object$tname%in%vlist) qd.xy[,object$tname] <- quad.pt
        if (length(x.list)) xx <- object$x.pt[,x.list,drop=FALSE]
        else xx <- NULL
        nphi <- object$terms[[label]]$nphi
        nrk <- object$terms[[label]]$nrk
        if (nphi) {
            phi <- object$terms[[label]]$phi
            for (i in 1:nphi) {
                n0.wk <- n0.wk + 1
                if (label=="1") {
                    d <- object$d[n0.wk]
                    nu <- nu + 1
                    qd.wk <- matrix(1,nqd,nx)
                    qd.s <- array(c(qd.s,qd.wk),c(nqd,nx,nu))
                    next
                }
                if (!any(label==include)) next
                d <- c(d,object$d[n0.wk])
                nu <- nu + 1
                if (is.null(xx))
                    qd.wk <- matrix(phi$fun(qd.xy[,,drop=TRUE],nu=i,env=phi$env),nqd,nx)
                else {
                    qd.wk <- NULL
                    for (j in 1:nx) {
                        qd.xy[,x.list] <- xx[rep(j,nqd),]
                        for (k in x.list)
                            if (is.factor(xx[,k])) qd.xy[,k] <- as.factor(qd.xy[,k])
                        qd.wk <- cbind(qd.wk,phi$fun(qd.xy[,,drop=TRUE],i,phi$env))
                    }
                }
                qd.s <- array(c(qd.s,qd.wk),c(nqd,nx,nu))
            }
        }
        if (nrk) {
            rk <- object$terms[[label]]$rk
            for (i in 1:nrk) {
                nq.wk <- nq.wk + 1
                if (!any(label==include)) next
                nq <- nq + 1
                theta <- c(theta,object$theta[nq.wk])
                q <- cbind(q,rk$fun(xy.basis,xy.basis,i,rk$env,out=FALSE))
                if (is.null(xx))
                    qd.r[[nq]] <- rk$fun(qd.xy[,,drop=TRUE],xy.basis,i,rk$env,out=TRUE)
                else {
                    qd.wk <- NULL
                    for (j in 1:nx) {
                        qd.xy[,x.list] <- xx[rep(j,nqd),]
                        for (k in x.list)
                            if (is.factor(xx[,k])) qd.xy[,k] <- as.factor(qd.xy[,k])
                        qd.wk <- array(c(qd.wk,rk$fun(qd.xy[,,drop=TRUE],xy.basis,i,rk$env,TRUE)),
                                       c(nqd,nbasis,j))
                    }
                    qd.r[[nq]] <- qd.wk
                }
            }
        }
    }
    if (!is.null(qd.s)) nnull <- dim(qd.s)[3]
    else nnull <- 0
    nn <- nxi + nnull
    ## calculate projection
    rkl <- function(theta1=NULL) {
        theta.wk <- 1:nq
        theta.wk[fix] <- theta[fix]
        if (nq-1) theta.wk[-fix] <- theta1
        qd.r.wk <- array(0,c(nqd,nxi,nx))
        for (i in 1:nq) {
            if (length(dim(qd.r[[i]]))==3) qd.r.wk <- qd.r.wk + 10^theta.wk[i]*qd.r[[i]]
            else qd.r.wk <- qd.r.wk + as.vector(10^theta.wk[i]*qd.r[[i]])
        }
        qd.r.wk <- aperm(qd.r.wk,c(1,3,2))
        qd.r.wk <- array(c(qd.r.wk,qd.s),c(nqd,nx,nn))
        qd.r.wk <- aperm(qd.r.wk,c(1,3,2))
        z <- .Fortran("hrkl",
                      cd=as.double(cd), as.integer(nn),
                      as.double(qd.r.wk), as.integer(nqd), as.integer(nx),
                      as.double(object$qd.wt), mesh=as.double(object$qd.wt*mesh0),
                      as.double(.Machine$double.eps), double(nqd*nx),
                      double(nn), double(nn), double(nn*nn), integer(nn), double(nn),
                      double(nn), double(nqd*nx), as.double(1e-6), as.integer(30),
                      info=integer(1), PACKAGE="gss")
        if (z$info==1)
            stop("gss error in project.sshzd: Newton iteration diverges")
        if (z$info==2)
            warning("gss warning in project.sshzd: Newton iteration fails to converge")
        assign("cd",z$cd,inherit=TRUE)
        assign("mesh1",z$mesh,inherit=TRUE)
        sum(object$qd.wt*(log(mesh0/mesh1)*mesh0-mesh0+mesh1))
    }
    ## initialization
    qd.r.wk <- array(0,c(nqd,nxi,nx))
    for (i in 1:nq) {
        if (length(dim(qd.r[[i]]))==3) qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[[i]]
        else qd.r.wk <- qd.r.wk + as.vector(10^theta[i]*qd.r[[i]])
    }
    v.s <- v.r <- 0
    for (i in 1:nx) {
        if (nnull) v.s <- v.s + apply(object$qd.wt[,i]*qd.s[,i,,drop=FALSE]^2,2,sum)
        v.r <- v.r + apply(object$qd.wt[,i]*qd.r.wk[,,i,drop=FALSE]^2,2,sum)
    }
    if (nnull) theta.wk <- log10(sum(v.s)/nnull/sum(v.r)*nxi) / 2
    else theta.wk <- 0
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
                warning("gss warning in project.sshzd: theta iteration fails to converge")
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
    kl0 <- sum(object$qd.wt*(log(mesh0/object$cfit)*mesh0-mesh0+object$cfit))
    obj <- list(ratio=kl/kl0,kl=kl)
    if (mesh) obj$mesh <- mesh1
    obj
}

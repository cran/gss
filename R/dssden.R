dssden <- ## Evaluate density estimate
function (object,x) {
    if (class(object)!="ssden") stop("gss error in dssden: not a ssden object")
    if (dim(object$mf)[2]==1&is.vector(x)) {
        x <- data.frame(x)
        colnames(x) <- colnames(object$mf)
    }
    s <- NULL
    r <- matrix(0,dim(x)[1],length(object$id.basis))
    nq <- 0
    for (label in object$terms$labels) {
        xx <- object$mf[object$id.basis,object$terms[[label]]$vlist]
        x.new <- x[,object$terms[[label]]$vlist]
        nphi <- object$terms[[label]]$nphi
        nrk <-  object$terms[[label]]$nrk
        if (nphi) {
            phi <-  object$terms[[label]]$phi
            for (i in 1:nphi) {
                s <- cbind(s,phi$fun(x.new,nu=i,env=phi$env))
            }
        }
        if (nrk) {
            rk <- object$terms[[label]]$rk
            for (i in 1:nrk) {
                nq <- nq + 1
                r <- r + 10^object$theta[nq]*rk$fun(x.new,xx,nu=i,env=rk$env,out=TRUE)
            }
        }
    }
    as.vector(exp(cbind(s,r)%*%c(object$d,object$c))/object$int)
}

pssden <- ## Compute cdf for univariate density estimate
function(object,q) {
    if (class(object)!="ssden") stop("gss error in pssden: not a ssden object")
    if (dim(object$mf)[2]!=1) stop("gss error in pssden: not a 1-D density")
    mn <- min(object$domain)
    mx <- max(object$domain)
    order.q <- rank(q)
    p <- q <- sort(q)
    q.dup <- duplicated(q)
    p[q<=mn] <- 0
    p[q>=mx] <- 1
    kk <- (1:length(q))[q>mn&q<mx]
    for (i in kk) {
        if (q.dup[i]) {
            p[i] <- p.dup
            next
        }
        nqd.l <- max(20,ceiling((q[i]-mn)/(mx-mn)*200))
        qd.l <- gauss.quad(nqd.l,c(mn,q[i]))
        p.l <- sum(dssden(object,qd.l$pt)*qd.l$wt)
        nqd.u <- max(20,ceiling((mx-q[i])/(mx-mn)*200))
        qd.u <- gauss.quad(nqd.u,c(q[i],mx))
        p.u <- sum(dssden(object,qd.u$pt)*qd.u$wt)
        p[i] <- p.dup <- p.l/(p.l+p.u)
    }
    p[order.q]
}

qssden <- ## Compute quantiles for univariate density estimate
function(object,p) {
    if (class(object)!="ssden") stop("gss error in qssden: not a ssden object")
    if (dim(object$mf)[2]!=1) stop("gss error in qssden: not a 1-D density")
    mn <- min(object$domain)
    mx <- max(object$domain)
    order.p <- rank(p)
    q <- p <- sort(p)
    p.dup <- duplicated(p)
    q[p<=0] <- mn
    q[p>=1] <- mx
    kk <- (1:length(p))[p>0&p<1]
    q.wk <- object$quad$pt[,1]
    p.wk <- cumsum(object$quad$wt*dssden(object,q.wk))
    for (i in kk) {
        if (p.dup[i]) {
            q[i] <- q.dup
            next
        }
        j <- which.min(abs(p[i]-p.wk))
        q0 <- q.wk[j]
        p0 <- pssden(object,q0)
        if (p0==p[i]) {
            q[i] <- q0
            next
        }
        if (p0<p[i]) {
            q.l <- q0
            p.l <- p0
            while (p0<p[i]) {
                j <- j+1
                q0 <- ifelse(is.null(q.wk[j]),mx,q.wk[j])
                p0 <- pssden(object,q0)
            }
            q.u <- q0
            p.u <- p0
        }
        else {
            q.u <- q0
            p.u <- p0
            while (p0>p[i]) {
                j <- j-1
                q0 <- ifelse(is.null(q.wk[j]),mn,q.wk[j])
                p0 <- pssden(object,q0)
            }
            q.l <- q0
            p.l <- p0
        }
        while (abs(p0-p[i])>1e-10) {
            q0 <- q.l+(p[i]-p.l)/(p.u-p.l)*(q.u-q.l)
            p0 <- pssden(object,q0)
            if (p0>p[i]) {
                q.u <- q0
                p.u <- p0
            }
            else {
                q.l <- q0
                p.l <- p0
            }
        }
        q[i] <- q.dup <- q0
    }
    q[order.p]
}

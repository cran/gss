cdssden <- ## Evaluate conditional density
function (object,x,cond,int=NULL) {
    if (class(object)!="ssden")
        stop("gss error in cdssden: not a ssden object")
    if (nrow(cond)!=1)
        stop("gss error in cdssden: condition has to be a single point")
    xnames <- NULL
    for (i in colnames(object$mf))
        if (all(i!=colnames(cond))) xnames <- c(xnames,i)
    if (any(length(xnames)==c(0,ncol(object$mf))))
        stop("gss error in cdssden: not a conditional density")
    if (length(xnames)==1&is.vector(x)) {
        x <- data.frame(x)
        colnames(x) <- xnames
    }
    if (!all(sort(xnames)==sort(colnames(x))))
        stop("gss error in cdssden: mismatched variable names")
    ## Calculate normalizing constant
    while (is.null(int)) {
        fac.list <- NULL
        num.list <- NULL
        for (xlab in xnames) {
            if (is.factor(x.wk <- x[[xlab]])) fac.list <- c(fac.list,xlab)
            else {
                if (!is.vector(x.wk)|is.null(object$domain[[xlab]])) {
                    warning("gss warning in cdssden: int set to 1")
                    int <- 1
                    next
                }
                else num.list <- c(num.list,xlab)
            }
        }
        ## Generate quadrature for numerical variables
        if (!is.null(num.list)) {
            if (length(num.list)==1) {
                ## Gauss-Legendre quadrature
                mn <- min(object$domain[,num.list])
                mx <- max(object$domain[,num.list])
                quad <- gauss.quad(200,c(mn,mx))
                quad$pt <- data.frame(quad$pt)
                colnames(quad$pt) <- num.list
            }
            else {
                ## Smolyak cubature
                domain.wk <- object$domain[,num.list]
                code <- c(15,14,13)
                quad <- smolyak.quad(ncol(domain.wk),code[ncol(domain.wk)-1])
                for (i in 1:ncol(domain.wk)) {
                    xlab <- colnames(domain.wk)[i]
                    wk <- object$mf[[xlab]]
                    jk <- ssden(~wk,domain=data.frame(wk=domain.wk[,i]),alpha=2,
                                id.basis=object$id.basis)
                    quad$pt[,i] <- qssden(jk,quad$pt[,i])
                    quad$wt <- quad$wt/dssden(jk,quad$pt[,i])
                }
                jk <- wk <- NULL
                quad$pt <- data.frame(quad$pt)
                colnames(quad.pt) <- colnames(domain.wk)
            }
        }
        else quad <- list(pt=data.frame(dum=1),wt=1)
        ## Incorporate factors in quadrature
        if (!is.null(fac.list)) {
            for (i in 1:length(fac.list)) {
                wk <- expand.grid(levels(object$mf[[fac.list[i]]]),1:length(quad$wt))
                quad$wt <- quad$wt[wk[,2]]
                col.names <- c(fac.list[i],colnames(quad$pt))
                quad$pt <- data.frame(wk[,1],quad$pt[wk[,2],])
                colnames(quad$pt) <- col.names
            }
        }
        xmesh <- quad$pt[,xnames,drop=FALSE]
        xx <- cond[rep(1,nrow(xmesh)),,drop=FALSE]
        int <- sum(dssden(object,cbind(xmesh,xx))*quad$wt)
    }
    ## Return value
    xx <- cond[rep(1,nrow(x)),,drop=FALSE]
    list(pdf=dssden(object,cbind(x,xx))/int,int=int)
}

cpssden <- ## Compute cdf for univariate conditional density
function(object,q,cond,int=NULL) {
    if (class(object)!="ssden")
        stop("gss error in cpssden: not a ssden object")
    xnames <- NULL
    for (i in colnames(object$mf))
        if (all(i!=colnames(cond))) xnames <- c(xnames,i)
    if ((length(xnames)!=1)|!is.vector(object$mf[,xnames]))
        stop("gss error in cpssden: not a 1-D conditional density")
    mn <- min(object$domain[,xnames])
    mx <- max(object$domain[,xnames])
    if (is.null(int)) int <- cdssden(object,mn,cond)$int
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
        p.l <- sum(cdssden(object,qd.l$pt,cond,int)$pdf*qd.l$wt)
        nqd.u <- max(20,ceiling((mx-q[i])/(mx-mn)*200))
        qd.u <- gauss.quad(nqd.u,c(q[i],mx))
        p.u <- sum(cdssden(object,qd.u$pt,cond,int)$pdf*qd.u$wt)
        p[i] <- p.dup <- p.l/(p.l+p.u)
    }
    p[order.q]
}

cqssden <- ## Compute quantiles for univariate conditional density
function(object,p,cond,int=NULL) {
    if (class(object)!="ssden")
        stop("gss error in cqssden: not a ssden object")
    xnames <- NULL
    for (i in colnames(object$mf))
        if (all(i!=colnames(cond))) xnames <- c(xnames,i)
    if ((length(xnames)!=1)|!is.vector(object$mf[,xnames]))
        stop("gss error in cqssden: not a 1-D conditional density")
    mn <- min(object$domain[,xnames])
    mx <- max(object$domain[,xnames])
    if (is.null(int)) int <- cdssden(object,mn,cond)$int
    order.p <- rank(p)
    q <- p <- sort(p)
    p.dup <- duplicated(p)
    q[p<=0] <- mn
    q[p>=1] <- mx
    kk <- (1:length(p))[p>0&p<1]
    quad <- gauss.quad(200,c(mn,mx))
    q.wk <- quad$pt
    p.wk <- cumsum(quad$wt*cdssden(object,q.wk,cond,int)$pdf)
    for (i in kk) {
        if (p.dup[i]) {
            q[i] <- q.dup
            next
        }
        j <- which.min(abs(p[i]-p.wk))
        q0 <- q.wk[j]
        p0 <- cpssden(object,q0,cond,int)
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
                p0 <- cpssden(object,q0,cond,int)
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
                p0 <- cpssden(object,q0,cond,int)
            }
            q.l <- q0
            p.l <- p0
        }
        while (abs(p0-p[i])>1e-10) {
            q0 <- q.l+(p[i]-p.l)/(p.u-p.l)*(q.u-q.l)
            p0 <- cpssden(object,q0,cond,int)
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

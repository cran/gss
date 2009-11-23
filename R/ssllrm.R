## Fit log-linear regression model
ssllrm <- function(formula,response,type=NULL,data=list(),weights,
                   subset,na.action=na.omit,partial=NULL,alpha=1.4,
                   id.basis=NULL,nbasis=NULL,seed=NULL,
                   prec=1e-7,maxiter=30,skip.iter=FALSE)
{
    ## Obtain model frame and model terms
    mf <- match.call()
    mf$response <- mf$type <- mf$partial <- NULL
    mf$alpha <- mf$id.basis <- mf$nbasis <- mf$seed <- NULL
    mf$prec <- mf$maxiter <- mf$skip.iter <- NULL
    term.wk <- terms.formula(mf$formula)
    ynames <- as.character(attr(terms(response),"variables"))[-1]
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf,parent.frame())
    cnt <- model.weights(mf)
    mf$"(weights)" <- NULL
    ## Generate sub-basis
    nobs <- nrow(mf)
    if (is.null(id.basis)) {
        if (is.null(nbasis))  nbasis <- max(30,ceiling(10*nobs^(2/9)))
        if (nbasis>=nobs)  nbasis <- nobs
        if (!is.null(seed))  set.seed(seed)
        id.basis <- sample(nobs,nbasis,prob=cnt)
    }
    else {
        if (max(id.basis)>nobs|min(id.basis)<1)
            stop("gss error in ssllrm: id.basis out of range")
        nbasis <- length(id.basis)
    }
    ## Check inputs
    mt <- attr(mf,"terms")
    vars <- as.character(attr(mt,"variables"))[-1]
    if(!all(ynames%in%vars)) stop("gss error in ssllrm: response missing in model")
    for (ylab in ynames) {
        if (!is.factor(mf[,ylab])) stop("gss error in ssllrm: response not a factor")
    }
    xnames <- vars[!(vars%in%ynames)]
    if (is.null(xnames)) stop("gss error in ssllrm: missing covariate")
    ## Generate terms    
    term <- mkterm(mf,type)
    term.labels <- labels(mt)
    facs <- attr(mt,"factors")
    ind.wk <- NULL
    for (lab in term.labels)
        ind.wk <- c(ind.wk,any(facs[ynames,lab]))
    term$labels <- term.labels[ind.wk]
    ## Generate quadrature
    qd.pt <- data.frame(levels(mf[,ynames[1]]))
    if (length(ynames)>1) {
        for (ylab in ynames[-1]) {
            wk <- expand.grid(levels(mf[,ylab]),1:dim(qd.pt)[1])
            qd.pt <- data.frame(qd.pt[wk[,2],],wk[,1])
        }
    }
    colnames(qd.pt) <- ynames
    x <- mf[,xnames,drop=FALSE]
    ## Partial
    if (!is.null(partial)) {
        if (is.vector(partial)) partial <- as.matrix(partial)
        partial <- scale(partial)
        if (dim(partial)[1]!=dim(mf)[1])
            stop("gss error in ssllrm: partial data are of wrong size")
        yterms <- labels(terms(response))
        mf$partial <- partial
    }
    else yterms <- NULL
    ## obtain unique covariate observations
    xx <- mf[,xnames,drop=FALSE]
    if (!is.null(partial)) xx <- cbind(xx,partial)
    xx <- as.matrix(xx)
    x.pt <- unique(xx)
    nx <- dim(x.pt)[1]
    x.dup.ind <- duplicated(xx)
    x.dup <- xx[x.dup.ind,,drop=FALSE]
    ## xx[i,]==x.pt[id.x[i],]
    id.x <- 1:nobs
    id.x[!x.dup.ind] <- 1:nx
    if (nobs-nx) {
        id.x.wk <- 1:(nobs-nx)
        for (i in 1:nx) {
            for (j in 1:(nobs-nx)) {
                if (sum(duplicated(rbind(x.pt[i,,drop=FALSE],x.dup[j,,drop=FALSE]))))
                    id.x.wk[j] <- i
            }
        }
        id.x[x.dup.ind] <- id.x.wk
    }
    ## integration weights at x.pt[i,]
    qd.wt <- rep(0,nx)
    if (is.null(cnt)) {
        for (i in 1:nobs) qd.wt[id.x[i]] <- qd.wt[id.x[i]]+1
    }
    else {
        for (i in 1:nobs) qd.wt[id.x[i]] <- qd.wt[id.x[i]]+cnt[i]
    }
    qd.wt <- qd.wt/sum(qd.wt)
    ## Generate s, r, qd.s, and qd.r
    nmesh <- dim(qd.pt)[1]
    s <- r <- qd.s <- NULL
    qd.r <- as.list(NULL)
    nu <- nq <- 0
    for (label in term$labels) {
        part <- label%in%yterms  ## y only: no main effect, no interaction with partial
        vlist <- term[[label]]$vlist
        x.list <- xnames[xnames%in%vlist]
        y.list <- ynames[ynames%in%vlist]
        xy <- mf[,vlist]
        xy.basis <- mf[id.basis,vlist]
        qd.xy <- data.frame(matrix(0,nmesh,length(vlist)))
        names(qd.xy) <- vlist
        qd.xy[,y.list] <- qd.pt[,y.list]
        if (length(x.list)) xx <- x[!x.dup.ind,x.list,drop=FALSE]
        else xx <- NULL
        nphi <- term[[label]]$nphi
        nrk <- term[[label]]$nrk
        if (nphi) {
            phi <- term[[label]]$phi
            for (i in 1:nphi) {
                nu <- nu+1
                s.wk <- phi$fun(xy,nu=i,env=phi$env)
                s <- cbind(s,s.wk)
                if (is.null(xx)) {
                    qd.s.wk <- phi$fun(qd.xy[,,drop=TRUE],nu=i,env=phi$env)
                    qd.wk <- matrix(qd.s.wk,nmesh,nx)
                }
                else {
                    qd.wk <- NULL
                    for (j in 1:nx) {
                        qd.xy[,x.list] <- xx[rep(j,nmesh),]
                        qd.wk <- cbind(qd.wk,phi$fun(qd.xy,i,phi$env))
                    }
                }
                qd.s <- array(c(qd.s,qd.wk),c(nmesh,nx,nu))
                if (part) {
                    for (j in 1:dim(partial)[2]) {
                        nu <- nu+1
                        s <- cbind(s,s.wk*partial[,j])
                        qd.s <- array(c(qd.s,outer(qd.s.wk,partial[!x.dup.ind,j])),
                                      c(nmesh,nx,nu))
                    }
                }
            }
        }
        if (nrk) {
            rk <- term[[label]]$rk
            for (i in 1:nrk) {
                nq <- nq+1
                r.wk <- rk$fun(xy,xy.basis,nu=i,env=rk$env,out=TRUE)
                r <- array(c(r,r.wk),c(nobs,nbasis,nq))
                if (is.null(xx)) {
                    qd.r.wk <- rk$fun(qd.xy[,,drop=TRUE],xy.basis,nu=i,env=rk$env,out=TRUE)
                    qd.r[[nq]] <- qd.r.wk
                }
                else {
                    qd.wk <- NULL
                    for (j in 1:nx) {
                        qd.xy[,x.list] <- xx[rep(j,nmesh),]
                        qd.wk <- array(c(qd.wk,rk$fun(qd.xy,xy.basis,i,rk$env,TRUE)),
                                       c(nmesh,nbasis,j))
                    }
                    qd.r[[nq]] <- qd.wk                    
                }
                if (part) {
                    nq <- nq+1
                    r <- array(c(r,r.wk*(partial%*%t(partial[id.basis,,drop=FALSE]))),
                                 c(nobs,nbasis,nq))
                    qd.wk <- NULL
                    part.wk <- partial[!x.dup.ind,,drop=FALSE]
                    for (j in 1:nx) {
                        ww <- t(qd.r.wk)*as.vector(partial[id.basis,,drop=FALSE]%*%part.wk[j,])
                        qd.wk <- array(c(qd.wk,t(ww)),c(nmesh,nbasis,j))
                    }
                    qd.r[[nq]] <- qd.wk
                }
            }
        }
    }
    ## Check s rank
    if (!is.null(s)) {
        nnull <- dim(s)[2]
        if (qr(s)$rank<nnull)
            stop("gss error in ssllrm: fixed effect MLE is not unique")
    }
    ## Fit the model
    z <- mspllrm(s,r,id.basis,cnt,qd.s,qd.r,qd.wt,id.x,prec,maxiter,alpha,skip.iter)
    ## independent constant model
    cfit <- 1
    for (lab in ynames) {
        wk <- table(mf[,lab])
        wk <- wk/sum(wk)
        cfit <- as.vector(outer(wk,cfit))
    }
    ## Brief description of model terms
    desc <- NULL
    for (label in term$labels)
        desc <- rbind(desc,as.numeric(c(term[[label]][c("nphi","nrk")])))
    desc <- rbind(desc,apply(desc,2,sum))
    rownames(desc) <- c(term$labels,"total")
    colnames(desc) <- c("Unpenalized","Penalized")
    ## Return the results
    obj <- c(list(call=match.call(),mf=mf,cnt=cnt,terms=term,desc=desc,qd.pt=qd.pt,
                  qd.wt=qd.wt,x.dup.ind=x.dup.ind,alpha=alpha,ynames=ynames,
                  xnames=xnames,yterms=yterms,cfit=cfit,id.basis=id.basis,
                  skip.iter=skip.iter),z)
    if (is.null(cnt)) obj$se.aux$v <- sqrt(nobs)*obj$se.aux$v
    else obj$se.aux$v <- sqrt(sum(cnt))*obj$se.aux$v
    class(obj) <- c("ssllrm")
    obj
}

## Fit (multiple smoothing parameter) log-linear regression model
mspllrm <- function(s,r,id.basis,cnt,qd.s,qd.r,qd.wt,id.x,prec,maxiter,alpha,skip.iter)
{
    nobs <- dim(r)[1]
    nxi <- dim(r)[2]
    nqd <- dim(qd.r[[1]])[1]
    nx <- length(qd.wt)
    if (!is.null(s)) nnull <- dim(s)[2]
    else nnull <- 0
    nxis <- nxi+nnull
    if (is.null(cnt)) cnt <- 0
    ## cv functions
    cv.s <- function(lambda) {
        fit <- .Fortran("llrmnewton",
                        cd=as.double(cd), as.integer(nxis),
                        as.double(10^lambda*q.wk), as.integer(nxi),
                        as.double(t(cbind(r.wk,s))), as.integer(nobs),
                        as.integer(sum(cnt)), as.integer(cnt),
                        as.double(qd.r.wk), as.integer(nqd), as.integer(nx),
                        as.double(qd.wt), as.integer(id.x),
                        as.double(prec), as.integer(maxiter),
                        as.double(.Machine$double.eps),
                        wk=double(2*(nqd+1)*nx+2*nobs+nxis*(2*nxis+6)),
                        info=integer(1),PACKAGE="gss")
        if (fit$info==1) stop("gss error in ssllrm: Newton iteration diverges")
        if (fit$info==2) warning("gss warning in ssllrm: Newton iteration fails to converge")
        assign("eta",fit$wk[1:(nqd*nx)],inherit=TRUE)
        assign("cd",fit$cd,inherit=TRUE)
        cv <- alpha*fit$wk[nqd*nx+2]-fit$wk[nqd*nx+1]
        alpha.wk <- max(0,log.la0-lambda-5)*(3-alpha) + alpha
        alpha.wk <- min(alpha.wk,3)
        adj <- ifelse (alpha.wk>alpha,(alpha.wk-alpha)*fit$wk[nqd*nx+2],0)
        cv+adj
    }
    cv.m <- function(theta) {
        ind.wk <- theta!=theta.old
        if (sum(ind.wk)==nq) {
            r.wk0 <- 0
            qd.r.wk0 <- array(0,c(nqd,nxi,nx))
            for (i in 1:nq) {
                r.wk0 <- r.wk0 + 10^theta[i]*r[,,i]
                if (length(dim(qd.r[[i]]))==3) qd.r.wk0 <- qd.r.wk0 + 10^theta[i]*qd.r[[i]]
                else qd.r.wk0 <- qd.r.wk0 + as.vector(10^theta[i]*qd.r[[i]])
            }
            assign("r.wk",r.wk0+0,inherit=TRUE)
            assign("qd.r.wk",qd.r.wk0+0,inherit=TRUE)
            assign("theta.old",theta+0,inherit=TRUE)
        }
        else {
            r.wk0 <- r.wk
            qd.r.wk0 <- qd.r.wk
            for (i in (1:nq)[ind.wk]) {
                theta.wk <- (10^(theta[i]-theta.old[i])-1)*10^theta.old[i]
                r.wk0 <- r.wk0 + theta.wk*r[,,i]
                if (length(dim(qd.r[[i]]))==3) qd.r.wk0 <- qd.r.wk0 + theta.wk*qd.r[[i]]
                else qd.r.wk0 <- qd.r.wk0 + as.vector(theta.wk*qd.r[[i]])
            }
        }
        q.wk <- r.wk0[id.basis,]
        qd.r.wk0 <- aperm(qd.r.wk0,c(1,3,2))
        qd.r.wk0 <- array(c(qd.r.wk0,qd.s),c(nqd,nx,nxis))
        qd.r.wk0 <- aperm(qd.r.wk0,c(1,3,2))
        fit <- .Fortran("llrmnewton",
                        cd=as.double(cd), as.integer(nxis),
                        as.double(10^lambda*q.wk), as.integer(nxi),
                        as.double(t(cbind(r.wk0,s))), as.integer(nobs),
                        as.integer(sum(cnt)), as.integer(cnt),
                        as.double(qd.r.wk0), as.integer(nqd), as.integer(nx),
                        as.double(qd.wt), as.integer(id.x),
                        as.double(prec), as.integer(maxiter),
                        as.double(.Machine$double.eps),
                        wk=double(2*(nqd+1)*nx+2*nobs+nxis*(2*nxis+6)),
                        info=integer(1),PACKAGE="gss")
        if (fit$info==1) stop("gss error in ssllrm: Newton iteration diverges")
        if (fit$info==2) warning("gss warning in ssllrm: Newton iteration fails to converge")
        assign("eta",fit$wk[1:(nqd*nx)],inherit=TRUE)
        assign("cd",fit$cd,inherit=TRUE)
        cv <- alpha*fit$wk[nqd*nx+2]-fit$wk[nqd*nx+1]
        alpha.wk <- max(0,theta-log.th0-5)*(3-alpha) + alpha
        alpha.wk <- min(alpha.wk,3)
        adj <- ifelse (alpha.wk>alpha,(alpha.wk-alpha)*fit$wk[nqd*nx+2],0)
        cv+adj
    }
    cv.wk <- function(theta) cv.scale*cv.m(theta)+cv.shift
    ## Initialization
    theta <- -log10(apply(r[id.basis,,,drop=FALSE],3,function(x)sum(diag(x))))
    nq <- length(theta)
    qd.r.wk <- array(0,c(nqd,nxi,nx))
    for (i in 1:nq) {
        if (length(dim(qd.r[[i]]))==3) qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[[i]]
        else qd.r.wk <- qd.r.wk + as.vector(10^theta[i]*qd.r[[i]])
    }
    if (!nnull) {
        vv.r <- 0
        for (i in 1:nx) {
            mu.r <- apply(qd.r.wk[,,i,drop=FALSE],2,sum)/nqd
            v.r <- apply(qd.r.wk[,,i,drop=FALSE]^2,2,sum)/nqd
            v.r <- v.r - mu.r^2
            vv.r <- vv.r + qd.wt[i]*v.r
        }
        theta.wk <- 0
    }
    else {
        vv.s <- vv.r <- 0
        for (i in 1:nx) {
            mu.s <- apply(qd.s[,i,,drop=FALSE],2,sum)/nqd
            v.s <- apply(qd.s[,i,,drop=FALSE]^2,2,sum)/nqd
            v.s <- v.s - mu.s^2
            mu.r <- apply(qd.r.wk[,,i,drop=FALSE],2,sum)/nqd
            v.r <- apply(qd.r.wk[,,i,drop=FALSE]^2,2,sum)/nqd
            v.r <- v.r - mu.r^2
            vv.s <- vv.s + qd.wt[i]*v.s
            vv.r <- vv.r + qd.wt[i]*v.r
        }
        theta.wk <- log10(sum(vv.s)/nnull/sum(vv.r)*nxi) / 2
    }
    theta <- theta + theta.wk
    qd.r.wk <- aperm(10^theta.wk*qd.r.wk,c(1,3,2))
    qd.r.wk <- array(c(qd.r.wk,qd.s),c(nqd,nx,nxis))
    qd.r.wk <- aperm(qd.r.wk,c(1,3,2))
    r.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
    }
    q.wk <- r.wk[id.basis,]
    log.la0 <- log10(sum(vv.r)/sum(diag(q.wk))) + 2*theta.wk
    ## fixed theta iteration
    eta <- NULL
    cd <- rep(0,nxi+nnull)
    la <- log.la0
    repeat {
        mn <- la-1
        mx <- la+1
        zz <- nlm0(cv.s,c(mn,mx))
        if (min(zz$est-mn,mx-zz$est)>=1e-3) break
        else la <- zz$est
    }
    if (nq==1) {
        lambda <- zz$est
        se.aux <- .Fortran("llrmaux",
                           as.double(cd), as.integer(nxis),
                           as.double(10^lambda*q.wk), as.integer(nxi),
                           as.double(qd.r.wk), as.integer(nqd),
                           as.integer(nx), as.double(qd.wt),
                           as.double(.Machine$double.eps), double(nqd*nx),
                           double(nx), double(nxis),
                           v=double(nxis*nxis), double(nxis*nxis),
                           jpvt=integer(nxis), PACKAGE="gss")[c("v","jpvt")]
        c <- cd[1:nxi]
        if (nnull) d <- cd[nxi+(1:nnull)]
        else d <- NULL
        eta <- matrix(eta,nqd,nx)
        for (i in 1:nx) eta[,i] <- eta[,i]/sum(eta[,i])
        return(list(lambda=zz$est,theta=theta,c=c,d=d,cv=zz$min,
                    fit=t(eta),se.aux=se.aux))
    }
    ## theta adjustment
    qd.r.wk <- array(0,c(nqd,nxi,nx))
    for (i in 1:nq) {
        theta[i] <- 2*theta[i] + log10(t(cd[1:nxi])%*%r[id.basis,,i]%*%cd[1:nxi])
        if (length(dim(qd.r[[i]]))==3) qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[[i]]
        else qd.r.wk <- qd.r.wk + as.vector(10^theta[i]*qd.r[[i]])
    }
    if (!nnull) {
        vv.r <- 0
        for (i in 1:nx) {
            mu.r <- apply(qd.r.wk[,,i,drop=FALSE],2,sum)/nqd
            v.r <- apply(qd.r.wk[,,i,drop=FALSE]^2,2,sum)/nqd
            v.r <- v.r - mu.r^2
            vv.r <- vv.r + qd.wt[i]*v.r
        }
        theta.wk <- 0
    }
    else {
        vv.s <- vv.r <- 0
        for (i in 1:nx) {
            mu.s <- apply(qd.s[,i,,drop=FALSE],2,sum)/nqd
            v.s <- apply(qd.s[,i,,drop=FALSE]^2,2,sum)/nqd
            v.s <- v.s - mu.s^2
            mu.r <- apply(qd.r.wk[,,i,drop=FALSE],2,sum)/nqd
            v.r <- apply(qd.r.wk[,,i,drop=FALSE]^2,2,sum)/nqd
            v.r <- v.r - mu.r^2
            vv.s <- vv.s + qd.wt[i]*v.s
            vv.r <- vv.r + qd.wt[i]*v.r
        }
        theta.wk <- log10(sum(vv.s)/nnull/sum(vv.r)*nxi) / 2
    }
    theta <- theta + theta.wk
    qd.r.wk <- aperm(10^theta.wk*qd.r.wk,c(1,3,2))
    qd.r.wk <- array(c(qd.r.wk,qd.s),c(nqd,nx,nxis))
    qd.r.wk <- aperm(qd.r.wk,c(1,3,2))
    r.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
    }
    q.wk <- r.wk[id.basis,]
    log.la0 <- log10(sum(vv.r)/sum(diag(q.wk))) + 2*theta.wk
    log.th0 <- theta-log.la0
    ## fixed theta iteration
    cd <- rep(0,nxi+nnull)
    la <- log.la0
    repeat {
        mn <- la-1
        mx <- la+1
        zz <- nlm0(cv.s,c(mn,mx))
        if (min(zz$est-mn,mx-zz$est)>=1e-3) break
        else la <- zz$est
    }
    lambda <- zz$est
    ## early return
    if (skip.iter) {
        q.wk <- 0
        qd.r.wk <- array(0,c(nqd,nxi,nx))
        for (i in 1:nq) {
            q.wk <- q.wk + 10^theta[i]*r[id.basis,,i]
            if (length(dim(qd.r[[i]]))==3) qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[[i]]
            else qd.r.wk <- qd.r.wk + as.vector(10^theta[i]*qd.r[[i]])
        }
        qd.r.wk <- aperm(qd.r.wk,c(1,3,2))
        qd.r.wk <- array(c(qd.r.wk,qd.s),c(nqd,nx,nxis))
        qd.r.wk <- aperm(qd.r.wk,c(1,3,2))
        se.aux <- .Fortran("llrmaux",
                           as.double(cd), as.integer(nxis),
                           as.double(10^lambda*q.wk), as.integer(nxi),
                           as.double(qd.r.wk), as.integer(nqd),
                           as.integer(nx), as.double(qd.wt),
                           as.double(.Machine$double.eps), double(nqd*nx),
                           double(nx), double(nxis),
                           v=double(nxis*nxis), double(nxis*nxis),
                           jpvt=integer(nxis), PACKAGE="gss")[c("v","jpvt")]
        c <- cd[1:nxi]
        if (nnull) d <- cd[nxi+(1:nnull)]
        else d <- NULL
        eta <- matrix(eta,nqd,nx)
        for (i in 1:nx) eta[,i] <- eta[,i]/sum(eta[,i])
        return(list(lambda=lambda,theta=theta,c=c,d=d,cv=zz$min,
                    fit=t(eta),se.aux=se.aux))
    }
    ## theta search
    counter <- 0
    r.wk <- 0
    qd.r.wk <- array(0,c(nqd,nxi,nx))
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
        if (length(dim(qd.r[[i]]))==3) qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[[i]]
        else qd.r.wk <- qd.r.wk + as.vector(10^theta[i]*qd.r[[i]])
    }
    theta.old <- theta
    tmp <- abs(cv.m(theta))
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
            warning("gss warning in ssllrm: CV iteration fails to converge")
            break
        }
    }
    ## return
    q.wk <- 0
    qd.r.wk <- array(0,c(nqd,nxi,nx))
    for (i in 1:nq) {
        q.wk <- q.wk + 10^theta[i]*r[id.basis,,i]
        if (length(dim(qd.r[[i]]))==3) qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[[i]]
        else qd.r.wk <- qd.r.wk + as.vector(10^theta[i]*qd.r[[i]])
    }
    qd.r.wk <- aperm(qd.r.wk,c(1,3,2))
    qd.r.wk <- array(c(qd.r.wk,qd.s),c(nqd,nx,nxis))
    qd.r.wk <- aperm(qd.r.wk,c(1,3,2))
    se.aux <- .Fortran("llrmaux",
                       as.double(cd), as.integer(nxis),
                       as.double(10^lambda*q.wk), as.integer(nxi),
                       as.double(qd.r.wk), as.integer(nqd),
                       as.integer(nx), as.double(qd.wt),
                       as.double(.Machine$double.eps), double(nqd*nx),
                       double(nx), double(nxis),
                       v=double(nxis*nxis), double(nxis*nxis),
                       jpvt=integer(nxis), PACKAGE="gss")[c("v","jpvt")]
    c <- cd[1:nxi]
    if (nnull) d <- cd[nxi+(1:nnull)]
    else d <- NULL
    cv <- (zz$min-cv.shift)/cv.scale
    eta <- matrix(eta,nqd,nx)
    for (i in 1:nx) eta[,i] <- eta[,i]/sum(eta[,i])
    list(lambda=lambda,theta=zz$est,c=c,d=d,cv=cv,
         fit=t(eta),se.aux=se.aux)
}

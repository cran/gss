## Fit hazard model
sshzd <- function(formula,type="cubic",data=list(),alpha=1.4,
                  weights=NULL,subset,na.action=na.omit,
                  id.basis=NULL,nbasis=NULL,seed=NULL,
                  ext=.05,order=2,prec=1e-7,maxiter=30)
{
    ## Local functions handling formula
    Surv <- function(time,status,start=0) {
        tname <- as.character(as.list(match.call())$time)
        if (!is.numeric(time)|!is.vector(time))
            stop("gss error in sshzd: time should be a numerical vector")
        if ((nobs <- length(time))-length(status))
            stop("gss error in sshzd: time and status mismatch in size")
        if ((length(start)-nobs)&(length(start)-1))
            stop("gss error in sshzd: time and start mismatch in size")
        if (any(start>time))
            stop("gss error in sshzd: start after follow-up time")
        if (min(start)<0)
            warning("gss warning in sshzd: start before time 0")
        time <- cbind(start,time)
        list(tname=tname,start=time[,1],end=time[,2],status=as.logical(status))
    }
    ## Obtain model frame and model terms
    mf <- match.call()
    mf$type <- mf$alpha <- NULL
    mf$id.basis <- mf$nbasis <- mf$seed <- NULL
    mf$ext <- mf$order <- mf$prec <- mf$maxiter <- NULL
    term.wk <- terms.formula(mf$formula)
    ## response
    resp <- attr(term.wk,"variable")[[2]]
    ind.wk <- length(strsplit(deparse(resp),'')[[1]])
    if ((substr(deparse(resp),1,5)!='Surv(')
        |(substr(deparse(resp),ind.wk,ind.wk)!=')'))
        stop("gss error in sshzd: response should be Surv(...)")
    attach(data)
    yy <- eval(resp)
    tname <- yy$tname
    ## model frame
    term.labels <- attr(term.wk,"term.labels")
    if (!(tname%in%term.labels))
        stop("gss error in sshzd: time main effect missing in model")
    mf[[1]] <- as.name("model.frame")
    mf[[2]] <- eval(parse(text=paste("~",paste(term.labels,collapse="+"))))
    mf <- eval(mf,sys.frame(sys.parent()))
    ## set domain
    xnames <- names(mf)
    xnames <- xnames[!xnames%in%tname]
    domain <- as.list(NULL)
    mn <- min(yy$start)
    mx <- max(yy$end)
    domain[[tname]] <- c(mn,mx)
    for (i in xnames) {
        if (is.factor(mf[[i]])) domain[[i]] <- levels(mf[[i]])[1:2]
        else {
            mn <- min(mf[[i]])
            mx <- max(mf[[i]])
            range <- mx-mn
            mn <- mn - ext*range
            mx <- mx + ext*range
            domain[[i]] <- c(mn,mx)
        }
    }
    domain <- as.data.frame(domain)
    ## Generate sub-basis
    cnt <- model.weights(mf)
    nobs <- nrow(mf)
    if (is.null(id.basis)) {
        if (is.null(nbasis)) nbasis <- max(30,ceiling(10*nobs^(2/9)))
        if (nbasis>sum(yy$status)) nbasis <- sum(yy$status)
        if (!is.null(seed)) set.seed(seed)
        id.basis <- sample((1:nobs)[yy$status],nbasis,prob=cnt[yy$status])
    }
    else {
        if (!all(id.basis%in%(1:nobs)[yy$status]))
            stop("gss error in sshzd: id.basis not all at failure cases")
        nbasis <- length(id.basis)
    }
    ## Generate terms    
    if (type=="cubic") term <- mkterm.cubic1(mf,domain)
    if (type=="linear") term <- mkterm.linear1(mf,domain)
    if (type=="tp") term <- mkterm.tp(mf,order,mf[id.basis,],1)
    if (is.null(term)) stop("gss error in sshzd: unknown type")
    ## Generate Gauss-Legendre quadrature
    nmesh <- 200
    quad <- gauss.quad(nmesh,domain[[tname]])
    ## obtain unique covariate observations
    if (length(xnames)) {
        xx <- mf[,xnames,drop=FALSE]
        x.pt <- unique(xx)
        nx <- dim(x.pt)[1]
        x.dup.ind <- duplicated(xx)
        x.dup <- xx[x.dup.ind,,drop=FALSE]
        ## xx[i,]==x.pt[x.ind[i],]
        x.ind <- 1:nobs
        x.ind[!x.dup.ind] <- 1:nx
        if (nobs-nx) {
            x.ind.wk <- 1:(nobs-nx)
            for (i in 1:nx) {
                for (j in 1:(nobs-nx)) {
                    if (sum(duplicated(rbind(x.pt[i,],x.dup[j,]))))
                        x.ind.wk[j] <- i
                }
            }
        }
        if (nobs-nx) x.ind[x.dup.ind] <- x.ind.wk
    }
    else {
        nx <- 1
        x.ind <- rep(1,nobs)
        x.pt <- NULL
    }
    ## integration weights at x.pt[i,]
    qd.wt <- matrix(0,nmesh,nx)
    for (i in 1:nobs) {
        wk <- (quad$pt<=yy$end[i])&(quad$pt>yy$start[i])
        if (is.null(cnt)) qd.wt[,x.ind[i]] <- qd.wt[,x.ind[i]]+wk
        else qd.wt[,x.ind[i]] <- qd.wt[,x.ind[i]]+cnt[i]*wk
    }
    if (is.null(cnt)) qd.wt <- quad$wt*qd.wt/nobs
    else qd.wt <- quad$wt*qd.wt/sum(cnt)
    ## Generate s, r, and q
    s <- r <- q <- qd.s <- NULL
    qd.r <- as.list(NULL)
    nT <- sum(yy$status)
    nq <- nu <- 0
    for (label in term$labels) {
        if (label=="1") {
            nu <- nu+1
            s <- cbind(s,rep(1,len=nT))
            qd.wk <- matrix(1,nmesh,nx)
            qd.s <- array(c(qd.s,qd.wk),c(nmesh,nx,nu))
            next
        }
        vlist <- term[[label]]$vlist
        x.list <- xnames[xnames%in%vlist]
        xy <- mf[yy$status,vlist]
        xy.basis <- mf[id.basis,vlist]
        qd.xy <- data.frame(matrix(0,nmesh,length(vlist)))
        names(qd.xy) <- vlist
        if (tname%in%vlist) qd.xy[,tname] <- quad$pt
        if (length(x.list)) xx <- x.pt[,x.list,drop=FALSE]
        else xx <- NULL
        nphi <- term[[label]]$nphi
        nrk <- term[[label]]$nrk
        if (nphi) {
            phi <- term[[label]]$phi
            for (i in 1:nphi) {
                nu <- nu+1
                s <- cbind(s,phi$fun(xy,nu=i,env=phi$env))
                if (is.null(xx))
                    qd.wk <- matrix(phi$fun(qd.xy[,,drop=TRUE],nu=i,env=phi$env),nmesh,nx)
                else {
                    qd.wk <- NULL
                    for (j in 1:nx) {
                        qd.xy[,x.list] <- xx[rep(j,nmesh),]
                        qd.wk <- cbind(qd.wk,phi$fun(qd.xy[,,drop=TRUE],i,phi$env))
                    }
                }
                qd.s <- array(c(qd.s,qd.wk),c(nmesh,nx,nu))
            }
        }
        if (nrk) {
            rk <- term[[label]]$rk
            for (i in 1:nrk) {
                nq <- nq+1
                r <- array(c(r,rk$fun(xy,xy.basis,nu=i,env=rk$env,out=TRUE)),c(nT,nbasis,nq))
                q <- array(c(q,rk$fun(xy.basis,xy.basis,nu=i,env=rk$env,out=TRUE)),
                           c(nbasis,nbasis,nq))
                if (is.null(xx))
                    qd.r[[nq]] <- rk$fun(qd.xy[,,drop=TRUE],xy.basis,i,rk$env,out=TRUE)
                else {
                    qd.wk <- NULL
                    for (j in 1:nx) {
                        qd.xy[,x.list] <- xx[rep(j,nmesh),]
                        qd.wk <- array(c(qd.wk,rk$fun(qd.xy[,,drop=TRUE],xy.basis,i,rk$env,TRUE)),
                                       c(nmesh,nbasis,j))
                    }
                    qd.r[[nq]] <- qd.wk
                }
            }
        }
    }
    if (!is.null(s)) {
        nnull <- dim(s)[2]
        ## Check s rank
        if (qr(s)$rank<nnull)
            stop("gss error in sscden: fixed effect MLE is not unique")
    }
    ## Fit the model
    Nobs <- ifelse(is.null(cnt),nobs,sum(cnt))
    if (!is.null(cnt)) cntt <- cnt[yy$status]
    else cntt <- NULL
    z <- msphzd(s,r,q,Nobs,cntt,qd.s,qd.r,qd.wt,prec,maxiter,alpha)
    cfit <- sum(yy$status)/Nobs/sum(qd.wt)
    ## Brief description of model terms
    desc <- NULL
    for (label in term$labels)
        desc <- rbind(desc,as.numeric(c(term[[label]][c("nphi","nrk")])))
    desc <- rbind(desc,apply(desc,2,sum))
    rownames(desc) <- c(term$labels,"total")
    colnames(desc) <- c("Unpenalized","Penalized")
    ## Return the results
    obj <- c(list(call=match.call(),mf=mf,tname=tname,xnames=xnames,
                  terms=term,desc=desc,alpha=alpha,domain=domain,cfit=cfit,
                  quad=quad,x.pt=x.pt,qd.wt=qd.wt,id.basis=id.basis),z)
    if (is.null(cnt)) obj$se.aux$v <- sqrt(nobs)*obj$se.aux$v
    else obj$se.aux$v <- sqrt(sum(cnt))*obj$se.aux$v
    class(obj) <- c("sshzd")
    obj
}

## Fit (multiple smoothing parameter) hazard function
msphzd <- function(s,r,q,Nobs,cnt,qd.s,qd.r,qd.wt,prec,maxiter,alpha)
{
    nT <- dim(r)[1]
    nxi <- dim(r)[2]
    nqd <- dim(qd.wt)[1]
    nx <- dim(qd.wt)[2]
    if (!is.null(s)) nnull <- dim(s)[2]
    else nnull <- 0
    nxis <- nxi+nnull
    if (is.null(cnt)) cnt <- 0
    ## cv functions
    cv.s <- function(lambda) {
        fit <- .Fortran("hzdnewton",
                        cd=as.double(cd), as.integer(nxis),
                        as.double(10^lambda*q.wk), as.integer(nxi),
                        as.double(t(cbind(r.wk,s))), as.integer(nT),
                        as.integer(Nobs), as.integer(sum(cnt)), as.integer(cnt),
                        as.double(qd.r.wk), as.integer(nqd),
                        as.double(qd.wt), as.integer(nx),
                        as.double(prec), as.integer(maxiter),
                        as.double(.Machine$double.eps),
                        wk=double(2*(nqd*nx+nT)+nxis*(2*nxis+5)+max(nxis,2)),
                        info=integer(1),PACKAGE="gss")
        if (fit$info==1) stop("gss error in sshzd: Newton iteration diverges")
        if (fit$info==2) warning("gss warning in sshzd: Newton iteration fails to converge")
        assign("cd",fit$cd,inherit=TRUE)
        assign("mesh0",matrix(fit$wk[max(nxis,2)+(1:(nqd*nx))],nqd,nx),inherit=TRUE)
        cv <- alpha*fit$wk[2]-fit$wk[1]
        alpha.wk <- max(0,log.la0-lambda-5)*(3-alpha) + alpha
        alpha.wk <- min(alpha.wk,3)
        adj <- ifelse (alpha.wk>alpha,(alpha.wk-alpha)*fit$wk[2],0)
        cv+adj
    }
    cv.m <- function(theta) {
        r.wk <- q.wk <- 0
        qd.r.wk <- array(0,c(nqd,nxi,nx))
        for (i in 1:nq) {
            r.wk <- r.wk + 10^theta[i]*r[,,i]
            q.wk <- q.wk + 10^theta[i]*q[,,i]
            if (length(dim(qd.r[[i]]))==3) qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[[i]]
            else qd.r.wk <- qd.r.wk + as.vector(10^theta[i]*qd.r[[i]])
        }
        qd.r.wk <- aperm(qd.r.wk,c(1,3,2))
        qd.r.wk <- array(c(qd.r.wk,qd.s),c(nqd,nx,nxis))
        qd.r.wk <- aperm(qd.r.wk,c(1,3,2))
        fit <- .Fortran("hzdnewton",
                        cd=as.double(cd), as.integer(nxis),
                        as.double(10^lambda*q.wk), as.integer(nxi),
                        as.double(t(cbind(r.wk,s))), as.integer(nT),
                        as.integer(Nobs), as.integer(sum(cnt)), as.integer(cnt),
                        as.double(qd.r.wk), as.integer(nqd),
                        as.double(qd.wt), as.integer(nx),
                        as.double(prec), as.integer(maxiter),
                        as.double(.Machine$double.eps),
                        wk=double(2*(nqd*nx+nT)+nxis*(2*nxis+5)+max(nxis,2)),
                        info=integer(1),PACKAGE="gss")
        if (fit$info==1) stop("gss error in sshzd: Newton iteration diverges")
        if (fit$info==2) warning("gss warning in sshzd: Newton iteration fails to converge")
        assign("cd",fit$cd,inherit=TRUE)
        assign("mesh0",matrix(fit$wk[max(nxis,2)+(1:(nqd*nx))],nqd,nx),inherit=TRUE)
        cv <- alpha*fit$wk[2]-fit$wk[1]
        alpha.wk <- max(0,theta-log.th0-5)*(3-alpha) + alpha
        alpha.wk <- min(alpha.wk,3)
        adj <- ifelse (alpha.wk>alpha,(alpha.wk-alpha)*fit$wk[2],0)
        cv+adj
    }
    cv.wk <- function(theta) cv.scale*cv.m(theta)+cv.shift
    ## Initialization
    theta <- -log10(apply(q,3,function(x)sum(diag(x))))
    nq <- length(theta)
    qd.r.wk <- array(0,c(nqd,nxi,nx))
    for (i in 1:nq) {
        if (length(dim(qd.r[[i]]))==3) qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[[i]]
        else qd.r.wk <- qd.r.wk + as.vector(10^theta[i]*qd.r[[i]])
    }
    v.s <- v.r <- 0
    for (i in 1:nx) {
        if (nnull) v.s <- v.s + apply(qd.wt[,i]*qd.s[,i,,drop=FALSE]^2,2,sum)
        v.r <- v.r + apply(qd.wt[,i]*qd.r.wk[,,i,drop=FALSE]^2,2,sum)
    }
    if (nnull) theta.wk <- log10(sum(v.s)/nnull/sum(v.r)*nxi) / 2
    else theta.wk <- 0
    theta <- theta + theta.wk
    qd.r.wk <- aperm(10^theta.wk*qd.r.wk,c(1,3,2))
    qd.r.wk <- array(c(qd.r.wk,qd.s),c(nqd,nx,nxis))
    qd.r.wk <- aperm(qd.r.wk,c(1,3,2))
    r.wk <- q.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
        q.wk <- q.wk + 10^theta[i]*q[,,i]
    }
    log.la0 <- log10(sum(v.r)/sum(diag(q.wk))) + 2*theta.wk
    ## fixed theta iteration
    mesh0 <- NULL
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
        se.aux <- .Fortran("hzdaux1",
                           as.double(cd), as.integer(nxis),
                           as.double(10^lambda*q.wk), as.integer(nxi),
                           as.double(qd.r.wk), as.integer(nqd),
                           as.double(qd.wt), as.integer(nx),
                           as.double(.Machine$double.eps), double(nqd*nx),
                           v=double(nxis*nxis), double(nxis*nxis),
                           jpvt=integer(nxis), PACKAGE="gss")[c("v","jpvt")]
        c <- cd[1:nxi]
        if (nnull) d <- cd[nxi+(1:nnull)]
        return(list(lambda=zz$est,theta=theta,c=c,d=d,cv=zz$min,mesh0=mesh0,se.aux=se.aux))
    }
    ## theta adjustment
    qd.r.wk <- array(0,c(nqd,nxi,nx))
    for (i in 1:nq) {
        theta[i] <- 2*theta[i] + log10(t(cd[1:nxi])%*%q[,,i]%*%cd[1:nxi])
        if (length(dim(qd.r[[i]]))==3) qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[[i]]
        else qd.r.wk <- qd.r.wk + as.vector(10^theta[i]*qd.r[[i]])
    }
    v.s <- v.r <- 0
    for (i in 1:nx) {
        if (nnull) v.s <- v.s + apply(qd.wt[,i]*qd.s[,i,,drop=FALSE]^2,2,sum)
        v.r <- v.r + apply(qd.wt[,i]*qd.r.wk[,,i,drop=FALSE]^2,2,sum)
    }
    if (nnull) theta.wk <- log10(sum(v.s)/nnull/sum(v.r)*nxi) / 2
    else theta.wk <- 0
    theta <- theta + theta.wk
    qd.r.wk <- aperm(10^theta.wk*qd.r.wk,c(1,3,2))
    qd.r.wk <- array(c(qd.r.wk,qd.s),c(nqd,nx,nxis))
    qd.r.wk <- aperm(qd.r.wk,c(1,3,2))
    r.wk <- q.wk <- 0
    for (i in 1:nq) {
        r.wk <- r.wk + 10^theta[i]*r[,,i]
        q.wk <- q.wk + 10^theta[i]*q[,,i]
    }
    log.la0 <- log10(sum(v.r)/sum(diag(q.wk))) + 2*theta.wk
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
    ## theta search
    lambda <- zz$est
    counter <- 0
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
            warning("gss warning in sscden: CV iteration fails to converge")
            break
        }
    }
    ## return
    theta <- zz$est
    cv <- (zz$min-cv.shift)/cv.scale
    q.wk <- 0
    qd.r.wk <- array(0,c(nqd,nxi,nx))
    for (i in 1:nq) {
        q.wk <- q.wk + 10^theta[i]*q[,,i]
        if (length(dim(qd.r[[i]]))==3) qd.r.wk <- qd.r.wk + 10^theta[i]*qd.r[[i]]
        else qd.r.wk <- qd.r.wk + as.vector(10^theta[i]*qd.r[[i]])
    }
    qd.r.wk <- aperm(qd.r.wk,c(1,3,2))
    qd.r.wk <- array(c(qd.r.wk,qd.s),c(nqd,nx,nxis))
    qd.r.wk <- aperm(qd.r.wk,c(1,3,2))
    se.aux <- .Fortran("hzdaux1",
                       as.double(cd), as.integer(nxis),
                       as.double(10^lambda*q.wk), as.integer(nxi),
                       as.double(qd.r.wk), as.integer(nqd),
                       as.double(qd.wt), as.integer(nx),
                       as.double(.Machine$double.eps), double(nqd*nx),
                       v=double(nxis*nxis), double(nxis*nxis),
                       jpvt=integer(nxis), PACKAGE="gss")[c("v","jpvt")]
    c <- cd[1:nxi]
    if (nnull) d <- cd[nxi+(1:nnull)]
    list(lambda=lambda,theta=theta,c=c,d=d,cv=cv,mesh0=mesh0,se.aux=se.aux)
}
## Fit density model
ssden9 <- function(formula,data=list(),alpha=1.4,weights=NULL,
                   subset,na.action=na.omit,
                   id.basis=NULL,nbasis=NULL,seed=NULL,
                   domain=as.list(NULL),prec=1e-7,maxiter=30)
{
    ## Obtain model frame and model terms
    mf <- match.call()
    mf$alpha <- mf$id.basis <- mf$nbasis <- mf$seed <- NULL
    mf$domain <- mf$prec <- mf$maxiter <- mf$skip.iter <- NULL
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf,parent.frame())
    cnt <- model.weights(mf)
    mf$"(weights)" <- NULL
    ## Generate sub-basis
    nobs <- dim(mf)[1]
    if (is.null(id.basis)) {
        if (is.null(nbasis))  nbasis <- max(30,ceiling(10*nobs^(2/9)))
        if (nbasis>=nobs)  nbasis <- nobs
        if (!is.null(seed))  set.seed(seed)
        id.basis <- sample(nobs,nbasis,prob=cnt)
    }
    else {
        if (max(id.basis)>nobs|min(id.basis)<1)
            stop("gss error in ssden9: id.basis out of range")
        nbasis <- length(id.basis)
    }
    ## Set domain, generate quadrature, and make phi, rk
    varlist <- names(mf)
    quad <- phi <- rk <- list(NULL)
    for (xlab in varlist) {
        x <- mf[[xlab]]
        if (!is.vector(x))
            stop("gss error in ssden9: variables must be 1-D")
        if (is.factor(x)) {
            domain[[xlab]] <- NULL
            pt <- levels(x)
            wt <- rep(1,nlevels(x))
            quad[[xlab]] <- list(pt=pt,wt=wt)
            fun.wk <- mkrk.nominal(levels(x))
            fun <- function(x,nu,env) {
                wk <- as.factor(names(env$env$code)[nu])
                env$fun(x,wk,env$env)
            }
            phi[[xlab]] <- list(fun=fun,env=fun.wk)
            rk[[xlab]] <- NULL
        }
        else {
            if (is.null(domain[[xlab]])) {
                mn <- min(x)
                mx <- max(x)
                domain[[xlab]] <- c(mn,mx)+c(-1,1)*(mx-mn)*.05
            }
            else domain[[xlab]] <- c(min(domain[[xlab]]),max(domain[[xlab]]))
            quad[[xlab]] <- gauss.quad(200,domain[[xlab]])
            phi[[xlab]] <- mkphi.cubic(domain[[xlab]])
            rk[[xlab]] <- mkrk.cubic(domain[[xlab]])
        }
    }
    ## Generate s and r
    sr <- qd.sr <- qd.wt <- list(NULL)
    for (xlab in varlist) {
        x <- mf[[xlab]]
        if (is.factor(x)) {
            x.qd <- quad[[xlab]]$pt
            phi.wk <- phi[[xlab]]
            s <- qd.s <- NULL
            for (nu in 1:(nlevels(x)-1)) {
                s <- cbind(s,phi.wk$fun(x,nu,phi.wk$env))
                qd.s <- cbind(qd.s,phi.wk$fun(x.qd,nu,phi.wk$env))
            }
            sr[[xlab]] <- s
            qd.sr[[xlab]] <- qd.s
            qd.wt[[xlab]] <- quad[[xlab]]$wt
        }
        else {
            x.basis <- x[id.basis]
            x.qd <- quad[[xlab]]$pt
            phi.wk <- phi[[xlab]]
            rk.wk <- rk[[xlab]]
            s <- phi.wk$fun(x,1,phi.wk$env)
            qd.s <- phi.wk$fun(x.qd,1,phi.wk$env)
            r <- rk.wk$fun(x,x.basis,rk.wk$env,TRUE)
            qd.r <- rk.wk$fun(x.qd,x.basis,rk.wk$env,TRUE)
            sr[[xlab]] <- cbind(s,150*r)
            qd.sr[[xlab]] <- cbind(qd.s,150*qd.r)
            qd.wt[[xlab]] <- quad[[xlab]]$wt
        }
    }
    ## Fit the model
#    if (length(varlist)) stop("gss error in ssden9: use ssden for 1-D estimation")
    z <- mspdsty9(4,mf,sr,id.basis,cnt,qd.sr,qd.wt,prec,maxiter,alpha)
    ## Return the results
    obj <- c(list(call=match.call(),alpha=alpha,domain=domain,
                  id.basis=id.basis,phi=phi,rk=rk,mf=mf[id.basis,]),z)
    class(obj) <- c("ssden9","ssden")
    obj
}

## Fit multiple smoothing parameter density
mspdsty9 <- function(ncomp,mf,sr,id.basis,cnt,qd.sr,qd.wt,prec,maxiter,alpha)
{
    ## Initialization
    varlist <- names(mf)
    nobs <- dim(mf)[1]
    if (ncomp>2) {
        wk <- apply(matrix(runif((ncomp-1)*nobs),ncomp-1,nobs),2,sort)
        ez <- t(rbind(wk[1,],apply(rbind(wk,1),2,diff)))
    }
    else {
        wk <- runif(nobs)
        ez <- cbind(wk,1-wk)
    }
    ff <- eznew <- ez
    iter <- 0
    disc <- 1
    la.wk <- dc.wk <- int.wk <- list(NULL)
    ## EM iteration
    repeat {
        iter <- iter+1
        for (j in 1:ncomp) {
            ez.wk <- ez[,j]
            if (!is.null(cnt)) ez.wk <- ez.wk*cnt
            ff[,j] <- 1
            for (xlab in varlist) {
                if (is.factor(mf[[xlab]])) {
                    x <- mf[[xlab]]
                    sr.wk <- sr[[xlab]]
                    qd.sr.wk <- qd.sr[[xlab]]
                    if (is.null(dc.wk[[xlab]])) {
                        dc.wk[[xlab]] <- matrix(0,dim(sr.wk)[2],ncomp)
                        int.wk[[xlab]] <- rep(0,ncomp)
                    }
                    f.wk <- rep(0,nlevels(x))
                    names(f.wk) <- levels(x)
                    for (i in 1:nobs) f.wk[x[i]] <- f.wk[x[i]]+ez.wk[i]
                    dc.wk[[xlab]][,j] <-
                      solve(t(qd.sr.wk)%*%qd.sr.wk,t(qd.sr.wk)%*%log(f.wk))
                    int.wk[[xlab]][j] <- sum(exp(qd.sr.wk%*%dc.wk[[xlab]][,j]))
                    ff[,j] <- ff[,j]*exp(sr.wk%*%dc.wk[[xlab]][,j])/int.wk[[xlab]][j]
                }
                else {
                    sr.wk <- sr[[xlab]]
                    q <- sr.wk[id.basis,-1]
                    qd.sr.wk <- qd.sr[[xlab]]
                    if (is.null(la.wk[[xlab]])) {
                        la.wk[[xlab]] <- int.wk[[xlab]] <- rep(0,ncomp)
                        dc.wk[[xlab]] <- matrix(0,dim(sr.wk)[2],ncomp)
                    }
                    if (disc>2e-1) init <- NULL
                    else init <- list(dc=dc.wk[[xlab]][,j],la=la.wk[[xlab]][j])
                    z <- sspdsty9(sr.wk,ez.wk,q,qd.sr.wk,qd.wt[[xlab]],prec,maxiter,alpha,init)
                    ff[,j] <- ff[,j]*exp(sr.wk%*%z$dc)/z$int
                    la.wk[[xlab]][j] <- z$lambda
                    dc.wk[[xlab]][,j] <- z$dc
                    int.wk[[xlab]][j] <- z$int
                }
            }
        }
        if (is.null(cnt)) ppnew <- apply(ez,2,mean)
        else ppnew <- apply(ez*cnt,2,sum)/sum(cnt)
        if (iter==1) pp <- ppnew
        for (i in 1:nobs) {
            ez.wk <- ff[i,]*pp
            eznew[i,] <- ez.wk/sum(ez.wk)
        }
        disc <- max(abs(ez-eznew)/(1+abs(ez)))
        if (disc<1e-2) break
        ez <- eznew
        lkhd <- -mean(log(ff%*%pp))
        lkhdnew <- -mean(log(ff%*%ppnew))
print(c(disc,lkhd,lkhdnew))
        if (lkhdnew<lkhd) pp <- ppnew
    }
    ## penalized likelihood of observed data
    if (is.null(cnt)) cntt <- 0
    else cntt <- cnt
    ez <- t(ff)*pp
    repeat {
        for (j in rev(order(pp))) {
            pp.wk <- pp
            pp.wk[j] <- 0
            f0 <- ff%*%pp.wk
            for (xlab in varlist) {
                f1 <- pp[j]*ff[,j]
                sr.wk <- sr[[xlab]]
                f1 <- f1/exp(sr.wk%*%dc.wk[[xlab]][,j])*int.wk[[xlab]][j]
                qd.sr.wk <- qd.sr[[xlab]]
                nsr <- dim(sr.wk)[2]
                nqd <- length(qd.wt[[xlab]])
                if (is.factor(mf[[xlab]])) q <- matrix(0,nsr,nsr)
                else q <- 10^la.wk[[xlab]][j]*sr.wk[id.basis,-1]*pp[j]
                fit <- .Fortran("dnewton99",
                                dc=as.double(dc.wk[[xlab]][,j]),
                                as.integer(nsr), as.double(t(sr.wk)), as.double(q),
                                as.integer(nobs), as.integer(sum(cntt)), as.integer(cntt),
                                as.double(cbind(f0,f1)), as.double(qd.sr.wk),
                                as.integer(nqd), as.double(qd.wt[[xlab]]),
                                as.double(prec), as.integer(maxiter),
                                as.double(.Machine$double.eps), integer(nsr),
                                double(2*(nqd+nobs)+nsr*(2*nsr+4)),
                                info=integer(1), PACKAGE="gss")
                if (fit$info==1) stop("gss error in ssden9: Newton iteration diverges")
                if (fit$info==2)
                    warning("gss warning in ssden9: Newton iteration fails to converge")
                dc.wk[[xlab]][,j] <- fit$dc
                int.wk[[xlab]][j] <- sum(qd.wt[[xlab]]*exp(qd.sr.wk%*%fit$dc))
                ff[,j] <- f1*exp(sr.wk%*%fit$dc)/int.wk[[xlab]][j]/pp[j]
            }
        }
        pp <- reweight9(pp,ff,cnt)
print(-mean(log(ff%*%pp)))
        eznew <- t(ff)*pp
        disc <- max(abs(ez-eznew)/(1+abs(ez)))
        if (disc<1e-7) break
        ez <- eznew
    }
    list(pp=pp,dc=dc.wk,int=int.wk)
}

sspdsty9 <- function(sr,ez,q,qd.sr,qdwt,prec,maxiter,alpha,init)
{
    nsr <- dim(sr)[2]
    nobs <- dim(sr)[1]
    nqd <- length(qdwt)
    ## cv function
    cv <- function(lambda) {
        fit <- .Fortran("dnewton9",
                        dc=as.double(dc), as.integer(nsr),
                        as.double(t(sr)), as.double(10^lambda*q),
                        as.integer(nobs), as.double(ez),
                        as.double(qd.sr), as.integer(nqd),
                        as.double(qdwt), as.double(prec),
                        as.integer(maxiter),
                        as.double(.Machine$double.eps),
                        integer(nsr),
                        wk=double(2*(nqd+nobs)+nsr*(nsr+4)),
                        info=integer(1), PACKAGE="gss")
        if (fit$info==1) stop("gss error in ssden9: Newton iteration diverges")
        if (fit$info==2) warning("gss warning in ssden9: Newton iteration fails to converge")
        assign("dc",fit$dc,inherits=TRUE)
        cv <- alpha*fit$wk[2]-fit$wk[1]
        alpha.wk <- max(0,log.la0-lambda-5)*(3-alpha) + alpha
        alpha.wk <- min(alpha.wk,3)
        adj <- ifelse (alpha.wk>alpha,(alpha.wk-alpha)*fit$wk[2],0)
        cv+adj
    }
    ## initialization
    mu.r <- apply(qdwt*qd.sr[,-1],2,sum)/sum(qdwt)
    v.r <- apply(qdwt*qd.sr[,-1]^2,2,sum)/sum(qdwt)
    log.la0 <- log10(sum(v.r-mu.r^2)/sum(diag(q)))
    ## lambda search
    if (is.null(init)) {
        dc <- rep(0,nsr)
        la <- log.la0
    }
    else {
        dc <- init$dc
        la <- init$la
    }
    mn0 <- log.la0-6
    mx0 <- log.la0+6
    repeat {
        mn <- max(la-1,mn0)
        mx <- min(la+1,mx0)
        zz <- nlm0(cv,c(mn,mx))
        if ((min(zz$est-mn,mx-zz$est)>=1e-1)||
            (min(zz$est-mn0,mx0-zz$est)<1e-1)) break
        else la <- zz$est
    }
    ## return
    jk1 <- cv(zz$est)
    int <- sum(qdwt*exp(qd.sr%*%dc))
    list(lambda=zz$est,dc=dc,int=int)
}

reweight9 <- function(pp,ff,cnt)
{
    m <- length(pp)
    if (m==2) {
        fun <- function(x) {
            if (is.null(cnt)) z <- -mean(log((1-x)*ff[,1]+x*ff[,2]))
            else z <- -sum(log((1-x)*ff[,1]+x*ff[,2])*cnt)/sum(cnt)
            z
        }
        z <- nlm0(fun,c(0,1))
        return(c(1-z$est,z$est))
    }
    nobs <- dim(ff)[1]
    fix <- rev(order(pp))[1]
    beta <- log(pp[-fix]/pp[fix])
    ppnew <- pp
    if (is.null(cnt)) lkhd <- -mean(log(ff%*%pp))
    else lkhd <- -sum(log(ff%*%pp)*cnt)/sum(cnt)
    repeat {
        ## gradient and hessian
        f0 <- NULL
        for (i in 1:nobs) f0 <- rbind(f0,ff[i,]*pp/sum(ff[i,]*pp))
        if (is.null(cnt)) {
            mu <- pp-apply(f0,2,mean)
            v <- 0
            for (i in 1:nobs) v <- v+diag(f0[i,])-outer(f0[i,],f0[i,])
            v <- diag(pp)-outer(pp,pp)-v/nobs
        }
        else {
            mu <- pp-apply(f0*cnt,2,sum)/sum(cnt)
            v <- 0
            for (i in 1:nobs) v <- v+cnt[i]*(diag(f0[i,])-outer(f0[i,],f0[i,]))
            v <- diag(pp)-outer(pp,pp)-v/sum(cnt)
        }
        mu <- mu[-fix]
        v <- v[-fix,-fix]
        ## modify hessian if necessary
        z <- .Fortran("dmcdc",
                      as.double(v), as.integer(m-1), as.integer(m-1),
                      ee=double(m-1),
                      pivot=integer(m-1),
                      integer(1), PACKAGE="gss")
        if (max(z$ee)) {
            z$ee[z$pivot] <- z$ee
            v <- v+diag(z$ee)
        }
        ## update beta
        mumax <- max(abs(mu))
        repeat {
            betanew <- beta-solve(v,mu)
            ppnew[-fix] <- exp(betanew)
            ppnew[fix] <- 1
            ppnew <- ppnew/sum(ppnew)
            if (is.null(cnt)) lkhdnew <- -mean(log(ff%*%ppnew))
            else lkhdnew <- -sum(log(ff%*%ppnew)*cnt)/sum(cnt)
            if (lkhdnew-lkhd<(1+abs(lkhd)*10*.Machine$double.eps)) break
            mu <- mu/2
            if (max(mu)/mumax<10*.Machine$double.eps) break
        }
        disc <- abs(lkhdnew-lkhd)/(1+abs(lkhd))
        disc <- max(disc,max(abs(pp-ppnew)/(1+abs(pp))))
        disc0 <- (mumax/(1+abs(lkhd)))^2
        pp <- ppnew
        beta <- betanew
        lkhd <- lkhdnew
        if (min(disc,disc0)<1e-7) break
    }
    pp
}  

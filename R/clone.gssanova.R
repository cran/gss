## Fit ssanova model
clone.gssanova <- function(object,nrep=1000,type=1,seed=NULL,eta=NULL,fix=FALSE,...)
{
    ## Extract relevant infomation from object
    mf <- object$mf
    nobs <- dim(mf)[1]
    term <- object$terms
    offset <- model.offset(mf)
    if (is.null(offset)) offset <- 0
    wt <- model.weights(mf)
    if (is.null(wt)) wt <- rep(1,nobs)
    ## Generate s and r
    s <- r <- NULL
    nxi <- length(object$id.basis)
    nq <- 0
    for (label in term$labels) {
        if (label=="1") {
            s <- cbind(s,rep(1,len=nobs))
            next
        }
        x <- mf[,term[[label]]$vlist]
        x.basis <- mf[object$id.basis,term[[label]]$vlist]
        nphi <- term[[label]]$nphi
        nrk <- term[[label]]$nrk
        if (nphi) {
            phi <- term[[label]]$phi
            for (i in 1:nphi)
                s <- cbind(s,phi$fun(x,nu=i,env=phi$env))
        }
        if (nrk) {
            rk <- term[[label]]$rk
            for (i in 1:nrk) {
                nq <- nq+1
                r <- array(c(r,rk$fun(x,x.basis,nu=i,env=rk$env,out=TRUE)),c(nobs,nxi,nq))
            }
        }
    }
    ## Partial term in s
    if (!is.null(object$partial)) {
        matx.p <- model.matrix(object$partial$mt.p)[,-1,drop=FALSE]
        matx.p <- scale(matx.p)
        s <- cbind(s,matx.p)
    }
    if (!is.null(s)) {
        if (is.vector(s)) nnull <- 1
        else nnull <- ncol(s)
    }
    else nnull <- 0
    ## Aggregate r
    rwk <- 0
    for (i in 1:nq) rwk <- rwk+10^object$theta[i]*r[,,i]
    ## Fitted fixed effects
    if (is.null(eta)) eta <- as.vector(s%*%object$d+rwk%*%object$c)+offset
    ## Random term
    if (!is.null(random<-object$random)) {
        nz <- ncol(as.matrix(random$z))
        sig.z <- random$sigma$fun(object$zeta,random$sigma$env)
        q.z <- 10^(2*object$ran.scal)*sig.z
        rwk.z <- 10^object$ran.scal*random$z
        sig.z <- t(chol(sig.z))
    }
    else nz <- 0
    ##
    nxiz <- nxi + nz
    nn <- nxiz + nnull
    if (fix) q <- rwk[object$id.basis,]
    ## Component roughness and weights
    qq <- r[object$id.basis,,,drop=FALSE]
    rho0 <- pwt <- NULL
    for (i in 1:nq) {
        rho0 <- c(rho0,log10(sum(object$c*as.vector(qq[,,i]%*%object$c))))
        pwt <- c(pwt,var(r[,,i]%*%object$c)*10^(2*object$theta[i]))
    }
    pwt <- pwt/sum(pwt)
    ## fixed theta
    fix1 <- function(nlam) {
        if (is.null(random)) qwk <- 10^(nlam)*q
        else {
            rwk <- cbind(rwk,rwk.z)
            qwk <- matrix(0,nxiz,nxiz)
            qwk[1:nxi,1:nxi] <- 10^(nlam)*q
            qwk[(nxi+1):nxiz,(nxi+1):nxiz] <- q.z
        }
        z <- ngreg(dc.wk,family,cbind(s,rwk),qwk,y,wt,offset,nu0,1.4)
        c <- z$dc[nnull+(1:nxi)]
        d <- z$dc[1:nnull]
        assign("dc.wk",z$dc,inherits=TRUE)
        assign("eta.wk",z$eta,inherits=TRUE)
        if (family=="polr") assign("nu0",z$nu,inherits=TRUE)
        if (type-1) {
            switch(family,
                   binomial=kl.binomial(eta,z$eta,y0$wt),
                   poisson=kl.poisson(eta,z$eta,wt),
                   Gamma=kl.Gamma(eta,z$eta,wt),
                   inverse.gaussian=kl.inverse.gaussian(eta,z$eta,wt),
                   nbinomial=kl.nbinomial(eta,z$eta,wt,y0$nu),
                   polr=kl.polr(list(eta=eta,nu=nu),list(eta=z$eta,nu=nu),wt),
                   weibull=kl.weibull(eta,z$eta,wt,nu,y0$int),
                   lognorm=kl.lognorm(eta,z$eta,wt,nu,y0),
                   loglogis=kl.loglogis(eta,z$eta,wt,nu,y0))
        }
        else {
            rho <- NULL
            for (i in 1:nq) rho <- c(rho,log10(sum(c*as.vector(qq[,,i]%*%c))))
            sum(pwt*(rho-rho0)^2)
        }
    }
    ## varying theta
    fix0 <- function(theta) {
        ## Aggregate r
        nq <- length(theta)
        rwk <- 0
        for (i in 1:nq) rwk <- rwk+10^theta[i]*r[,,i]
        q <- rwk[object$id.basis,]
        if (is.null(random)) qwk <- 10^(nlam)*q
        else {
            rwk <- cbind(rwk,rwk.z)
            qwk <- matrix(0,nxiz,nxiz)
            qwk[1:nxi,1:nxi] <- 10^(nlam)*q
            qwk[(nxi+1):nxiz,(nxi+1):nxiz] <- q.z
        }
        z <- ngreg(dc.wk,family,cbind(s,rwk),qwk,y,wt,offset,nu0,1.4)
        c <- z$dc[nnull+(1:nxi)]
        d <- z$dc[1:nnull]
        assign("dc.wk",z$dc,inherits=TRUE)
        assign("eta.wk",z$eta,inherits=TRUE)
        if (family=="polr") assign("nu0",z$nu,inherits=TRUE)
        if (type-1) {
            switch(family,
                   binomial=kl.binomial(eta,z$eta,y0$wt),
                   poisson=kl.poisson(eta,z$eta,wt),
                   Gamma=kl.Gamma(eta,z$eta,wt),
                   inverse.gaussian=kl.inverse.gaussian(eta,z$eta,wt),
                   nbinomial=kl.nbinomial(eta,z$eta,wt,y0$nu),
                   polr=kl.polr(list(eta=eta,nu=nu),list(eta=z$eta,nu=nu),wt),
                   weibull=kl.weibull(eta,z$eta,wt,nu,y0$int),
                   lognorm=kl.lognorm(eta,z$eta,wt,nu,y0),
                   loglogis=kl.loglogis(eta,z$eta,wt,nu,y0))
        }
        else {
            rho <- NULL
            for (i in 1:nq) rho <- c(rho,log10(sum(c*as.vector(qq[,,i]%*%c))))
            sum(pwt*(rho-rho0+2*(theta-object$theta))^2)
        }
    }
    ## sampling environment
    family <- object$family
    nu <- object$nu
    nu0 <- list(NULL,FALSE)
    if (family!="polr") y <- model.response(mf,"numeric")
    env <- NULL
    if (family=="binomial") {
        if (!is.vector(y)) env <- apply(y[,1:2],1,sum)
        else env <- rep(1,nobs)
        env <- wt*env
        object$nu <- 1
    }
    if (family%in%c("Gamma","inverse.gaussian")) env <- object$varht
    if (family=="nbinomial") {
        if (!is.vector(y)) env <- y[,2]
        else env <- rep(nu,length(y))
    }
    if (family=="polr") {
        y <- model.response(mf)
        env <- list(lvls=levels(y),cut=cumsum(c(0,nu)))
        nu0[[1]] <- nu
    }
    if (family%in%c("weibull","lognorm","loglogis")) {
        env <- list(nu=object$nu,y=y)
        nu0[[1]] <- nu
    }
    ## cloning    
    dc <- y9 <- eta9 <- kl <- theta <- NULL
    if (is.null(seed)) {
        seed <- abs(object$c[1])+abs(rev(object$c)[1])
        seed <- round((seed/ceiling(seed))^2,6)*1000000
        seed <- seed%%10000
    }
    set.seed(seed)
    for (i in 1:nrep) {
        ## sampling
        if (!is.null(random)) eta <- eta + random$z%*%(sig.z%*%rnorm(nz))
        if (family=="binomial") {
            y <- rbinom(eta,env,plogis(eta))/env
            wt <- env
        }
        if (family=="poisson") y <- rpois(eta,exp(eta))
        if (family=="Gamma") y <- rgamma(eta,shape=1/env,scale=env*exp(eta))
        if (family=="inverse.gaussian") y <- rinvgauss(eta,exp(eta),,env)
        if (family=="nbinomial") {
            ywk <- rnbinom(eta,env,plogis(eta))
            y <- cbind(ywk,env)
        }
        if (family=="polr") {
            ywk <- outer(-eta+rlogis(eta),env$cut,">")
            y <- ordered(env$lvls[apply(ywk,1,sum)+1],levels=env$lvls)
            y <- outer(y,env$lvls,"==")
        }
        if (family%in%c("weibull","lognorm","loglogis")) {
            y <- env$y
            ywk <- switch(family,
                          weibull=rweibull(eta,env$nu,exp(eta)),
                          lognorm=exp(rnorm(eta)/env$nu+eta),
                          loglogis=exp(rlogis(eta)/env$nu+eta))
            if (dim(y)[2]==3) {
                ok <- ywk>y[,3]
                while(m<-sum(!ok)) {
                    ywk[!ok] <- switch(family,
                                       weibull=rweibull(m,env$nu,exp(eta[!ok])),
                                       lognorm=exp(rnorm(m)/env$nu+eta[!ok]),
                                       loglogis=exp(rlogis(m)/env$nu+eta[!ok]))
                    ok <- ywk>y[,3]
                }
            }
            idx <- (1:dim(y)[1])[y[,2]==0]
            idx <- idx[ywk[idx]>y[idx,1]]
            y[-idx,1] <- ywk[-idx]
            y[-idx,2] <- 1
        }
        ## fitting
        eta.wk <- NULL
        dc.wk <- c(object$d,object$c,object$b)
        nlam <- object$nlambda
        y0 <- switch(family,
                     binomial=y0.binomial(y,eta,wt),
                     poisson=y0.poisson(eta),
                     Gamma=y0.Gamma(eta),
                     inverse.gaussian=y0.inverse.gaussian(eta),
                     nbinomial=y0.nbinomial(y,eta,nu),
                     polr=y0.polr(list(eta=eta,nu=nu)),
                     weibull=y0.weibull(y,eta,nu),
                     lognorm=y0.lognorm(y,eta,nu),
                     loglogis=y0.loglogis(y,eta,nu))
        if (fix) z <- nlm(fix1,nlam,stepmax=0.5)
        else z <- nlm(fix0,object$theta,stepmax=0.5)
        if (type-1) kl.wk <- as.vector(z[[1]])
        else {
            kl.wk <- switch(family,
                            binomial=kl.binomial(eta,eta.wk,y0$wt),
                            poisson=kl.poisson(eta,eta.wk,wt),
                            Gamma=kl.Gamma(eta,eta.wk,wt),
                            inverse.gaussian=kl.inverse.gaussian(eta,eta.wk,wt),
                            nbinomial=kl.nbinomial(eta,eta.wk,wt,y0$nu),
                            polr=kl.polr(list(eta=eta,nu=nu),list(eta=eta.wk,nu=nu),wt),
                            weibull=kl.weibull(eta,eta.wk,wt,nu,y0$int),
                            lognorm=kl.lognorm(eta,eta.wk,wt,nu,y0),
                            loglogis=kl.loglogis(eta,eta.wk,wt,nu,y0))
        }
        if (!fix) theta <- cbind(theta,as.vector(z[[2]]))
        kl <- c(kl,kl.wk)
        dc <- cbind(dc,dc.wk)
        y9 <- c(y9,list(y))
        eta9 <- c(eta9,list(eta.wk))
    }
    c <- dc[nnull+(1:nxi),]
    if (nnull) d <- dc[1:nnull,,drop=FALSE]
    else d <- NULL
    if (nz) b <- 10^(object$ran.scal)*dc[nnull+nxi+(1:nz),,drop=FALSE]
    else b <- NULL
    obj <- list(fit=object,nrep=nrep,c=c,d=d,b=b,y=y9,eta=eta9,kl=kl,type=type,theta=theta)
    class(obj) <- "clone"
    obj
}
## routine stolen from package statmod
rinvgauss <- function(n, mean=1, shape=NULL, dispersion=1)
#	Random variates from inverse Gaussian distribution
#	Gordon Smyth
#	Created 15 Jan 1998.  Last revised 27 Feb 2017.
{
#	Dispersion is reciprocal of shape
	if(!is.null(shape)) dispersion <- 1/shape

#	Check n
	if(length(n)>1L) n <- length(n) else n <- as.integer(n)
	if(n<0L) stop("n can't be negative")
	if(n==0L || length(mean)==0L || length(dispersion)==0L) return(numeric(0L))

#	Make arguments same length
	mu <- rep_len(mean,n)
	phi <- rep_len(dispersion,n)

#	Setup output vector
	r <- rep_len(0,n)

#	Non-positive parameters give NA
	mu.ok <- (mu > 0 & is.finite(mu))
	phi.ok <- (phi > 0 & is.finite(phi))
	i <- (mu.ok & phi.ok)
	if(!all(i)) {
		j <- !i
#		Infinite mu is special case
		invchisq <- (mu[j]==Inf & phi.ok[j])
		invchisq[is.na(invchisq)] <- FALSE
		if(any(invchisq)) {
			m <- sum(invchisq)
			r[j][invchisq] <- rnorm(m)^(-2) / phi[j][invchisq]
			j[j][invchisq] <- FALSE
		}
		infdisp <- (phi[j]==Inf)
		infdisp[is.na(infdisp)] <- FALSE
		if(any(infdisp)) {
			r[j][infdisp] <- 0
			j[j][infdisp] <- FALSE
		}
		r[j] <- NA
		n <- sum(i)
		if(n==0L) return(r)
	}

#	Generate chisquare on 1 df
	Y <- rnorm(n)^2

#	Divide out mu	
	Yphi <- Y*phi[i]*mu[i]

#	Taylor series is more accurate when Y*phi is large
	bigphi <- (Yphi > 5e5)
	if(any(bigphi)) {
		X1 <- Y
		X1[bigphi] <- 1 / Yphi[bigphi]
		X1[!bigphi] <- 1 + Yphi[!bigphi]/2 * (1 - sqrt(1 + 4/Yphi[!bigphi]))
	} else {
		X1 <- 1 + Yphi/2 * (1 - sqrt(1 + 4/Yphi))
	}
	firstroot <- (runif(n) < 1/(1+X1))
	r[i][firstroot] <- X1[firstroot]
	r[i][!firstroot] <- 1/X1[!firstroot]

#	Add mu back in again
	r[i] <- mu[i]*r[i]

	r
}

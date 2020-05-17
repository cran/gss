##%%%%%%%%%%  Binomial Family %%%%%%%%%%

## Make pseudo data for logistic regression
mkdata.binomial <- function(y,eta,wt,offset)
{
    if (is.vector(y)) y <- as.matrix(y)
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    if (is.null(offset)) offset <- rep(0,dim(y)[1])
    if (dim(y)[2]==1) {
        if ((max(y)>1)|(min(y)<0))
            stop("gss error: binomial responses should be between 0 and 1")
    }
    else {
        if (min(y)<0)
            stop("gss error: paired binomial response should be nonnegative")
        wt <- wt * (y[,1]+y[,2])
        y <- y[,1]/(y[,1]+y[,2])
    }
    odds <- exp(eta)
    p <- odds/(1+odds)
    u <- p - y
    w <- p/(1+odds)
    ywk <- eta-u/w-offset
    wt <- w*wt
    list(ywk=ywk,wt=wt,u=u*wt)
}

## Calculate deviance residuals for logistic regression
dev.resid.binomial <- function(y,eta,wt)
{
    if (is.vector(y)) y <- as.matrix(y)
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    if (dim(y)[2]>1) {
        wt <- wt * (y[,1]+y[,2])
        y <- y[,1]/(y[,1]+y[,2])
    }
    odds <- exp(eta)
    as.vector(2*wt*(y*log(ifelse(y==0,1,y*(1+odds)/odds))
                    +(1-y)*log(ifelse(y==1,1,(1-y)*(1+odds)))))
}

## Calculate null deviance for logistic regression
dev.null.binomial <- function(y,wt,offset)
{
    if (is.vector(y)) y <- as.matrix(y)
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    if (dim(y)[2]>1) {
        wt <- wt * (y[,1]+y[,2])
        y <- y[,1]/(y[,1]+y[,2])
    }
    p <- sum(wt*y)/sum(wt)
    odds <- p/(1-p)
    if (!is.null(offset)) {
        eta <- log(odds) - mean(offset)
        repeat {
            odds <- exp(eta+offset)
            p <- odds/(1+odds)
            u <- p - y
            w <- p/(1+odds)
            eta.new <- eta-sum(wt*u)/sum(wt*w)
            if (abs(eta-eta.new)/(1+abs(eta))<1e-7) break
            eta <- eta.new    
        }
    }
    sum(2*wt*(y*log(ifelse(y==0,1,y*(1+odds)/odds))
              +(1-y)*log(ifelse(y==1,1,(1-y)*(1+odds)))))
}


##%%%%%%%%%%  Poisson Family %%%%%%%%%%

## Make pseudo data for Poisson regression
mkdata.poisson <- function(y,eta,wt,offset)
{
    if (is.null(wt)) wt <- rep(1,length(y))
    if (is.null(offset)) offset <- rep(0,length(y))
    if (min(y)<0)
        stop("gss error: Poisson response should be nonnegative")
    lambda <- exp(eta)
    u <- lambda - y
    w <- lambda
    ywk <- eta-u/w-offset
    wt <- w*wt
    list(ywk=ywk,wt=wt,u=u*wt)
}

## Calculate deviance residuals for Poisson regression
dev.resid.poisson <- function(y,eta,wt)
{
    if (is.null(wt)) wt <- rep(1,length(y))
    lambda <- exp(eta)
    as.vector(2*wt*(y*log(ifelse(y==0,1,y/lambda))-(y-lambda)))
}

## Calculate null deviance for Poisson regression
dev.null.poisson <- function(y,wt,offset)
{
    if (is.null(wt)) wt <- rep(1,length(y))
    lambda <- sum(wt*y)/sum(wt)
    if (!is.null(offset)) {
        eta <- log(lambda) - mean(offset)
        repeat {
            lambda <- exp(eta+offset)
            u <- lambda - y
            w <- lambda
            eta.new <- eta-sum(wt*u)/sum(wt*w)
            if (abs(eta-eta.new)/(1+abs(eta))<1e-7) break
            eta <- eta.new    
        }
    }
    sum(2*wt*(y*log(ifelse(y==0,1,y/lambda))-(y-lambda)))
}


##%%%%%%%%%%  Gamma Family %%%%%%%%%%

## Make pseudo data for Gamma regression
mkdata.Gamma <- function(y,eta,wt,offset)
{
    if (is.null(wt)) wt <- rep(1,length(y))
    if (is.null(offset)) offset <- rep(0,length(y))
    if (min(y)<=0)
        stop("gss error: gamma responses should be positive")
    mu <- exp(eta)
    u <- 1-y/mu
    ywk <- eta-u-offset
    list(ywk=ywk,wt=wt,u=u*wt)
}

## Calculate deviance residuals for Gamma regression
dev.resid.Gamma <- function(y,eta,wt)
{
    if (is.null(wt)) wt <- rep(1,length(y))
    mu <- exp(eta)
    as.vector(2*wt*(-log(y/mu)+(y-mu)/mu))
}

## Calculate null deviance for Gamma regression
dev.null.Gamma <-
function(y,wt,offset) {
  if (is.null(wt)) wt <- rep(1,length(y))
  mu <- sum(wt*y)/sum(wt)
  if (!is.null(offset)) {
    eta <- log(mu)-mean(offset)
    repeat {
      mu <- exp(eta+offset)
      u <- 1-y/mu
      eta.new <- eta-sum(wt*u)/sum(wt)
      if (abs(eta-eta.new)/(1+abs(eta))<1e-7) break
      eta <- eta.new    
    }
  }
  sum(2*wt*(-log(y/mu)+(y-mu)/mu))
}


##%%%%%%%%%%  Inverse Gaussian Family %%%%%%%%%%

## Make pseudo data for IG regression
mkdata.inverse.gaussian <- function(y,eta,wt,offset)
{
    if (is.null(wt)) wt <- rep(1,length(y))
    if (is.null(offset)) offset <- rep(0,length(y))
    if (min(y)<=0)
        stop("gss error: inverse gaussian responses should be positive")
    mu <- exp(eta)
    u <- (1-y/mu)/mu
    w <- 1/mu
    ywk <- eta-u/w-offset
    wt <- w*wt
    list(ywk=ywk,wt=wt,u=u*wt)
}

## Calculate deviance residuals for IG regression
dev.resid.inverse.gaussian <- function(y,eta,wt)
{
    if (is.null(wt)) wt <- rep(1,length(y))
    mu <- exp(eta)
    as.vector(wt*((y-mu)^2/(y*mu^2)))
}

## Calculate null deviance for IG regression
dev.null.inverse.gaussian <- function(y,wt,offset)
{
    if (is.null(wt)) wt <- rep(1,length(y))
    mu <- sum(wt*y)/sum(wt)
    if (!is.null(offset)) {
        eta <- log(mu)-mean(offset)
        repeat {
            mu <- exp(eta+offset)
            u <- (1-y/mu)/mu
            w <- 1/mu
            eta.new <- eta-sum(wt*u)/sum(wt*w)
            if (abs(eta-eta.new)/(1+abs(eta))<1e-7) break
            eta <- eta.new    
        }
    }
    sum(wt*((y-mu)^2/(y*mu^2)))
}


##%%%%%%%%%%  Negative Binomial Family %%%%%%%%%%

## Make pseudo data for NB regression
mkdata.nbinomial <- function(y,eta,wt,offset,nu)
{
    if (is.vector(y)) y <- as.matrix(y)
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    if (is.null(offset)) offset <- rep(0,dim(y)[1])
    if (dim(y)[2]==2) {
        if (min(y[,1])<0)
            stop("gss error: negative binomial response should be nonnegative")
        if (min(y[,2])<=0)
            stop("gss error: negative binomial size should be positive")
        odds <- exp(eta)
        p <- odds/(1+odds)
        q <- 1/(1+odds)
        u <- y[,1]*p-y[,2]*q
        w <- y[,2]*q
        ywk <- eta-u/w-offset
        wt <- w*wt
        list(ywk=ywk,wt=wt,nu=1,u=u*wt)
    }
    else {
        if (min(y)<0)
            stop("gss error: negative binomial response should be nonnegative")
        odds <- exp(eta)
        p <- odds/(1+odds)
        q <- 1/(1+odds)
        if (is.null(nu)) log.nu <- log(mean(y*odds))
        else log.nu <- log(nu)
        lkhd <- function(log.nu) {
            nu <- exp(log.nu)
            lgamma(nu)-sum(wt*lgamma(nu+y))/sum(wt)-nu*sum(wt*log(p))/sum(wt)
        }
        nu <- exp(nlm(lkhd,log.nu,stepmax=.5)$est)
        u <- y*p-nu*q
        w <- nu*q
        ywk <- eta-u/w-offset
        wt <- w*wt
        list(ywk=ywk,wt=wt,nu=nu,u=u*wt)
    }
}

## Calculate deviance residuals for NB regression
dev.resid.nbinomial <- function(y,eta,wt)
{
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    odds <- exp(eta)
    p <- odds/(1+odds)
    q <- 1/(1+odds)
    as.vector(2*wt*(y[,1]*log(ifelse(y[,1]==0,1,y[,1]/(y[,1]+y[,2])/q))
                    +y[,2]*log(y[,2]/(y[,1]+y[,2])/p)))
}

## Calculate null deviance for NB regression
dev.null.nbinomial <- function(y,wt,offset)
{
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    p <- sum(wt*y[,2])/sum(wt*y)
    if (!is.null(offset)) {
        eta <- log(p/(1-p)) - mean(offset)
        repeat {
            odds <- exp(eta+offset)
            p <- odds/(1+odds)
            q <- 1/(1+odds)
            u <- y[,1]*p-y[,2]*q
            w <- y[,2]*q
            eta.new <- eta-sum(wt*u)/sum(wt*w)
            if (abs(eta-eta.new)/(1+abs(eta))<1e-7) break
            eta <- eta.new    
        }
    }
    sum(2*wt*(y[,1]*log(ifelse(y[,1]==0,1,y[,1]/(y[,1]+y[,2])/q))
              +y[,2]*log(y[,2]/(y[,1]+y[,2])/p)))
}


##%%%%%%%%%%  Proportional Odds Logistic Regression %%%%%%%%%%

## Make pseudo data for PO logistic regression
mkdata.polr <- function(y,eta,wt,offset,nu)
{
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    if (is.null(offset)) offset <- rep(0,dim(y)[1])
    nnu <- length(nu)
    hess <- matrix(0,nnu,nnu)
    G <- c(0,cumsum(nu))
    P <- exp(outer(eta,G,"+"))
    lkhd <- 0
    for (i in 1:(nnu+1))
        lkhd <- lkhd+sum(wt*(y[,i]+y[,i+1])*log(1+P[,i]))/sum(wt)
    for (i in 1:nnu) lkhd <- lkhd-sum(wt*y[,i+1])/sum(wt)*log(exp(nu[i])-1)
    if (nnu>1) {
        for (i in 1:(nnu-1)) {
            tmp <- 0
            for (j in (i+1):nnu) tmp <- tmp+sum(wt*y[,j+1])/sum(wt)
            lkhd <- lkhd-tmp*nu[i]
        }
    }
    dd <- log(nu)
    repeat {
        ## gradient and hessian
        nu <- exp(dd)
        G <- c(0,cumsum(nu))
        P <- exp(outer(eta,G,"+"))
        grad <- hess.wk <- NULL
        for (i in 1:nnu) {
            g.wk <- h.wk <- 0
            for (j in (i+1):(nnu+1)) {
                g.wk <- g.wk+sum(wt*(y[,j]+y[,j+1])*P[,j]/(1+P[,j]))/sum(wt)
                if (j<=nnu) g.wk <- g.wk-sum(wt*y[,j+1])/sum(wt)
                h.wk <- h.wk+sum(wt*(y[,j]+y[,j+1])*P[,j]/(1+P[,j])^2)/sum(wt)
            }
            g.wk <- g.wk-exp(nu[i])/(exp(nu[i])-1)*sum(wt*y[,i+1])/sum(wt)
            grad <- c(grad,g.wk)
            hess.wk <- c(hess.wk,h.wk)
        }
        for (i in 1:nnu) {
            hess[1:i,i] <- hess.wk[i]
            hess[i,i] <- hess[i,i]+exp(nu[i])/(exp(nu[i])-1)^2*sum(wt*y[,i+1])/sum(wt)
        }
        grad <- grad*nu
        hess <- hess*outer(nu,nu)
        diag(hess) <- diag(hess)+grad
        ## modify hessian if necessary
        if (nnu>1) {
            z <- .Fortran("dmcdc",
                          as.double(hess), as.integer(nnu), as.integer(nnu),
                          ee=double(nnu),
                          pivot=integer(nnu),
                          integer(1), PACKAGE="gss")
            if (max(z$ee)) {
                z$ee[z$pivot] <- z$ee
                hess <- hess+diag(z$ee)
            }
        }
        else hess <- abs(hess)
        ## update nu
        mumax <- max(abs(grad))
        dd.diff <- solve(hess,grad)
        repeat {
            lkhd.line <- function(x) {
                ddnew <- dd-c(x)*dd.diff
                nu <- exp(ddnew)
                G <- c(0,cumsum(nu))
                P <- exp(outer(eta,G,"+"))
                lkhd <- 0
                for (i in 1:(nnu+1))
                    lkhd <- lkhd+sum(wt*(y[,i]+y[,i+1])*log(1+P[,i]))/sum(wt)
                for (i in 1:nnu) lkhd <- lkhd-sum(wt*y[,i+1])/sum(wt)*log(exp(nu[i])-1)
                if (nnu>1) {
                    for (i in 1:(nnu-1)) {
                        tmp <- 0
                        for (j in (i+1):nnu) tmp <- tmp+sum(wt*y[,j+1])/sum(wt)
                        lkhd <- lkhd-tmp*nu[i]
                    }
                }
                lkhd
            }
            if (!is.finite(lkhdnew <- lkhd.line(1))) {
                dd.diff <- dd.diff/2
                next
            }
            ddnew <- dd-dd.diff
            if (lkhdnew-lkhd<(1+abs(lkhd))*10*.Machine$double.eps) break
            z <- nlm0(lkhd.line,c(0,1))
            ddnew <- dd-z$est*dd.diff
            lkhdnew <- z$min
            break
        }
        disc <- abs(lkhdnew-lkhd)/(1+abs(lkhd))
        disc <- max(disc,max(abs(dd-ddnew)/(1+abs(dd))))
        disc0 <- (mumax/(1+abs(lkhd)))^2
        dd <- ddnew
        lkhd <- lkhdnew
        if (min(disc,disc0)<1e-7) break
    }
    u <- -1+y[,nnu+2]
    for (i in 1:(nnu+1)) u <- u+(y[,i]+y[,i+1])*P[,i]/(1+P[,i])
    w <- P[,2]/(1+P[,2])*P[,1]/(1+P[,1])^2
    w <- w+1/(1+P[,nnu])*P[,nnu+1]/(1+P[,nnu+1])^2
    if (nnu>1) {
        for (i in 2:nnu)
            w <- w+(P[,i+1]-P[,i-1])/(1+P[,i+1])/(1+P[,i-1])*P[,i]/(1+P[,i])^2
    }
    ywk <- eta-u/w-offset
    wt <- w*wt
    list(ywk=ywk,wt=wt,nu=nu,u=u*wt)
}

## Calculate deviance residuals for PO logistic regression
dev.resid.polr <- function(y,eta,wt,nu)
{
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    nnu <- length(nu)
    G <- c(0,cumsum(nu))
    P <- plogis(outer(eta,G,"+"))
    wk <- NULL
    for (i in 1:length(wt)) {
        idx <- (1:(nnu+2))[y[i,]]
        if (idx==1) wk <- c(wk,P[i,1])
        if (idx==nnu+2) wk <- c(wk,1-P[i,nnu+1])
        if ((idx>1)&(idx<nnu+2)) wk <- c(wk,P[i,idx]-P[i,idx-1])
    }
    as.vector(-2*wt*log(wk))
}

## Calculate null deviance for PO logistic regression
dev.null.polr <- function(y,wt,offset)
{
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    nobs <- dim(y)[1]
    P <- apply(y*wt,2,sum)
    P <- P/sum(P)
    wk <- NULL
    for (i in 1:nobs) {
        idx <- (1:length(P))[y[i,]]
        wk <- c(wk,P[idx])
    }
    if (!is.null(offset)) {
        P <- cumsum(P)
        J <- length(P)
        eta0 <- qlogis(P[-J])-mean(offset)
        eta0[-1] <- log(diff(P))
        lkhd <- function(eta) {
            eta <- cumsum(c(eta[1],exp(eta[-1])))
            tmp <- 0
            for (i in 1:nobs) {
              idx <- (1:J)[y[i,]]
              if (idx==1) wk <- wk-wt[i]*log(plogis(eta[1]+offset[i]))
              if (idx==J) wk <- wk-wt[i]*log(1-plogis(eta[J-1]+offset[i]))
              if ((idx>1)&(idx<J))
                  wk <- wk-wt[i]*log(plogis(eta[idx]+offset[i])-plogis(eta[idx-1]+offset[i]))
            }
        }
        eta <- nlm(lkhd,eta0,stepmax=1)$est
        eta <- cumsum(c(eta[1],exp(eta[-1])))
        wk <- NULL
        for (i in 1:nobs) {
            idx <- (1:length(P))[y[i,]]
            if (idx==1) wk <- c(wk,plogis(eta[1]+offset[i]))
            if (idx==J) wk <- c(wk,1-plogis(eta[J-1]+offset[i]))
            if ((idx>1)&(idx<J))
                wk <- c(wk,log(plogis(eta[idx]+offset[i])-plogis(eta[idx-1]+offset[i])))
        }
    }
    sum(-2*wt*log(wk))
}

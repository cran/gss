##%%%%%%%%%%  Binomial Family %%%%%%%%%%
y0.binomial <- function(y,eta0,wt)
{
    if (is.matrix(y)) wt <- wt * (y[,1]+y[,2])
    p <- plogis(eta0)
    list(p=p,eta=eta0,wt=wt)
}
proj0.binomial <- function(y0,eta,offset)
{
    if (is.null(offset)) offset <- rep(0,length(eta))
    p <- plogis(eta)
    u <- p - y0$p
    w <- p*(1-p)
    ywk <- eta-u/w-offset
    wt <- w*y0$wt
    kl <- mean(y0$wt*(y0$p*(y0$eta-eta)+log((1-y0$p)/(1-p))))
    list(ywk=ywk,wt=wt,kl=kl,u=wt*u)
}
kl.binomial <- function(eta0,eta1,wt)
{
    p0 <- plogis(eta0)
    p1 <- plogis(eta1)
    mean(wt*(p0*(eta0-eta1)+log((1-p0)/(1-p1))))
}
cfit.binomial <- function(y,wt,offset)
{
    if (is.vector(y)) y <- as.matrix(y)
    if (dim(y)[2]>1) {
        wt <- wt * (y[,1]+y[,2])
        y <- y[,1]/(y[,1]+y[,2])
    }
    p <- sum(wt*y)/sum(wt)
    if (is.null(offset)) eta <- rep(qlogis(p),length(y))
    else {
        eta <- qlogis(p)-mean(offset)
        repeat {
            p <- plogis(eta+offset)
            u <- p - y
            w <- p*(1-p)
            eta.new <- eta-sum(wt*u)/sum(wt*w)
            if (abs(eta-eta.new)/(1+abs(eta))<1e-7) break
            eta <- eta.new    
        }
        eta <- eta + offset
    }
    eta
}


##%%%%%%%%%%  Poisson Family %%%%%%%%%%
y0.poisson <- function(eta0)
{
    lambda <- exp(eta0)
    list(lambda=lambda,eta=eta0)
}
proj0.poisson <- function(y0,eta,wt,offset)
{
    if (is.null(offset)) offset <- rep(0,length(eta))
    lambda <- exp(eta)
    u <- lambda - y0$lambda
    w <- lambda
    ywk <- eta-u/w-offset
    kl <- mean(wt*(y0$lambda*(y0$eta-eta)-y0$lambda+lambda))
    wt <- w*wt
    list(ywk=ywk,wt=wt,kl=kl,u=wt*u)
}
kl.poisson <- function(eta0,eta1,wt)
{
    lambda0 <- exp(eta0)
    lambda1 <- exp(eta1)
    mean(wt*(lambda0*(eta0-eta1)-lambda0+lambda1))
}
cfit.poisson <- function(y,wt,offset)
{
    lambda <- sum(wt*y)/sum(wt)
    if (is.null(offset)) eta <- rep(log(lambda),length(y))
    else {
        eta <- log(lambda) - mean(offset)
        repeat {
            lambda <- exp(eta+offset)
            u <- lambda - y
            w <- lambda
            eta.new <- eta-sum(wt*u)/sum(wt*w)
            if (abs(eta-eta.new)/(1+abs(eta))<1e-7) break
            eta <- eta.new    
        }
        eta <- eta + offset
    }
    eta
}


##%%%%%%%%%%  Gamma Family %%%%%%%%%%
y0.Gamma <- function(eta0)
{
    mu <- exp(eta0)
    list(mu=mu)
}
proj0.Gamma <- function(y0,eta,wt,offset)
{
    if (is.null(offset)) offset <- rep(0,length(y))
    mu <- exp(eta)
    u <- 1-y0$mu/mu
    w <- y0$mu/mu
    ywk <- eta-u/w-offset
    kl <- mean(wt*(y0$mu*(-1/y0$mu+1/mu)+log(mu/y0$mu)))
    wt <- w*wt
    list(ywk=ywk,wt=wt,kl=kl,u=wt*u)
}
kl.Gamma <- function(eta0,eta1,wt)
{
    mu0 <- exp(eta0)
    mu1 <- exp(eta1)
    mean(wt*(mu0*(-1/mu0+1/mu1)+log(mu1/mu0)))
}
cfit.Gamma <- function(y,wt,offset)
{
    mu <- sum(wt*y)/sum(wt)
    if (is.null(offset)) eta <- rep(log(mu),length(y))
    else {
        eta <- log(mu)-mean(offset)
        repeat {
            mu <- exp(eta+offset)
            u <- 1-y/mu
            w <- y/mu
            eta.new <- eta-sum(wt*u)/sum(wt*w)
            if (abs(eta-eta.new)/(1+abs(eta))<1e-7) break
            eta <- eta.new    
        }
        eta <- eta + offset
    }
    eta
}


##%%%%%%%%%%  Negative Binomial Family %%%%%%%%%%
y0.nbinomial <- function(y,eta0,nu)
{
    if (!is.vector(y)) {
        nu <- y[,2]
        y <- y[,1]
    }
    mu <- nu*exp(-eta0)
    list(y=y,nu=nu,mu=mu,eta=eta0)
}
proj0.nbinomial <- function(y0,eta,wt,offset)
{
    if (is.null(offset)) offset <- rep(0,length(eta))
    p <- plogis(eta)
    u <- (y0$mu+y0$nu)*p-y0$nu
    w <- (y0$mu+y0$nu)*p*(1-p)
    ywk <- eta-u/w-offset
    kl <- mean(wt*((y0$nu+y0$mu)*log((1+exp(eta))/(1+exp(y0$eta)))
                   +y0$nu*(y0$eta-eta)))
    wt <- w*wt
    list(ywk=ywk,wt=wt,kl=kl,u=wt*u)
}
kl.nbinomial <- function(eta0,eta1,wt,nu)
{
    mu0 <- nu*exp(-eta0)
    mean(wt*((nu+mu0)*log((1+exp(eta1))/(1+exp(eta0)))+nu*(eta0-eta1)))
}
cfit.nbinomial <- function(y,wt,offset,nu)
{
    if (!is.vector(y)) {
        nu <- y[,2]
        y <- y[,1]
    }
    p <- sum(wt*nu)/sum(wt*(y+nu))
    if (is.null(offset)) eta <- rep(qlogis(p),length(y))
    else {
        eta <- qlogis(p)-mean(offset)
        repeat {
            p <- 1-1/(1+exp(eta+offset))
            u <- (y+nu)*p-nu
            w <- (y+nu)*p*(1-p)
            eta.new <- eta-sum(wt*u)/sum(wt*w)
            if (abs(eta-eta.new)/(1+abs(eta))<1e-7) break
            eta <- eta.new    
        }
        eta <- eta + offset
    }
    eta
}


##%%%%%%%%%%  Weibull Family %%%%%%%%%%
y0.weibull <- function(y,eta0,nu)
{
    xx <- y[,1]
    if (dim(y)[2]>=3) zz <- y[,3]
    else zz <- rep(0,length(xx))
    lam <- exp(-nu*eta0)
    list(lam=lam,eta=eta0,int=(xx^nu-zz^nu))
}
proj0.weibull <- function(y0,eta,wt,offset,nu)
{
    if (is.null(offset)) offset <- rep(0,length(eta))
    u <- nu*(y0$lam-exp(-nu*eta))
    w <- nu*nu*exp(-nu*eta)
    ywk <- eta-u/w-offset
    kl <- mean(wt*y0$int*(y0$lam*nu*(eta-y0$eta)+exp(-nu*eta)-y0$lam))
    u <- y0$int*u
    w <- y0$int*w
    wt <- w*wt
    list(ywk=ywk,wt=wt,kl=kl,u=wt*u)
}
kl.weibull <- function(eta0,eta1,wt,nu,int)
{
    lam0 <- exp(-nu*eta0)
    lam1 <- exp(-nu*eta1)
    mean(wt*int*(lam0*nu*(eta1-eta0)+lam1-lam0))
}
cfit.weibull <- function(y,wt,offset,nu)
{
    xx <- y[,1]
    delta <- as.logical(y[,2])
    if (dim(y)[2]>=3) zz <- y[,3]
    else zz <- rep(0,length(xx))
    if (is.null(offset)) offset <- rep(0,length(xx))
    eta <- log(sum(wt*(xx^nu-zz^nu)*exp(-nu*offset))/sum(wt*delta))/nu
    eta + offset    
}


##%%%%%%%%%%  Lognorm Family %%%%%%%%%%
y0.lognorm <- function(y,eta0,nu)
{
    xx <- y[,1]
    if (dim(y)[2]>=3) zz <- y[,3]
    else zz <- rep(0,length(xx))
    quad <- gauss.quad(50,c(0,1))
    list(eta=eta0,xx=xx,zz=zz,q.pt=quad$pt,q.wt=quad$wt)
}
proj0.lognorm <- function(y0,eta,wt,offset,nu)
{
    if (is.null(offset)) offset <- rep(0,length(eta))
    u <- NULL
    kl <- 0
    for (i in 1:length(eta)) {
        q.pt <- y0$q.pt*(y0$xx[i]-y0$zz[i])+y0$zz[i]
        q.wt <- y0$q.wt*(y0$xx[i]-y0$zz[i])
        z0 <- nu*(log(q.pt)-y0$eta[i])
        z1 <- nu*(log(q.pt)-eta[i])
        lam0 <- ifelse(z0<7,dnorm(z0)/(1-pnorm(z0)),z0+.15)
        lam1 <- ifelse(z1<7,dnorm(z1)/(1-pnorm(z1)),z1+.15)
        u <- c(u,nu*nu*sum(q.wt*(lam0-lam1)*(lam1-z1)/q.pt))
        kl <- kl + nu*sum(q.wt*(lam0*log(lam0/lam1)+lam1-lam0)/q.pt)
    }
    xx <- nu*(log(y0$xx)-eta)
    zz <- nu*(log(y0$zz)-eta)
    s.xx <- ifelse(xx<7,dnorm(xx)/(1-pnorm(xx)),xx+.15)
    s.zz <- ifelse(zz<7,dnorm(zz)/(1-pnorm(zz)),zz+.15)
    s.xx <- pmax(s.xx,s.zz)
    w <- (s.xx^2/2-xx*s.xx+xx^2/2+log(s.xx)+log(2*pi)/2)
    w <- nu^2*(w-ifelse(s.zz==0,0,(s.zz^2/2-zz*s.zz+zz^2/2+log(s.zz)+log(2*pi)/2)))
    w <- pmax(w,1e-6)
    ywk <- eta-u/w-offset
    wt <- w*wt
    list(ywk=ywk,wt=wt,kl=kl/length(eta),u=wt*u)
}
kl.lognorm <- function(eta0,eta1,wt,nu,y0)
{
    kl <- 0
    for (i in 1:length(eta0)) {
        q.pt <- y0$q.pt*(y0$xx[i]-y0$zz[i])+y0$zz[i]
        q.wt <- y0$q.wt*(y0$xx[i]-y0$zz[i])
        z0 <- nu*(log(q.pt)-eta0[i])
        z1 <- nu*(log(q.pt)-eta1[i])
        lam0 <- ifelse(z0<7,dnorm(z0)/(1-pnorm(z0)),z0+.15)
        lam1 <- ifelse(z1<7,dnorm(z1)/(1-pnorm(z1)),z1+.15)
        kl <- kl + nu*sum(q.wt*(lam0*log(lam0/lam1)+lam1-lam0)/q.pt)
    }
    kl/length(eta0)
}
cfit.lognorm <- function(y,wt,offset,nu)
{
    xx <- y[,1]
    delta <- as.logical(y[,2])
    if (dim(y)[2]>=3) zz <- y[,3]
    else zz <- rep(0,length(xx))
    if (is.null(offset)) offset <- rep(0,length(xx))
    lkhd <- function(eta) {
        eta <- eta + offset
        xx.wk <- nu*(log(xx)-eta)
        zz.wk <- nu*(log(zz)-eta)
        -sum(wt*(delta*(-xx.wk^2/2-log(1-pnorm(xx.wk)))
                 +log((1-pnorm(xx.wk))/(1-pnorm(zz.wk)))))
    }
    nlm(lkhd,mean(log(xx)-offset),stepmax=1)$est + offset
}


##%%%%%%%%%%  Loglogis Family %%%%%%%%%%
y0.loglogis <- function(y,eta0,nu)
{
    xx <- y[,1]
    if (dim(y)[2]>=3) zz <- y[,3]
    else zz <- rep(0,length(xx))
    quad <- gauss.quad(50,c(0,1))
    list(eta=eta0,xx=xx,zz=zz,q.pt=quad$pt,q.wt=quad$wt)
}
proj0.loglogis <- function(y0,eta,wt,offset,nu)
{
    if (is.null(offset)) offset <- rep(0,length(eta))
    e0 <- exp(-nu*y0$eta)
    e1 <- exp(-nu*eta)
    kl <- sum(log((1+y0$xx^nu*e1)*(1+y0$zz^nu*e0)/(1+y0$zz^nu*e1)/(1+y0$xx^nu*e0))
              +nu*(eta-y0$eta)*log((1+y0$xx^nu*e0)/(1+y0$zz^nu*e0)))
    xx <- 1/(1+y0$xx^nu*e1)
    zz <- 1/(1+y0$zz^nu*e1)
    u <- -nu*(zz-xx)
    w <- nu^2/2*(zz^2-xx^2)
    for (i in 1:length(eta)) {
        q.pt <- y0$q.pt*(y0$xx[i]-y0$zz[i])+y0$zz[i]
        q.wt <- y0$q.wt*(y0$xx[i]-y0$zz[i])
        u[i] <- u[i]+nu^2*sum(q.wt*q.pt^(nu-1)*e0[i]
                              /(1+q.pt^nu*e0[i])/(1+q.pt^nu*e1[i]))
        kl <- kl + nu*sum(q.wt*q.pt^(nu-1)*e0[i]/(1+q.pt^nu*e0[i])
                          *log((1+q.pt^nu*e1[i])/(1+q.pt^nu*e0[i])))
    }
    w <- pmax(w,1e-6)
    ywk <- eta-u/w-offset
    wt <- w*wt
    list(ywk=ywk,wt=wt,kl=kl/length(eta),u=wt*u)
}
kl.loglogis <- function(eta0,eta1,wt,nu,y0)
{
    e0 <- exp(-nu*eta0)
    e1 <- exp(-nu*eta1)
    kl <- sum(log((1+y0$xx^nu*e1)*(1+y0$zz^nu*e0)/(1+y0$zz^nu*e1)/(1+y0$xx^nu*e0))
              +nu*(eta1-eta0)*log((1+y0$xx^nu*e0)/(1+y0$zz^nu*e0)))
    for (i in 1:length(eta0)) {
        q.pt <- y0$q.pt*(y0$xx[i]-y0$zz[i])+y0$zz[i]
        q.wt <- y0$q.wt*(y0$xx[i]-y0$zz[i])
        kl <- kl + nu*sum(q.wt*q.pt^(nu-1)*e0[i]/(1+q.pt^nu*e0[i])
                          *log((1+q.pt^nu*e1[i])/(1+q.pt^nu*e0[i])))
    }
    kl/length(eta0)
}
cfit.loglogis <- function(y,wt,offset,nu)
{
    xx <- y[,1]
    delta <- as.logical(y[,2])
    if (dim(y)[2]>=3) zz <- y[,3]
    else zz <- rep(0,length(xx))
    if (is.null(offset)) offset <- rep(0,length(xx))
    lkhd <- function(eta) {
        eta <- eta + offset
        xx.wk <- nu*(log(xx)-eta)
        zz.wk <- nu*(log(zz)-eta)
        -sum(wt*(delta*(xx.wk-log(1+exp(xx.wk)))
                 -log((1+exp(xx.wk))/(1+exp(zz.wk)))))
    }
    nlm(lkhd,mean(log(xx)-offset),stepmax=1)$est + offset
}

##%%%%%%%%%%  Weibull Family %%%%%%%%%%

## Make pseudo data for Weibull regression
mkdata.weibull <- function(y,eta,wt,offset,alpha)
{
    if (is.vector(y)) stop("gss error: missing censoring indicator")
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    if (is.null(offset)) offset <- rep(0,dim(y)[1])
    xx <- y[,1]
    delta <- as.logical(y[,2])
    if (dim(y)[2]>=3) zz <- y[,3]
    else zz <- rep(0,length(xx))
    if (any(zz<0)|any(zz>=xx))
        stop("gss error: inconsistent life time data")
    if (alpha[[2]]) {
        lkhd <- function(log.alpha) {
            alpha <- exp(log.alpha)
            -sum(wt*(delta*(alpha*(log(xx)-eta)+log.alpha)
                 -(xx^alpha-zz^alpha)*exp(-alpha*eta)))
        }
        if (is.null(alpha[[1]])) alpha[[1]] <- 1
        alpha[[1]] <- exp(nlm(lkhd,log(alpha[[1]]),stepmax=.5)$est)
    }
    u <- alpha[[1]]*(delta-(xx^alpha[[1]]-zz^alpha[[1]])*exp(-alpha[[1]]*eta))
    w <- alpha[[1]]^2*(xx^alpha[[1]]-zz^alpha[[1]])*exp(-alpha[[1]]*eta)
    ywk <- eta-u/w-offset
    wt <- w*wt
    list(ywk=ywk,wt=wt,alpha=alpha)
}

## Calculate deviance residuals for Weibull regression
dev.resid.weibull <- function(y,eta,wt,alpha)
{
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    xx <- y[,1]
    delta <- as.logical(y[,2])
    if (dim(y)[2]>=3) zz <- y[,3]
    else zz <- rep(0,length(xx))
    z <- -delta*alpha*(log(xx)-eta)+(xx^alpha-zz^alpha)*exp(-alpha*eta)
    as.numeric(2*wt*(z+delta*(log(xx^alpha)-log(xx^alpha-zz^alpha)-1)))
}

## Calculate null deviance for Weibull regression
dev.null.weibull <- function(y,wt,offset,alpha)
{
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    xx <- y[,1]
    delta <- as.logical(y[,2])
    if (dim(y)[2]>=3) zz <- y[,3]
    else zz <- rep(0,length(xx))
    if (is.null(offset)) offset <- rep(0,length(xx))
    eta <- log(sum(wt*(xx^alpha-zz^alpha)*exp(-alpha*offset))/sum(wt*delta))/alpha
    eta <- eta + offset    
    z <- -delta*alpha*(log(xx)-eta)+(xx^alpha-zz^alpha)*exp(-alpha*eta)
    sum(2*wt*(z+delta*(log(xx^alpha)-log(xx^alpha-zz^alpha)-1)))
}


##%%%%%%%%%%  Log Normal Family %%%%%%%%%%

## Make pseudo data for log normal regression
mkdata.lognorm <- function(y,eta,wt,offset,alpha)
{
    if (is.vector(y)) stop("gss error: missing censoring indicator")
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    if (is.null(offset)) offset <- rep(0,dim(y)[1])
    xx <- y[,1]
    delta <- as.logical(y[,2])
    if (dim(y)[2]>=3) zz <- y[,3]
    else zz <- rep(0,length(xx))
    if (any(zz<0)|any(zz>=xx))
        stop("gss error: inconsistent life time data")
    if (alpha[[2]]) {
        lkhd <- function(log.alpha) {
            alpha <- exp(log.alpha)
            xx.wk <- alpha*(log(xx)-eta)
            zz.wk <- alpha*(log(zz)-eta)
            -sum(wt*(delta*(-xx.wk^2/2-log(1-pnorm(xx.wk))+log.alpha)
                 +log((1-pnorm(xx.wk))/(1-pnorm(zz.wk)))))
        }
        if (is.null(alpha[[1]])) alpha[[1]] <- 1
        alpha[[1]] <- exp(nlm(lkhd,log(alpha[[1]]),stepmax=.5)$est)
    }
    xx <- alpha[[1]]*(log(xx)-eta)
    zz <- alpha[[1]]*(log(zz)-eta)
    s.xx <- dnorm(xx)/(1-pnorm(xx))
    s.zz <- dnorm(zz)/(1-pnorm(zz))
    u <- alpha[[1]]*(delta*(s.xx-xx)-(s.xx-s.zz))
    w <- (s.xx^2/2-xx*s.xx-log(1-pnorm(xx)))
    w <- alpha[[1]]^2*(w-ifelse(s.zz==0,0,(s.zz^2/2-zz*s.zz-log(1-pnorm(zz)))))
    ywk <- eta-u/w-offset
    wt <- w*wt
    list(ywk=ywk,wt=wt,alpha=alpha)
}

## Calculate deviance residuals for log normal regression
dev.resid.lognorm <- function(y,eta,wt,alpha)
{
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    xx <- y[,1]
    delta <- as.logical(y[,2])
    if (dim(y)[2]>=3) zz <- y[,3]
    else zz <- rep(0,length(xx))
    xx <- alpha[[1]]*(log(xx)-eta)
    zz <- alpha[[1]]*(log(zz)-eta)
    s.zz <- dnorm(zz)/(1-pnorm(zz))
    z <- -delta*(-xx^2/2-log(1-pnorm(xx)))-log((1-pnorm(xx))/(1-pnorm(zz)))
    as.numeric(2*wt*(z+delta*(-s.zz^2/2-log(1-pnorm(zz)))))
}

## Calculate null deviance for log normal regression
dev.null.lognorm <- function(y,wt,offset,alpha)
{
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    xx <- y[,1]
    delta <- as.logical(y[,2])
    if (dim(y)[2]>=3) zz <- y[,3]
    else zz <- rep(0,length(xx))
    if (is.null(offset)) offset <- rep(0,length(xx))
    lkhd <- function(eta) {
        eta <- eta + offset
        xx.wk <- alpha*(log(xx)-eta)
        zz.wk <- alpha*(log(zz)-eta)
        -sum(wt*(delta*(-xx.wk^2/2-log(1-pnorm(xx.wk)))
                 +log((1-pnorm(xx.wk))/(1-pnorm(zz.wk)))))
    }
    eta <- nlm(lkhd,mean(log(xx)-offset))$est + offset
    xx <- alpha[[1]]*(log(xx)-eta)
    zz <- alpha[[1]]*(log(zz)-eta)
    s.zz <- dnorm(zz)/(1-pnorm(zz))
    z <- -delta*(-xx^2/2-log(1-pnorm(xx)))-log((1-pnorm(xx))/(1-pnorm(zz)))
    sum(2*wt*(z+delta*(-s.zz^2/2-log(1-pnorm(zz)))))
}


##%%%%%%%%%%  Log Logistic Family %%%%%%%%%%

## Make pseudo data for log logistic regression
mkdata.loglogis <- function(y,eta,wt,offset,alpha)
{
    if (is.vector(y)) stop("gss error: missing censoring indicator")
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    if (is.null(offset)) offset <- rep(0,dim(y)[1])
    xx <- y[,1]
    delta <- as.logical(y[,2])
    if (dim(y)[2]>=3) zz <- y[,3]
    else zz <- rep(0,length(xx))
    if (any(zz<0)|any(zz>=xx))
        stop("gss error: inconsistent life time data")
    if (alpha[[2]]) {
        lkhd <- function(log.alpha) {
            alpha <- exp(log.alpha)
            xx.wk <- alpha*(log(xx)-eta)
            zz.wk <- alpha*(log(zz)-eta)
            -sum(wt*(delta*(xx.wk-log(1+exp(xx.wk))+log.alpha)
                 -log((1+exp(xx.wk))/(1+exp(zz.wk)))))
        }
        if (is.null(alpha[[1]])) alpha[[1]] <- 1
        alpha[[1]] <- exp(nlm(lkhd,log(alpha[[1]]),stepmax=.5)$est)
    }
    xx <- 1/(1+exp(alpha[[1]]*(log(xx)-eta)))
    zz <- 1/(1+exp(alpha[[1]]*(log(zz)-eta)))
    u <- alpha[[1]]*(delta*xx-(zz-xx))
    w <- alpha[[1]]^2/2*(zz^2-xx^2)
    ywk <- eta-u/w-offset
    wt <- w*wt
    list(ywk=ywk,wt=wt,alpha=alpha)
}

## Calculate deviance residuals for log logistic regression
dev.resid.loglogis <- function(y,eta,wt,alpha)
{
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    xx <- y[,1]
    delta <- as.logical(y[,2])
    if (dim(y)[2]>=3) zz <- y[,3]
    else zz <- rep(0,length(xx))
    xx <- alpha*(log(xx)-eta)
    zz <- alpha*(log(zz)-eta)
    z <- -delta*(xx-log(1+exp(xx)))+log((1+exp(xx))/(1+exp(zz)))
    as.numeric(2*wt*(z+delta*(log(1+2*exp(zz))-log(1+exp(zz))-2*log(2))))
}

## Calculate null deviance for log logistic regression
dev.null.loglogis <- function(y,wt,offset,alpha)
{
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    xx <- y[,1]
    delta <- as.logical(y[,2])
    if (dim(y)[2]>=3) zz <- y[,3]
    else zz <- rep(0,length(xx))
    if (is.null(offset)) offset <- rep(0,length(xx))
    lkhd <- function(eta) {
        eta <- eta + offset
        xx.wk <- alpha*(log(xx)-eta)
        zz.wk <- alpha*(log(zz)-eta)
        -sum(wt*(delta*(xx.wk-log(1+exp(xx.wk)))
                 -log((1+exp(xx.wk))/(1+exp(zz.wk)))))
    }
    eta <- nlm(lkhd,mean(log(xx)-offset))$est + offset
    xx <- alpha*(log(xx)-eta)
    zz <- alpha*(log(zz)-eta)
    z <- -delta*(xx-log(1+exp(xx)))+log((1+exp(xx))/(1+exp(zz)))
    sum(2*wt*(z+delta*(log(1+2*exp(zz))-log(1+exp(zz))-2*log(2))))
}

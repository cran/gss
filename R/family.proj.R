##%%%%%%%%%%  Binomial Family %%%%%%%%%%

## Make pseudo data for logistic regression
proj0.binomial <- function(y,eta,wt,offset)
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
    p <- 1-1/(1+exp(eta))
    u <- p - y
    w <- p*(1-p)
    ywk <- eta-u/w-offset
    wt <- w*wt
    list(ywk=ywk,wt=wt,mu=p,theta=eta,b=-log(1-p),u=wt*u)
}


##%%%%%%%%%%  Poisson Family %%%%%%%%%%

## Make pseudo data for Poisson regression
proj0.poisson <- function(y,eta,wt,offset)
{
    if (is.null(wt)) wt <- rep(1,length(y))
    if (is.null(offset)) offset <- rep(0,length(y))
    if (min(y)<0)
        stop("gss error: paired binomial response should be nonnegative")
    lambda <- exp(eta)
    u <- lambda - y
    w <- lambda
    ywk <- eta-u/w-offset
    wt <- w*wt
    list(ywk=ywk,wt=wt,mu=lambda,theta=eta,b=lambda,u=wt*u)
}


##%%%%%%%%%%  Gamma Family %%%%%%%%%%

## Make pseudo data for Gamma regression
proj0.Gamma <- function(y,eta,wt,offset)
{
    if (is.null(wt)) wt <- rep(1,length(y))
    if (is.null(offset)) offset <- rep(0,length(y))
    if (min(y)<=0)
        stop("gss error: gamma responses should be positive")
    mu <- exp(eta)
    u <- 1-y/mu
    w <- y/mu
    ywk <- eta-u/w-offset
    wt <- w*wt
    list(ywk=ywk,wt=wt,mu=mu,theta=-1/mu,b=log(mu),u=wt*u)
}

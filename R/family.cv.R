##%%%%%%%%%%  Binomial Family %%%%%%%%%%

## Calculate CV score for binomial regression
cv.binomial <- function(y,eta,wt,hat,alpha)
{
    if (is.vector(y)) y <- as.matrix(y)
    if (is.null(wt)) wt <- rep(1,dim(y)[1])
    if (dim(y)[2]==1) {
        if ((max(y)>1)|(min(y)<0))
            stop("gss error: binomial responses should be between 0 and 1")
        m <- rep(1,dim(y)[1])
    }
    else {
        if (min(y)<0)
            stop("gss error: paired binomial response should be nonnegative")
        m <- y[,1]+y[,2]
        y <- y[,1]/m
    }
    wtt <- wt * m
    p <- 1-1/(1+exp(eta))
    w <- p*(1-p)
    lkhd <- -sum(wtt*(y*eta+log(1-p)))/sum(wtt)
    aux1 <- sum(hat/w)/(sum(wtt)-sum(hat))
    aux2 <- sum(wtt*y*(1-p))/sum(wtt)
alpha <- 1    
    list(score=lkhd+abs(alpha)*aux1*aux2,varht=1,w=as.vector(wtt*w))
}


##%%%%%%%%%%  Poisson Family %%%%%%%%%%

## Calculate CV score for Poisson regression
cv.poisson <- function(y,eta,wt,hat,alpha,sr,q)
{
    if (is.null(wt)) wt <- rep(1,length(y))
    if (min(y)<0)
        stop("gss error: Poisson response should be nonnegative")
    nxi <- ncol(q)
    nn <- ncol(sr)
    nnull <- nn-nxi
    lambda <- exp(eta)
    w <- as.vector(lambda)
    lkhd <- -sum(wt*(y*eta-lambda))/sum(wt*y)
    ## matrix H
    mu <- apply(wt*w*sr,2,sum)/sum(wt*w)
    v <- t(sr)%*%(wt*w*sr)/sum(wt*w)-outer(mu,mu)
    v[(nnull+1):nn,(nnull+1):nn] <- v[(nnull+1):nn,(nnull+1):nn]+q/sum(wt*y)
    ## Cholesky decomposition of H
    z <- .Fortran("dchdc",
                  v=as.double(v), as.integer(nn), as.integer(nn),
                  double(nn), jpvt=as.integer(rep(0,nn)),
                  as.integer(1), rkv=integer(1),
                  PACKAGE="base")[c("v","jpvt","rkv")]
    v <- matrix(z$v,nn,nn)
    rkv <- z$rkv
    while (v[rkv,rkv]<v[1,1]*sqrt(.Machine$double.eps)) rkv <- rkv-1
    if (rkv<nn) v[(rkv+1):nn,(rkv+1):nn] <- diag(v[1,1],nn-rkv)
    ## trace
    mu <- apply(wt*y*sr,2,sum)/sum(wt*y)
    sr <- sqrt(wt*y)*t(t(sr)-mu)
    sr <- backsolve(v,t(sr[,z$jpvt]),tran=TRUE)
    aux1 <- sum(sr^2)
    aux2 <- 1/sum(wt*y)/(sum(wt*y)-1)
    list(score=lkhd+abs(alpha)*aux1*aux2,varht=1,w=as.vector(wt*w))
}


##%%%%%%%%%%  Gamma Family %%%%%%%%%%

## Calculate CV score for Gamma regression
cv.Gamma <- function(y,eta,wt,hat,rss,alpha)
{
    if (is.null(wt)) wt <- rep(1,length(y))
    if (min(y)<=0)
        stop("gss error: gamma responses should be positive")
    mu <- exp(eta)
    u <- 1-y/mu
    w <- y/mu
    lkhd <- sum(wt*(y/mu+eta))/sum(wt)
    aux1 <- sum(hat/w)/(sum(wt)-sum(hat))
    aux2 <- -sum(wt*y*u/mu)/sum(wt)
    list(score=lkhd+alpha*aux1*aux2,varht=rss/(1-mean(hat)),w=as.vector(wt*w))
}

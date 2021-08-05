## Composition estimation with one sample
sscomp <- function(x,wt=rep(1,length(x)),alpha=1.4)
{
    ## Check inputs
    if ((nlvl <- length(x))<3)
        stop("gss error in sscomp: length of x should be 3 or more")
    if (length(x)!=length(wt))
        stop("gss error in sscomp: x and wt mismatch in lengths")
    ## Generate terms
    cnt <- x
    x <- as.factor(1:nlvl)
    mf <- model.frame(~x)
    term <- mkterm(mf,NULL)
    rk <- term$x$rk
    ## get basis functions
    id.basis <- 1:nlvl
    if (max(abs(wt-mean(wt)))/mean(wt)<.Machine$double.eps)
        id.basis <- id.basis[cnt>0]
    if (length(id.basis)==nlvl) id.basis <- id.basis[-nlvl]
    ## generate matrices
    r <- rk$fun(x[id.basis],x,nu=1,env=rk$env,out=TRUE)
    q <- r[,id.basis]
    qd.wt <- as.vector(wt)
    ## Fit the model
    nt <- b.wt <- 1
    t.wt <- matrix(1,nlvl,1)
    bias0 <- list(nt=nt,wt=b.wt,qd.wt=t.wt)
    z <- sspdsty(NULL,r,q,cnt,NULL,r,qd.wt,1e-7,30,alpha,bias0)
    ## return fitted probabilities
    fit <- exp(t(r)%*%z$c)*qd.wt
    rownames(fit) <- rownames(x)
    fit/sum(fit)
}
## Composition estimation with a matrix input
sscomp2 <- function(x,alpha=1.4)
{
    if (!is.matrix(x))
        stop("gss error in sscomp2: x should be a matrix")
    if (min(x)<0)
        stop("gss error in sscomp2: x should have non-negative entries")
    if (any(apply(x,2,sum)==0))
        stop("gss error in sscomp2: column totals of x must be positive")
    nlvl <- dim(x)[1]
    yy <- apply(x,1,sum)
    p0 <- sscomp(yy)
    fit <- NULL
    for (i in 1:dim(x)[2]) {
        fit <- cbind(fit,sscomp(x[,i],p0,alpha))
    }
    rownames(fit) <- rownames(x)
    colnames(fit) <- colnames(x)
    fit
}

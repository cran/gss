## Fit gssanova model
gssanova1 <- function(formula,family,type="cubic",data=list(),
                      weights,subset,offset,na.action=na.omit,
                      partial=NULL,alpha=NULL,nu=NULL,
                      id.basis=NULL,nbasis=NULL,seed=NULL,random=NULL,
                      ext=.05,order=2)
{
    if (!(family%in%c("binomial","poisson","Gamma")))
        stop("gss error in gssanova1: family not implemented")
    if (is.null(alpha)) {
        alpha <- switch(family,
                        binomial=1,
                        poisson=1.4,
                        Gamma=1.4)
    }
    ## Obtain model frame and model terms
    mf <- match.call()
    mf$family <- mf$type <- mf$partial <- NULL
    mf$method <- mf$varht <- mf$nu <- NULL
    mf$alpha <- mf$id.basis <- mf$nbasis <- mf$seed <- NULL
    mf$random <- mf$ext <- mf$order <- NULL
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf,sys.frame(sys.parent()))
    wt <- model.weights(mf)
    ## Generate sub-basis
    nobs <- dim(mf)[1]
    if (is.null(id.basis)) {
        if (is.null(nbasis))  nbasis <- max(30,ceiling(10*nobs^(2/9)))
        if (nbasis>=nobs)  nbasis <- nobs
        if (!is.null(seed))  set.seed(seed)
        id.basis <- sample(nobs,nbasis,prob=wt)
    }
    else {
        if (max(id.basis)>nobs|min(id.basis)<1)
            stop("gss error in gssanova1: id.basis out of range")
        nbasis <- length(id.basis)
    }
    ## Generate terms
    if (type=="cubic") term <- mkterm.cubic(mf,ext)
    if (type=="linear") term <- mkterm.linear(mf,ext)
    if (type=="tp") term <- mkterm.tp(mf,order,mf[id.basis,],1)
    if (is.null(term)) stop("gss error in gssanova1: unknown type")
    ## Generate random
    if (!is.null(random)) {
        if (class(random)=="formula") random <- mkran(random,data)
    }
    ## Generate s, r, q, and y
    s <- r <- NULL
    nq <- 0
    for (label in term$labels) {
        if (label=="1") {
            s <- cbind(s,rep(1,len=nobs))
            next
        }
        x <- mf[,term[[label]]$vlist]
        x.basis <- mf[id.basis,term[[label]]$vlist]
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
                r <- array(c(r,rk$fun(x,x.basis,nu=i,env=rk$env,out=TRUE)),c(nobs,nbasis,nq))
            }
        }
    }
    if (is.null(r))
        stop("gss error in gssanova1: use glm for models with only fixed effects")
    else q <- r[id.basis,,,drop=FALSE]
    ## Add the partial term
    if (!is.null(partial)) {
        if (is.vector(partial)) partial <- as.matrix(partial)
        if (dim(partial)[1]!=dim(mf)[1])
            stop("gss error in gssanova1: partial data are of wrong size")
        term$labels <- c(term$labels,"partial")
        term$partial <- list(nphi=dim(partial)[2],nrk=0,
                             iphi=ifelse(is.null(s),0,dim(s)[2])+1)
        s <- cbind(s,partial)
        mf$partial <- partial
    }
    if (qr(s)$rank<dim(s)[2])
        stop("gss error in gssanova1: fixed effects are linearly dependent")
    ## Prepare the data
    y <- model.response(mf,"numeric")
    offset <- model.offset(mf)
    if (!is.null(offset)) {
        term$labels <- c(term$labels,"offset")
        term$offset <- list(nphi=0,nrk=0)
    }
    ## Fit the model
    if (nq==1) {
        r <- r[,,1]
        q <- q[,,1]
        z <- sspngreg1(family,s,r,q,y,wt,offset,alpha,nu,random)
    }
    else z <- mspngreg1(family,s,r,q,y,wt,offset,alpha,nu,random)
    ## Brief description of model terms
    desc <- NULL
    for (label in term$labels)
        desc <- rbind(desc,as.numeric(c(term[[label]][c("nphi","nrk")])))
    desc <- rbind(desc,apply(desc,2,sum))
    rownames(desc) <- c(term$labels,"total")
    colnames(desc) <- c("Unpenalized","Penalized")
    ## Return the results
    obj <- c(list(call=match.call(),family=family,mf=mf,terms=term,desc=desc,
                  alpha=alpha,id.basis=id.basis,random=random),z)
    class(obj) <- c("gssanova1","ssanova1","gssanova","ssanova")
    obj
}

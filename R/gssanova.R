## Fit gssanova model
gssanova <- function(formula,family,type="cubic",data=list(),
                     weights,subset,offset,na.action=na.omit,
                     partial=NULL,method=NULL,varht=1,alpha=NULL,
                     prec=1e-7,maxiter=30,ext=.05,order=2)
{
    ## Obtain model frame and model terms
    mf <- match.call()
    mf$family <- mf$type <- mf$partial <- NULL
    mf$method <- mf$varht <- mf$alpha <- NULL
    mf$prec <- mf$maxiter <- mf$ext <- mf$order <- NULL
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf,sys.frame(sys.parent()))
    if (type=="cubic") term <- mkterm.cubic(mf,ext)
    if (type=="linear") term <- mkterm.linear(mf,ext)
    if (type=="tp") term <- mkterm.tp(mf,order,mf,1)
    ## Specify default method
    if (is.null(method)) {
        method <- switch(family,
                         binomial="u",
                         nbinomial="u",
                         poisson="u",
                         inverse.gaussian="v",
                         Gamma="v",
                         weibull="u",
                         lognorm="u",
                         loglogis="u")
    }
    ## Generate s, q, and y
    nobs <- dim(mf)[1]
    s <- q <- NULL
    nq <- 0
    for (label in term$labels) {
        if (label=="1") {
            s <- cbind(s,rep(1,len=nobs))
            next
        }
        x <- mf[,term[[label]]$vlist]
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
                q <- array(c(q,rk$fun(x,x,nu=i,env=rk$env,out=TRUE)),c(nobs,nobs,nq))
            }
        }
    }
    ## Add the partial term
    if (!is.null(partial)) {
        if (is.vector(partial)) partial <- as.matrix(partial)
        if (dim(partial)[1]!=dim(mf)[1])
            stop("gss error: partial data are of wrong size")
        term$labels <- c(term$labels,"partial")
        term$partial <- list(nphi=dim(partial)[2],nrk=0,
                             iphi=ifelse(is.null(s),0,dim(s)[2])+1)
        s <- cbind(s,partial)
        mf$partial <- partial
    }
    if (qr(s)$rank<dim(s)[2])
        stop("gss error: fixed effects are linearly dependent")
    y <- model.response(mf,"numeric")
    wt <- model.weights(mf)
    offset <- model.offset(mf)
    if (!is.null(offset)) {
        term$labels <- c(term$labels,"offset")
        term$offset <- list(nphi=0,nrk=0)
    }
    if (!nq) stop("use glm for models with only fixed effects")
    ## Fit the model
    if (nq==1) {
        q <- q[,,1]
        z <- sspregpoi(family,s,q,y,wt,offset,method,varht,alpha,prec,maxiter)
    }
    else z <- mspregpoi(family,s,q,y,wt,offset,method,varht,alpha,prec,maxiter)
    ## Brief description of model terms
    desc <- NULL
    for (label in term$labels)
        desc <- rbind(desc,as.numeric(c(term[[label]][c("nphi","nrk")])))
    desc <- rbind(desc,apply(desc,2,sum))
    rownames(desc) <- c(term$labels,"total")
    colnames(desc) <- c("Fixed","Random")
    ## Return the results
    obj <- c(list(call=match.call(),family=family,mf=mf,terms=term,desc=desc),z)
    class(obj) <- c("gssanova","ssanova")
    obj
}

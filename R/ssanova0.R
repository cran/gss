## Fit ssanova model
ssanova0 <- function(formula,type=NULL,data=list(),weights,subset,
                     offset,na.action=na.omit,partial=NULL,
                     method="v",varht=1,prec=1e-7,maxiter=30)
{
    ## Obtain model frame and model terms
    mf <- match.call()
    mf$type <- mf$method <- mf$varht <- mf$partial <- NULL
    mf$prec <- mf$maxiter <- NULL
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf,parent.frame())
    ## Generate terms
    term <- mkterm(mf,type)
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
            stop("gss error in ssanova0: partial data are of wrong size")
        term$labels <- c(term$labels,"partial")
        term$partial <- list(nphi=dim(partial)[2],nrk=0,
                             iphi=ifelse(is.null(s),0,dim(s)[2])+1)
        s <- cbind(s,partial)
        mf$partial <- partial
    }
    ## Prepare the data
    y <- model.response(mf,"numeric")
    w <- model.weights(mf)
    offset <- model.offset(mf)
    if (!is.null(offset)) {
        term$labels <- c(term$labels,"offset")
        term$offset <- list(nphi=0,nrk=0)
        y <- y - offset
    }
    if (!is.null(w)) {
        w <- sqrt(w)
        y <- w*y
        s <- w*s
        for (i in 1:nq) q[,,i] <- w*t(w*q[,,i])
    }
    if (qr(s)$rank<dim(s)[2])
        stop("gss error in ssanova0: fixed effects are linearly dependent")
    if (!nq) stop("gss error in ssanova0: use lm for models with only fixed effects")
    ## Fit the model
    if (nq==1) {
        q <- q[,,1]
        z <- sspreg0(s,q,y,method,varht)
    }
    else z <- mspreg0(s,q,y,method,varht,prec,maxiter)
    ## Brief description of model terms
    desc <- NULL
    for (label in term$labels)
        desc <- rbind(desc,as.numeric(c(term[[label]][c("nphi","nrk")])))
    desc <- rbind(desc,apply(desc,2,sum))
    rownames(desc) <- c(term$labels,"total")
    colnames(desc) <- c("Unpenalized","Penalized")
    ## Return the results
    obj <- c(list(call=match.call(),mf=mf,terms=term,desc=desc),z)
    class(obj) <- c("ssanova0","ssanova")
    obj
}

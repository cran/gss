## Fit density model
ssden <- function(formula,type="cubic",data=list(),alpha=1.4,
                  weights=NULL,subset,na.action=na.omit,
                  id.basis=NULL,nbasis=NULL,seed=NULL,
                  domain=as.list(NULL),quadrature=NULL,ext=.05,order=2,
                  prec=1e-7,maxiter=30)
{
    ## Obtain model frame and model terms
    mf <- match.call()
    mf$type <- mf$alpha <- NULL
    mf$id.basis <- mf$nbasis <- mf$seed <- NULL
    mf$domain <- mf$quadrature <- mf$ext  <- NULL
    mf$prec <- mf$maxiter <- mf$order <- NULL
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf,sys.frame(sys.parent()))
    cnt <- model.weights(mf)
    mf$"(weights)" <- NULL
    ## set domain
    for (i in names(mf)) {
        if (is.factor(mf[[i]])) domain[[i]] <- levels(mf[[i]])[1:2]
        else {
            if (is.null(domain[[i]])) {
                mn <- min(mf[[i]])
                mx <- max(mf[[i]])
                range <- mx-mn
                mn <- mn - ext*range
                mx <- mx + ext*range
                domain[[i]] <- c(mn,mx)
            }
        }
    }
    domain <- as.data.frame(domain)
    ## Generate sub-basis
    nobs <- dim(mf)[1]
    if (is.null(id.basis)) {
        if (is.null(nbasis))  nbasis <- max(30,ceiling(10*nobs^(2/9)))
        if (nbasis>=nobs)  nbasis <- nobs
        if (!is.null(seed))  set.seed(seed)
        id.basis <- sample(nobs,nbasis,prob=cnt)
    }
    else {
        if (max(id.basis)>nobs|min(id.basis)<1)
            stop("gss error in ssden: id.basis out of range")
        nbasis <- length(id.basis)
    }
    ## Generate terms
    if (type=="cubic") term <- mkterm.cubic1(mf,domain)
    if (type=="linear") term <- mkterm.linear1(mf,domain)
    if (type=="tp") {
        if (is.null(quadrature))
            stop("gss error in ssden: quadrature needed for type tp")
        term <- mkterm.tp(mf,order,mf[id.basis,],1)
    }
    if (is.null(term)) stop("gss error in ssden: unknown type")
    term$labels <- term$labels[term$labels!="1"]
    ## Generate default quadrature
    if (is.null(quadrature)) {
        ## TO DO: HANDLING OF FACTORS
        domain <- domain[,colnames(mf),drop=FALSE]
        mn <- apply(domain,2,min)
        mx <- apply(domain,2,max)
        if (ncol(mf)==1) {
            ## Gauss-Legendre quadrature
            quad <- gauss.quad(200,c(mn,mx))
            quad$pt <- data.frame(quad$pt)
            colnames(quad$pt) <- colnames(mf)
        }
        else {
            ## Smolyak cubature
            if (ncol(mf)>4)
                stop("gss error in ssden: dimension higher than 4 unsupported")
            code <- c(15,14,13)
            quad <- smolyak.quad(ncol(mf),code[ncol(mf)-1])
            for (i in 1:ncol(mf)) {
                wk <- mf[,i]
                jk <- ssden(~wk,domain=data.frame(wk=domain[,i]),alpha=2,
                            id.basis=id.basis)
                quad$pt[,i] <- qssden(jk,quad$pt[,i])
                quad$wt <- quad$wt/dssden(jk,quad$pt[,i])
            }
            jk <- wk <- NULL
            quad$pt <- data.frame(quad$pt)
            colnames(quad$pt) <- colnames(mf)
        }
        quadrature <- list(pt=quad$pt,wt=quad$wt)
    }
    ## Generate s, r, and q
    qd.pt <- quadrature$pt
    qd.wt <- quadrature$wt
    nmesh <- length(qd.wt)
    s <- qd.s <- r <- qd.r <- q <- NULL
    nq <- 0
    for (label in term$labels) {
        x <- mf[,term[[label]]$vlist]
        x.basis <- mf[id.basis,term[[label]]$vlist]
        qd.x <- qd.pt[,term[[label]]$vlist]
        nphi <- term[[label]]$nphi
        nrk <- term[[label]]$nrk
        if (nphi) {
            phi <- term[[label]]$phi
            for (i in 1:nphi) {
                s <- cbind(s,phi$fun(x,nu=i,env=phi$env))
                qd.s <- cbind(qd.s,phi$fun(qd.x,nu=i,env=phi$env))
            }
        }
        if (nrk) {
            rk <- term[[label]]$rk
            for (i in 1:nrk) {
                nq <- nq+1
                r <- array(c(r,rk$fun(x.basis,x,nu=i,env=rk$env,out=TRUE)),c(nbasis,nobs,nq))
                qd.r <- array(c(qd.r,rk$fun(x.basis,qd.x,nu=i,env=rk$env,out=TRUE)),
                              c(nbasis,nmesh,nq))
                q <- array(c(q,rk$fun(x.basis,x.basis,nu=i,env=rk$env,out=TRUE)),
                           c(nbasis,nbasis,nq))
            }
        }
    }
    if (!is.null(s)) {
        nnull <- dim(s)[2]
        ## Check s rank
        if (qr(s)$rank<nnull)
            stop("gss error in ssden: fixed effect MLE is not unique")
        s <- t(s)
        qd.s <- t(qd.s)
    }
    ## Fit the model
    if (nq==1) {
        r <- r[,,1]
        qd.r <- qd.r[,,1]
        q <- q[,,1]
        z <- sspdsty(s,r,q,cnt,qd.s,qd.r,qd.wt,prec,maxiter,alpha)
    }
    else z <- mspdsty(s,r,q,cnt,qd.s,qd.r,qd.wt,prec,maxiter,alpha)
    ## Brief description of model terms
    desc <- NULL
    for (label in term$labels)
        desc <- rbind(desc,as.numeric(c(term[[label]][c("nphi","nrk")])))
    desc <- rbind(desc,apply(desc,2,sum))
    rownames(desc) <- c(term$labels,"total")
    colnames(desc) <- c("Fixed","Random")
    ## Return the results
    obj <- c(list(call=match.call(),mf=mf,cnt=cnt,terms=term,desc=desc,alpha=alpha,
                  domain=domain,quad=quadrature,id.basis=id.basis),z)
    class(obj) <- "ssden"
    obj
}

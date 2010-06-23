## Fit density model
ssden <- function(formula,type=NULL,data=list(),alpha=1.4,
                  weights=NULL,subset,na.action=na.omit,
                  id.basis=NULL,nbasis=NULL,seed=NULL,
                  domain=as.list(NULL),quad=NULL,qdsz.depth=NULL,
                  prec=1e-7,maxiter=30,skip.iter=FALSE)
{
    ## Obtain model frame and model terms
    mf <- match.call()
    mf$type <- mf$alpha <- NULL
    mf$id.basis <- mf$nbasis <- mf$seed <- NULL
    mf$domain <- mf$quad <- mf$qdsz.depth <- NULL
    mf$prec <- mf$maxiter <- mf$skip.iter <- NULL
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf,parent.frame())
    cnt <- model.weights(mf)
    mf$"(weights)" <- NULL
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
    ## Set domain and/or generate quadrature
    if (is.null(quad)) {
        ## Set domain and type
        fac.list <- NULL
        for (xlab in names(mf)) {
            x <- mf[[xlab]]
            if (is.factor(x)) {
                fac.list <- c(fac.list,xlab)
                domain[[xlab]] <- NULL
            }
            else {
                if (!is.vector(x))
                    stop("gss error in ssden: no default quadrature")
                if (is.null(domain[[xlab]])) {
                    mn <- min(x)
                    mx <- max(x)
                    domain[[xlab]] <- c(mn,mx)+c(-1,1)*(mx-mn)*.05
                }
                else domain[[xlab]] <- c(min(domain[[xlab]]),max(domain[[xlab]]))
                if (is.null(type[[xlab]]))
                    type[[xlab]] <- list("cubic",domain[[xlab]])
                else {
                    if (length(type[[xlab]])==1)
                        type[[xlab]] <- list(type[[xlab]][[1]],domain[[xlab]])
                }
            }
        }
        ## Generate numerical quadrature
        domain <- data.frame(domain)
        mn <- domain[1,]
        mx <- domain[2,]
        dm <- ncol(domain)
        if (dm==1) {
            ## Gauss-Legendre quadrature
            quad <- gauss.quad(200,c(mn,mx))
            quad$pt <- data.frame(quad$pt)
            colnames(quad$pt) <- colnames(domain)
        }
        else {
            ## Smolyak cubature
            if (is.null(qdsz.depth)) qdsz.depth <- switch(min(dm,6)-1,18,14,10,9,7)
            quad <- smolyak.quad(dm,qdsz.depth)
            for (i in 1:ncol(domain)) {
                xlab <- colnames(domain)[i]
                wk <- mf[[xlab]]
                jk <- ssden(~wk,domain=data.frame(wk=domain[,i]),alpha=2,
                            id.basis=id.basis,weights=cnt)
                quad$pt[,i] <- qssden(jk,quad$pt[,i])
                quad$wt <- quad$wt/dssden(jk,quad$pt[,i])
            }
            jk <- wk <- NULL
            quad$pt <- data.frame(quad$pt)
            colnames(quad$pt) <- colnames(domain)
        }
        ## Incorporate factors in quadrature
        if (!is.null(fac.list)) {
            for (i in 1:length(fac.list)) {
                wk <-
                  expand.grid(levels(mf[[fac.list[i]]]),1:length(quad$wt))
                quad$wt <- quad$wt[wk[,2]]
                col.names <- c(fac.list[i],colnames(quad$pt))
                quad$pt <- data.frame(wk[,1],quad$pt[wk[,2],])
                colnames(quad$pt) <- col.names
            }
        }
        quad <- list(pt=quad$pt,wt=quad$wt)
    }
    else {
        for (xlab in names(mf)) {
            x <- mf[[xlab]]
            if (is.vector(x)&!is.factor(x)) {
                mn <- min(x,quad$pt[[xlab]])
                mx <- max(x,quad$pt[[xlab]])
                range <- c(mn,mx)+c(-1,1)*(mx-mn)*.05
                if (is.null(type[[xlab]]))
                    type[[xlab]] <- list("cubic",range)
                else {
                    if (length(type[[xlab]])==1)
                        type[[xlab]] <- list(type[[xlab]][[1]],range)
                    else {
                        mn0 <- min(type[[xlab]][[2]])
                        mx0 <- max(type[[xlab]][[2]])
                        if ((mn0>mn)|(mx0<mx))
                            stop("gss error in ssden: range not covering domain")
                    }
                }
            }
        }
    }
    ## Generate terms
    term <- mkterm(mf,type)
    term$labels <- term$labels[term$labels!="1"]
    ## Generate s and r
    qd.pt <- quad$pt
    qd.wt <- quad$wt
    nmesh <- length(qd.wt)
    s <- qd.s <- r <- qd.r <- NULL
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
        z <- sspdsty(s,r,r[,id.basis],cnt,qd.s,qd.r,qd.wt,prec,maxiter,alpha)
    }
    else z <- mspdsty(s,r,id.basis,cnt,qd.s,qd.r,qd.wt,prec,maxiter,alpha,skip.iter)
    ## Brief description of model terms
    desc <- NULL
    for (label in term$labels)
        desc <- rbind(desc,as.numeric(c(term[[label]][c("nphi","nrk")])))
    desc <- rbind(desc,apply(desc,2,sum))
    rownames(desc) <- c(term$labels,"total")
    colnames(desc) <- c("Unpenalized","Penalized")
    ## Return the results
    obj <- c(list(call=match.call(),mf=mf,cnt=cnt,terms=term,desc=desc,
                  alpha=alpha,domain=domain,quad=quad,id.basis=id.basis,
                  qdsz.depth=qdsz.depth,skip.iter=skip.iter),z)
    class(obj) <- "ssden"
    obj
}

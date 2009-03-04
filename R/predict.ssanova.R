## Calculate prediction and Bayesian SE from ssanova objects
predict.ssanova <- function(object,newdata,se.fit=FALSE,
                            include=object$terms$labels,...)
{
    nnew <- nrow(newdata)
    nbasis <- length(object$id.basis)
    nnull <- length(object$d)
    nz <- length(object$b)
    nn <- nbasis + nnull + nz
    ## Extract included terms
    term <- object$terms
    philist <- rklist <- NULL
    s <- r <- NULL
    nq <- 0
    for (label in include) {
        if (label=="1") {
            philist <- c(philist,term[[label]]$iphi)
            s <- cbind(s,rep(1,len=nnew))
            next
        }
        if (label=="partial") next
        if (label=="offset") next
        xnew <- newdata[,term[[label]]$vlist]
        x <- object$mf[object$id.basis,term[[label]]$vlist]
        nphi <- term[[label]]$nphi
        nrk <- term[[label]]$nrk
        if (nphi) {
            iphi <- term[[label]]$iphi
            phi <- term[[label]]$phi
            for (i in 1:nphi) {
                philist <- c(philist,iphi+(i-1))
                s <- cbind(s,phi$fun(xnew,nu=i,env=phi$env))
            }
        }
        if (nrk) {
            irk <- term[[label]]$irk
            rk <- term[[label]]$rk
            for (i in 1:nrk) {
                rklist <- c(rklist,irk+(i-1))
                nq <- nq+1
                r <- array(c(r,rk$fun(xnew,x,nu=i,env=rk$env,out=TRUE)),c(nnew,nbasis,nq))
            }
        }
    }
    if (any(include=="partial")) {
        nphi <- term$partial$nphi
        iphi <- term$partial$iphi
        for (i in 1:nphi) philist <- c(philist,iphi+(i-1))
        s <- cbind(s,newdata$partial)
    }
    r.wk <- matrix(0,nnew,nbasis)
    nq <- 0
    for (i in rklist) {
        nq <- nq + 1
        r.wk <- r.wk + 10^object$theta[i]*r[,,nq]
    }
    ## random effects
    if (nz) {
        if (is.null(newdata$random)) z.wk <- matrix(0,nnew,nz)
        else z.wk <- newdata$random
        r.wk <- cbind(r.wk,z.wk)
    }
    ## Compute posterior mean
    nphi <- length(philist)
    pmean <- as.vector(r.wk%*%c(object$c,object$b))
    if (nphi) pmean <- pmean + as.vector(s%*%object$d[philist])
    if (any(include=="offset")) {
        if (is.null(model.offset(object$mf)))
            stop("gss error: no offset in the fit")
        offset <- newdata$offset
        if (is.null(offset)) offset <- newdata$"(offset)"
        if (is.null(offset)) stop("gss error: missing offset")
        pmean <- pmean + offset
    }
    if (se.fit) {
        b <- object$varht/10^object$nlambda
        ## Get cr, dr, and sms
        z <- .Fortran("regaux",
                      as.double(object$chol), as.integer(nn),
                      as.integer(object$jpvt), as.integer(object$rkv),
                      drcr=as.double(object$se.aux[[1]]%*%t(r.wk)), as.integer(nnew),
                      sms=double(nnull^2), as.integer(nnull), double(nn*nnull),
                      PACKAGE="gss")[c("drcr","sms")]
        drcr <- matrix(z$drcr,nn,nnew)
        dr <- drcr[1:nnull,,drop=FALSE][philist,,drop=FALSE]
        sms <- 10^object$nlambda*matrix(z$sms,nnull,nnull)[philist,philist]
        ## Compute posterior variance
        rr <- r.wk%*%object$qinv
        cr <- r.wk%*%t(object$se.aux[[2]])
        fn2 <- function(x,n) x[1:n]%*%x[n+(1:n)]
        pvar <- apply(t(cbind(r.wk,rr)),2,fn2,nbasis+nz)
        pvar <- pvar - apply(cbind(cr,cr),1,fn2,dim(cr)[2])
        if (nphi) {
            fn1 <- function(x,sms) t(x)%*%sms%*%x
            pvar <- pvar + apply(s,1,fn1,sms)
            pvar <- pvar - 2*apply(rbind(t(s),dr),2,fn2,nphi)
        }
        pse <- as.numeric(sqrt(b*pvar))
        list(fit=pmean,se.fit=pse)
    }
    else pmean
}

## Calculate prediction and Bayesian SE from ssanova objects
predict.ssanova <- function(obj,newdata,se.fit=FALSE,
                            include=obj$terms$labels)
{
    nnew <- dim(newdata)[1]
    nobs <- length(obj$c)
    ## Extract included terms
    term <- obj$terms
    philist <- rklist <- NULL
    s <- q <- NULL
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
        x <- obj$mf[,term[[label]]$vlist]
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
                q <- array(c(q,rk$fun(xnew,x,nu=i,env=rk$env,out=TRUE)),c(nnew,nobs,nq))
            }
        }
    }
    if (any(include=="partial")) {
        nphi <- term$partial$nphi
        iphi <- term$partial$iphi
        for (i in 1:nphi) philist <- c(philist,iphi+(i-1))
        s <- cbind(s,newdata$partial)
    }
    qq <- matrix(0,nnew,nobs)
    nq <- 0
    for (i in rklist) {
        nq <- nq + 1
        qq <- qq + 10^obj$theta[i]*q[,,nq]
    }
    if (!is.null(obj$w)) w <- obj$w
    else w <- model.weights(obj$mf)
    if (!is.null(w)) qq <- t(sqrt(w)*t(qq))
    ## Compute posterior mean
    nphi <- length(philist)
    pmean <- as.vector(qq%*%obj$c)
    if (nphi) pmean <- pmean + as.vector(s%*%obj$d[philist])
    if (any(include=="offset")) {
        if (is.null(model.offset(obj$mf)))
            stop("gss error: no offset in the fit")
        offset <- newdata$offset
        if (is.null(offset)) offset <- newdata$"(offset)"
        if (is.null(offset)) stop("gss error: missing offset")
        pmean <- pmean + offset
    }
    if (se.fit) {
        b <- obj$varht/10^obj$nlambda
        ## Get cr, dr, and sms
        crdr <- getcrdr(obj,t(qq))
        cr <- crdr$cr
        dr <- crdr$dr[philist,,drop=FALSE]
        sms <- getsms(obj)[philist,philist]
        ## Compute posterior variance
        r <- 0
        for (label in include) {
            if (label=="1") next
            xnew <- newdata[,term[[label]]$vlist]
            nrk <- term[[label]]$nrk
            if (nrk) {
                irk <- term[[label]]$irk
                rk <- term[[label]]$rk
                for (i in 1:nrk) {
                    ind <- irk+(i-1)
                    r <- r + 10^obj$theta[ind]*rk$fun(xnew,xnew,nu=i,env=rk$env)
                }
            }
        }
        fn2 <- function(x,n) x[1:n]%*%x[n+(1:n)]
        pvar <- r - apply(rbind(t(qq),cr),2,fn2,nobs)
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

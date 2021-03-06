## S3 method
predict1 <- function (object,...) UseMethod("predict1")
## Calculate prediction and Bayesian SE from ssanova objects
predict1.ssanova <- function(object,contr=c(1,-1),newdata,se.fit=TRUE,
                             include=c(object$terms$labels,object$lab.p),...)
{
    if (class(object)[1]%in%c("ssanova0","gssanova0"))
        stop("gss error: predict1 is not implemented for ssanova0/gssanova0")
    ncontr <- length(contr)
    nnew <- nrow(newdata[[1]])
    nbasis <- length(object$id.basis)
    nnull <- length(object$d)
    nz <- length(object$b)
    nn <- nbasis + nnull + nz
    labels.p <- object$lab.p
    ## Extract included terms
    term <- object$terms
    pmean.c <- ss.c <- rr.c <- 0
    for (jj in 1:ncontr) {
        philist <- rklist <- NULL
        s <- r <- NULL
        nq <- 0
        for (label in include) {
            if (label=="1") {
                philist <- c(philist,term[[label]]$iphi)
                s <- cbind(s,rep(1,len=nnew))
                next
            }
            if (label%in%labels.p) next
            if (label=="offset") next
            xnew <- newdata[[jj]][,term[[label]]$vlist]
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
        if (!is.null(object$partial)) {
            vars.p <- as.character(attr(object$partial$mt,"variables"))[-1]
            facs.p <- attr(object$partial$mt,"factors")
            vlist <- vars.p[as.logical(apply(facs.p,1,sum))]
            for (lab in labels.p) {
                if (lab%in%include) {
                    vlist.wk <- vars.p[as.logical(facs.p[,lab])]
                    vlist <- vlist[!(vlist%in%vlist.wk)]
                }
            }
            if (length(vlist)) {
                for (lab in vlist) newdata[[jj]][[lab]] <- 0
            }
            matx.p <- model.matrix(object$partial$mt,newdata[[jj]])[,-1,drop=FALSE]
            matx.p <- sweep(matx.p,2,object$partial$center)
            matx.p <- sweep(matx.p,2,object$partial$scale,"/")
            nu <- nnull-dim(matx.p)[2]
            for (label in labels.p) {
                nu <- nu+1
                if (label%in%include) {
                    philist <- c(philist,nu)
                    s <- cbind(s,matx.p[,label])
                }
            }
        }
        r.wk <- matrix(0,nnew,nbasis)
        nq <- 0
        for (i in rklist) {
            nq <- nq + 1
            r.wk <- r.wk + 10^object$theta[i]*r[,,nq]
        }
        ## random effects
        if (nz) {
            if (is.null(newdata[[jj]]$random)) z.wk <- matrix(0,nnew,nz)
            else z.wk <- newdata[[jj]]$random
            r.wk <- cbind(r.wk,z.wk)
        }
        ## Compute posterior mean
        nphi <- length(philist)
        pmean <- as.vector(r.wk%*%c(object$c,object$b))
        if (nphi) pmean <- pmean + as.vector(s%*%object$d[philist])
        if (any(include=="offset")) {
            if (is.null(model.offset(object$mf)))
                stop("gss error: no offset in the fit")
            offset <- newdata[[jj]]$offset
            if (is.null(offset)) offset <- newdata[[jj]]$"(offset)"
            if (is.null(offset)) stop("gss error: missing offset")
            pmean <- pmean + offset
        }
        pmean.c <- pmean.c + contr[jj]*pmean
        if (se.fit) {
            ss <- matrix(0,nnull,nnew)
            if (!is.null(philist)) ss[philist,] <- t(s)
            rr <- t(r.wk%*%object$se.aux$vec)
        }
        ss.c <- ss.c + contr[jj]*ss
        rr.c <- rr.c + contr[jj]*rr
    }
    if (se.fit) {
        b <- object$varht/10^object$nlambda
        ## Compute posterior variance
        wk <- object$se.aux$hfac%*%rbind(ss.c,rr.c)
        pse <- sqrt(b*apply(wk^2,2,sum))
        list(fit=pmean.c,se.fit=pse)
    }
    else pmean.c
}

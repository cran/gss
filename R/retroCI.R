retroCI <- function(object,newdata,limits=c(0.025,0.975),
                    include=c(object$fit$terms$labels,object$fit$lab.p)) {
    if (!(class(object)%in%"clone")) stop("gss error in retroCI: not a clone object")
    ## 
    nnew <- nrow(newdata)
    obj0 <- object$fit
    if (any(include=="offset")) {
        if (is.null(model.offset(obj0$mf)))
            stop("gss error: no offset in the fit")
        offset <- newdata$offset
        if (is.null(offset)) offset <- newdata$"(offset)"
        if (is.null(offset)) stop("gss error: missing offset")
    }
    nbasis <- length(obj0$id.basis)
    nnull <- length(obj0$d)
    nz <- length(obj0$b)
    nn <- nbasis + nnull + nz
    labels.p <- obj0$lab.p
    ## Extract included terms
    term <- obj0$terms
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
        xnew <- newdata[,term[[label]]$vlist]
        x <- obj0$mf[obj0$id.basis,term[[label]]$vlist]
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
    if (!is.null(obj0$partial)) {
        vars.p <- as.character(attr(obj0$partial$mt,"variables"))[-1]
        facs.p <- attr(obj0$partial$mt,"factors")
        vlist <- vars.p[as.logical(apply(facs.p,1,sum))]
        for (lab in labels.p) {
            if (lab%in%include) {
                vlist.wk <- vars.p[as.logical(facs.p[,lab])]
                vlist <- vlist[!(vlist%in%vlist.wk)]
            }
        }
        if (length(vlist)) {
            for (lab in vlist) newdata[[lab]] <- 0
        }
        matx.p <- model.matrix(obj0$partial$mt,newdata)[,-1,drop=FALSE]
        matx.p <- sweep(matx.p,2,obj0$partial$center)
        matx.p <- sweep(matx.p,2,obj0$partial$scale,"/")
        nu <- nnull-dim(matx.p)[2]
        for (label in labels.p) {
            nu <- nu+1
            if (label%in%include) {
                philist <- c(philist,nu)
                s <- cbind(s,matx.p[,label])
            }
        }
    }
    if (nz) {
        if (is.null(newdata$random)) z.wk <- matrix(0,nnew,nz)
        else z.wk <- newdata$random
    }
    ##
    nphi <- length(philist)
    nrep <- dim(object$c)[2]
    eta <- NULL
    for (i in 1:nrep) {
        r.wk <- matrix(0,nnew,nbasis)
        nq <- 0
        for (j in rklist) {
            nq <- nq + 1
            r.wk <- r.wk + 10^object$theta[j,i]*r[,,nq]
        }
        pmean <- as.vector(r.wk%*%object$c[,i])
        if (nz) pmean <- pmean + as.vector(z.wk%*%object$b[,i])
        if (nphi) pmean <- pmean + as.vector(s%*%object$d[philist,i])
        if (any(include=="offset")) pmean <- pmean + offset
        eta <- cbind(eta,pmean)
    }
    t(apply(eta,1,quantile,limits))
}

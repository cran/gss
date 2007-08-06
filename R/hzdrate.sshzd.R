hzdrate.sshzd <- ## Evaluate hazard estimate
function (object,x,se=FALSE) {
    if (class(object)!="sshzd")
        stop("gss error in hzdrate.sshzd: not a sshzd object")
    if (dim(object$mf)[2]==1&is.vector(x)) {
        x <- data.frame(x)
        colnames(x) <- colnames(object$mf)
    }
    s <- NULL
    r <- matrix(0,dim(x)[1],length(object$id.basis))
    nq <- 0
    for (label in object$terms$labels) {
        if (label=="1") {
            s <- cbind(s,rep(1,dim(x)[1]))
            next
        }
        xx <- object$mf[object$id.basis,object$terms[[label]]$vlist]
        x.new <- x[,object$terms[[label]]$vlist]
        nphi <- object$terms[[label]]$nphi
        nrk <-  object$terms[[label]]$nrk
        if (nphi) {
            phi <-  object$terms[[label]]$phi
            for (i in 1:nphi) {
                s <- cbind(s,phi$fun(x.new,nu=i,env=phi$env))
            }
        }
        if (nrk) {
            rk <- object$terms[[label]]$rk
            for (i in 1:nrk) {
                nq <- nq + 1
                r <- r + 10^object$theta[nq]*rk$fun(x.new,xx,nu=i,env=rk$env,out=TRUE)
            }
        }
    }
    rs <- cbind(r,s)
    if (!se) as.vector(exp(rs%*%c(object$c,object$d)))
    else {
        fit <- as.vector(exp(rs%*%c(object$c,object$d)))
        se.fit <- .Fortran("hzdaux2",
                           as.double(object$se.aux$v), as.integer(dim(rs)[2]),
                           as.integer(object$se.aux$jpvt),
                           as.double(t(rs)), as.integer(dim(rs)[1]),
                           se=double(dim(rs)[1]), PACKAGE="gss")[["se"]]
        list(fit=fit,se.fit=se.fit)
    }
}

hzdcurve.sshzd <- ## Evaluate hazard curve for plotting
function (object,time,covariates=NULL,se=FALSE) {
    tname <- object$tname
    xnames <- object$xnames
    if (class(object)!="sshzd")
        stop("gss error in hzdcurve.sshzd: not a sshzd object")
    if (length(xnames)&&(!all(xnames%in%names(covariates))))
        stop("gss error in survexp.sshzd: missing covariates")
    mn <- min(object$tdomain)
    mx <- max(object$tdomain)
    if ((min(time)<mn)|(max(time)>mx))
        stop("gss error in hzdcurve.sshzd: time range over the domain")
    if (length(xnames)) {
        xx <- covariates[,xnames,drop=FALSE]
        xy <- data.frame(matrix(0,length(time),length(xnames)+1))
        names(xy) <- c(tname,xnames)
        xy[,tname] <- time
    }
    else xx <- NULL
    if (!se) {
        if (is.null(xx))
            zz <- hzdrate.sshzd(object,time)
        else {
            zz <- NULL
            for (i in 1:dim(xx)[1]) {
                xy[,xnames] <- xx[rep(i,length(time)),]
                zz <- cbind(zz,hzdrate.sshzd(object,xy))
            }
            zz <- zz[,,drop=TRUE]
        }
        zz
    }
    else {
        if (is.null(xx))
            zz <- hzdrate.sshzd(object,time,TRUE)
        else {
            fit <- se.fit <- NULL
            for (i in 1:dim(xx)[1]) {
                xy[,xnames] <- xx[rep(i,length(time)),]
                wk <- hzdrate.sshzd(object,xy,TRUE)
                fit <- cbind(fit,wk$fit)
                se.fit <- cbind(se.fit,wk$se.fit)
            }
            zz <- list(fit=fit[,,drop=TRUE],se.fit=se.fit[,,drop=TRUE])
        }
        zz
    }
}

survexp.sshzd <- ## Compute expected survival
function(object,time,covariates=NULL,start=0) {
    tname <- object$tname
    xnames <- object$xnames
    ## Check inputs
    if (class(object)!="sshzd")
        stop("gss error in survexp.sshzd: not a sshzd object")
    if (length(xnames)&&(!all(xnames%in%names(covariates))))
        stop("gss error in survexp.sshzd: missing covariates")
    lmt <- cbind(start,time)
    if (any(lmt[,1]>lmt[,2]))
        stop("gss error in survexp.sshzd: start after follow-up time")
    nt <- dim(lmt)[1]
    if (is.null(covariates)) ncov <- 1
    else ncov <- dim(covariates)[1]
    if (length(xnames)&&(nt-1)&&(ncov-1)&&(nt-ncov))
        stop("gss error in survexp.sshzd: size mismatch")
    mn <- min(object$tdomain)
    mx <- max(object$tdomain)
    if ((min(start)<mn)|(max(time)>mx))
        stop("gss error in survexp.sshzd: time range over the domain")
    ## Calculate
    if (is.null(covariates)) {
        zz <- NULL
        for (i in 1:nt) {
            nqd <- max(20,ceiling((lmt[i,2]-lmt[i,1])/(mx-mn)*200))
            quad <- gauss.quad(nqd,lmt[i,])
            zz <- c(zz,sum(quad$wt*hzdrate.sshzd(object,quad$pt)))
        }
    }
    else {
        if (ncov>nt)
            lmt <- matrix(lmt,ncov,2,byrow=TRUE)
        if (ncov<nt)
            covariates <- covariates[rep(1,nt),,drop=FALSE]
        zz <- NULL
        for (i in 1:max(ncov,nt)) {
            nqd <- max(20,ceiling((lmt[i,2]-lmt[i,1])/(mx-mn)*200))
            quad <- gauss.quad(nqd,lmt[i,])
            wk <- covariates[rep(i,nqd),,drop=FALSE]
            wk[[tname]] <- quad$pt
            zz <- c(zz,sum(quad$wt*hzdrate.sshzd(object,wk)))
        }
    }
    exp(-zz)
}

## Make phi and rk for thin-plate spline model terms
mkterm.tp <- function(mf,order,mesh,weight)
{
    order <- max(order,1)
    ## Obtain model terms
    mt <- attr(mf,"terms")
    xvars <- as.character(attr(mt,"variables"))[-1]
    xfacs <- attr(mt,"factors")
    term.labels <- labels(mt)
    if (attr(attr(mf,"terms"),"intercept"))
        term.labels <- c("1",term.labels)
    ## Create the phi and rk functions
    term <- list(labels=term.labels)
    iphi.wk <- 1
    irk.wk <- 1
    for (label in term.labels) {
        iphi <- irk <- phi <- rk <- NULL
        if (label=="1") {
            ## the constant term
            iphi <- iphi.wk
            iphi.wk <- iphi.wk + 1
            term[[label]] <- list(iphi=iphi,nphi=1,nrk=0)
            next
        }
        vlist <- xvars[as.logical(xfacs[,label])]
        x <- mf[,vlist]
        xmesh <- mesh[,vlist]
        dm <- length(vlist)
        if (dm==1) {
            if (!is.factor(x)) {
                ## numeric variable
                if (is.vector(x)) xdim <- 1
                else xdim <- dim(x)[2]
                ## phi
                phi.env <- mkphi.tp(xdim,order,xmesh,weight)
                phi.fun <- function(x,nu,env) {
                    env$fun(x,nu+1,env$env)
                }
                nphi <- choose(xdim+order-1,xdim)-1
                iphi <- iphi.wk
                iphi.wk <- iphi.wk + nphi
                phi <- list(fun=phi.fun,env=phi.env)
                ## rk
                rk.env <- mkrk.tp(xdim,order,xmesh,weight)
                rk.fun <- function(x,y,nu=1,env,outer.prod=FALSE) {
                    env$fun(x,y,env$env,outer.prod)
                }
                nrk <- 1
                irk <- irk.wk
                irk.wk <- irk.wk + nrk
                rk <- list(fun=rk.fun,env=rk.env)
            }
            else {
                ## factor variable
                if (!is.ordered(x)) fun.env <- mkrk.nominal(levels(x))
                else fun.env <- mkrk.ordinal(levels(x))
                if (nlevels(x)>2) {
                    ## phi
                    nphi <- 0
                    ## rk
                    rk.fun <- function(x,y,nu=1,env,outer.prod=FALSE) {
                        env$fun(x,y,env$env,outer.prod)
                    }
                    nrk <- 1
                    irk <- irk.wk
                    irk.wk <- irk.wk + nrk
                    rk <- list(fun=rk.fun,env=fun.env)
                }
                else {
                    ## phi
                    phi.fun <- function(x,nu=1,env) {
                        wk <- as.factor(names(env$env$code)[1])
                        env$fun(x,wk,env$env)
                    }
                    nphi <- 1
                    iphi <- iphi.wk
                    iphi.wk <- iphi.wk + nphi
                    phi <- list(fun=phi.fun,env=fun.env)
                    ## rk
                    nrk <- 0
                }
            }
        }
        else {
            bin.fac <- xdim <- phi.list <- rk.list <- NULL
            for (i in 1:dm) {
                if (!is.factor(x[[i]])) {
                    ## numeric variable
                    if (is.vector(x[[i]])) xdim <- c(xdim,1)
                    else xdim <- c(xdim,dim(x[[i]])[2])
                    phi.wk <- mkphi.tp(xdim[i],order,xmesh[[i]],weight)
                    rk.wk <- mkrk.tp(xdim[i],order,xmesh[[i]],weight)
                    bin.fac <- c(bin.fac,0)
                }
                else {
                    ## factor variable
                    xdim <- c(xdim,0)
                    if (!is.ordered(x[[i]]))
                        rk.wk <- mkrk.nominal(levels(x[[i]]))
                    else rk.wk <- mkrk.ordinal(levels(x[[i]]))
                    phi.wk <- rk.wk
                    bin.fac <- c(bin.fac,!(nlevels(x[[i]])>2))
                }
                phi.list <- c(phi.list,list(phi.wk))
                rk.list <- c(rk.list,list(rk.wk))
            }
            n.phi <- choose(xdim+order-1,xdim)-1
            ## phi
            if (!all(n.phi+bin.fac)) nphi <- 0
            else {
                phi.env <- list(dim=dm,phi=phi.list,n.phi=n.phi,bin.fac=bin.fac)
                phi.fun <- function(x,nu,env) {
                    ind <- nu - 1
                    z <- 1
                    for (i in 1:env$dim) {
                        if (env$bin.fac[i]) {
                            wk <- as.factor(names(env$phi[[i]]$env$code)[1])
                            z <- z * env$phi[[i]]$fun(x[[i]],wk,env$phi[[i]]$env)
                        }
                        else {
                            code <- ind%%env$n.phi[i] + 1
                            ind <- ind%/%env$n.phi[i]
                            z <- z * env$phi[[i]]$fun(x[[i]],code+1,env$phi[[i]]$env)
                        }
                    }
                    z
                }
                nphi <- prod(n.phi+bin.fac)
                iphi <- iphi.wk
                iphi.wk <- iphi.wk + nphi
                phi <- list(fun=phi.fun,env=phi.env)
            }
            ## rk
            rk.env <- list(dim=dm,n.phi=n.phi,nphi=nphi,
                           phi=phi.list,rk=rk.list)
            rk.fun <- function(x,y,nu,env,outer.prod=FALSE) {
                n.rk <- ifelse(env$n.phi,2,1)
                ind <- nu - !env$nphi
                z <- 1
                for (i in 1:env$dim) {
                    code <- ind%%n.rk[i] + 1
                    ind <- ind%/%n.rk[i]
                    if (code==n.rk[i]) {
                        z <- z * env$rk[[i]]$fun(x[[i]],y[[i]],
                                                 env$rk[[i]]$env,outer.prod)
                    }
                    else {
                        z.wk <- 0
                        for (j in 1:env$n.phi[i]) {
                            phix <- env$phi[[i]]$fun(x[[i]],j+1,env$phi[[i]]$env)
                            phiy <- env$phi[[i]]$fun(y[[i]],j+1,env$phi[[i]]$env)
                            if (outer.prod) z.wk <- z.wk + outer(phix,phiy)
                            else z.wk <- z.wk + phix * phiy
                        }
                        z <- z * z.wk
                    }
                }
                z
            }
            n.rk <- ifelse(n.phi,2,1)
            nrk <- prod(n.rk) - as.logical(nphi)
            irk <- irk.wk
            irk.wk <- irk.wk + nrk
            rk <- list(fun=rk.fun,env=rk.env)
        }
        term[[label]] <- list(vlist=vlist,
                              iphi=iphi,nphi=nphi,phi=phi,
                              irk=irk,nrk=nrk,rk=rk)
    }
    term
}

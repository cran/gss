## Make phi and rk for cubic spline model terms
mkterm.cubic1 <- function(mf,range)
{
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
        lmt <- range[,vlist]
        dm <- length(vlist)
        if (dm==1) {
            if (!is.factor(x)) {
                ## numeric variable
                mn <- min(lmt)
                mx <- max(lmt)
                ## phi
                phi.env <- mkphi.cubic(c(mn,mx))
                phi.fun <- function(x,nu=1,env) env$fun(x,env$env)
                nphi <- 1
                iphi <- iphi.wk
                iphi.wk <- iphi.wk + nphi
                phi <- list(fun=phi.fun,env=phi.env)
                ## rk
                rk.env <- mkrk.cubic(c(mn,mx))
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
            bin.fac <- n.phi <- phi.list <- rk.list <- NULL
            for (i in 1:dm) {
                if (!is.factor(x[[i]])) {
                    ## numeric variable
                    mn <- min(lmt[[i]])
                    mx <- max(lmt[[i]])
                    phi.wk <- mkphi.cubic(c(mn,mx))
                    rk.wk <- mkrk.cubic(c(mn,mx))
                    n.phi <- c(n.phi,1)
                    bin.fac <- c(bin.fac,0)
                }
                else {
                    ## factor variable
                    if (!is.ordered(x[[i]]))
                        rk.wk <- mkrk.nominal(levels(x[[i]]))
                    else rk.wk <- mkrk.ordinal(levels(x[[i]]))
                    phi.wk <- rk.wk
                    n.phi <- c(n.phi,0)
                    bin.fac <- c(bin.fac,!(nlevels(x[[i]])>2))
                }
                phi.list <- c(phi.list,list(phi.wk))
                rk.list <- c(rk.list,list(rk.wk))
            }
            ## phi
            if (sum(n.phi+bin.fac)<dm) nphi <- 0
            else {
                phi.env <- list(dim=dm,n.phi=n.phi,phi=phi.list)
                phi.fun <- function(x,nu=1,env) {
                    z <- 1
                    for (i in 1:env$dim) {
                        if (env$n.phi[i])
                            z <- z * env$phi[[i]]$fun(x[[i]],env$phi[[i]]$env)
                        else {
                            wk <- as.factor(names(env$phi[[i]]$env$code)[1])
                            z <- z * env$phi[[i]]$fun(x[[i]],wk,env$phi[[i]]$env)
                        }
                    }
                    z
                }
                nphi <- 1
                iphi <- iphi.wk
                iphi.wk <- iphi.wk + nphi
                phi <- list(fun=phi.fun,env=phi.env)
            }
            ## rk
            rk.env <- list(dim=dm,n.phi=n.phi,nphi=nphi,phi=phi.list,rk=rk.list)
            rk.fun <- function(x,y,nu,env,outer.prod=FALSE) {
                div <- env$n.phi + 1
                ind <- nu - 1 + env$nphi
                z <- 1
                for (i in 1:env$dim) {
                    code <- ind%%div[i] + 1
                    ind <- ind%/%div[i]
                    if (code==div[i])
                        z <- z * env$rk[[i]]$fun(x[[i]],y[[i]],
                                                 env$rk[[i]]$env,outer.prod)
                    else {
                        phix <- env$phi[[i]]$fun(x[[i]],env$phi[[i]]$env)
                        phiy <- env$phi[[i]]$fun(y[[i]],env$phi[[i]]$env)
                        if (outer.prod) z <- z * outer(phix,phiy)
                        else z <- z * phix * phiy
                    }
                }
                z
            }
            nrk <- prod(n.phi+1) - nphi
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

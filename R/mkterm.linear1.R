## Make phi and rk for linear spline model terms
mkterm.linear1 <- function(mf,range)
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
    iphi.wk <- irk.wk <- 1
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
                mx <- max(lmt)
                mn <- min(lmt)
                ## phi
                nphi <- 0
                ## rk
                rk.env <- mkrk.linear(c(mn,mx))
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
            bin.fac <- rk.list <- NULL
            for (i in 1:dm) {
                if (!is.factor(x[[i]])) {
                    ## numeric variable
                    mx <- max(lmt[[i]])
                    mn <- min(lmt[[i]])
                    rk.wk <- mkrk.linear(c(mn,mx))
                    bin.fac <- c(bin.fac,0)
                }
                else {
                    ## factor variable
                    if (!is.ordered(x[[i]]))
                        rk.wk <- mkrk.nominal(levels(x[[i]]))
                    else rk.wk <- mkrk.ordinal(levels(x[[i]]))
                    bin.fac <- c(bin.fac,!(nlevels(x[[i]])>2))
                }
                rk.list <- c(rk.list,list(rk.wk))
            }
            rk.env <- list(dim=dm,rk=rk.list)
            if (sum(bin.fac)==dm) {
                ## phi
                phi.fun <- function(x,nu=1,env) {
                    z <- 1
                    for (i in 1:env$dim) {
                        wk <- as.factor(names(env$rk[[i]]$env$code)[1])
                        z <- z * env$rk[[i]]$fun(x[[i]],wk,env$rk[[i]]$env)
                    }
                    z
                }
                nphi <- 1
                iphi <- iphi.wk
                iphi.wk <- iphi.wk + nphi
                phi <- list(fun=phi.fun,env=rk.env)
                ## rk
                nrk <- 0
            }
            else {             
                ## phi
                nphi <- 0
                ## rk
                rk.fun <- function(x,y,nu,env,outer.prod=FALSE) {
                    z <- 1
                    for (i in 1:env$dim)
                        z <- z * env$rk[[i]]$fun(x[[i]],y[[i]],
                                                 env$rk[[i]]$env,outer.prod)
                    z
                }
                nrk <- 1
                irk <- irk.wk
                irk.wk <- irk.wk + nrk
                rk <- list(fun=rk.fun,env=rk.env)
            }
        }
        term[[label]] <- list(vlist=vlist,
                              iphi=iphi,nphi=nphi,phi=phi,
                              irk=irk,nrk=nrk,rk=rk)
    }
    term
}

mkterm.cubic <- ## Make phi and rk for cubic spline model terms
function(mf,ext) {
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
    if (label=="1") {                   # the constant term
      iphi <- iphi.wk
      iphi.wk <- iphi.wk + 1
      term[[label]] <- list(iphi=iphi,nphi=1,nrk=0)
      next
    }
    vlist <- xvars[as.logical(xfacs[,label])]
    x <- mf[,vlist]
    dm <- length(vlist)
    if (dm==1) {
      mx <- max(x)
      mn <- min(x)
      range <- mx - mn
      ## phi
      phi.env <- list(phi=mkphi.cubic()$fun,
                      mn=mn-ext*range,mx=mx+ext*range)
      phi.fun <- function(x,nu=1,env) {
        env$phi((x-env$mn)/(env$mx-env$mn),2)
      }
      nphi <- 1
      iphi <- iphi.wk
      iphi.wk <- iphi.wk + nphi
      phi <- list(fun=phi.fun,env=phi.env)
      ## rk
      rk.env <- list(rk=mkrk.cubic()$fun,
                     mn=mn-ext*range,mx=mx+ext*range)
      rk.fun <- function(x,y,nu=1,env,outer.prod=FALSE) {
        x <- (x-env$mn)/(env$mx-env$mn)
        y <- (y-env$mn)/(env$mx-env$mn)
        env$rk(x,y,out=outer.prod)
      }
      nrk <- 1
      irk <- irk.wk
      irk.wk <- irk.wk + nrk
      rk <- list(fun=rk.fun,env=rk.env)
    }
    else {
      mx <- apply(x,2,max)
      mn <- apply(x,2,min)
      range <- mx - mn
      ## phi
      phi.env <- list(phi=mkphi.cubic()$fun,
                      dim=dm,mn=mn-ext*range,mx=mx+ext*range)
      phi.fun <- function(x,nu=1,env) {
        if (is.vector(x)) x <- t(as.matrix(x))
        if (env$dim!=dim(x)[2]) {
          stop("gss error in phi: inputs are of wrong dimensions")
        }
        x <- t((t(x)-env$mn)/(env$mx-env$mn))
        z <- env$phi(x[,1],2)
        for (i in 2:env$dim) z <- z * env$phi(x[,i],2)
        z
      }
      nphi <- 1
      iphi <- iphi.wk
      iphi.wk <- iphi.wk + nphi
      phi <- list(fun=phi.fun,env=phi.env)
      ## rk
      rk.env <- list(rk=mkrk.cubic()$fun,phi=mkphi.cubic()$fun,
                     dim=dm,mn=mn-ext*range,mx=mx+ext*range)
      rk.fun <- function(x,y,nu,env,outer.prod=FALSE) {
        if (is.vector(x)) x <- t(as.matrix(x))
        if (env$dim!=dim(x)[2]) {
          stop("gss error in rk: inputs are of wrong dimensions")
        }
        x <- t((t(x)-env$mn)/(env$mx-env$mn))
        if (is.vector(y)) y <- t(as.matrix(y))
        if (env$dim!=dim(y)[2]) {
          stop("gss error in rk: inputs are of wrong dimensions")
        }
        y <- t((t(y)-env$mn)/(env$mx-env$mn))
        z <- 1
        ind <- nu
        for (i in 1:env$dim) {
          code <- ind%%2
          ind <- ind%/%2
          if (code) z <- z * env$rk(x[,i],y[,i],out=outer.prod)
          else {
            phix <- env$phi(x[,i],2)
            phiy <- env$phi(y[,i],2)
            if (outer.prod) z <- z * outer(phix,phiy)
            else z <- z * phix * phiy
          }
        }
        z
      }
      nrk <- 2^dm-1
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

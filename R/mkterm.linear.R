mkterm.linear <- ## Make phi and rk for cubic spline model terms
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
      nphi <- 0
      ## rk
      rk.env <- list(rk=mkrk.linear()$fun,
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
      nphi <- 0
      ## rk
      rk.env <- list(rk=mkrk.linear()$fun,
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
        for (i in 1:env$dim)
          z <- z * env$rk(x[,i],y[,i],out=outer.prod)
        z
      }
      nrk <- 1
      irk <- irk.wk
      irk.wk <- irk.wk + nrk
      rk <- list(fun=rk.fun,env=rk.env)
    }
    term[[label]] <- list(vlist=vlist,nphi=nphi,
                          irk=irk,nrk=nrk,rk=rk)
  }
  term
}

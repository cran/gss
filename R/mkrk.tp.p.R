mkrk.tp.p <- ## Make pseudo RK for thin-plate splines
function(dm,order) {
  ## Check inputs
  if (!((2*order>dm)&(dm>=1))) {
    stop("gss error: thin-plate spline undefined for the parameters")
  }
  ## Create the environment
  if (dm%%2) {                    
    theta <- gamma(dm/2-order)/2^(2*order)/pi^(dm/2)/gamma(order)
  }
  else {
    theta <- (-1)^(dm/2+order+1)/2^(2*order-1)/pi^(dm/2)/
      gamma(order)/gamma(order-dm/2+1)
  }
  env <- list(dim=dm,order=order,theta=theta)
  ## Create the rk.p function
  fun <- function(x,y,env,outer.prod=FALSE) {
    ## Check inputs
    if (env$dim==1) {
      if (!(is.vector(x)&is.vector(y))) {
        stop("gss error in rk: inputs are of wrong types")
      }
    }
    else {
      if (is.vector(x)) x <- t(as.matrix(x))
      if (env$dim!=dim(x)[2]) {
        stop("gss error in rk: inputs are of wrong dimensions")
      }
      if (is.vector(y)) y <- t(as.matrix(y))
      if (env$dim!=dim(y)[2]) {
        stop("gss error in rk: inputs are of wrong dimensions")
      }
    }
    ## Return the results
    if (outer.prod) {               
      if (env$dim==1) {
        fn1 <- function(x,y) abs(x-y)
        d <- outer(x,y,fn1)
      }
      else {
        fn2 <- function(x,y) sqrt(sum((x-y)^2))
        d <- NULL
        for (i in 1:dim(y)[1]) d <- cbind(d,apply(x,1,fn2,y[i,]))
      }
    }
    else {
      if (env$dim==1) d <- abs(x-y)
      else {
        N <- max(dim(x)[1],dim(y)[1])
        x <- t(matrix(t(x),env$dim,N))
        y <- t(matrix(t(y),env$dim,N))
        fn <- function(x) sqrt(sum(x^2))
        d <- apply(x-y,1,fn)
      }
    }
    power <- 2*env$order-env$dim
    switch(1+env$dim%%2,
           env$theta*d^power*log(ifelse(d>0,d,1)),
           env$theta*d^power)
  }
  ## Return the function and the environment
  list(fun=fun,env=env)
}

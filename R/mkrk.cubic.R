mkrk.cubic <- ## Make RK for cubic splines
function() {
  ## Create the rk function
  fun <- function(x,y,env=NULL,outer.prod=FALSE) {
    ## Check inputs
    if (!(is.vector(x)&is.vector(y))) {
      stop("gss error in rk: inputs are of wrong types")
    }
    if ((min(x,y)<0)|(max(x,y)>1)) {
      stop("gss error in rk: inputs are out of range")
    }
    ## Return the results
    rk <- function(x,y) {
      k2 <- function(x) ((x-.5)^2-1/12)/2
      k4 <- function(x) ((x-.5)^4-(x-.5)^2/2+7/240)/24
      k2(x)*k2(y)-k4(abs(x-y))
    }
    if (outer.prod) outer(x,y,rk)
    else rk(x,y)
  }
  ## Return the function and the (null) environment
  list(fun=fun,env=NULL)
}

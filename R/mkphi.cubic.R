mkphi.cubic <- ## Make phi function for cubic splines
function() {
  ## Create the phi function
  fun <- function(x,nu,env=NULL) {
    ## Check inputs
    if (!is.vector(x)) {
      stop("gss error in phi: inputs are of wrong types")
    }
    if ((min(x)<0)|(max(x)>1)) {
      stop("gss error in phi: inputs are out of range")
    }
    ## Return the results
    phi1 <- function(x) rep(1,len=length(x))
    phi2 <- function(x) x-.5
    switch(nu, phi1(x), phi2(x))
  }
  ## Return the function and the (null) environment
  list(fun=fun,env=NULL)
}

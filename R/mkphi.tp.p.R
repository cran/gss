mkphi.tp.p <-  ## Make pseudo phi function for thin-plate splines
function(dm,order) {
  ## Check inputs
  if (!((2*order>dm)&(dm>=1))) {
    stop("gss error: thin-plate spline undefined for the parameters")
  }
  ## Create the environment
  pol.code <- NULL
  for (i in 0:(order^dm-1)) {
    ind <- i; code <- NULL
    for (j in 1:dm) {
      code <- c(code,ind%%order)
      ind <- ind%/%order
    }
    if (sum(code)<order) pol.code <- cbind(pol.code,code)
  }
  env <- list(dim=dm,pol.code=pol.code)
  ## Create the phi function  
  fun <- function(x,nu,env) {
    if (env$dim==1) x <- as.matrix(x)
    if (env$dim!=dim(x)[2]) {
      stop("gss error in phi: inputs are of wrong dimensions")
    }
    apply(t(x)^env$pol.code[,nu],2,prod)
  }
  ## Return the function and the environment
  list(fun=fun,env=env)
}

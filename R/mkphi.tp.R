mkphi.tp <-  ## Make phi function for thin-plate splines
function(dm,order,mesh,weight) {
  ## Check inputs
  if (!((2*order>dm)&(dm>=1))) {
    stop("gss error: thin-plate spline undefined for the parameters")
  }
  if (xor(is.vector(mesh),dm==1)
      |xor(is.matrix(mesh),dm>=2)) {
    stop("gss error in mkphi: mismatched inputs")
  }
  if ((min(weight)<0)|(max(weight)<=0)) {
    stop("gss error in mkphi: negative weights")
  }
  ## Set weights
  if (is.vector(mesh)) N <- length(mesh)
  else N <- dim(mesh)[1]
  weight <- rep(weight,len=N)
  weight <- sqrt(weight/sum(weight))
  ## Create the environment
  phi.p <- mkphi.tp.p(dm,order)
  nnull <- choose(dm+order-1,dm)
  s <- NULL
  for (nu in 1:nnull) s <- cbind(s,phi.p$fun(mesh,nu,phi.p$env))
  s <- qr(weight*s)
  if (s$rank<nnull) {
    stop("gss error in mkphi: insufficient normalizing mesh for thin-plate spline")
  }
  r <- qr.R(s)
  env <- list(dim=dm,order=order,phi.p=phi.p,r=r)
  ## Create the phi function
  fun <- function(x,nu,env) {
    nnull <- choose(env$dim+env$order-1,env$dim)
    phix <- NULL
    for(i in 1:nnull)
      phix <- rbind(phix,env$phi.p$fun(x,i,env$phi.p$env))
    t(backsolve(env$r,phix,tr=TRUE))[,nu]
  }
  ## Return the function and the environment
  list(fun=fun,env=env)
}

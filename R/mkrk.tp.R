mkrk.tp <- ## Make RK for thin-plate splines
function(dm,order,mesh,weight=1) {
  ## Check inputs
  if (!((2*order>dm)&(dm>=1))) {
    stop("gss error: thin-plate spline undefined for the parameters")
  }
  if (xor(is.vector(mesh),dm==1)
      |xor(is.matrix(mesh),dm>=2)) {
    stop("gss error in mkrk: mismatched inputs")
  }
  if ((min(weight)<0)|(max(weight)<=0)) {
    stop("gss error in mkrk: negative weights")
  }
  ## Set weights
  if (is.vector(mesh)) N <- length(mesh)
  else N <- dim(mesh)[1]
  weight <- rep(weight,len=N)
  weight <- sqrt(weight/sum(weight))
  ## Obtain orthonormal basis
  phi.p <- mkphi.tp.p(dm,order)
  nnull <- choose(dm+order-1,dm)
  s <- NULL
  for (nu in 1:nnull) s <- cbind(s,phi.p$fun(mesh,nu,phi.p$env))
  s <- qr(weight*s)
  if (s$rank<nnull) {
    stop("gss error in mkrk: insufficient normalizing mesh for thin-plate spline")
  }
  q <- qr.Q(s)
  r <- qr.R(s)
  ## Set Q^{T}E(|u_{i}-u_{j}|)Q
  rk.p <- mkrk.tp.p(dm,order)
  pep <- weight*t(weight*rk.p$fun(mesh,mesh,rk.p$env,out=TRUE))
  pep <- t(q)%*%pep%*%q
  ## Create the environment
  env <- list(dim=dm,order=order,weight=weight,
              phi.p=phi.p,rk.p=rk.p,q=q,r=r,mesh=mesh,pep=pep)
  ## Create the rk function
  fun <- function(x,y,env,outer.prod=FALSE) {
    ## Check inputs
    if (env$dim==1) {
      if (!(is.vector(x)&is.vector(y))) {
        stop("gss error in rk: inputs are of wrong types")
      }
      nx <- length(x)
      ny <- length(y)
    }
    else {
      if (is.vector(x)) x <- t(as.matrix(x))
      if (env$dim!=dim(x)[2]) {
        stop("gss error in rk: inputs are of wrong dimensions")
      }
      nx <- dim(x)[1]
      if (is.vector(y)) y <- t(as.matrix(y))
      if (env$dim!=dim(y)[2]) {
        stop("gss error in rk: inputs are of wrong dimensions")
      }
      ny <- dim(y)[1]
    }
    ## Return the results
    nnull <- choose(env$dim+env$order-1,env$dim)
    if (outer.prod) {
      phix <- phiy <- NULL
      for (nu in 1:nnull) {
        phix <- rbind(phix,env$phi.p$fun(x,nu,env$phi.p$env))
        phiy <- rbind(phiy,env$phi.p$fun(y,nu,env$phi.p$env))
      }
      phix <- backsolve(env$r,phix,tr=TRUE)
      phiy <- backsolve(env$r,phiy,tr=TRUE)
      ex <- env$rk.p$fun(env$mesh,x,env$rk.p$env,out=TRUE)
      ex <- env$weight*ex
      ex <- t(env$q)%*%ex
      ey <- env$rk.p$fun(env$mesh,y,env$rk.p$env,out=TRUE)
      ey <- env$weight*ey
      ey <- t(env$q)%*%ey
      env$rk.p$fun(x,y,env$rk.p$env,out=TRUE)-t(phix)%*%ey-
        t(ex)%*%phiy+t(phix)%*%env$pep%*%phiy
    }
    else {
      N <- max(nx,ny)
      phix <- phiy <- NULL
      for (nu in 1:nnull) {
        phix <- rbind(phix,env$phi.p$fun(x,nu,env$phi.p$env))
        phiy <- rbind(phiy,env$phi.p$fun(y,nu,env$phi.p$env))
      }
      phix <- backsolve(env$r,phix,tr=TRUE)
      phix <- matrix(phix,nnull,N)
      phiy <- backsolve(env$r,phiy,tr=TRUE)
      phiy <- matrix(phiy,nnull,N)
      ex <- env$rk.p$fun(env$mesh,x,env$rk.p$env,out=TRUE)
      ex <- env$weight*ex
      ex <- t(env$q)%*%ex
      ex <- matrix(ex,nnull,N)
      ey <- env$rk.p$fun(env$mesh,y,env$rk.p$env,out=TRUE)
      ey <- env$weight*ey
      ey <- t(env$q)%*%ey
      ey <- matrix(ey,nnull,N)
      fn1 <- function(x,n) x[1:n]%*%x[n+(1:n)]
      fn2 <- function(x,pep,n) t(x[1:n])%*%pep%*%x[n+(1:n)]
      env$rk.p$fun(x,y,env$rk.p$env)-apply(rbind(phix,ey),2,fn1,nnull)-
        apply(rbind(phiy,ex),2,fn1,nnull)+
          apply(rbind(phix,phiy),2,fn2,env$pep,nnull)
    }
  }
  ## Return the function and the environment
  list(fun=fun,env=env)
}

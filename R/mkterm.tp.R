mkterm.tp <- ## Make phi and rk for thin-plate spline model terms
function(mf,order,mesh,weight) {
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
    xmesh <- mesh[,vlist]
    dm <- length(vlist)
    if (dm==1) {
      if (is.vector(x)) xdim <- 1
      else xdim <- dim(x)[2]
      ## phi
      phi.env <- list(phi=mkphi.tp(xdim,order,xmesh,weight))
      phi.fun <- function(x,nu,env) {
        env$phi$fun(x,nu+1,env=env$phi$env)
      }
      nphi <- choose(xdim+order-1,xdim)-1
      iphi <- iphi.wk
      iphi.wk <- iphi.wk + nphi
      phi <- list(fun=phi.fun,env=phi.env)
      ## rk
      rk.env <- list(rk=mkrk.tp(xdim,order,xmesh,weight))
      rk.fun <- function(x,y,nu=1,env,outer.prod=FALSE) {
        env$rk$fun(x,y,env=env$rk$env,out=outer.prod)
      }
      nrk <- 1
      irk <- irk.wk
      irk.wk <- irk.wk + nrk
      rk <- list(fun=rk.fun,env=rk.env)
    }
    else {
      xdim <- phi.list <- rk.list <- NULL
      for (i in 1:dm) {
        if (is.vector(x[[i]])) xdim <- c(xdim,1)
        else xdim <- c(xdim,dim(x[[i]])[2])
        phi <- mkphi.tp(xdim[i],order,xmesh[[i]],weight)
        phi.list <- c(phi.list,list(phi))
        rk <- mkrk.tp(xdim[i],order,xmesh[[i]],weight)
        rk.list <- c(rk.list,list(rk))
      }
      ## phi
      nnphi <- choose(xdim+order-1,xdim)-1
      phi.env <- list(dim=dm,phi=phi.list,nnphi=nnphi)
      phi.fun <- function(x,nu,env) {
        nu.wk <- nu - 1
        code <- NULL
        for (i in 1:env$dim) {
          code <- c(code,nu.wk%%env$nnphi[i]+1)
          nu.wk <- nu.wk%/%env$nnphi[i]
        }
        z <- 1
        for (i in 1:env$dim) {
          z <- z * env$phi[[i]]$fun(x[[i]],code[i]+1,env=env$phi[[i]]$env)
        }
        z
      }
      nphi <- prod(nnphi)
      iphi <- iphi.wk
      iphi.wk <- iphi.wk + nphi
      phi <- list(fun=phi.fun,env=phi.env)
      ## rk
      rk.env <- list(dim=dm,phi=phi.list,rk=rk.list,nnphi=nnphi)
      rk.fun <- function(x,y,nu,env,outer.prod=FALSE) {
        nnrk <- ifelse(env$nnphi,2,1)
        ind <- nu - 1 + ifelse(all(nnrk==2),1,0)
        z <- 1
        for (i in 1:env$dim) {
          code <- ind%%nnrk[i] + 1
          ind <- ind%/%nnrk[i]
          if (code==nnrk[i]) {
            z <- z * env$rk[[i]]$fun(x[[i]],y[[i]],env=env$rk[[i]]$env,out=outer.prod)
          }
          else {
            z.wk <- 0
            for (j in 1:env$nnphi[i]) {
              phix <- env$phi[[i]]$fun(x[[i]],j+1,env=env$phi[[i]]$env)
              phiy <- env$phi[[i]]$fun(y[[i]],j+1,env=env$phi[[i]]$env)
              if (outer.prod) z.wk <- z.wk + outer(phix,phiy)
              else z.wk <- z.wk + phix * phiy
            }
            z <- z * z.wk
          }
        }
        z
      }
      nnrk <- ifelse(nnphi,2,1)
      nrk <- prod(nnrk) - ifelse(all(nnrk==2),1,0)
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

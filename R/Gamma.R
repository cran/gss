mkdata.Gamma <- ## Make pseudo data for gamma regression
function(y,eta,wt,offset) {
  if (is.null(wt)) wt <- rep(1,length(y))
  if (is.null(offset)) offset <- rep(0,length(y))
  if (min(y)<=0)
    stop("gss error: gamma responses should be positive")
  mu <- exp(eta)
  u <- 1-y/mu
  w <- y/mu
  ywk <- eta-u/w-offset
  wt <- w*wt
  list(ywk=ywk,wt=wt)
}

dev.resid.Gamma <- ## Calculate the deviance residuals of gamma fit
function(y,eta,wt) {
  if (is.null(wt)) wt <- rep(1,length(y))
  mu <- exp(eta)
  as.vector(2*wt*(-log(y/mu)+(y-mu)/mu))
}

dev.null.Gamma <-
function(y,wt,offset) {
  if (is.null(wt)) wt <- rep(1,length(y))
  mu <- sum(wt*y)/sum(wt)
  if (!is.null(offset)) {
    eta <- log(mu)-mean(offset)
    repeat {
      mu <- exp(eta+offset)
      u <- 1-y/mu
      w <- y/mu
      eta.new <- eta-sum(wt*u)/sum(wt*w)
      if (abs(eta-eta.new)/(1+abs(eta))<1e-7) break
      eta <- eta.new    
    }
  }
  sum(2*wt*(-log(y/mu)+(y-mu)/mu))
}

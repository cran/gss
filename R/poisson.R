mkdata.poisson <- ## Make pseudo data for Poisson regression
function(y,eta,wt,offset) {
  if (is.null(wt)) wt <- rep(1,length(y))
  if (is.null(offset)) offset <- rep(0,length(y))
  if (min(y)<0)
    stop("gss error: paired binomial response should be nonnegative")
  lambda <- exp(eta)
  u <- lambda - y
  w <- lambda
  ywk <- eta-u/w-offset
  wt <- w*wt
  list(ywk=ywk,wt=wt)
}

dev.resid.poisson <- ## Calculate deviance residuals for Poisson regression
function(y,eta,wt) {
  if (is.null(wt)) wt <- rep(1,length(y))
  lambda <- exp(eta)
  as.vector(2*wt*(y*log(ifelse(y==0,1,y/lambda))-(y-lambda)))
}

dev.null.poisson <-
function(y,wt,offset) {
  if (is.null(wt)) wt <- rep(1,length(y))
  lambda <- mean(y)
  if (!is.null(offset)) {
    eta <- log(lambda) - mean(offset)
    repeat {
      lambda <- exp(eta+offset)
      u <- lambda - y
      w <- lambda
      eta.new <- eta-sum(wt*u)/sum(wt*w)
      if (abs(eta-eta.new)/(1+abs(eta))<1e-7) break
      eta <- eta.new    
    }
  }
  sum(2*wt*(y*log(ifelse(y==0,1,y/lambda))-(y-lambda)))
}

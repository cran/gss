mkdata.nbinomial <- ## Make pseudo data for negative binomial regression
function(y,eta,wt,offset,alpha) {
  if (is.vector(y)) y <- as.matrix(y)
  if (is.null(wt)) wt <- rep(1,dim(y)[1])
  if (is.null(offset)) offset <- rep(0,dim(y)[1])
  if (dim(y)[2]==2) {
    if (min(y[,1])<0)
      stop("gss error: negative binomial response should be nonnegative")
    if (min(y[,2])<=0)
      stop("gss error: negative binomial size should be positive")
    p <- 1-1/(1+exp(eta))
    u <- (y[,1]+y[,2])*p-y[,2]
    w <- (y[,1]+y[,2])*p*(1-p)
    ywk <- eta-u/w-offset
    wt <- w*wt
    list(ywk=ywk,wt=wt)
  }
  else {
    if (min(y)<0)
      stop("gss error: negative binomial response should be nonnegative")
    p <- 1-1/(1+exp(eta))
    if (is.null(alpha)) log.alpha <- log(mean(y*exp(eta)))
    else log.alpha <- log(alpha)
    repeat {
      alpha <- exp(log.alpha)
      ua <- sum(digamma(y+alpha)-digamma(alpha)+log(p))*alpha
      wa <- sum(trigamma(y+alpha)-trigamma(alpha))*alpha*alpha+ua
      log.alpha.new <- log.alpha - ua/wa
      if (abs(log.alpha-log.alpha.new)/(1+abs(log.alpha))<1e-7) break
      log.alpha <- log.alpha.new
    }
    u <- (y+alpha)*p-alpha
    w <- (y+alpha)*p*(1-p)
    ywk <- eta-u/w-offset
    wt <- w*wt
    list(ywk=ywk,wt=wt,alpha=alpha)
  }
}

dev.resid.nbinomial <-
function(y,eta,wt) {
  if (is.null(wt)) wt <- rep(1,dim(y)[1])
  p <- 1-1/(1+exp(eta))
  as.vector(2*wt*(y[,1]*log(ifelse(y[,1]==0,1,y[,1]/(y[,1]+y[,2])/(1-p)))
                  +y[,2]*log(y[,2]/(y[,1]+y[,2])/p)))
}

dev.null.nbinomial <-
function(y,wt,offset) {
  if (is.null(wt)) wt <- rep(1,dim(y)[1])
  p <- sum(y[,2])/sum(y)
  if (!is.null(offset)) {
    eta <- log(p/(1-p)) - mean(offset)
    repeat {
      p <- 1-1/(1+exp(eta+offset))
      u <- (y[,1]+y[,2])*p-y[,2]
      w <- (y[,1]+y[,2])*p*(1-p)
      eta.new <- eta-sum(wt*u)/sum(wt*w)
      if (abs(eta-eta.new)/(1+abs(eta))<1e-7) break
      eta <- eta.new    
    }
  }
  sum(2*wt*(y[,1]*log(ifelse(y[,1]==0,1,y[,1]/(y[,1]+y[,2])/(1-p)))
            +y[,2]*log(y[,2]/(y[,1]+y[,2])/p)))
}

mkdata.binomial <- ## Make pseudo data for logistic regression
function(y,eta,wt,offset) {
  if (is.vector(y)) y <- as.matrix(y)
  if (is.null(wt)) wt <- rep(1,dim(y)[1])
  if (is.null(offset)) offset <- rep(0,dim(y)[1])
  if (dim(y)[2]==1) {
    if ((max(y)>1)|(min(y)<0))
      stop("gss error: binomial responses should be between 0 and 1")
  }
  else {
    if (min(y)<0)
      stop("gss error: paired binomial response should be nonnegative")
    wt <- wt * (y[,1]+y[,2])
    y <- y[,1]/(y[,1]+y[,2])
  }
  p <- 1-1/(1+exp(eta))
  u <- p - y
  w <- p*(1-p)
  ywk <- eta-u/w-offset
  wt <- w*wt
  list(ywk=ywk,wt=wt)
}

dev.resid.binomial <- ## Calculate the deviance residuals of logistic fit
function(y,eta,wt) {
  if (is.vector(y)) y <- as.matrix(y)
  if (is.null(wt)) wt <- rep(1,dim(y)[1])
  if (dim(y)[2]>1) {
    wt <- wt * (y[,1]+y[,2])
    y <- y[,1]/(y[,1]+y[,2])
  }
  p <- 1-1/(1+exp(eta))
  as.vector(2*wt*(y*log(ifelse(y==0,1,y/p))
                  +(1-y)*log(ifelse(y==1,1,(1-y)/(1-p)))))
}

dev.null.binomial <-
function(y,wt,offset) {
  if (is.vector(y)) y <- as.matrix(y)
  if (is.null(wt)) wt <- rep(1,dim(y)[1])
  if (dim(y)[2]>1) {
    wt <- wt * (y[,1]+y[,2])
    y <- y[,1]/(y[,1]+y[,2])
  }
  p <- sum(wt*y)/sum(wt)
  if (!is.null(offset)) {
    eta <- log(p/(1-p)) - mean(offset)
    repeat {
      p <- 1-1/(1+exp(eta+offset))
      u <- p - y
      w <- p*(1-p)
      eta.new <- eta-sum(wt*u)/sum(wt*w)
      if (abs(eta-eta.new)/(1+abs(eta))<1e-7) break
      eta <- eta.new    
    }
  }
  sum(2*wt*(y*log(ifelse(y==0,1,y/p))
            +(1-y)*log(ifelse(y==1,1,(1-y)/(1-p)))))
}

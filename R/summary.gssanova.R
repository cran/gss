summary.gssanova <- ## Summarize gssanova objects
function(obj,diagnostics=FALSE) {
  y <- model.response(obj$mf,"numeric")
  wt <- model.weights(obj$mf)
  offset <- model.offset(obj$mf)
  if (!is.null(obj$alpha)) y <- cbind(y,obj$alpha)
  dev.resid <- switch(obj$family,
                      binomial=dev.resid.binomial(y,obj$eta,wt),
                      nbinomial=dev.resid.nbinomial(y,obj$eta,wt),
                      poisson=dev.resid.poisson(y,obj$eta,wt),
                      inverse.gaussian=dev.resid.inverse.gaussian(y,obj$eta,wt),
                      Gamma=dev.resid.Gamma(y,obj$eta,wt))
  dev.null <- switch(obj$family,
                     binomial=dev.null.binomial(y,wt,offset),
                     nbinomial=dev.null.nbinomial(y,wt,offset),
                     poisson=dev.null.poisson(y,wt,offset),
                     inverse.gaussian=dev.null.inverse.gaussian(y,wt,offset),
                     Gamma=dev.null.Gamma(y,wt,offset))
  w <- obj$w
  if (is.null(offset)) offset <- rep(0,length(obj$eta))
  ## Calculate the summaries
  res <- 10^obj$nlambda*obj$c           # Residuals
  fitted <- obj$eta                     # Fitted values
  fitted.off <- fitted-offset
  sigma2 <- obj$varht                   # dispersion
  rss <- sum(res^2)                     # Residual sum of squares
  dev <- sum(dev.resid)
  penalty <- sum(obj$c*fitted.off*sqrt(w))
  penalty <-                            # Penalty associated with the fit
    as.vector(10^obj$nlambda*penalty)
  ## Calculate the diagnostics
  if (diagnostics) {
    comp <- NULL
    for (label in obj$terms$labels) {
      if (label=="1") next
      if (label=="offset") next
      comp <- cbind(comp,predict(obj,obj$mf,inc=label))
    }
    comp <- cbind(comp,yhat=fitted.off,y=fitted.off+res/sqrt(w),e=res/sqrt(w))
    term.label <- obj$terms$labels[obj$terms$labels!="1"]
    term.label <- term.label[term.label!="offset"]
    colnames(comp) <- c(term.label,"yhat","y","e")
    comp <- sqrt(w)*comp - outer(sqrt(w),apply(w*comp,2,sum))/sum(w)
    corr <- t(comp)%*%comp
    corr <- t(corr/sqrt(diag(corr)))/sqrt(diag(corr))
    norm <- apply(comp,2,function(x){sqrt(sum(x^2))})
    cosines <- rbind(corr[c("y","e"),],norm)
    rownames(cosines) <- c("cos.y","cos.e","norm")
    corr <- corr[term.label,term.label,drop=FALSE]
    if (qr(corr)$rank<dim(corr)[2]) kappa <- rep(Inf,len=dim(corr)[2])
    else kappa <- as.numeric(diag(solve(corr)))
    rough <- as.vector(10^obj$nlambda*t(comp[,term.label])%*%obj$c/penalty)
    names(kappa) <- names(rough) <- term.label
  }
  else kappa <- cosines <- rough <- NULL
  ## Return the summaries
  z <- list(call=obj$call,family=obj$family,method=obj$method,iter=obj$iter,
            fitted=fitted,dispersion=sigma2,residuals=res,rss=rss,
            deviance=dev,dev.resid=sqrt(dev.resid)*sign(res),
            dev.null=dev.null,penalty=penalty,
            kappa=kappa,cosines=cosines,roughness=rough)
  class(z) <- "summary.gssanova"
  z
}

summary.ssanova <- ## Summarize ssanova objects
function(obj,diagnostics=FALSE) {
  y <- model.response(obj$mf,"numeric")
  w <- model.weights(obj$mf)
  offset <- model.offset(obj$mf)
  if (is.null(offset)) offset <- rep(0,length(obj$c))
  ## Calculate the summaries
  res <- 10^obj$nlambda*obj$c           # Residuals
  if (!is.null(w)) res <- res/sqrt(w)
  fitted <- as.numeric(y-res)           # Fitted values
  fitted.off <- fitted-offset
  sigma <- sqrt(obj$varht)              # (estimated) sigma
  if (!is.null(w)) {                    # R^2
    r.squared <- sum(w*(fitted-sum(w*fitted)/sum(w))^2)
    r.squared <- r.squared/sum(w*(y-sum(w*y)/sum(w))^2)
  }
  else r.squared <- var(fitted)/var(y)       
  if (is.null(w)) rss <- sum(res^2)     # Residual sum of squares
  else rss <- sum(w*res^2)
  if (is.null(w))                       # Penalty associated with the fit
    penalty <- sum(obj$c*fitted.off)
  else penalty <- sum(obj$c*fitted.off*sqrt(w))
  penalty <- as.vector(10^obj$nlambda*penalty)
  ## Calculate the diagnostics
  if (diagnostics) {
    comp <- NULL
    for (label in obj$terms$labels) {
      if (label=="1") next
      if (label=="offset") next
      comp <- cbind(comp,predict(obj,obj$mf,inc=label))
    }
    comp <- cbind(comp,yhat=fitted.off,y=fitted.off+res,e=res)
    term.label <- obj$terms$labels[obj$terms$labels!="1"]
    term.label <- term.label[term.label!="offset"]
    colnames(comp) <- c(term.label,"yhat","y","e")
    if (!is.null(w))
      comp <- sqrt(w)*comp - outer(sqrt(w),apply(w*comp,2,sum))/sum(w)
    else comp <- sweep(comp,2,apply(comp,2,mean))
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
  z <- list(call=obj$call,method=obj$method,fitted=fitted,residuals=res,
            sigma=sigma,r.squared=r.squared,rss=rss,penalty=penalty,
            kappa=kappa,cosines=cosines,roughness=rough)
  class(z) <- "summary.ssanova"
  z
}

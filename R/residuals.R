residuals.ssanova <-
function(obj) {
  w <- model.weights(obj$mf)
  res <- 10^obj$nlambda*obj$c
  if (!is.null(w)) res <- res/sqrt(w)
  res
}

residuals.gssanova <-
function(obj,type="working") {
  res <- 10^obj$nlambda*obj$c
  if (!is.na(charmatch(type,"deviance"))) {
    y <- model.response(obj$mf,"numeric")
    wt <- model.weights(obj$mf)
    dev.resid <- switch(obj$family,
                        binomial=dev.resid.binomial(y,obj$eta,wt),
                        poisson=dev.resid.poisson(y,obj$eta,wt))
    res <- sqrt(dev.resid)*sign(res)
  }
  res
}

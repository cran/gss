fitted.ssanova <-
function(obj) {
  y <- model.response(obj$mf,"numeric")
  w <- model.weights(obj$mf)
  res <- 10^obj$nlambda*obj$c
  if (!is.null(w)) res <- res/sqrt(w)
  as.numeric(y-res)
}

fitted.gssanova <-
function(obj) obj$eta

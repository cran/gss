## Obtain fitted values from ssanova objects
fitted.ssanova <- function(obj)
{
    y <- model.response(obj$mf,"numeric")
    w <- model.weights(obj$mf)
    res <- 10^obj$nlambda*obj$c
    if (!is.null(w)) res <- res/sqrt(w)
    as.numeric(y-res)
}

## Obtain fitted values in working scale from gssanova objects
fitted.gssanova <- function(obj) obj$eta

## Obtain residuals from ssanova objects
residuals.ssanova <- function(obj)
{
    w <- model.weights(obj$mf)
    res <- 10^obj$nlambda*obj$c
    if (!is.null(w)) res <- res/sqrt(w)
    res
}

## Obtain residuals from gssanova objects
residuals.gssanova <- function(obj,type="working")
{
    res <- 10^obj$nlambda*obj$c/sqrt(obj$w)
    if (!is.na(charmatch(type,"deviance"))) {
        y <- model.response(obj$mf,"numeric")
        wt <- model.weights(obj$mf)
        dev.resid <- switch(obj$family,
                            binomial=dev.resid.binomial(y,obj$eta,wt),
                            poisson=dev.resid.poisson(y,obj$eta,wt),
                            poisson=dev.resid.poisson(y,obj$eta,wt),
                            inverse.gaussian=dev.resid.inverse.gaussian(y,obj$eta,wt),
                            Gamma=dev.resid.Gamma(y,obj$eta,wt),
                            weibull=dev.resid.weibull(y,obj$eta,wt,obj$alpha),
                            lognorm=dev.resid.lognorm(y,obj$eta,wt,obj$alpha),
                            loglogis=dev.resid.loglogis(y,obj$eta,wt,obj$alpha))
        res <- sqrt(dev.resid)*sign(res)
    }
    res
}

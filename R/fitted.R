## Obtain fitted values from ssanova objects
fitted.ssanova <- function(object,...)
{
    predict(object,object$mf)
}

## Obtain residuals from ssanova objects
residuals.ssanova <- function(object,...)
{
    y <- model.response(object$mf,"numeric")
    as.numeric(y-predict(object,object$mf))
}

## Obtain fitted values in working scale from gssanova objects
fitted.gssanova <- function(object,...) object$eta

## Obtain residuals from gssanova objects
residuals.gssanova <- function(object,type="working",...)
{
    res <- 10^object$nlambda*object$c/sqrt(object$w)
    if (!is.na(charmatch(type,"deviance"))) {
        y <- model.response(object$mf,"numeric")
        wt <- model.weights(object$mf)
        dev.resid <- switch(object$family,
                            binomial=dev.resid.binomial(y,object$eta,wt),
                            poisson=dev.resid.poisson(y,object$eta,wt),
                            poisson=dev.resid.poisson(y,object$eta,wt),
                            inverse.gaussian=dev.resid.inverse.gaussian(y,object$eta,wt),
                            Gamma=dev.resid.Gamma(y,object$eta,wt),
                            weibull=dev.resid.weibull(y,object$eta,wt,object$alpha),
                            lognorm=dev.resid.lognorm(y,object$eta,wt,object$alpha),
                            loglogis=dev.resid.loglogis(y,object$eta,wt,object$alpha))
        res <- sqrt(dev.resid)*sign(res)
    }
    res
}

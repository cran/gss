\name{mkdata.binomial}
\alias{mkdata.binomial}
\alias{dev.resid.binomial}
\alias{dev.null.binomial}
\alias{mkdata.poisson}
\alias{dev.resid.poisson}
\alias{dev.null.poisson}
\alias{mkdata.Gamma}
\alias{dev.resid.Gamma}
\alias{dev.null.Gamma}
\alias{mkdata.inverse.gaussian}
\alias{dev.resid.inverse.gaussian}
\alias{dev.null.inverse.gaussian}
\alias{mkdata.nbinomial}
\alias{dev.resid.nbinomial}
\alias{dev.null.nbinomial}
\alias{mkdata.weibull}
\alias{dev.resid.weibull}
\alias{dev.null.weibull}
\alias{mkdata.lognorm}
\alias{dev.resid.lognorm}
\alias{dev.null.lognorm}
\alias{mkdata.loglogis}
\alias{dev.resid.loglogis}
\alias{dev.null.loglogis}
\title{Utility Functions for Error Families}
\description{
    Utility functions for fitting Smoothing Spline ANOVA models with
    non-Gaussian responses.
}
\usage{
mkdata.binomial(y, eta, wt, offset)
dev.resid.binomial(y, eta, wt)
dev.null.binomial(y, wt, offset)

mkdata.poisson(y, eta, wt, offset)
dev.resid.poisson(y, eta, wt)
dev.null.poisson(y, wt, offset)

mkdata.Gamma(y, eta, wt, offset)
dev.resid.Gamma(y, eta, wt)
dev.null.Gamma(y, wt, offset)

mkdata.inverse.gaussian(y, eta, wt, offset)
dev.resid.inverse.gaussian(y, eta, wt)
dev.null.inverse.gaussian(y, wt, offset)

mkdata.nbinomial(y, eta, wt, offset, alpha)
dev.resid.nbinomial(y, eta, wt)
dev.null.nbinomial(y, wt, offset)

mkdata.weibull(y, eta, wt, offset, alpha)
dev.resid.weibull(y, eta, wt, alpha)
dev.null.weibull(y, wt, offset, alpha)

mkdata.lognorm(y, eta, wt, offset, alpha)
dev.resid.lognorm(y, eta, wt, alpha)
dev.null.lognorm(y, wt, offset, alpha)

mkdata.loglogis(y, eta, wt, offset, alpha)
dev.resid.loglogis(y, eta, wt, alpha)
dev.null.loglogis(y, wt, offset, alpha)
}
\arguments{
    \item{y}{Model response.}
    \item{eta}{Fitted values on link scale.}
    \item{wt}{Model weights.}
    \item{offset}{Model offset.}
    \item{alpha}{Size for nbinomial.  Inverse scale for log life time.}
}
\details{
    These are not to be called by the user.

    \code{mkdata.x} create the pseudo data to be used in iterated
    penalized least squares fitting.  \code{dev.resid.x} calculate the
    deviance residuals.  \code{dev.null.x} calculate the deviance of the
    constant null model.
}
\seealso{
    \code{\link{gssanova}}.
}
\keyword{internal}

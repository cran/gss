\name{family}
\alias{mkdata.binomial}
\alias{dev.resid.binomial}
\alias{dev.null.binomial}
\alias{cv.binomial}
\alias{y0.binomial}
\alias{proj0.binomial}
\alias{kl.binomial}
\alias{cfit.binomial}
\alias{mkdata.poisson}
\alias{dev.resid.poisson}
\alias{dev.null.poisson}
\alias{cv.poisson}
\alias{y0.poisson}
\alias{proj0.poisson}
\alias{kl.poisson}
\alias{cfit.poisson}
\alias{mkdata.Gamma}
\alias{dev.resid.Gamma}
\alias{dev.null.Gamma}
\alias{cv.Gamma}
\alias{y0.Gamma}
\alias{proj0.Gamma}
\alias{kl.Gamma}
\alias{cfit.Gamma}
\alias{mkdata.inverse.gaussian}
\alias{dev.resid.inverse.gaussian}
\alias{dev.null.inverse.gaussian}
\alias{cv.inverse.gaussian}
\alias{y0.inverse.gaussian}
\alias{proj0.inverse.gaussian}
\alias{kl.inverse.gaussian}
\alias{cfit.inverse.gaussian}
\alias{mkdata.nbinomial}
\alias{dev.resid.nbinomial}
\alias{dev.null.nbinomial}
\alias{cv.nbinomial}
\alias{y0.nbinomial}
\alias{proj0.nbinomial}
\alias{kl.nbinomial}
\alias{cfit.nbinomial}
\alias{mkdata.weibull}
\alias{dev.resid.weibull}
\alias{dev.null.weibull}
\alias{cv.weibull}
\alias{y0.weibull}
\alias{proj0.weibull}
\alias{kl.weibull}
\alias{cfit.weibull}
\alias{mkdata.lognorm}
\alias{dev.resid.lognorm}
\alias{dev0.resid.lognorm}
\alias{dev.null.lognorm}
\alias{cv.lognorm}
\alias{y0.lognorm}
\alias{proj0.lognorm}
\alias{kl.lognorm}
\alias{cfit.lognorm}
\alias{mkdata.loglogis}
\alias{dev.resid.loglogis}
\alias{dev0.resid.loglogis}
\alias{dev.null.loglogis}
\alias{cv.loglogis}
\alias{y0.loglogis}
\alias{proj0.loglogis}
\alias{kl.loglogis}
\alias{cfit.loglogis}
\title{Utility Functions for Error Families}
\description{
    Utility functions for fitting Smoothing Spline ANOVA models with
    non-Gaussian responses.
}
\usage{
mkdata.binomial(y, eta, wt, offset)
dev.resid.binomial(y, eta, wt)
dev.null.binomial(y, wt, offset)
cv.binomial(y, eta, wt, hat, alpha)
y0.binomial(y, eta0, wt)
proj0.binomial(y0, eta, offset)
kl.binomial(eta0, eta1, wt)
cfit.binomial(y, wt, offset)

mkdata.poisson(y, eta, wt, offset)
dev.resid.poisson(y, eta, wt)
dev.null.poisson(y, wt, offset)
cv.poisson(y, eta, wt, hat, alpha, sr, q)
y0.poisson(eta0)
proj0.poisson(y0, eta, wt, offset)
kl.poisson(eta0, eta1, wt)
cfit.poisson(y, wt, offset)

mkdata.Gamma(y, eta, wt, offset)
dev.resid.Gamma(y, eta, wt)
dev.null.Gamma(y, wt, offset)
cv.Gamma(y, eta, wt, hat, rss, alpha)
y0.Gamma(eta0)
proj0.Gamma(y0, eta, wt, offset)
kl.Gamma(eta0, eta1, wt)
cfit.Gamma(y, wt, offset)

mkdata.inverse.gaussian(y, eta, wt, offset)
dev.resid.inverse.gaussian(y, eta, wt)
dev.null.inverse.gaussian(y, wt, offset)
cv.inverse.gaussian(y, eta, wt, hat, rss, alpha)
y0.inverse.gaussian(eta0)
proj0.inverse.gaussian(y0, eta, wt, offset)
kl.inverse.gaussian(eta0, eta1, wt)
cfit.inverse.gaussian(y, wt, offset)

mkdata.nbinomial(y, eta, wt, offset, nu)
dev.resid.nbinomial(y, eta, wt)
dev.null.nbinomial(y, wt, offset)
cv.nbinomial(y, eta, wt, hat, alpha)
y0.nbinomial(y,eta0,nu)
proj0.nbinomial(y0, eta, wt, offset)
kl.nbinomial(eta0, eta1, wt, nu)
cfit.nbinomial(y, wt, offset, nu)

mkdata.weibull(y, eta, wt, offset, nu)
dev.resid.weibull(y, eta, wt, nu)
dev.null.weibull(y, wt, offset, nu)
cv.weibull(y, eta, wt, hat, nu, alpha)
y0.weibull(y, eta0, nu)
proj0.weibull(y0, eta, wt, offset, nu)
kl.weibull(eta0, eta1, wt, nu, int)
cfit.weibull(y, wt, offset, nu)

mkdata.lognorm(y, eta, wt, offset, nu)
dev.resid.lognorm(y, eta, wt, nu)
dev0.resid.lognorm(y, eta, wt, nu)
dev.null.lognorm(y, wt, offset, nu)
cv.lognorm(y, eta, wt, hat, nu, alpha)
y0.lognorm(y, eta0, nu)
proj0.lognorm(y0, eta, wt, offset, nu)
kl.lognorm(eta0, eta1, wt, nu, y0)
cfit.lognorm(y, wt, offset, nu)

mkdata.loglogis(y, eta, wt, offset, nu)
dev.resid.loglogis(y, eta, wt, nu)
dev0.resid.loglogis(y, eta, wt, nu)
dev.null.loglogis(y, wt, offset, nu)
cv.loglogis(y, eta, wt, hat, nu, alpha)
y0.loglogis(y, eta0, nu)
proj0.loglogis(y0, eta, wt, offset, nu)
kl.loglogis(eta0, eta1, wt, nu, y0)
cfit.loglogis(y, wt, offset, nu)
}
\arguments{
    \item{y}{Model response.}
    \item{eta}{Fitted values on link scale.}
    \item{wt}{Model weights.}
    \item{offset}{Model offset.}
    \item{nu}{Size for nbinomial.  Inverse scale for log life time.}
}
\note{
    \code{\link{gssanova0}} uses \code{mkdata.x}, \code{dev.resid.x},
    and \code{dev.null.x}.  \code{\link{gssanova}} uses the above plus
    \code{dev0.resid.x} and \code{cv.x}.

    \code{y0.x}, \code{proj0.x}, \code{kl.x}, and \code{cfit.x} are used
    by \code{\link{project.gssanova}}.
}
\keyword{internal}

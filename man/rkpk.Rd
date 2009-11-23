\name{rkpk}
\alias{sspreg1}
\alias{mspreg1}
\alias{sspngreg}
\alias{mspngreg}
\alias{ngreg}
\alias{regaux}
\alias{ngreg.proj}
\title{Numerical Engine for ssanova and gssanova}
\description{
    Perform numerical calculations for the \code{\link{ssanova}} and
    \code{\link{gssanova}} suites.
}
\usage{
sspreg1(s, r, q, y, method, alpha, varht, random)
mspreg1(s, r, id.basis, y, method, alpha, varht, random, skip.iter)

sspngreg(family, s, r, q, y, wt, offset, alpha, nu, random)
mspngreg(family, s, r, id.basis, y, wt, offset, alpha, nu, random, skip.iter)
ngreg(dc, family, sr, q, y, wt, offset, nu, alpha)

regaux(s, r, q, nlambda, fit)

ngreg.proj(dc, family, sr, q, y0, wt, offset, nu)
}
\details{
    \code{sspreg1} is used by \code{\link{ssanova}} to compute
    regression estimates with a single smoothing parameter.
    \code{mspreg1} is used by \code{\link{ssanova}} to compute
    regression estimates with multiple smoothing parameters.

    \code{ssngpreg} is used by \code{\link{gssanova}} to compute
    non-Gaussian regression estimates with a single smoothing
    parameter.  \code{mspngreg} is used by \code{\link{gssanova}} to
    compute non-Gaussian regression estimates with multiple smoothing
    parameters.  \code{ngreg} is used by \code{ssngpreg} and
    \code{mspngreg} to perform Newton iteration with fixed smoothing
    parameters and to calculate cross-validation scores on return.

    \code{regaux} is used by \code{sspreg1}, \code{mspreg1},
    \code{ssngpreg}, and \code{mspngreg} to obtain auxiliary information
    needed by \code{predict.ssanova} for standard error calculation.

    \code{ngreg.proj} is used by \code{\link{project.gssanova}} to
    calculate the Kullback-Leibler projection for non-Gaussian
    regression.
}
\arguments{
    \item{family}{Description of the error distribution.  Supported
	are exponential families \code{"binomial"}, \code{"poisson"},
	\code{"Gamma"}, and \code{"nbinomial"}.  Also supported are
	accelerated life model families \code{"weibull"},
	\code{"lognorm"}, and \code{"loglogis"}.}
    \item{s}{Unpenalized terms evaluated at data points.}
    \item{r}{Basis of penalized terms evaluated at data points.}
    \item{q}{Penalty matrix.}
    \item{id.basis}{Index of observations to be used as "knots."}
    \item{y}{Response vector.}
    \item{wt}{Model weights.}
    \item{offset}{Model offset.}
    \item{method}{\code{"v"} for GCV, \code{"m"} for GML, or \code{"u"}
        for Mallows' CL.}
    \item{alpha}{Parameter modifying GCV or Mallows' CL scores for
        smoothing parameter selection.}
    \item{nu}{Optional argument for future support of nbinomial,
      weibull, lognorm, and loglogis families.}
    \item{varht}{External variance estimate needed for \code{method="u"}.}
    \item{random}{Input for parametric random effects in nonparametric
        mixed-effect models.}
    \item{skip.iter}{Flag indicating whether to use initial values of
      theta and skip theta iteration.}

    \item{nlambda}{Smoothing parameter in effect.}
    \item{fit}{Fitted model.}

    \item{dc}{Coefficients of fits.}
    \item{sr}{\code{cbind(s,r)}.}
    \item{y0}{Components of the fit to be projected.}
}
\references{
    Kim, Y.-J. and Gu, C. (2004), Smoothing spline Gaussian regression:
    more scalable computation via efficient approximation.
    \emph{Journal of the Royal Statistical Society, Ser. B}, \bold{66},
    337--356.
}
\keyword{internal}

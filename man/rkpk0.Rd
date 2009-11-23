\name{rkpk0}
\alias{sspreg0}
\alias{mspreg0}
\alias{getcrdr}
\alias{getsms}
\alias{sspregpoi}
\alias{mspregpoi}
\title{Interface to RKPACK}
\description{
    Call RKPACK routines for numerical calculations concerning the
    \code{\link{ssanova0}} and \code{\link{gssanova0}} suites.
}
\usage{
sspreg0(s, q, y, method="v", varht=1)
mspreg0(s, q, y, method="v", varht=1, prec=1e-7, maxiter=30)
sspregpoi(family, s, q, y, wt, offset, method="u", varht=1, nu, prec=1e-7, maxiter=30)
mspregpoi(family, s, q, y, wt, offset, method="u", varht=1, nu, prec=1e-7, maxiter=30)
getcrdr(obj, r)
getsms(obj)
}
\details{
    \code{sspreg0} is used by \code{\link{ssanova0}} to fit Gaussian
    models with a single smoothing parameter.  \code{mspreg0} is used to
    fit Gaussian models with multiple smoothing parameters.

    \code{sspregpoi} is used by \code{\link{gssanova0}} to fit non
    Gaussian models with a single smoothing parameter.  \code{mspregpoi}
    is used to fit non Gaussian models with multiple smoothing
    parameters.

    \code{getcrdr} and \code{getsms} are used by
    \code{\link{predict.ssanova0}} to calculate standard errors of the
    fitted terms.
}
\arguments{
    \item{s}{Design matrix of unpenalized terms.}
    \item{q}{Penalty matrices of penalized terms.}
    \item{y}{Model response.}
    \item{method}{Method for smoothing parameter selection.}
    \item{varht}{Assumed dispersion parameter, needed only for
	\code{method="u"}.}
    \item{prec}{Precision requirement for iterations.}
    \item{maxiter}{Maximum number of iterations allowed.}
    \item{family}{Error family.}	
    \item{wt}{Model weights.}
    \item{offset}{Model offset.}
    \item{obj}{Object returned from a call to \code{sspreg},
      \code{mspreg}, \code{sspregpoi}, or \code{mspregpoi}.}
    \item{nu}{Optional argument for nbinomial, weibull, lognorm, and
      loglogis families.}
    \item{r}{Inputs for standard error calculation.}
}
\references{
    Gu, C. (1989), RKPACK and its applications: Fitting smoothing spline
    models.  In \emph{ASA Proceedings of Statistical Computing Section},
    pp. 42--51.

    Gu, C. (1992), Cross validating non Gaussian data.  \emph{Journal of
    Computational and Graphical Statistics}, \bold{1}, 169--179. 
}
\keyword{internal}

\name{rkpk}
\alias{sspreg}
\alias{mspreg}
\alias{getcrdr}
\alias{getsms}
\alias{sspregpoi}
\alias{mspregpoi}
\title{Interface to RKPACK}
\description{
    Call RKPACK routines for numerical calculations for fitting and
    predicting from Smoothing Spline ANOVA models.
}
\usage{
sspreg(s, q, y, method="v", varht=1)
mspreg(s, q, y, method="v", varht=1, prec=1e-7, maxiter=30)
sspregpoi(family, s, q, y, wt, offset, method="u",
          varht=1, prec=1e-7, maxiter=30)
mspregpoi(family, s, q, y, wt, offset, method="u",
          varht=1, prec=1e-7, maxiter=30)
getcrdr(obj, r)
getsms(obj)
}
\details{
    \code{sspreg} is used by \code{\link{ssanova}} to fit Gaussian
    models with a single smoothing parameter.  \code{mspreg} is used to
    fit Gaussian models with multiple single smoothing parameters.

    \code{sspregpoi} is used by \code{\link{gssanova}} to fit non
    Gaussian models with a single smoothing parameter.  \code{mspregpoi}
    is used to fit non Gaussian models with multiple single smoothing
    parameters.

    \code{getcrdr} and \code{getsms} are used by
    \code{\link{predict.ssanova}} to calculate standard errors of the
    fitted terms.
}
\arguments{
    \item{s}{Design matrix of fixed effects.}
    \item{q}{Inverse variance-covariance matrices of random effects.}
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
}
\references{
    Gu, C. (1989), RKPACK and its applications: Fitting smoothing spline
    models.  In \emph{ASA Proceedings of Statistical Computing Section},
    pp. 42--51.
}
\keyword{gss-internal}

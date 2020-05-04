\name{hzdrate.sshzd}
\alias{hzdrate.sshzd}
\alias{hzdcurve.sshzd}
\alias{survexp.sshzd}
\title{Evaluating Smoothing Spline Hazard Estimates}
\description{
    Evaluate smoothing spline hazard estimates by \code{sshzd}.
}
\usage{
hzdrate.sshzd(object, x, se=FALSE, include=c(object$terms$labels,object$lab.p))
hzdcurve.sshzd(object, time, covariates=NULL, se=FALSE)
survexp.sshzd(object, time, covariates=NULL, start=0)
}
\arguments{
    \item{object}{Object of class \code{"sshzd"}.}
    \item{x}{Data frame or vector of points on which hazard is to be
      evaluated.}
    \item{se}{Flag indicating if standard errors are required.}
    \item{include}{List of model terms to be included in the
	evaluation.}
    \item{time}{Vector of time points.}
    \item{covariates}{Vector of covariate values.}
    \item{start}{Optional starting times of the intervals.}
}
\value{
    For \code{se=FALSE}, \code{hzdrate.sshzd} returns a vector of hazard
    evaluations, and \code{hzdcurve.sshzd} returns a vector or columns of
    hazard curve(s) evaluated on \code{time} points at the
    \code{covariates} values.  For \code{se=TRUE}, \code{hzdrate.sshzd}
    and \code{hzdcurve.sshzd} return a list consisting of the following
    elements.
    \item{fit}{Vector or columns of hazard.}
    \item{se.fit}{Vector or columns of standard errors for log hazard.}

    \code{survexp.sshzd} returns a vector or columns of expected
    survivals based on the cumulative hazards over (\code{start},
    \code{time}) at the \code{covariates} values, which in fact are the
    (conditional) survival probabilities \eqn{S(time)/S(start)}.
}
\note{
    For left-truncated data, \code{start} must be at or after the
    earliest truncation point.
}
\seealso{
    Fitting function \code{\link{sshzd}}.
}
\keyword{models}
\keyword{smooth}
\keyword{survival}

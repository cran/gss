\name{hzdrate.sshzd}
\alias{hzdrate.sshzd}
\alias{hzdcurve.sshzd}
\alias{survexp.sshzd}
\title{Evaluating Smoothing Spline Hazard Estimates}
\description{
    Evaluate smoothing spline hazard estimates by \code{sshzd}.
}
\usage{
hzdrate.sshzd(object, x, se=FALSE)
hzdcurve.sshzd(object, time, covariates=NULL, se=FALSE)
survexp.sshzd(object, time, covariates=NULL, start=0)
}
\arguments{
    \item{object}{Object of class \code{"sshzd"}.}
    \item{x}{Data frame or vector of points on which hazard is to be
      evaluated.}
    \item{se}{Flag indicating if standard errors are required.}
    \item{time}{Vector of time points.}
    \item{covariates}{Vector of covariate values.}
    \item{start}{Optional starting times of the intervals.}
}
\value{
    For \code{se=FALSE}, \code{hzdrate.sshzd} returns a vector of hazard
    evaluations.  For \code{se=TRUE}, \code{hzdrate.sshzd} returns a
    list consisting of the following components.
    \item{fit}{Vector of hazard evaluations.}
    \item{se.fit}{Vector of standard errors for log hazard.}
    
    For \code{se=FALSE}, \code{hzdcurve.sshzd} returns a vector or
    columns of hazard curve(s) evaluated on \code{time} points at the
    \code{covariates} values.  For \code{se=TRUE}, \code{hzdrate.sshzd}
    returns a list consisting of the following components.
    \item{fit}{Vector or columns of hazard curve(s).}
    \item{se.fit}{Vector or columns of standard errors for log hazard
      curve(s).}

    \code{survexp.sshzd} returns a vector of expected survivals, which
    are based on the cumulative hazards over the specified intervals.
}
\seealso{
    \code{\link{sshzd}}.
}
\keyword{models}
\keyword{smooth}
\keyword{survival}

\name{hzdrate.sshzd}
\alias{hzdrate.sshzd}
\alias{hzdcurve.sshzd}
\alias{survexp.sshzd}
\title{Evaluating Smoothing Spline Hazard Estimates}
\description{
    Evaluate smoothing spline hazard estimates by \code{sshzd}.
}
\usage{
hzdrate.sshzd(object, x)
hzdcurve.sshzd(object,time,covariates=NULL)
survexp.sshzd(object,time,covariates=NULL,start=0)
}
\arguments{
    \item{object}{Object of class \code{"sshzd"}.}
    \item{x}{Data frame or vector of points on which hazard is to be
        evaluated.}
    \item{time}{Vector of time points.}
    \item{covariates}{Vector of covariate values.}
    \item{start}{Optional starting times of the intervals.}
}
\value{
    \code{hzdrate.sshzd} returns a vector of hazard evaluations.

    \code{hzdcurve.sshzd} returns a vector or columns of hazard curve(s)
    evaluated on \code{time} points at the \code{covariates} values.

    \code{survexp.sshzd} returns a vector of expected survivals, which
    are based on the cumulative hazards over the specified intervals.
}
\seealso{
    \code{\link{sshzd}}.
}
\keyword{models}
\keyword{smooth}
\keyword{survival}

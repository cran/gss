\name{hzdrate.sshzd2d}
\alias{hzdrate.sshzd2d}
\alias{survexp.sshzd2d}
\title{Evaluating 2-D Smoothing Spline Hazard Estimates}
\description{
    Evaluate 2-D smoothing spline hazard estimates by \code{sshzd2d}.
}
\usage{
hzdrate.sshzd2d(object, time, covariates=NULL)
survexp.sshzd2d(object, time, covariates=NULL, job=3) 
}
\arguments{
    \item{object}{Object of class \code{"sshzd2d"}.}
    \item{time}{Matrix or vector of time points on which hazard or
        survival function is to be evaluated.}
    \item{covariates}{Data frame of covariate values.}
    \item{job}{Flag indicating which survival function to evaluate.}
}
\value{
    A vector of hazard or survival values.
}
\note{
    For \code{job=1,2}, \code{survexp.sshzd2d} returns marginal survival
    \eqn{S1(t)} or \eqn{S2(t)}.  For \code{job=3},
    \code{survexp.sshzd2d} returns the 2-D survival \eqn{S(t1,t2)}.

    For \code{hzdrate.sshzd2d} and \code{survexp.sshzd2d} with
    \code{job=3}, \code{time} should be a matrix of two columns.  For
    \code{survexp.sshzd2d} with \code{job=1,2}, \code{time} should be a
    vector.

    When \code{covariates} is present, its length should be either 1 or
    that of \code{time}.
}
\seealso{
    Fitting function \code{\link{sshzd2d}}.
}
\keyword{models}
\keyword{smooth}
\keyword{survival}

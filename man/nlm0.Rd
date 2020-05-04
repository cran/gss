\name{nlm0}
\alias{nlm0}
\title{Minimizing Univariate Functions on Finite Intervals}
\description{
   Minimize univariate functions on finite intervals using 3-point
   quadratic fit, with golden-section safe-guard.
}
\usage{
nlm0(fun, range, prec=1e-7)
}
\arguments{
    \item{fun}{Function to be minimized.}
    \item{range}{Interval on which the function to be minimized.}
    \item{prec}{Desired precision of the solution.}
}
\value{
    \code{nlm0} returns a list object with the following elements.
    \item{estimate}{Minimizer.}
    \item{minimum}{Minimum.}
    \item{evaluations}{Number of function evaluations.}
}
\keyword{math}

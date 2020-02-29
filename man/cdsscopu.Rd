\name{cdsscopu}
\alias{cdsscopu}
\alias{cpsscopu}
\alias{cqsscopu}
\title{Evaluating 1-D Conditional PDF, CDF, and Quantiles of Copula
    Density Estimates}
\description{
    Evaluate conditional pdf, cdf, and quantiles of copula density
    estimates.
}
\usage{
cdsscopu(object, x, cond, pos=1, int=NULL)
cpsscopu(object, q, cond, pos=1)
cqsscopu(object, p, cond, pos=1)
}
\arguments{
    \item{object}{Object of class \code{"sscopu"}.}
    \item{x}{Vector of points on which conditional pdf is to be
        evaluated.}
    \item{cond}{Value of conditioning variables.}
    \item{pos}{Position of variable of interest.}
    \item{int}{Normalizing constant.}
    \item{q}{Vector of points on which conditional cdf is to be
        evaluated.}
    \item{p}{Vector of probabilities for which conditional quantiles are
        to be calculated.}
}
\value{
    A vector of conditional pdf, cdf, or quantiles.
}
\seealso{
    Fitting functions \code{\link{sscopu}} and \code{\link{sscopu2}},
    and \code{\link{dsscopu}}.
}
\keyword{models}
\keyword{distribution}
\keyword{smooth}

\name{dssden}
\alias{dssden}
\alias{pssden}
\alias{qssden}
\alias{d.ssden}
\alias{d.ssden1}
\title{Evaluating PDF, CDF, and Quantiles of Smoothing Spline Density Estimates}
\description{
    Evaluate pdf, cdf, and quantiles for smoothing spline density
    estimates.
}
\usage{
dssden(object, x)
pssden(object, q)
qssden(object, p)
d.ssden(object, x)
d.ssden1(object, x)
}
\arguments{
    \item{object}{Object of class \code{"ssden"}.}
    \item{x}{Data frame or vector of points on which density is to be
        evaluated.}
    \item{q}{Vector of points on which cdf is to be evaluated.}
    \item{p}{Vector of probabilities for which quantiles are to be
        calculated.}
}
\value{
    A vector of pdf, cdf, or quantiles.
}
\details{
    The argument \code{x} in \code{dssden} is of the same form as the
    argument \code{newdata} in \code{\link{predict.lm}}, but can take a
    vector for 1-D densities.

    \code{pssden} and \code{qssden} naturally only work for 1-D
    densities.
}
\seealso{
    Fitting function \code{\link{ssden}} and \code{\link{cdssden}}.
}
\keyword{models}
\keyword{distribution}
\keyword{smooth}

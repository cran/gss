\name{dsscden}
\alias{dsscden}
\alias{psscden}
\alias{qsscden}
\alias{d.sscden}
\alias{d.sscden1}
\title{Evaluating PDF, CDF, and Quantiles of Smoothing Spline
    Conditional Density Estimates}
\description{
    Evaluate pdf, cdf, and quantiles for smoothing spline conditional
    density estimates.
}
\usage{
dsscden(object, y, x)
psscden(object, q, x)
qsscden(object, p, x)
d.sscden(object, x, y)
d.sscden1(object, x, y, scale=TRUE)
}
\arguments{
    \item{object}{Object of class \code{"sscden"} or \code{"sscden1"}.}
    \item{x}{Data frame of x values on which conditional density f(y|x)
        is to be evaluated.}
    \item{y}{Data frame or vector of points on which conditional density
        f(y|x) is to be evaluated.}
    \item{q}{Vector of points on which cdf is to be evaluated.}
    \item{p}{Vector of probabilities for which quantiles are to be
        calculated.}
    \item{scale}{Flag indicating whether to use approximate scaling
        without quadrature.}
}
\value{
    A matrix or vector of pdf, cdf, or quantiles of f(y|x), with each
    column corresponding to a distinct x value.
}
\details{
    The arguments \code{x} and \code{y} are of the same form as the
    argument \code{newdata} in \code{\link{predict.lm}}, but \code{y} in
    \code{dsscden} can take a vector for 1-D responses.

    \code{psscden} and \code{qsscden} naturally only work for 1-D
    responses.
}
\seealso{
    Fitting function \code{\link{sscden}} and \code{\link{cdsscden}}.
}
\keyword{models}
\keyword{distribution}
\keyword{smooth}

\name{cdsscden}
\alias{cdsscden}
\alias{cpsscden}
\alias{cqsscden}
\title{Evaluating Conditional PDF, CDF, and Quantiles of Smoothing
   Spline Conditional Density Estimates}
\description{
    Evaluate conditional pdf, cdf, and quantiles of f(y1|x,y2) for
    smoothing spline conditional density estimates f(y|x).
}
\usage{
cdsscden(object, y, x, cond, int=NULL)
cpsscden(object, q, x, cond)
cqsscden(object, p, x, cond)
}
\arguments{
    \item{object}{Object of class \code{"sscden"} or \code{"sscden1"}.}
    \item{x}{Data frame of x values on which conditional density
        f(y1|x,y2) is to be evaluated.}
    \item{y}{Data frame or vector of y1 points on which conditional
      density f(y1|x,y2) is to be evaluated.}
    \item{cond}{One row data frame of conditioning variables y2.}
    \item{q}{Vector of points on which cdf is to be evaluated.}
    \item{p}{Vector of probabilities for which quantiles are to be
        calculated.}
    \item{int}{Vector of normalizing constants.}
}
\value{
    \code{cdsscden} returns a list object with the following
    components.
    \item{pdf}{Matrix or vector of conditional pdf f(y1|x,y2), with each
        column corresponding to a distinct x value.}
    \item{int}{Vector of normalizing constants.}
      
    \code{cpsscden} and \code{cqsscden} return a matrix or vector of
    conditional cdf or quantiles of f(y1|x,y2).
}
\details{
    The arguments \code{x} and \code{y} are of the same form as the
    argument \code{newdata} in \code{\link{predict.lm}}, but \code{y} in
    \code{cdsscden} can take a vector for 1-D y1.

    \code{cpsscden} and \code{cqsscden} naturally only work for 1-D y1.
}
\note{
    If variables other than factors or numerical vectors are involved in
    \code{y1}, the normalizing constants can not be computed.
}
\seealso{
    Fitting function \code{\link{sscden}} and \code{\link{dsscden}}.
}
\keyword{models}
\keyword{distribution}
\keyword{smooth}

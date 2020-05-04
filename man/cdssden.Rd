\name{cdssden}
\alias{cdssden}
\alias{cpssden}
\alias{cqssden}
\title{Evaluating Conditional PDF, CDF, and Quantiles of Smoothing
    Spline Density Estimates}
\description{
    Evaluate conditional pdf, cdf, and quantiles for smoothing spline
    density estimates.
}
\usage{
cdssden(object, x, cond, int=NULL)
cpssden(object, q, cond)
cqssden(object, p, cond)
}
\arguments{
    \item{object}{Object of class \code{"ssden"}.}
    \item{x}{Data frame or vector of points on which conditional density
        is to be evaluated.}
    \item{cond}{One row data frame of conditioning variables.}
    \item{int}{Normalizing constant.}
    \item{q}{Vector of points on which conditional cdf is to be
        evaluated.}
    \item{p}{Vector of probabilities for which conditional quantiles are
        to be calculated.}
}
\value{
    \code{cdssden} returns a list object with the following elements.
    \item{pdf}{Vector of conditional pdf.}
    \item{int}{Normalizing constant.}
      
    \code{cpssden} and \code{cqssden} return a vector of conditional cdf
    or quantiles.
}
\details{
    The argument \code{x} in \code{cdssden} is of the same form as the
    argument \code{newdata} in \code{\link{predict.lm}}, but can take a
    vector for 1-D conditional densities.

    \code{cpssden} and \code{cqssden} naturally only work for 1-D
    conditional densities of a numerical variable.
}
\note{
    If variables other than factors or numerical vectors are involved in
    \code{x}, the normalizing constant can not be computed.
}
\seealso{
    Fitting function \code{\link{ssden}} and \code{\link{dssden}}.
}
\keyword{models}
\keyword{distribution}
\keyword{smooth}

\name{mkterm.ssanova}
\alias{mkterm.ssanova}
\alias{mkterm.linear}
\alias{mkterm.cubic}
\alias{mkterm.tp}
\title{
    Assembling Model Terms for Smoothing Spline ANOVA Models
}
\description{
    Assemble numerical functions for calculating model terms in a
    Smoothing Spline ANOVA Model.
}
\usage{
mkterm.linear(mf, ext)
mkterm.cubic(mf, ext)
mkterm.tp(mf, order, mesh, weight)
}
\arguments{
    \item{mf}{Model frame of the model formula.}
    \item{ext}{Size of the "buffer zone" beyond the data range.}
    \item{order}{Order of the differential operator.}
    \item{mesh}{Normalizing mesh.}
    \item{weight}{Normalizing weights.}
}
\details{
    These are not to be called by the user.

    For polynomial splines, \code{ext} specifies how far to go beyond
    the data range percentage wise.  For example, if the minimum and
    maximum values of a variable in \code{mf} is 0 and 1 and
    \code{ext=.05}, then the marginal domain on which the model is
    defined would be \eqn{[-.95,1.05]}.

    See \code{\link{mkfun.tp}} for \code{order}, \code{mesh},
    \code{weight}.
}
\value{
    A list object with a component \code{labels} containing the labels
    of all model terms.  For each of the model terms, there is a
    component holding the numerical functions for calculating the fixed
    and random effects within the term.
}
\note{
    The numerical functions are assembled using building blocks crafted
    by \code{\link{mkfun.poly}}, \code{\link{mkfun.tp}},
    \code{\link{mkrk.factor}}.
}
\keyword{internal}

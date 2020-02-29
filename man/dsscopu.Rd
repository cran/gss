\name{dsscopu}
\alias{dsscopu}
\title{Evaluating Copula Density Estimates}
\description{
    Evaluate copula density estimates.
}
\usage{
dsscopu(object, x, copu=TRUE)
}
\arguments{
    \item{object}{Object of class \code{"sscopu"}.}
    \item{x}{Vector or matrix of point(s) on which copula density is to
        be evaluated.}
    \item{copu}{Flag indicating whether to apply copularization.}
}
\value{
    A vector of copula density values.
}
\seealso{
    Fitting functions \code{\link{sscopu}} and \code{\link{sscopu2}}.
}
\keyword{models}
\keyword{distribution}
\keyword{smooth}

\name{rkpk1}
\alias{sspreg1}
\alias{mspreg1}
\title{Numerical Engine for ssanova1}
\description{
    Calculate penalized least squares regression estimates via the
    normal equation and evaluate the GCV, GML, or Mallows' CL scores, as
    implemented in the RATFOR routine \code{reg.r}, and minimize the
    cross-validation score using \code{\link{nlm}}.
}
\usage{
sspreg1(s,r,q,y,method,alpha,varht,random)
mspreg1(s,r,q,y,method,alpha,varht,random)
}
\details{
    \code{sspreg1} is used by \code{\link{ssanova1}} to compute
    regression estimates with a single smoothing parameter.
    \code{mspreg2} is used by \code{\link{ssanova1}} to compute
    regression estimates with multiple smoothing parameters.
}
\arguments{
    \item{s}{Fixed effects evaluated at data points.}
    \item{r}{Basis of random effects evaluated at data points.}
    \item{q}{Penalty matrix.}
    \item{y}{Response vector.}
    \item{method}{\code{"v"} for GCV, \code{"m"} for GML, or \code{"u"}
        for Mallows' CL.}
    \item{alpha}{Parameter modifying GCV or Mallows' CL scores for
        smoothing parameter selection.}
    \item{varht}{External variance estimate needed for \code{method="u"}.}
    \item{random}{Input for parametric random effects in nonparametric
        mixed-effect models.}
}
\references{
    Kim, Y.-J. and Gu, C. (2002) \emph{Penalized Least Squares
    Regression: Fast Computation via Efficient Approximation}.
    Available at \url{http://stat.purdue.edu/~chong/manu.html}.
}
\keyword{smooth}

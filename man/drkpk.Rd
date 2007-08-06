\name{drkpk}
\alias{sspdsty}
\alias{mspdsty}
\alias{msphzd}
\title{Numerical Engine for ssden and sshzd}
\description{
    Perform numerical calculations for the \code{\link{ssden}} and
    \code{\link{sshzd}} suites.
}
\usage{
sspdsty(s, r, q, cnt, qd.s, qd.r, qd.wt, prec, maxiter, alpha)
mspdsty(s, r, q, cnt, qd.s, qd.r, qd.wt, prec, maxiter, alpha)

msphzd(s, r, q, Nobs, cnt, qd.s, qd.r, qd.wt, prec, maxiter, alpha)
}
\details{
    \code{sspdsty} is used by \code{\link{ssden}} to compute
    cross-validated density estimate with a single smoothing
    parameter. \code{mspdsty} is used by \code{\link{ssden}} to compute
    cross-validated density estimate with multiple smoothing
    parameters.

    \code{msphzd} is used by \code{\link{sshzd}} to compute
    cross-validated hazard estimate with single or multiple smoothing
    parameters.
}
\arguments{
    \item{s}{Unpenalized terms evaluated at data points.}
    \item{r}{Basis of penalized terms evaluated at data points.}
    \item{q}{Penalty matrix.}
    \item{Nobs}{Total number of lifetime observations.}
    \item{cnt}{Bin-counts for histogram data.}
    \item{qd.s}{Unpenalized terms evaluated at quadrature nodes.}
    \item{qd.r}{Basis of penalized terms evaluated at quadrature nodes.}
    \item{qd.wt}{Quadrature weights.}
    \item{prec}{Precision requirement for internal iterations.}
    \item{maxiter}{Maximum number of iterations allowed for
        internal iterations.}
    \item{alpha}{Parameter defining cross-validation score for smoothing
        parameter selection.}
}
\references{
    Gu, C. (2002), \emph{Smoothing Spline ANOVA Models}.  New York:
    Springer-Verlag.

    Gu, C. and Wang, J. (2003), Penalized likelihood density
    estimation: Direct cross-validation and scalable approximation.
    \emph{Statistica Sinica}, \bold{13}, 811--826.
}
\keyword{internal}

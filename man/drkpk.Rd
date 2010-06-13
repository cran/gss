\name{drkpk}
\alias{sspdsty}
\alias{mspdsty}
\alias{msphzd}
\alias{msphzd1}
\alias{mspllrm}
\title{Numerical Engine for ssden, sshzd, and sshzd1}
\description{
    Perform numerical calculations for the \code{\link{ssden}} and
    \code{\link{sshzd}} suites.
}
\usage{
sspdsty(s, r, q, cnt, qd.s, qd.r, qd.wt, prec, maxiter, alpha)
mspdsty(s, r, id.basis, cnt, qd.s, qd.r, qd.wt, prec, maxiter, alpha, skip.iter)

msphzd(s, r, id.wk, Nobs, cnt, qd.s, qd.r, qd.wt, prec, maxiter, alpha, skip.iter)

msphzd1(s, r, id.wk, Nobs, cnt, int.s, int.r, rho, prec, maxiter, alpha, skip.iter)

mspllrm(s, r, id.basis, cnt, qd.s, qd.r, xx.wt, random, prec, maxiter, alpha, skip.iter)
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
    \item{id.basis}{Index of observations to be used as "knots."}
    \item{id.wk}{Index of observations to be used as "knots."}
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
    \item{skip.iter}{Flag indicating whether to use initial values of
        theta and skip theta iteration.}
    \item{int.s}{Integrals of unpenalized terms.}
    \item{int.r}{Integrals of basis of penalized terms.}
    \item{rho}{rho function value on failure times.}
    \item{xx.wt}{Weights at unique x.}
    \item{random}{Input for parametric random effects in nonparametric
        mixed-effect models.}
}
\references{
    Gu, C. (2002), \emph{Smoothing Spline ANOVA Models}.  New York:
    Springer-Verlag.

    Gu, C. and Wang, J. (2003), Penalized likelihood density
    estimation: Direct cross-validation and scalable approximation.
    \emph{Statistica Sinica}, \bold{13}, 811--826.

    Du, P. and Gu, C. (2006), Penalized likelihood hazard estimation:
    efficient approximation and Bayesian confidence intervals.
    \emph{Statistics and Probability Letters}, \bold{76}, 244--254.

    Du, P. and Gu, C. (2008), Penalized Pseudo-Likelihood Hazard
    Estimation: A Fast Alternative to Penalized Likelihood.
    \emph{Journal of Statistical Planning and Inference}, \bold{00},
    000--000.
}
\keyword{internal}

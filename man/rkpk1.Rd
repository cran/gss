\name{rkpk1}
\alias{sspreg1}
\alias{mspreg1}
\alias{sspngreg1}
\alias{mspngreg1}
\alias{ngreg1}
\alias{ngreg.proj}
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

sspngreg1(family,s,r,q,y,wt,offset,alpha,nu,random)
mspngreg1(family,s,r,q,y,wt,offset,alpha,nu,random)
ngreg1(dc,family,sr,q,y,wt,offset,nu,alpha)

ngreg.proj(dc,family,sr,q,fit0,wt,offset,nu)
}
\details{
    \code{sspreg1} is used by \code{\link{ssanova1}} to compute
    regression estimates with a single smoothing parameter.
    \code{mspreg1} is used by \code{\link{ssanova1}} to compute
    regression estimates with multiple smoothing parameters.

    \code{ssngpreg1} is used by \code{\link{gssanova1}} to compute
    non-Gaussian regression estimates with a single smoothing
    parameter.  \code{mspngreg1} is used by \code{\link{gssanova1}} to
    compute non-Gaussian regression estimates with multiple smoothing
    parameters.  \code{ngreg1} is used by \code{ssngpreg1} and
    \code{mspngreg1} to perform Newton iteration with fixed smoothing
    parameters and to calculate cross-validation scores on return.

    \code{ngreg.proj} is used by \code{project.gssanova1} to calculate
    Kullback-Leibler projection for non-Gaussian regression.
}
\arguments{
    \item{family}{Description of the error distribution.  Supported
	are exponential families \code{"binomial"}, \code{"poisson"},
	and \code{"Gamma"}}.
    \item{s}{Fixed effects evaluated at data points.}
    \item{r}{Basis of random effects evaluated at data points.}
    \item{q}{Penalty matrix.}
    \item{y}{Response vector.}
    \item{wt}{Model weights.}
    \item{offset}{Model offset.}
    \item{method}{\code{"v"} for GCV, \code{"m"} for GML, or \code{"u"}
        for Mallows' CL.}
    \item{alpha}{Parameter modifying GCV or Mallows' CL scores for
        smoothing parameter selection.}
    \item{nu}{Optional argument for future support of nbinomial,
      weibull, lognorm, and loglogis families.}
    \item{varht}{External variance estimate needed for \code{method="u"}.}
    \item{random}{Input for parametric random effects in nonparametric
        mixed-effect models.}

    \item{dc}{Coefficients of fits.}
    \item{sr}{\code{cbind(s,r)}.}
    \item{fit0}{The fit to be projected.}
}
\references{
    Kim, Y.-J. and Gu, C. (2002) \emph{Penalized Least Squares
    Regression: Fast Computation via Efficient Approximation}.
    Available at \url{http://stat.purdue.edu/~chong/manu.html}.
}
\keyword{smooth}

\name{sscopu}
\alias{sscopu}
\alias{sscopu2}
\title{Estimating Copula Density Using Smoothing Splines}
\description{
    Estimate copula densities using tensor-product cubic splines.
}
\usage{
sscopu(x, symmetry=FALSE, alpha=1.4, order=NULL, exclude=NULL,
       weights=NULL, id.basis=NULL, nbasis=NULL, seed=NULL,
       qdsz.depth=NULL, prec=1e-7, maxiter=30, skip.iter=dim(x)[2]!=2)

sscopu2(x, censoring=NULL, truncation=NULL, symmetry=FALSE, alpha=1.4,
        weights=NULL, id.basis=NULL, nbasis=NULL, seed=NULL, prec=1e-7,
        maxiter=30)
}
\arguments{
    \item{x}{Matrix of observations on unit cubes.}
    \item{symmetry}{Flag indicating whether to enforce symmetry, or
        invariance under coordinate permutation.}
    \item{order}{Highest order of interaction terms in log density.
        When \code{NULL}, it is set to \code{dim(x)[2]} internally.}
    \item{exclude}{Pair(s) of marginals whose interactions to be
        excluded in log density.}
    \item{alpha}{Parameter defining cross-validation score for smoothing
        parameter selection.}
    \item{weights}{Optional vector of bin-counts for histogram data.}
    \item{id.basis}{Index of observations to be used as "knots."}
    \item{nbasis}{Number of "knots" to be used.  Ignored when
        \code{id.basis} is specified.}
    \item{seed}{Seed to be used for the random generation of "knots."
        Ignored when \code{id.basis} is specified.}
    \item{qdsz.depth}{Depth to be used in \code{\link{smolyak.quad}} for
        the generation of quadrature.}
    \item{prec}{Precision requirement for internal iterations.}
    \item{maxiter}{Maximum number of iterations allowed for
        internal iterations.}
    \item{skip.iter}{Flag indicating whether to use initial values of
        theta and skip theta iteration.  See \code{\link{ssanova}} for
	notes on skipping theta iteration.}
    \item{censoring}{Optional censoring indicator.}
    \item{truncation}{Optional truncation points.}
}
\details{
    \code{sscopu} is essentially \code{\link{ssden}} applied to
    observations on unit cubes.  Instead of variables in data frames,
    the data are entered as a numerical matrix, and model complexity is
    globally controlled by the highest order of interactions allowed in    log density.

    \code{sscopu2} further restricts the domain to the unit square, but
    allows for possible censoring and truncation.  With
    \code{censoring==0,1,2,3}, a data point \eqn{(x1,x2)} represents
    exact observation, \eqn{[0,x1]x{x2}}, \eqn{{x1}x[0,x2]}, or
    \eqn{[0,x1]x[0,x2]}.  With \code{truncation} point \eqn{(t1,t2)},
    the sample is taken from \eqn{[0,t1]x[0,t2]} instead of the unit
    square.

    With \code{symmetriy=TRUE}, one may enforce the interchangeability
    of coordinates so that \eqn{f(x1,x2)=f(x2,x1)}, say.

    When \code{(1,2)} is a row in \code{exclude}, interaction terms
    involving coordinates \code{1} and \code{2} are excluded.
}
\note{
    For reasonable execution time in higher dimensions, set
    \code{skip.iter=TRUE} in calls to \code{sscopu}.

    When \code{"Newton iteration diverges"} in \code{sscopu}, try to use
    a larger \code{qdsz.depth}; the default values for dimensions 2, 3,
    4, 5, 6+ are 24, 14, 12, 11, 10.  To be sure a larger
    \code{qdsz.depth} indeed makes difference, verify the cubature size
    using \code{\link{smolyak.size}}.

    The results may vary from run to run.  For consistency, specify
    \code{id.basis} or set \code{seed}.
}
\value{
    \code{sscopu} and \code{sscopu2} return a list object of class
    \code{"sscopu"}.  \code{\link{dsscopu}} can be used to evaluate the
    estimated copula density.  A "copularization" process is applied to
    the estimated density by default so the resulting marginal densities
    are guaranteed to be uniform.

    \code{\link{cdsscopu}}, \code{\link{cpsscopu}}, and
    \code{\link{cqsscopu}} can be used to evaluate 1-D conditional pdf,
    cdf, and quantiles.
}
\author{Chong Gu, \email{chong@stat.purdue.edu}}
\references{
    Gu, C. (2013), \emph{Smoothing Spline ANOVA Models (2nd Ed)}.  New
    York: Springer-Verlag.
    
    Gu, C. (2015), Hazard estimation with bivariate survival data and
    copula density estimation. \emph{Journal of Computational and
    Graphical Statistics}, \bold{24}, 1053-1073.
}
\examples{
## simulate 2-D data
x <- matrix(runif(200),100,2)
## fit copula density
fit <- sscopu(x)
## "same fit"
fit2 <- sscopu2(x,id=fit$id)
## symmetric fit
fit.s <- sscopu(x,sym=TRUE,id=fit$id)
## Kendall's tau and Spearman's rho
summary(fit); summary(fit2); summary(fit.s)
## clean up
\dontrun{rm(x,fit,fit2,fit.s)}
}
\keyword{smooth}
\keyword{models}
\keyword{distribution}

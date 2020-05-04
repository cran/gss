\name{sscden}
\alias{sscden}
\alias{sscden1}
\title{Estimating Conditional Probability Density Using Smoothing
    Splines}
\description{
    Estimate conditional probability densities using smoothing spline
    ANOVA models.  The symbolic model specification via \code{formula}
    follows the same rules as in \code{\link{lm}}.
}
\usage{
sscden(formula, response, type=NULL, data=list(), weights, subset,
       na.action=na.omit, alpha=1.4, id.basis=NULL, nbasis=NULL,
       seed=NULL, ydomain=as.list(NULL), yquad=NULL, prec=1e-7,
       maxiter=30, skip.iter=FALSE)

sscden1(formula, response, type=NULL, data=list(), weights, subset,
        na.action=na.omit, alpha=1.4, id.basis=NULL, nbasis=NULL,
        seed=NULL, rho=list("xy"), ydomain=as.list(NULL), yquad=NULL,
        prec=1e-7, maxiter=30, skip.iter=FALSE)
}
\arguments{
    \item{formula}{Symbolic description of the model to be fit.}
    \item{response}{Formula listing response variables.}
    \item{type}{List specifying the type of spline for each variable.
        See \code{\link{mkterm}} for details.}
    \item{data}{Optional data frame containing the variables in the
        model.}
    \item{weights}{Optional vector of counts for duplicated data.}
    \item{subset}{Optional vector specifying a subset of observations
	to be used in the fitting process.}
    \item{na.action}{Function which indicates what should happen when
        the data contain NAs.}
    \item{alpha}{Parameter defining cross-validation scores for
        smoothing parameter selection.}
    \item{id.basis}{Index of observations to be used as "knots."}
    \item{nbasis}{Number of "knots" to be used.  Ignored when
        \code{id.basis} is specified.}
    \item{seed}{Seed to be used for the random generation of "knots."
        Ignored when \code{id.basis} is specified.}
    \item{ydomain}{Data frame specifying marginal support of conditional
        density.}
    \item{yquad}{Quadrature for calculating integral on Y domain.
        Mandatory if response variables other than factors or numerical
	vectors are involved.}
    \item{prec}{Precision requirement for internal iterations.}
    \item{maxiter}{Maximum number of iterations allowed for
        internal iterations.}
    \item{skip.iter}{Flag indicating whether to use initial values of
        theta and skip theta iteration.  See \code{\link{ssanova}} for
	notes on skipping theta iteration.}
    \item{rho}{rho function needed for sscden1.}
}
\details{
    The model is specified via \code{formula} and \code{response}, where
    \code{response} lists the response variables.  For example,
    \code{sscden(~y*x,~y)} prescribe a model of the form
    \deqn{
        log f(y|x) = g_{y}(y) + g_{xy}(x,y) + C(x)
    }
    with the terms denoted by \code{"y"}, \code{"y:x"}; the term(s) not
    involving response(s) are removed and the constant \code{C(x)} is
    determined by the fact that a conditional density integrates to one
    on the \code{y} axis.  \code{sscden1} does keep terms not involving
    response(s) during estimation, although those terms cancel out when
    one evaluates the estimated conditional density.    

    The model terms are sums of unpenalized and penalized
    terms. Attached to every penalized term there is a smoothing
    parameter, and the model complexity is largely determined by the
    number of smoothing parameters.

    A subset of the observations are selected as "knots."  Unless
    specified via \code{id.basis} or \code{nbasis}, the number of
    "knots" \eqn{q} is determined by \eqn{max(30,10n^{2/9})}, which is
    appropriate for the default cubic splines for numerical vectors.
}
\note{
    Default quadrature on the Y domain will be constructed for numerical
    vectors on a hyper cube, then outer product with factor levels will
    be taken if factors are involved.  The sides of the hyper cube are
    specified by \code{ydomain}; for \code{ydomain$y} missing, the default
    is \code{c(min(y),max(y))+c(-1,1)*(max(y)-mimn(y))*.05}.

    On a 1-D interval, the quadrature is the 200-point Gauss-Legendre
    formula returned from \code{\link{gauss.quad}}.  For multiple
    numerical vectors, delayed Smolyak cubatures from
    \code{\link{smolyak.quad}} are used on cubes with the marginals
    properly transformed; see Gu and Wang (2003) for the marginal
    transformations.

    The results may vary from run to run.  For consistency, specify
    \code{id.basis} or set \code{seed}.

    For reasonable execution time in high dimensions, set
    \code{skip.iter=TRUE}.
}
\value{
    \code{sscden} returns a list object of class \code{"sscden"}.
    \code{sscden1} returns a list object of class
    \code{c("sscden1","sscden")}.

    \code{\link{dsscden}} and \code{\link{cdsscden}} can be used to
    evaluate the estimated conditional density \eqn{f(y|x)} and
    \eqn{f(y1|x,y2)}; \code{\link{psscden}}, \code{\link{qsscden}},
    \code{\link{cpsscden}}, and \code{\link{cqsscden}} can be used to
    evaluate conditional cdf and quantiles.  The methods
    \code{\link{project.sscden}} or \code{\link{project.sscden1}} can
    be used to calculate the Kullback-Leibler or square-error
    projections for model selection.
}
\references{
    Gu, C. (1995), Smoothing spline density estimation: Conditional
    distribution.  \emph{Statistica Sinica}, \bold{5}, 709--726.
    Springer-Verlag.

    Gu, C. (2014), Smoothing Spline ANOVA Models: R Package gss.
    \emph{Journal of Statistical Software}, 58(5), 1-25. URL
    http://www.jstatsoft.org/v58/i05/.
}
\examples{
data(penny); set.seed(5732)
fit <- sscden1(~year*mil,~mil,data=penny,
              ydomain=data.frame(mil=c(49,61)))
yy <- 1944+(0:92)/2
quan <- qsscden(fit,c(.05,.25,.5,.75,.95),
                data.frame(year=yy))
plot(penny$year+.1*rnorm(90),penny$mil,ylim=c(49,61))
for (i in 1:5) lines(yy,quan[i,])
## Clean up
\dontrun{rm(penny,yy,quan)}
}
\keyword{smooth}
\keyword{models}
\keyword{distribution}

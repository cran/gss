\name{ssanova0}
\alias{ssanova0}
\title{Fitting Smoothing Spline ANOVA Models}
\description{
    Fit smoothing spline ANOVA models in Gaussian regression.  The
    symbolic model specification via \code{formula} follows the same
    rules as in \code{\link{lm}}.
}
\usage{
ssanova0(formula, type=NULL, data=list(), weights, subset,
         offset, na.action=na.omit, partial=NULL, method="v",
         varht=1, prec=1e-7, maxiter=30)
}
\arguments{
    \item{formula}{Symbolic description of the model to be fit.}
    \item{type}{List specifying the type of spline for each variable.
        See \code{\link{mkterm}} for details.}
    \item{data}{Optional data frame containing the variables in the
	model.}
    \item{weights}{Optional vector of weights to be used in the
	fitting process.}
    \item{subset}{Optional vector specifying a subset of observations
	to be used in the fitting process.}
    \item{offset}{Optional offset term with known parameter 1.}
    \item{na.action}{Function which indicates what should happen when
	the data contain NAs.}
    \item{partial}{Optional symbolic description of parametric terms in
        partial spline models.}
    \item{method}{Method for smoothing parameter selection.  Supported
	are \code{method="v"} for GCV, \code{method="m"} for GML (REML),
	and \code{method="u"} for Mallow's CL.}
    \item{varht}{External variance estimate needed for
	\code{method="u"}.  Ignored when \code{method="v"} or
	\code{method="m"} are specified.}
    \item{prec}{Precision requirement in the iteration for multiple
	smoothing parameter selection.  Ignored when only one smoothing
	parameter is involved.}
    \item{maxiter}{Maximum number of iterations allowed for multiple
	smoothing parameter selection.  Ignored when only one smoothing
	parameter is involved.}
}
\details{
    The model specification via \code{formula} is intuitive.  For
    example, \code{y~x1*x2} yields a model of the form
    \deqn{
	y = C + f_{1}(x1) + f_{2}(x2) + f_{12}(x1,x2) + e
    }
    with the terms denoted by \code{"1"}, \code{"x1"}, \code{"x2"}, and
    \code{"x1:x2"}.

    The model terms are sums of unpenalized and penalized
    terms. Attached to every penalized term there is a smoothing
    parameter, and the model complexity is largely determined by the
    number of smoothing parameters.

    \code{ssanova0} and the affiliated methods provide a front end to
    RKPACK, a collection of RATFOR routines for nonparametric regression
    via the penalized least squares.  The algorithms implemented in
    RKPACK are of the order \eqn{O(n^{3})}.
}
\note{
    For complex models and large sample sizes, the approximate solution
    of \code{\link{ssanova}} can be faster.

    The method \code{\link{project}} is not implemented for
    \code{ssanova0}, nor is the mixed-effect model support through
    \code{\link{mkran}}.

    In \emph{gss} versions earlier than 1.0, \code{ssanova0} was under
    the name \code{ssanova}.
}
\value{
    \code{ssanova0} returns a list object of class
    \code{c("ssanova0","ssanova")}.

    The method \code{\link{summary.ssanova0}} can be used to obtain
    summaries of the fits.  The method \code{\link{predict.ssanova0}}
    can be used to evaluate the fits at arbitrary points along with
    standard errors.  The methods \code{\link{residuals.ssanova}} and
    \code{\link{fitted.ssanova}} extract the respective traits from the
    fits.
}
\author{Chong Gu, \email{chong@stat.purdue.edu}}
\references{
    Wahba, G. (1990), \emph{Spline Models for Observational Data}.
    Philadelphia: SIAM.

    Gu, C. (2013), \emph{Smoothing Spline ANOVA Models (2nd Ed)}.  New
    York: Springer-Verlag.

    Gu, C. (2014), Smoothing Spline ANOVA Models: R Package gss.
    \emph{Journal of Statistical Software}, 58(5), 1-25. URL
    http://www.jstatsoft.org/v58/i05/.
}
\examples{
## Fit a cubic spline
x <- runif(100); y <- 5 + 3*sin(2*pi*x) + rnorm(x)
cubic.fit <- ssanova0(y~x,method="m")
## Obtain estimates and standard errors on a grid
new <- data.frame(x=seq(min(x),max(x),len=50))
est <- predict(cubic.fit,new,se=TRUE)
## Plot the fit and the Bayesian confidence intervals
plot(x,y,col=1); lines(new$x,est$fit,col=2)
lines(new$x,est$fit+1.96*est$se,col=3)
lines(new$x,est$fit-1.96*est$se,col=3)
## Clean up
\dontrun{rm(x,y,cubic.fit,new,est)
dev.off()}

## Fit a tensor product cubic spline
data(nox)
nox.fit <- ssanova0(log10(nox)~comp*equi,data=nox)
## Fit a spline with cubic and nominal marginals
nox$comp<-as.factor(nox$comp)
nox.fit.n <- ssanova0(log10(nox)~comp*equi,data=nox)
## Fit a spline with cubic and ordinal marginals
nox$comp<-as.ordered(nox$comp)
nox.fit.o <- ssanova0(log10(nox)~comp*equi,data=nox)
## Clean up
\dontrun{rm(nox,nox.fit,nox.fit.n,nox.fit.o)}
}
\keyword{smooth}
\keyword{models}
\keyword{regression}

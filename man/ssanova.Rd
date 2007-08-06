\name{ssanova}
\alias{ssanova}
\title{Fitting Smoothing Spline ANOVA Models}
\description{
    Fit smoothing spline ANOVA models in Gaussian regression.  The
    symbolic model specification via \code{formula} follows the same
    rules as in \code{\link{lm}}.
}
\usage{
ssanova(formula, type=NULL, data=list(), weights, subset, offset,
        na.action=na.omit, partial=NULL, method="v", alpha=1.4,
        varht=1, id.basis=NULL, nbasis=NULL, seed=NULL, random=NULL)
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
    \item{partial}{Optional extra unpenalized terms in partial spline
        models.}
    \item{method}{Method for smoothing parameter selection.  Supported
	are \code{method="v"} for GCV, \code{method="m"} for GML (REML),
	and \code{method="u"} for Mallows' CL.}
    \item{alpha}{Parameter modifying GCV or Mallows' CL; larger absolute
        values yield smoother fits; negative value invokes a stable and
	more accurate GCV/CL evaluation algorithm but may take two to
	five times as long.  Ignored when \code{method="m"} are
	specified.}
    \item{varht}{External variance estimate needed for
	\code{method="u"}.  Ignored when \code{method="v"} or
	\code{method="m"} are specified.}
    \item{id.basis}{Index designating selected "knots".}
    \item{nbasis}{Number of "knots" to be selected.  Ignored when
	\code{id.basis} is supplied.}
    \item{seed}{Seed to be used for the random generation of "knots".
	Ignored when \code{id.basis} is supplied.}
    \item{random}{Input for parametric random effects in nonparametric
        mixed-effect models.  See \code{\link{mkran}} for details.}
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

    A subset of the observations are selected as "knots."  Unless
    specified via \code{id.basis} or \code{nbasis}, the number of
    "knots" \eqn{q} is determined by \eqn{max(30,10n^{2/9})}, which is
    appropriate for the default cubic splines for numerical vectors.

    Using \eqn{q} "knots," \code{ssanova} calculates an approximate
    solution to the penalized least squares problem using algorithms of
    the order \eqn{O(nq^{2})}, which for \eqn{q<<n} scale better than
    the \eqn{O(n^{3})} algorithms of \code{\link{ssanova0}}.  For the
    exact solution, one may set \eqn{q=n} in \code{ssanova}, but
    \code{\link{ssanova0}} would be much faster.
}
\note{
    To use GCV and Mallows' CL unmodified, set \code{alpha=1}.

    For simpler models and moderate sample sizes, the exact solution of
    \code{\link{ssanova0}} can be faster.

    The results may vary from run to run. For consistency, specify
    \code{id.basis} or set \code{seed}.

    In \emph{gss} versions earlier than 1.0, \code{ssanova} was under
    the name \code{ssanova1}.
}
\value{
    \code{ssanova} returns a list object of class \code{"ssanova"}.

    The method \code{\link{summary.ssanova}} can be used to obtain
    summaries of the fits.  The method \code{\link{predict.ssanova}} can
    be used to evaluate the fits at arbitrary points along with standard
    errors.  The method \code{\link{project.ssanova}} can be used to
    calculate the Kullback-Leibler projection for model selection.  The
    methods \code{\link{residuals.ssanova}} and
    \code{\link{fitted.ssanova}} extract the respective traits
    from the fits.
}
\author{Chong Gu, \email{chong@stat.purdue.edu}}
\references{
    Gu, C. (2002), \emph{Smoothing Spline ANOVA Models}.  New York:
    Springer-Verlag.

    Kim, Y.-J. and Gu, C. (2004), Smoothing spline Gaussian regression:
    more scalable computation via efficient approximation.
    \emph{Journal of the Royal Statistical Society, Ser. B}, \bold{66},
    337--356.
    
    Wahba, G. (1990), \emph{Spline Models for Observational Data}.
    Philadelphia: SIAM.
}
\examples{
## Fit a cubic spline
x <- runif(100); y <- 5 + 3*sin(2*pi*x) + rnorm(x)
cubic.fit <- ssanova(y~x)
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
nox.fit <- ssanova(log10(nox)~comp*equi,data=nox)
## Fit a spline with cubic and nominal marginals
nox$comp<-as.factor(nox$comp)
nox.fit.n <- ssanova(log10(nox)~comp*equi,data=nox)
## Fit a spline with cubic and ordinal marginals
nox$comp<-as.ordered(nox$comp)
nox.fit.o <- ssanova(log10(nox)~comp*equi,data=nox)
## Clean up
\dontrun{rm(nox,nox.fit,nox.fit.n,nox.fit.o)}
}
\keyword{smooth}
\keyword{models}
\keyword{regression}

\name{ssanova1}
\alias{ssanova1}
\title{Fitting Smoothing Spline ANOVA Models}
\description{
    Fit smoothing spline ANOVA models with cubic spline, linear spline,
    or thin-plate spline marginals for numerical variables.  Factors are
    also accepted.  The symbolic model specification via \code{formula}
    follows the same rules as in \code{\link{lm}}.
}
\usage{
ssanova1(formula, type="cubic", data=list(), weights, subset,
         offset, na.action=na.omit, partial=NULL, method="v",
         alpha=1.4, varht=1, id.basis=NULL, nbasis=NULL, seed=NULL,
         random=NULL, ext=.05, order=2)
}
\arguments{
    \item{formula}{Symbolic description of the model to be fit.}
    \item{type}{Type of numerical marginals to be used.  Supported are
	\code{type="cubic"} for cubic spline marginals,
	\code{type="linear"} for linear spline marginals, and
	\code{type="tp"} for thin-plate spline marginals.}
    \item{data}{Optional data frame containing the variables in the
	model.}
    \item{weights}{Optional vector of weights to be used in the
	fitting process.}
    \item{subset}{Optional vector specifying a subset of observations
	to be used in the fitting process.}
    \item{offset}{Optional offset term with known parameter 1.}
    \item{na.action}{Function which indicates what should happen when
	the data contain NAs.}
    \item{partial}{Optional extra fixed effect terms in partial spline
	models.}
    \item{method}{Method for smoothing parameter selection.  Supported
	are \code{method="v"} for GCV, \code{method="m"} for GML (REML),
	and \code{method="u"} for Mallows' CL.}
    \item{alpha}{Parameter modifying GCV or Mallows' CL; larger values
        yield smoother fits.  Ignored when \code{method="m"} are
	specified.}
    \item{varht}{External variance estimate needed for
	\code{method="u"}.  Ignored when \code{method="v"} or
	\code{method="m"} are specified.}
    \item{id.basis}{Index designating selected "knots".}
    \item{nbasis}{Number of "knots" to be selected.  Ignored when
	\code{id.basis} is supplied.}
    \item{seed}{Seed for reproducible random selection of "knots".
	Ignored when \code{id.basis} is supplied.}
    \item{random}{Input for parametric random effects in nonparametric
        mixed-effect models.  The utility is not yet ready for public
	use.}
    \item{ext}{For cubic spline and linear spline marginals, this option
	specifies how far to extend the domain beyond the minimum and
	the maximum as a percentage of the range.  The default
	\code{ext=.05} specifies marginal domains of lengths 110 percent
	of their respective ranges.  Prediction outside of the domain
	will result in an error.  Ignored if \code{type="tp"} is
	specified.}
    \item{order}{For thin-plate spline marginals, this option specifies
	the order of the marginal penalties.  Ignored if
	\code{type="cubic"} or \code{type="linear"} are specified.}
}
\details{
    \code{ssanova1} implements an alternative approach to the fitting of
    \code{ssanova} models.  It works in some \eqn{q}-dimensional
    function space determined by a set of "knots" typically through
    random selection; the default dimension (\code{nbasis}) is set to be
    \eqn{q=10n^(2/9)}, adequate for cubic splines.  The algorithms are
    of the order \eqn{O(nq^2)}, scaling much better than the
    \eqn{O(n^3)} of \code{ssanova}; the memory requirements are
    \eqn{O(nq)} versus \eqn{O(n^2)}.

    A simple modification of GCV and Mallows' CL is incorporated in the
    code, to help preventing the occasional severe undersmoothing by
    these methods.  To use the methods unmodified, set \code{alpha=1}.

    To "duplicate" \code{ssanova} fits, use \code{ssanova1} with
    \code{nbasis} set to the sample size and \code{alpha} set to 1;
    numerical discrepancies are to be found, but should be largely
    negligible for practical purposes.  The \code{ssanova1} "duplicates"
    of the \code{ssanova} fits would take much longer to calculate,
    however.

    The model specification via \code{formula} is intuitive.  For
    example, \code{y~x1*x2} yields a model of the form
    \deqn{
	y = c + f_{1}(x1) + f_{2}(x2) + f_{12}(x1,x2) + e
    }
    with the terms denoted by \code{"1"}, \code{"x1"}, \code{"x2"}, and
    \code{"x1:x2"}.  Through the specifications of the side conditions,
    these terms are uniquely defined.  In the current implementation,
    \eqn{f_{1}} and \eqn{f_{12}} integrate to \eqn{0} on the \code{x1}
    domain for cubic spline and linear spline marginals, and add to
    \eqn{0} over the \code{x1} (marginal) sampling points for thin-plate
    spline marginals.

    The penalized least squares problem is equivalent to a certain
    empirical Bayes model or a mixed effect model, and the model terms
    themselves are generally sums of finer terms of two types, the
    unpenalized \emph{fixed effects} and the penalized \emph{random}
    \emph{effects}.  Attached to every random effect there is a
    smoothing parameter, and the model complexity is largely determined
    by the number of smoothing parameters.

    The method \code{predict} can be used to evaluate the sum of
    selected or all model terms at arbitrary points within the domain,
    along with standard errors derived from a certain Bayesian
    calculation.  The method \code{summary} has a flag to request
    diagnostics for the practical identifiability and significance of
    the model terms.
}
\note{
    The independent variables appearing in \code{formula} can be
    multivariate themselves.  In particular,
    \code{ssanova(y~x,"tp",order=order)} can be used to fit ordinary
    thin-plate splines in any dimension, of any order permissible, and
    with standard errors available for Bayesian confidence intervals.
    Note that thin-plate splines reduce to polynomial splines in one
    dimension.

    For univariate marginals, the additive models using
    \code{type="cubic"} and \code{type="tp"} yield identical fit through
    different internal makes.  For example,
    \code{ssanova(y~x1+x2,"cubic")} and \code{ssanova(y~x1+x2,"tp")}
    yield the same fit.  The same is not true for models with
    interactions, however.
    
    Mathematically, the domain (through \code{ext} for
    \code{type="cubic"}) or the order (through \code{order} for
    \code{type="tp"}) could be specified individually for each of the
    variables.  Such flexibility is not provided in our implementation,
    however, as it would be more a source for confusion than a practical
    utility.
}
\section{Factors}{
    Factors are accepted as predictors.  When a factor has 3 or more
    levels, all terms involving it are treated as random effects, with
    the "level means" being shrunk towards each other.  The shrinking is
    done differently for nominal and ordinal factors; see
    \code{\link{mkrk.factor}} for details.
}
\value{
    \code{ssanova1} returns a list object of primary class
    \code{"ssanova1"} and secondary class \code{"ssanova"}.

    The method \code{\link{summary}} is used to obtain summaries of the
    fits.  The method \code{\link{predict}} can be used to evaluate the
    fits at arbitrary points, along with the standard errors to be used
    in Bayesian confidence intervals.  The methods
    \code{\link{residuals}} and \code{\link{fitted.values}} extract the
    respective traits from the fits.
}
\seealso{
    \code{\link{ssanova}} and methods \code{\link{predict.ssanova1}},
    \code{\link{summary.ssanova1}}, and \code{\link{fitted.ssanova}}.
}
\author{Chong Gu, \email{chong@stat.purdue.edu}}
\references{
    Kim, Y.-J. and Gu, C. (2002) \emph{Penalized Least Squares
    Regression: Fast Computation via Efficient Approximation}.
    Available at \url{http://stat.purdue.edu/~chong/manu.html}.
  
    Gu, C. (2002), \emph{Smoothing Spline ANOVA Models}.  New York:
    Springer-Verlag.
    
    Wahba, G. (1990), \emph{Spline Models for Observational Data}.
    Philadelphia: SIAM.
}
\examples{
## Fit a cubic spline
x <- runif(100); y <- 5 + 3*sin(2*pi*x) + rnorm(x)
cubic.fit <- ssanova1(y~x)
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
nox.fit <- ssanova1(log10(nox)~comp*equi,data=nox)
## Fit a spline with cubic and nominal marginals
nox$comp<-as.factor(nox$comp)
nox.fit.n <- ssanova1(log10(nox)~comp*equi,data=nox)
## Fit a spline with cubic and ordinal marginals
nox$comp<-as.ordered(nox$comp)
nox.fit.o <- ssanova1(log10(nox)~comp*equi,data=nox)
## Clean up
\dontrun{rm(nox,nox.fit,nox.fit.n,nox.fit.o)}
}
\keyword{smooth}
\keyword{models}
\keyword{regression}

\name{gssanova}
\title{Fitting Smoothing Spline ANOVA Models with Non Gaussian Responses}
\usage{
gssanova(formula, family, type="cubic", data=list(), weights, subset,
        offset, na.action=na.omit, partial=NULL, method=NULL,
        varht=1, prec=1e-7, maxiter=30, ext=.05, order=2)
}
\arguments{
    \item{formula}{a symbolic description of the model to be fit.}
    \item{family}{a description of the error distribution.  Supported
	are \code{"binomial"}, \code{"poisson"}, \code{"Gamma"},
	\code{"inverse.gaussian"}, and \code{"nbinomial"}.}
    \item{type}{the type of marginals to be used.  Supported currently
	are \code{type="cubic"} for cubic spline marginals,
	\code{type="linear"} for linear spline marginals, and
	\code{type="tp"} for thin-plate spline marginals.}
    \item{data}{an optional data frame containing the variables in the
	model.}
    \item{weights}{an optional vector of weights to be used in the
	fitting process.}
    \item{subset}{an optional vector specifying a subset of observations
	to be used in the fitting process.}
    \item{offset}{an optional offset term with known parameter 1.}
    \item{na.action}{a function which indicates what should happen when
	the data contain NAs.}
    \item{partial}{optional extra fixed effect terms in partial spline
	models.}
    \item{method}{the score used to drive the performance-oriented
	iteration.  Supported are \code{method="v"} for GCV, \code{method="m"} for
	type-II ML, and \code{method="u"} for Mallow's CL.}
    \item{varht}{an external variance estimate needed for
	\code{method="u"}.  It is ignored when \code{method="v"} or
	\code{method="m"} are specified.}
    \item{prec}{the precision in the fit required to stop the iteration
	for multiple smoothing parameter selection.  It is ignored when
	only one smoothing parameter is involved.}
    \item{maxiter}{the maximum number of iterations allowed for
	performance-oriented iteration, and for inner-loop multiple
	smoothing parameter selection when applicable.}
    \item{ext}{for cubic spline and linear spline marginals, this option
	specifies how far to extend the domain beyond the minimum and
	the maximum as a percentage of the range.  The default
	\code{ext=.05} specifies marginal domains of lengths 110 percent
	of their respective ranges.  Prediction outside of the domain
	will result in an error.  It is ignored if \code{type="tp"} is
	specified.}
    \item{order}{for thin-plate spline marginals, this option specifies
	the order of the marginal penalties.  It is ignored
	if \code{type="cubic"} or \code{type="linear"} are specified.}
}
\description{
    \code{gssanova} fits smoothing spline ANOVA models with cubic spline,
    linear spline, or thin-plate spline marginals to responses from
    selected exponential families.  The symbolic model specification via
    \code{formula} follows the same rule as in \code{\link{lm}} and
    \code{\link{glm}}.
}
\details{
    The models are fitted by penalized likelihood method through the
    performance-oriented iteration, as described in the reference cited
    below.

    Only one link is implemented for each family.  It is the logit link
    for \code{family="binomial"}, and the log link for \code{"poisson"},
    \code{"Gamma"}, and \code{"inverse.gaussian"}.  For
    \code{"nbinomial"}, the working parameter is the logit of the
    probability \eqn{p}, which is proportional to the reciprocal of the
    mean.

    For \code{family="binomial"}, \code{"poisson"}, and
    \code{"nbinomial"}, the score driving the performance-oriented
    iteration defaults to \code{method="u"} with \code{varht=1}.  For
    \code{family="Gamma"} and \code{"inverse.gaussian"}, the default is
    \code{method="v"}.

    See \code{\link{ssanova}} for details and notes concerning smoothing
    spline ANOVA models.
}
\note{
    For \code{family="binomial"}, the responses can be specified either
    as two columns of counts or as a column of sample proportion plus a
    column of weights, the same as in \code{\link{glm}}.

    For \code{family="nbinomial"}, the responses may be specified
    as two columns with the second being the known sizes, or simply a
    single column with the common unknown size to be estimated by ML.
}
\value{
    \code{gssanova} returns a list object of \code{\link{class}}
    \code{"gssanova"} which inherits from the class \code{"ssanova"}.

    The method \code{\link{summary}} is used to obtain summaries of the
    fits.  The method \code{\link{predict}} can be used to evaluate the
    fits at arbitrary points, along with the standard errors to be used
    in Bayesian confidence intervals, both on the scale of the link.
    The methods \code{\link{residuals}} and \code{\link{fitted.values}}
    extract the respective traits from the fits.
}
\seealso{
    \code{\link{predict.ssanova}} for predictions and
    \code{\link{summary.gssanova}} for summaries.
}
\author{Chong Gu, \email{chong@stat.purdue.edu}}
\references{
    Gu, C. (1992), "Cross-validating non Gaussian data," \emph{Journal
	of Computational and Graphical Statistics}, \bold{1}, 169-179.
}
\examples{
## Fit a cubic smoothing spline logistic regression model
test <- function(x)
        {.3*(1e6*(x^11*(1-x)^6)+1e4*(x^3*(1-x)^10))-2}
x <- (0:100)/100
p <- 1-1/(1+exp(test(x)))
y <- rbinom(x,3,p)
logit.fit <- gssanova(cbind(y,3-y)~x,family="binomial")

## The same fit
logit.fit1 <- gssanova(y/3~x,"binomial",weights=rep(3,101))

## Obtain estimates and standard errors on a grid
est <- predict(logit.fit,data.frame(x=x),se=TRUE)

## Plot the fit and the Bayesian confidence intervals
plot(x,y/3,ylab="p")
lines(x,p,col=1)
lines(x,1-1/(1+exp(est$fit)),col=2)
lines(x,1-1/(1+exp(est$fit+1.96*est$se)),col=3)
lines(x,1-1/(1+exp(est$fit-1.96*est$se)),col=3)
}
\keyword{gssanova}
\keyword{model}
\keyword{regression}
\keyword{smoothing}
\keyword{spline}

\name{gssanova}
\alias{gssanova}
\title{Fitting Smoothing Spline ANOVA Models with Non Gaussian Responses}
\description{
    Fit smoothing spline ANOVA models to responses from selected
    exponential families with cubic spline, linear spline, or thin-plate
    spline marginals for numerical variables.  Factors are also
    accepted.  The symbolic model specification via \code{formula}
    follows the same rules as in \code{\link{lm}} and \code{\link{glm}}.
}
\usage{
gssanova(formula, family, type="cubic", data=list(), weights, subset,
        offset, na.action=na.omit, partial=NULL, method=NULL,
        varht=1, alpha=NULL, prec=1e-7, maxiter=30, ext=.05, order=2)
}
\arguments{
    \item{formula}{Symbolic description of the model to be fit.}
    \item{family}{Description of the error distribution.  Supported
	are exponential families \code{"binomial"}, \code{"poisson"},
	\code{"Gamma"}, \code{"inverse.gaussian"}, and
	\code{"nbinomial"}.  Also supported are accelerated life model
	families \code{weibull}, \code{lognorm}, and \code{loglogis}.}
    \item{type}{Type of numerical marginals to be used.  Supported
	are \code{type="cubic"} for cubic spline marginals,
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
    \item{method}{Score used to drive the performance-oriented
	iteration.  Supported are \code{method="v"} for GCV,
	\code{method="m"} for GML, and \code{method="u"} for Mallow's CL.}
    \item{varht}{Dispersion parameter needed for \code{method="u"}.
	Ignored when \code{method="v"} or \code{method="m"} are
	specified.}
    \item{alpha}{Inverse scale parameter in accelerate life model
        families.  Ignored for exponential families.}
    \item{prec}{Precision requirement for the iterations.}
    \item{maxiter}{Maximum number of iterations allowed for
	performance-oriented iteration, and for inner-loop multiple
	smoothing parameter selection when applicable.}
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
    The models are fitted by penalized likelihood method through the
    performance-oriented iteration, as described in the reference cited
    below.

    Only one link is implemented for each \code{\link{family}}.  It is the
    logit link for \code{"binomial"}, and the log link for
    \code{"poisson"}, \code{"Gamma"}, and \code{"inverse.gaussian"}.
    For \code{"nbinomial"}, the working parameter is the logit of the
    probability \eqn{p}; see \code{\link{NegBinomial}}.

    For \code{family="binomial"}, \code{"poisson"}, \code{"nbinomial"},
    \code{"weibull"}, \code{"lognorm"}, and \code{"loglogis"}, the score
    driving the performance-oriented     iteration defaults to
    \code{method="u"} with \code{varht=1}.  For \code{family="Gamma"}
    and \code{"inverse.gaussian"}, the default is \code{method="v"}.

    See \code{\link{ssanova}} for details and notes concerning smoothing
    spline ANOVA models.
}
\note{
    For \code{family="binomial"}, the response can be specified either
    as two columns of counts or as a column of sample proportion plus a
    column of weights, as in \code{\link{glm}}.

    For \code{family="nbinomial"}, the response may be specified as two
    columns with the second being the known sizes, or simply as a single
    column with the common unknown size to be estimated through the
    maximum likelihood method.

    For \code{family="weibull"}, \code{"lognorm"}, or \code{"loglogis"},
    the response consists of three columns, with the first giving the
    follow-up time, the second the censoring status, and the third the
    left-truncation time.  For data with no truncation, the third column
    can be omitted.
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
    Methods \code{\link{predict.ssanova}},
    \code{\link{summary.gssanova}}, and \code{\link{fitted.gssanova}}.
}
\author{Chong Gu, \email{chong@stat.purdue.edu}}
\references{
    Gu, C. (1992), Cross-validating non Gaussian data. \emph{Journal
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
\keyword{models}
\keyword{regression}
\keyword{smooth}

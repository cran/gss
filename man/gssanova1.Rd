\name{gssanova1}
\alias{gssanova1}
\title{Fitting Smoothing Spline ANOVA Models with Non Gaussian Responses}
\description{
    Fit smoothing spline ANOVA models to responses from selected
    exponential families with cubic spline, linear spline, or thin-plate
    spline marginals for numerical variables.  Factors are also
    accepted.  The symbolic model specification via \code{formula}
    follows the same rules as in \code{\link{lm}} and \code{\link{glm}}.
}
\usage{
gssanova1(formula, family, type="cubic", data=list(), weights, subset,
          offset, na.action=na.omit, partial=NULL, alpha=NULL, nu=NULL,
          id.basis=NULL, nbasis=NULL, seed=NULL, random=NULL,
          ext=.05,order=2)
}
\arguments{
    \item{formula}{Symbolic description of the model to be fit.}
    \item{family}{Description of the error distribution.  Supported
	are exponential families \code{"binomial"}, \code{"poisson"},
	and \code{"Gamma"}.}
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
    \item{partial}{Optional extra unpenalized terms in partial spline
	models.}
    \item{alpha}{Tuning parameter defining cross-validation; larger
        values yield smoother fits.  Defaults are \code{alpha=1} for
	\code{family="binomial"} and \code{alpha=1.4} for
	\code{family="poisson"} and \code{family="Gamma"}.}
    \item{nu}{Input for future support of accelerated life model
        families.}
    \item{id.basis}{Index designating selected "knots".}
    \item{nbasis}{Number of "knots" to be selected.  Ignored when
	\code{id.basis} is supplied.}
    \item{seed}{Seed for reproducible random selection of "knots".
	Ignored when \code{id.basis} is supplied.}
    \item{random}{Input for parametric random effects in nonparametric
        mixed-effect models.  See \code{\link{mkran}} for details.}
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
    \code{gssanova1} implements an alternative approach to the fitting of
    \code{gssanova} models.  It works in some \eqn{q}-dimensional
    function space determined by a set of "knots" typically through
    random selection; the default dimension (\code{nbasis}) is set to be
    \eqn{q=10n^(2/9)}, adequate for cubic splines.  The algorithms are
    of the order \eqn{O(nq^2)}, scaling much better than the
    \eqn{O(n^3)} of \code{ssanova}; the memory requirements are
    \eqn{O(nq)} versus \eqn{O(n^2)}.

    \code{gssanova1} adopts direct cross-validation instead of the
    indirect cross-validation of \code{gssanova}.  Currently, only three
    families are supported, \code{"binomial"}, \code{"poisson"}, and
    \code{"Gamma"}.  Other families of \code{gssanova} shall be added in
    the future.

    Only one link is implemented for each \code{\link{family}}.  It is
    the logit link for \code{"binomial"}, and the log link for
    \code{"poisson"}, and \code{"Gamma"}.

    See \code{\link{ssanova}} for details and notes concerning smoothing
    spline ANOVA models.
}
\note{
    For \code{family="binomial"}, the response can be specified either
    as two columns of counts or as a column of sample proportions plus a
    column of total counts entered through the argument \code{weights},
    as in \code{\link{glm}}.
}
\value{
    \code{gssanova1} returns a list object of \code{\link{class}}
    \code{"gssanova1"} which inherits from the classes
    \code{"ssanova1"}, \code{"gssanova"}, and \code{"ssanova"}.

    The method \code{\link{summary}} is used to obtain summaries of the
    fits.  The method \code{\link{predict}} can be used to evaluate the
    fits at arbitrary points, along with the standard errors to be used
    in Bayesian confidence intervals, both on the scale of the link.
    The methods \code{\link{residuals}} and \code{\link{fitted.values}}
    extract the respective traits from the fits.
}
\seealso{
    \code{\link{gssanova}} and methods \code{\link{predict.ssanova1}},
    \code{\link{summary.gssanova1}}, and \code{\link{fitted.gssanova}}.
}
\author{Chong Gu, \email{chong@stat.purdue.edu}}
\references{
}
\examples{
## Fit a cubic smoothing spline logistic regression model
test <- function(x)
        {.3*(1e6*(x^11*(1-x)^6)+1e4*(x^3*(1-x)^10))-2}
x <- (0:100)/100
p <- 1-1/(1+exp(test(x)))
y <- rbinom(x,3,p)
logit.fit <- gssanova1(cbind(y,3-y)~x,family="binomial")
## The same fit
logit.fit1 <- gssanova1(y/3~x,"binomial",weights=rep(3,101),
                        id.basis=logit.fit$id.basis)
## Obtain estimates and standard errors on a grid
est <- predict(logit.fit,data.frame(x=x),se=TRUE)
## Plot the fit and the Bayesian confidence intervals
plot(x,y/3,ylab="p")
lines(x,p,col=1)
lines(x,1-1/(1+exp(est$fit)),col=2)
lines(x,1-1/(1+exp(est$fit+1.96*est$se)),col=3)
lines(x,1-1/(1+exp(est$fit-1.96*est$se)),col=3)

## Fit a mixed-effect logistic model
data(bacteriuria)
bact.fit <- gssanova1(infect~trt+time,family="binomial",data=bacteriuria,
                     id.basis=(1:820)[bacteriuria$id\%in\%c(3,38)],random=~1|id)
## Predict fixed effects
predict(bact.fit,data.frame(time=2:16,trt=as.factor(rep(1,15))),se=TRUE)
## Estimated random effects
bact.fit$b

## Clean up
\dontrun{rm(test,x,p,y,logit.fit,logit.fit1,est,bacteriuria,bact.fit)
dev.off()}
}
\keyword{models}
\keyword{regression}
\keyword{smooth}

\name{gssanova}
\alias{gssanova}
\title{Fitting Smoothing Spline ANOVA Models with Non-Gaussian Responses}
\description{
    Fit smoothing spline ANOVA models in non-Gaussian regression.  The
    symbolic model specification via \code{formula} follows the same
    rules as in \code{\link{lm}} and \code{\link{glm}}.
}
\usage{
gssanova(formula, family, type=NULL, data=list(), weights, subset,
         offset, na.action=na.omit, partial=NULL, alpha=NULL, nu=NULL,
         id.basis=NULL, nbasis=NULL, seed=NULL, random=NULL,
         skip.iter=FALSE)
}
\arguments{
    \item{formula}{Symbolic description of the model to be fit.}
    \item{family}{Description of the error distribution.  Supported
	are exponential families \code{"binomial"}, \code{"poisson"},
	\code{"Gamma"}, \code{"inverse.gaussian"}, and
	\code{"nbinomial"}.  Also supported are accelerated life model
	families \code{"weibull"}, \code{"lognorm"}, and
	\code{"loglogis"}.}
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
    \item{alpha}{Tuning parameter defining cross-validation; larger
        values yield smoother fits.  Defaults are \code{alpha=1} for
	\code{family="binomial"} and \code{alpha=1.4} otherwise.}
    \item{nu}{Inverse scale parameter in accelerated life model
        families.  Ignored for exponential families.}
    \item{id.basis}{Index designating selected "knots".}
    \item{nbasis}{Number of "knots" to be selected.  Ignored when
	\code{id.basis} is supplied.}
    \item{seed}{Seed for reproducible random selection of "knots".
	Ignored when \code{id.basis} is supplied.}
    \item{random}{Input for parametric random effects in nonparametric
        mixed-effect models.  See \code{\link{mkran}} for details.}
    \item{skip.iter}{Flag indicating whether to use initial values of
        theta and skip theta iteration.  See \code{\link{ssanova}} for
	notes on skipping theta iteration.}
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

    Only one link is implemented for each \code{family}.  It is the
    logit link for \code{"binomial"}, and the log link for
    \code{"poisson"}, and \code{"Gamma"}.  For \code{"nbinomial"}, the
    working parameter is the logit of the probability \eqn{p}; see
    \code{\link{NegBinomial}}.  For \code{"weibull"}, \code{"lognorm"},
    and \code{"loglogis"}, it is the location parameter for the log
    lifetime.

    The selection of smoothing parameters is through direct
    cross-validation.  The cross-validation score used for
    \code{family="poisson"} is taken from density estimation as in Gu
    and Wang (2003), and those used for other families are derived
    following the lines of Gu and Xiang (2001).

    A subset of the observations are selected as "knots."  Unless
    specified via \code{id.basis} or \code{nbasis}, the number of
    "knots" \eqn{q} is determined by \eqn{max(30,10n^{2/9})}, which is
    appropriate for the default cubic splines for numerical vectors.
}
\section{Responses}{
    For \code{family="binomial"}, the response can be specified either
    as two columns of counts or as a column of sample proportions plus a
    column of total counts entered through the argument \code{weights},
    as in \code{\link{glm}}.

    For \code{family="nbinomial"}, the response may be specified as two
    columns with the second being the known sizes, or simply as a single
    column with the common unknown size to be estimated through the
    maximum likelihood.

    For \code{family="weibull"}, \code{"lognorm"}, or \code{"loglogis"},
    the response consists of three columns, with the first giving the
    follow-up time, the second the censoring status, and the third the
    left-truncation time.  For data with no truncation, the third column
    can be omitted.
}
\note{
    For simpler models and moderate sample sizes, the exact solution of
    \code{\link{gssanova0}} can be faster.

    The results may vary from run to run. For consistency, specify
    \code{id.basis} or set \code{seed}.

    In \emph{gss} versions earlier than 1.0, \code{gssanova} was under
    the name \code{gssanova1}.
}
\value{
    \code{gssanova} returns a list object of class
    \code{c("gssanova","ssanova")}.

    The method \code{\link{summary.gssanova}} can be used to obtain
    summaries of the fits.  The method \code{\link{predict.ssanova}} can
    be used to evaluate the fits at arbitrary points along with standard
    errors, on the link scale.  The method
    \code{\link{project.gssanova}} can be used to calculate the
    Kullback-Leibler projection for model selection.  The methods
    \code{\link{residuals.gssanova}} and \code{\link{fitted.gssanova}}
    extract the respective traits from the fits.
}
\author{Chong Gu, \email{chong@stat.purdue.edu}}
\references{
    Gu, C. and Xiang, D. (2001), Cross validating non Gaussian data:
    generalized approximate cross validation revisited.  \emph{Journal
    of Computational and Graphical Statistics}, \bold{10}, 581--591.

    Gu, C. and Wang, J. (2003), Penalized likelihood density
    estimation: Direct cross-validation and scalable approximation.
    \emph{Statistica Sinica}, \bold{13}, 811--826.
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
logit.fit1 <- gssanova(y/3~x,"binomial",weights=rep(3,101),
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
bact.fit <- gssanova(infect~trt+time,family="binomial",data=bacteriuria,
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

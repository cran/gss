\name{ssllrm}
\alias{ssllrm}
\title{Fitting Smoothing Spline Log-Linear Regression Models}
\description{
    Fit smoothing spline log-linear regression models.  The symbolic
    model specification via \code{formula} follows the same rules as in
    \code{\link{lm}}.
}
\usage{
ssllrm(formula, response, type=NULL, data=list(), weights, subset,
       na.action=na.omit, alpha=1, id.basis=NULL, nbasis=NULL,
       seed=NULL, random=NULL, prec=1e-7, maxiter=30, skip.iter=FALSE)
}
\arguments{
    \item{formula}{Symbolic description of the model to be fit.}
    \item{response}{Formula listing response variables.}
    \item{type}{List specifying the type of spline for each variable.
        See \code{\link{mkterm}} for details.}
    \item{data}{Optional data frame containing the variables in the
	model.}
    \item{weights}{Optional vector of weights to be used in the
	fitting process.}
    \item{subset}{Optional vector specifying a subset of observations
	to be used in the fitting process.}
    \item{na.action}{Function which indicates what should happen when
	the data contain NAs.}
    \item{alpha}{Parameter modifying GCV or Mallows' CL; larger absolute
        values yield smoother fits; negative value invokes a stable and
	more accurate GCV/CL evaluation algorithm but may take two to
	five times as long.  Ignored when \code{method="m"} are
	specified.}
    \item{id.basis}{Index designating selected "knots".}
    \item{nbasis}{Number of "knots" to be selected.  Ignored when
	\code{id.basis} is supplied.}
    \item{seed}{Seed to be used for the random generation of "knots".
        Ignored when \code{id.basis} is supplied.}
    \item{random}{Input for parametric random effects in nonparametric
        mixed-effect models.  See \code{\link{mkran}} for details.}
    \item{prec}{Precision requirement for internal iterations.}
    \item{maxiter}{Maximum number of iterations allowed for
        internal iterations.}
    \item{skip.iter}{Flag indicating whether to use initial values of
        theta and skip theta iteration.  See \code{\link{ssanova}} for
	notes on skipping theta iteration.}
}
\details{
    The model is specified via \code{formula} and \code{response}, where
    \code{response} lists the response variables.  For example,
    \code{ssllrm(~y1*y2*x,~y1+y2)} prescribe a model of the form
    \deqn{
      log f(y1,y2|x) = g_{1}(y1) + g_{2}(y2) + g_{12}(y1,y2)
              + g_{x1}(x,y1) + g_{x2}(x,y2) + g_{x12}(x,y1,y2) + C(x)
    }
    with the terms denoted by \code{"y1"}, \code{"y2"}, \code{"y1:y2"},
    \code{"y1:x"}, \code{"y2:x"}, and \code{"y1:y2:x"}; the term(s) not
    involving response(s) are removed and the constant \code{C(x)} is
    determined by the fact that a conditional density integrates (adds)
    to one on the \code{y} axis.

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
    The responses, or y-variables, must be factors, and there must be at
    least one numerical x's.  For \code{response}, there is no difference
    between \code{~y1+y2} and \code{~y1*y2}.
    
    The results may vary from run to run. For consistency, specify
    \code{id.basis} or set \code{seed}.
}
\value{
    \code{ssllrm} returns a list object of class \code{"ssllrm"}.

    The method \code{\link{predict.ssllrm}} can be used to evaluate
    \code{f(y|x)} at arbitrary x, or contrasts of \code{log{f(y|x)}}
    such as the odds ratio along with standard errors.  The method
    \code{\link{project.ssllrm}} can be used to calculate the
    Kullback-Leibler projection for model selection.
}
\author{Chong Gu, \email{chong@stat.purdue.edu}}
\references{
    Gu, C. and Ma, P. (2009), Log-linear nonparametric regression
    models.  Manuscript.
}
\examples{
## Simulate data
test <- function(x)
        {.3*(1e6*(x^11*(1-x)^6)+1e4*(x^3*(1-x)^10))-2}
x <- (0:100)/100
p <- 1-1/(1+exp(test(x)))
y <- rbinom(x,3,p)
y1 <- as.ordered(y)
y2 <- as.factor(rbinom(x,1,p))
## Fit model
fit <- ssllrm(~y1*y2*x,~y1+y2)

## Evaluate f(y|x)
est <- predict(fit,data.frame(x=x),
               data.frame(y1=as.factor(0:3),y2=as.factor(rep(0,4))))
## f(y|x) at all y values (fit$qd.pt)
est <- predict(fit,data.frame(x=x))

## Evaluate contrast of log f(y|x)
est <- predict(fit,data.frame(x=x),odds=c(-1,.5,.5,0),
               data.frame(y1=as.factor(0:3),y2=as.factor(rep(0,4))),se=TRUE)
## Odds ratio log{f(0,0|x)/f(3,0|x)}
est <- predict(fit,data.frame(x=x),odds=c(1,-1),
               data.frame(y1=as.factor(c(0,3)),y2=as.factor(c(0,1))),se=TRUE)

## KL projection
kl <- project(fit,include=c("y2:x","y1:y2","y1:x","y2:x"))

## Clean up
\dontrun{rm(test,x,p,y,y1,y2,fit,est,kl)
dev.off()}
}
\keyword{smooth}
\keyword{models}
\keyword{regression}

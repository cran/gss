\name{sshzd}
\alias{sshzd}
\title{Estimating Hazard Function Using Smoothing Splines}
\description{
    Estimate hazard function using smoothing spline ANOVA models
    with cubic spline, linear spline, or thin-plate spline marginals for
    numerical variables.  The symbolic model specification via
    \code{formula} follows the same rules as in \code{\link{lm}}, but
    with the response of a special form.
}
\usage{
sshzd(formula, type="cubic", data=list(), alpha=1.4, weights=NULL,
      subset, na.action=na.omit, id.basis=NULL, nbasis=NULL, seed=NULL,
      ext=.05, order=2, prec=1e-7, maxiter=30)
}
\arguments{
    \item{formula}{Symbolic description of the model to be fit.  Details
        are given below.}
    \item{type}{Type of numerical marginals to be used.  Supported are
	\code{type="cubic"} for cubic spline marginals,
	\code{type="linear"} for linear spline marginals, and
	\code{type="tp"} for thin-plate spline marginals.}
    \item{data}{Optional data frame containing the variables in the
        model.}
    \item{alpha}{Parameter defining cross-validation score for smoothing
        parameter selection.}
    \item{weights}{Optional vector of bin-counts for histogram data.}
    \item{subset}{Optional vector specifying a subset of observations
	to be used in the fitting process.}
    \item{na.action}{Function which indicates what should happen when
        the data contain NAs.}
    \item{id.basis}{Index of observations to be used as "knots."}
    \item{nbasis}{Number of "knots" to be used.  Ignored when
        \code{id.basis} is specified.}
    \item{seed}{Seed to be used for the random generation of "knots."
        Ignored when \code{id.basis} is specified.}
    \item{ext}{For cubic spline and linear spline marginals, this option
	specifies how far to extend the domain beyond the minimum and
	the maximum as a percentage of the range.  The default
	\code{ext=.05} specifies marginal domains of lengths 110 percent
	of their respective ranges.  Evaluation outside of the domain
	will result in an error.  Ignored if \code{type="tp"} or
	\code{domain} are specified.}
    \item{order}{For thin-plate spline marginals, this option specifies
	the order of the marginal penalties.  Ignored if
	\code{type="cubic"} or \code{type="linear"} are specified.}
    \item{prec}{Precision requirement for internal iterations.}
    \item{maxiter}{Maximum number of iterations allowed for
	internal iterations.}
}
\details{
    The model specification via \code{formula} is for the log hazard.
    For example, \code{~x1*x2} prescribes a model of the form
    \deqn{
	log f(x1,x2) = C + g_{1}(x1) + g_{2}(x2) + g_{12}(x1,x2)
    }
    with the terms denoted by \code{"x1"}, \code{"x2"}, and
    \code{"x1:x2"}.

    \code{sshzd} takes standard right-censored lifetime data, with
    possible left-truncation and covariates.  The response in
    \code{formula} must be of the form
    \code{Surv(futime,status,start=0)}, where \code{futime} is the 
    follow-up time, \code{status} is the censoring indicator, and
    \code{start} is the optional left-truncation time.  The function
    \code{Surv} is defined and parsed inside \code{sshzd}, not quite the
    same as the one in the \code{survival} package.

    The main effect of \code{futime} must appear in the model terms.
    The absence of interactions between \code{futime} and covariates
    characterizes proportional hazard models.

    Parallel to those in a \code{\link{ssanova}} object, the model terms
    are sums of unpenalized and penalized terms.  Attached to every
    penalized term there is a smoothing parameter, and the model
    complexity is largely determined by the number of smoothing
    parameters.

    The selection of smoothing parameters is through a cross-validation
    mechanism described in the references, with a parameter
    \code{alpha}; \code{alpha=1} is "unbiased" for the minimization of
    Kullback-Leibler loss but may yield severe undersmoothing, whereas
    larger \code{alpha} yields smoother estimates.

    A subset of the observations are selected as "knots."  Unless
    specified via \code{id.basis} or \code{nbasis}, the subset size is
    determined by \eqn{max(30,10n^(2/9))}, which is appropriate for
    \code{type="cubic"} but not necessarily for \code{type="linear"} or
    \code{type="tp"}.
}
\note{
    Integration on the time axis is done by the 200-point Gauss-Legendre
    formula on [0,T], where T is the largest follow-up time.
}
\value{
    \code{sshzd} returns a list object of \code{\link{class} "sshzd"}.

    \code{\link{hzdrate.sshzd}} can be used to evaluate the estimated
    hazard function.  \code{\link{hzdcurve.sshzd}} can be used to
    evaluate hazard curves with fixed covariates.
    \code{\link{survexp.sshzd}} can be used to calculated estimated
    expected survival.
}
\seealso{
    \code{\link{hzdrate.sshzd}}, \code{\link{hzdcurve.sshzd}}, and
    \code{\link{survexp.sshzd}}.
}
\author{Chong Gu, \email{chong@stat.purdue.edu}}
\references{
    Gu, C. and Wang, J. (2002), \emph{Penalized Likelihood Density
    Estimation: Direct Cross-Validation and Scalable Approximation}.
    Available at \url{http://stat.purdue.edu/~chong/manu.html}.
  
    Gu, C. (2002), \emph{Smoothing Spline ANOVA Models}.  New York:
    Springer-Verlag.
}
\examples{
## Model with interaction
data(gastric)
gastric.fit <- sshzd(Surv(futime,status)~futime*trt,data=gastric)
## exp(-Lambda(600)), exp(-(Lambda(1200)-Lambda(600))), and exp(-Lambda(1200))
survexp.sshzd(gastric.fit,c(600,1200,1200),data.frame(trt=as.factor(1)),c(0,600,0))
## Clean up
\dontrun{rm(gastric,gastric.fit)
dev.off()}

## THE FOLLOWING EXAMPLE IS TIME-CONSUMING
## Proportional hazard model
\dontrun{
data(stan)
stan.fit <- sshzd(Surv(futime,status)~futime+age,data=stan)
## Evaluate fitted hazard
hzdrate.sshzd(stan.fit,data.frame(futime=c(10,20),age=c(20,30)))
## Plot lambda(t,age=20)
tt <- seq(0,60,leng=101)
hh <- hzdcurve.sshzd(stan.fit,tt,data.frame(age=20))
plot(tt,hh,type="l")
## Clean up
rm(stan,stan.fit,tt,hh)
dev.off()
}
}
\keyword{smooth}
\keyword{models}
\keyword{survival}

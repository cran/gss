\name{sshzd}
\alias{sshzd}
\alias{sshzd1}
\title{Estimating Hazard Function Using Smoothing Splines}
\description{
    Estimate hazard function using smoothing spline ANOVA models.  The
    symbolic model specification via \code{formula} follows the same
    rules as in \code{\link{lm}}, but with the response of a special
    form.
}
\usage{
sshzd(formula, type=NULL, data=list(), alpha=1.4, weights=NULL,
      subset, offset, na.action=na.omit, partial=NULL, id.basis=NULL,
      nbasis=NULL, seed=NULL, random=NULL, prec=1e-7, maxiter=30,
      skip.iter=FALSE)
sshzd1(formula, type=NULL, data=list(), alpha=1.4, weights=NULL,
       subset, na.action=na.omit, rho=list("marginal"), partial=NULL,
       id.basis=NULL, nbasis=NULL, seed=NULL, random=NULL, prec=1e-7,
       maxiter=30, skip.iter=FALSE)
}
\arguments{
    \item{formula}{Symbolic description of the model to be fit, where
        the response is of the form \code{Surv(futime,status,start=0)}.}
    \item{type}{List specifying the type of spline for each variable.
        See \code{\link{mkterm}} for details.}
    \item{data}{Optional data frame containing the variables in the
        model.}
    \item{alpha}{Parameter defining cross-validation score for smoothing
        parameter selection.}
    \item{weights}{Optional vector of counts for duplicated data.}
    \item{subset}{Optional vector specifying a subset of observations
	to be used in the fitting process.}
    \item{offset}{Optional offset term with known parameter 1.}
    \item{na.action}{Function which indicates what should happen when
        the data contain NAs.}
    \item{partial}{Optional symbolic description of parametric terms in
        partial spline models.}
    \item{id.basis}{Index of observations to be used as "knots."}
    \item{nbasis}{Number of "knots" to be used.  Ignored when
        \code{id.basis} is specified.}
    \item{seed}{Seed to be used for the random generation of "knots."
        Ignored when \code{id.basis} is specified.}
    \item{random}{Input for parametric random effects (frailty) in
        nonparametric mixed-effect models.  See \code{\link{mkran}} for
	details.}
    \item{prec}{Precision requirement for internal iterations.}
    \item{maxiter}{Maximum number of iterations allowed for
	internal iterations.}
    \item{skip.iter}{Flag indicating whether to use initial values of
        theta and skip theta iteration.  See \code{\link{ssanova}} for
	notes on skipping theta iteration.}
    \item{rho}{rho function needed for sshzd1.}
}
\details{
    The model specification via \code{formula} is for the log hazard.
    For example, \code{Suve(t,d)~t*u} prescribes a model of the form
    \deqn{
	log f(t,u) = C + g_{t}(t) + g_{u}(u) + g_{t,u}(t,u)
    }
    with the terms denoted by \code{"1"}, \code{"t"}, \code{"u"}, and
    \code{"t:u"}.  Replacing \code{t*u} by \code{t+u} in the
    \code{formula}, one gets a proportional hazard model with
    \eqn{g_{t,u}=0}.

    \code{sshzd} takes standard right-censored lifetime data, with
    possible left-truncation and covariates; in
    \code{Surv(futime,status,start=0)~...}, \code{futime} is the  
    follow-up time, \code{status} is the censoring indicator, and
    \code{start} is the optional left-truncation time.  The main effect
    of \code{futime} must appear in the model terms specified via
    \code{...}.

    Parallel to those in a \code{\link{ssanova}} object, the model terms
    are sums of unpenalized and penalized terms.  Attached to every
    penalized term there is a smoothing parameter, and the model
    complexity is largely determined by the number of smoothing
    parameters.

    The selection of smoothing parameters is through a cross-validation
    mechanism described in Gu (2002, Sec. 7.2), with a parameter
    \code{alpha}; \code{alpha=1} is "unbiased" for the minimization of
    Kullback-Leibler loss but may yield severe undersmoothing, whereas
    larger \code{alpha} yields smoother estimates.

    A subset of the observations are selected as "knots."  Unless
    specified via \code{id.basis} or \code{nbasis}, the number of
    "knots" \eqn{q} is determined by \eqn{max(30,10n^{2/9})}, which is
    appropriate for the default cubic splines for numerical vectors.
}
\note{
    The function \code{Surv(futime,status,start=0)} is defined and
    parsed inside \code{sshzd}, not quite the same as the one in the
    \code{survival} package.

    Integration on the time axis is done by the 200-point Gauss-Legendre
    formula on \code{c(min(start),max(futime))}, returned from
    \code{\link{gauss.quad}}.

    \code{sshzd1} can be up to 50 times faster than \code{sshzd}, at the
    cost of performance degradation.

    The results may vary from run to run.  For consistency, specify
    \code{id.basis} or set \code{seed}.
}
\value{
    \code{sshzd} returns a list object of class \code{"sshzd"}.
    \code{sshzd1} returns a list object of class
    \code{c("sshzd1","sshzd")}.

    \code{\link{hzdrate.sshzd}} can be used to evaluate the estimated
    hazard function.  \code{\link{hzdcurve.sshzd}} can be used to
    evaluate hazard curves with fixed covariates.
    \code{\link{survexp.sshzd}} can be used to calculated estimated
    expected survival.

    The method \code{\link{project.sshzd}} can be used to calculate the
    Kullback-Leibler projection of \code{"sshzd"} objects for model
    selection; \code{\link{project.sshzd1}} can be used to calculate the
    square error projection of \code{"sshzd1"} objects.
}
\author{Chong Gu, \email{chong@stat.purdue.edu}}
\references{
    Gu, C. (2002), \emph{Smoothing Spline ANOVA Models}.  New York:
    Springer-Verlag.

    Du, P. and Gu, C. (2006), Penalized likelihood hazard estimation:
    efficient approximation and Bayesian confidence intervals.
    \emph{Statistics and Probability Letters}, \bold{76}, 244--254.

    Du, P. and Gu, C. (2009), Penalized Pseudo-Likelihood Hazard
    Estimation: A Fast Alternative to Penalized Likelihood.
    \emph{Journal of Statistical Planning and Inference}, \bold{139},
    891--899.

    Du, P. and Ma, S. (2010), Frailty Model with Spline Estimated
    Nonparametric Hazard Function, \emph{Statistica Sinica}, \bold{20},
    561--580.
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

\name{sscox}
\alias{sscox}
\title{Estimating Relative Risk Using Smoothing Splines}
\description{
    Estimate relative risk using smoothing spline ANOVA models.  The
    symbolic model specification via \code{formula} follows the same
    rules as in \code{\link{lm}}, but with the response of a special
    form.
}
\usage{
sscox(formula, type=NULL, data=list(), weights=NULL, subset,
      na.action=na.omit, partial=NULL, alpha=1.4, id.basis=NULL,
      nbasis=NULL, seed=NULL, random=NULL, prec=1e-7, maxiter=30,
      skip.iter=FALSE)
}
\arguments{
    \item{formula}{Symbolic description of the model to be fit, where
        the response is of the form \code{Surv(futime,status,start=0)}.}
    \item{type}{List specifying the type of spline for each variable.
        See \code{\link{mkterm}} for details.}
    \item{data}{Optional data frame containing the variables in the
        model.}
    \item{weights}{Optional vector of counts for duplicated data.}
    \item{subset}{Optional vector specifying a subset of observations
	to be used in the fitting process.}
    \item{na.action}{Function which indicates what should happen when
        the data contain NAs.}
    \item{partial}{Optional symbolic description of parametric terms in
        partial spline models.}
    \item{alpha}{Parameter defining cross-validation score for smoothing
        parameter selection.}
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
}
\details{
    A proportional hazard model is assumed, and the relative risk is
    estimated via penalized partial likelihood.  The model specification
    via \code{formula} is for the log relative risk.  For example,
    \code{Suve(t,d)~u*v} prescribes a model of the form
    \deqn{
	log f(u,v) = g_{u}(u) + g_{v}(v) + g_{u,v}(u,v)
    }
    with the terms denoted by \code{"u"}, \code{"v"}, and \code{"u:v"};
    relative risk is defined only up to a multiplicative constant, so
    the constant term is not included in the model.

    \code{sscox} takes standard right-censored lifetime data, with
    possible left-truncation and covariates; in
    \code{Surv(futime,status,start=0)~...}, \code{futime} is the  
    follow-up time, \code{status} is the censoring indicator, and
    \code{start} is the optional left-truncation time.

    Parallel to those in a \code{\link{ssanova}} object, the model terms
    are sums of unpenalized and penalized terms.  Attached to every
    penalized term there is a smoothing parameter, and the model
    complexity is largely determined by the number of smoothing
    parameters.

    The selection of smoothing parameters is through a cross-validation
    mechanism designed for density estimation under biased sampling,
    with a fudge factor \code{alpha}; \code{alpha=1} is "unbiased" for
    the minimization of Kullback-Leibler loss but may yield severe
    undersmoothing, whereas larger \code{alpha} yields smoother
    estimates.

    A subset of the observations are selected as "knots."  Unless
    specified via \code{id.basis} or \code{nbasis}, the number of
    "knots" \eqn{q} is determined by \eqn{max(30,10n^{2/9})}, which is
    appropriate for the default cubic splines for numerical vectors.
}
\note{
    The function \code{Surv(futime,status,start=0)} is defined and
    parsed inside \code{sscox}, not quite the same as the one in the
    \code{survival} package.  The estimation is invariant of monotone
    transformations of time.

    The results may vary from run to run.  For consistency, specify
    \code{id.basis} or set \code{seed}.
}
\value{
    \code{sscox} returns a list object of class \code{"sscox"}.

    The method \code{\link{predict.sscox}} can be used to evaluate the
    fits at arbitrary points along with standard errors.  The method
    \code{\link{project.sscox}} can be used to calculate the
    Kullback-Leibler projection for model selection.
}
\author{Chong Gu, \email{chong@stat.purdue.edu}}
\references{
    Gu, C. (2013), \emph{Smoothing Spline ANOVA Models (2nd Ed)}.  New
    York: Springer-Verlag.

    Gu, C. (2014), Smoothing Spline ANOVA Models: R Package gss.
    \emph{Journal of Statistical Software}, 58(5), 1-25. URL
    http://www.jstatsoft.org/v58/i05/.
}
\examples{
## Relative Risk
data(stan)
fit.rr <- sscox(Surv(futime,status)~age,data=stan)
est.rr <- predict(fit.rr,data.frame(age=c(35,40)),se=TRUE)
## Base Hazard
risk <- predict(fit.rr,stan)
fit.bh <- sshzd(Surv(futime,status)~futime,data=stan,offset=log(risk))
tt <- seq(0,max(stan$futime),length=51)
est.bh <- hzdcurve.sshzd(fit.bh,tt,se=TRUE)
## Clean up
\dontrun{rm(stan,fit.rr,est.rr,risk,fit.bh,tt,est.bh)}
}
\keyword{smooth}
\keyword{models}
\keyword{survival}

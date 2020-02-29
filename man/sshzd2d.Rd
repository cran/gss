\name{sshzd2d}
\alias{sshzd2d}
\alias{sshzd2d1}
\title{Estimating 2-D Hazard Function Using Smoothing Splines}
\description{
    Estimate 2-D hazard function using smoothing spline ANOVA models.
}
\usage{
sshzd2d(formula1, formula2, symmetry=FALSE, data, alpha=1.4,
        weights=NULL, subset=NULL, id.basis=NULL, nbasis=NULL, seed=NULL,
        prec=1e-7, maxiter=30, skip.iter=FALSE)

sshzd2d1(formula1, formula2, symmetry=FALSE, data, alpha=1.4,
         weights=NULL, subset=NULL, rho="marginal",
         id.basis=NULL, nbasis=NULL, seed=NULL, prec=1e-7, maxiter=30,
         skip.iter=FALSE)
}
\arguments{
    \item{formula1}{Description of the hazard model to be fit on the
        first axis.}
    \item{formula2}{Description of the hazard model to be fit on the
        second axis.}
    \item{symmetry}{Flag indicating whether to enforce symmetry of the
        two axes.}
    \item{data}{Data frame containing the variables in the model.}
    \item{alpha}{Parameter defining cross-validation scores for smoothing
        parameter selection.}
    \item{weights}{Optional vector of counts for duplicated data.}
    \item{subset}{Optional vector specifying a subset of observations
	to be used in the fitting process.}
    \item{id.basis}{Index of observations to be used as "knots."}
    \item{nbasis}{Number of "knots" to be used.  Ignored when
        \code{id.basis} is specified.}
    \item{seed}{Seed to be used for the random generation of "knots."
        Ignored when \code{id.basis} is specified.}
    \item{prec}{Precision requirement for internal iterations.}
    \item{maxiter}{Maximum number of iterations allowed for
	internal iterations.}
    \item{skip.iter}{Flag indicating whether to use initial values of
        theta and skip theta iteration in marginal hazard estimation.}
    \item{rho}{Choice of rho function for sshzd2d1: \code{"marginal"} or
        \code{"weibull"}.}
}
\details{
    The 2-D survival function is expressed as
    \eqn{S(t1,t2)=C(S1(t1),S2(t2))}, where \eqn{S1(t1)}, \eqn{S2(t2)}
    are marginal survival functions and \eqn{C(u1,u2)} is a 2-D copula.
    The marginal survival functions are estimated via the marginal
    hazards as in \code{\link{sshzd}}, and the copula is estimated
    nonparametrically by calling \code{\link{sscopu2}}.

    When \code{symmetry=TRUE}, a common marginal survial function
    S1(t)=S2(t) is estimated, and a symmetric copula is estimated such
    that \eqn{C(u1,u2)=C(u2,u1)}.

    Covariates can be incorporated in the marginal hazard models as in
    \code{\link{sshzd}}, including parametric terms via \code{partial}
    and frailty terms via \code{random}.  Arguments \code{formula1} and
    \code{formula2} are typically model formulas of the same form as the
    argument \code{formula} in \code{\link{sshzd}}, but when
    \code{partial} or \code{random} are needed, \code{formula1} and
    \code{formula2} should be lists with model formulas as the first
    elements and \code{partial}/\code{random} as named elements; when
    necessary, variable configurations (that are done via argument
    \code{type} in \code{\link{sshzd}}) should also be entered as named
    elements of lists \code{formula1}/\code{formula2}.

    When \code{symmetry=TRUE}, parallel model formulas must be
    consistent of each other, such as
    \tabular{l}{
    \code{formula1=list(Surv(t1,d1)~t1*u1,partial=~z1,random=~1|id1)}\cr
    \code{formula2=list(Surv(t2,d2)~t2*u2,partial=~z2,random=~1|id2)}
    }
    where pairs \code{t1}-\code{t2}, \code{d2}-\code{d2} respectively
    are different elements in \code{data}, pairs \code{u1}-\code{u2},
    \code{z1}-\code{z2} respectively may or may not be different
    elements in \code{data}, and factors \code{id1} and \code{id2}
    are typically the same but at least should have the same levels.
}
\note{
    \code{sshzd2d1} executes faster than \code{sshzd2d}, but often at
    the cost of performance degradation.
   
    The results may vary from run to run.  For consistency, specify
    \code{id.basis} or set \code{seed}.
}
\value{
    \code{sshzd2d} and \code{sshzd2d1} return a list object of class
    \code{"sshzd2d"}.

    \code{\link{hzdrate.sshzd2d}} can be used to evaluate the estimated
    2-D hazard function.  \code{\link{survexp.sshzd2d}} can be used to
    calculate estimated survival functions.
}
\author{Chong Gu, \email{chong@stat.purdue.edu}}
\references{
    Gu, C. (2015), Hazard estimation with bivariate survival data and
    copula density estimation. \emph{Journal of Computational and
    Graphical Statistics}, \bold{24}, 1053-1073.
}
\examples{
## THE FOLLOWING EXAMPLE IS TIME-CONSUMING
\dontrun{
data(DiaRet)
## Common proportional hazard model on the margins
fit <- sshzd2d(Surv(time1,status1)~time1+trt1*type,
               Surv(time2,status2)~time2+trt2*type,
               data=DiaRet,symmetry=TRUE)
## Evaluate fitted survival and hazard functions
time <- cbind(c(50,70),c(70,70))
cova <- data.frame(trt1=as.factor(c(1,1)),trt2=as.factor(c(1,0)),
                   type=as.factor(c("juvenile","adult")))
survexp.sshzd2d(fit,time,cov=cova)
hzdrate.sshzd2d(fit,time,cov=cova)
## Association between margins: Kendall's tau and Spearman's rho
summary(fit$copu)
## Clean up
rm(DiaRet,fit,time,cova)
dev.off()
}
}
\keyword{smooth}
\keyword{models}
\keyword{survival}

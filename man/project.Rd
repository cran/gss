\name{project}
\alias{project}
\alias{project.ssanova1}
\alias{project.gssanova1}
\alias{project.ssden}
\alias{project.sshzd}
\title{Projecting Smoothing Spline ANOVA Fits for Model Diagnostics}
\description{
    Calculate Kullback-Leibler projection of smoothing spline ANOVA fits
    for Model Diagnostics.
}
\usage{
project(object, ...)
project.ssanova1(object, include, ...)
project.gssanova1(object, include, ...)
project.ssden(object, include, mesh=FALSE, ...)
project.sshzd(object, include, mesh=FALSE, ...)
}
\arguments{
    \item{object}{Object of class \code{"ssanova1"}, \code{"gssanova1"},
        \code{"ssden"}, or \code{"sshzd"}.}
    \item{...}{Additional arguments.  Ignored in \code{project.x}.}
    \item{include}{List of model terms to be included in the reduced
        model space.  The \code{partial} and \code{offset} terms, if
	present, are to be specified by \code{"partial"} and
	\code{"offset"}, respectively.}
    \item{mesh}{Flag indicating whether to return evaluations of the
        projection.}
}
\value{
    The functions return a list consisting of the following components.
    \item{ratio}{KL(fit0,fit1)/KL(fit0,null); the smaller the value, the
        more feasible the reduced model is.}
    \item{kl}{KL(fit0,fit1).}

    For regression fits, the list also contains the following component.
    \item{check}{KL(fit0,fit1)/KL(fit0,null)+KL(fit1,null)/KL(fit0,null);
        a value closer to 1 is preferred.}

    For density and hazard fits, the list may contain the following
    optional component.
    \item{mesh}{The evaluations of the projection.}
}
\details{
    The entropy KL(fit0,null) can be decomposed as the sum of
    KL(fit0,fit1) and KL(fit1,null), where fit0 is the fit to be
    projected, fit1 is the projection in the reduced model space, and
    null is the constant fit.  The ratio KL(fit0,fit1)/KL(fit0,null)
    serves as a diagnostic of the feasibility of the reduced model.

    For regression fits, smoothness safe-guard is used to prevent
    interpolation, and KL(fit0,fit1)+KL(fit1,null) may not match
    KL(fit0,null) perfectly.

    For mixed-effect models from \code{ssanova1} and \code{gssanova1},
    the estimated random effects are treated as offset.
}
\author{Chong Gu, \email{chong@stat.purdue.edu}}
\references{
    Gu, C. (2003), \emph{Model Diagnostics for Smoothing Spline ANOVA
    Models}.  Available at
    \url{http://stat.purdue.edu/~chong/manu.html}.
}
\seealso{
    Fitting functions \code{\link{ssanova1}}, \code{\link{gssanova1}},
    \code{\link{ssden}}, and \code{\link{sshzd}}.
}
\keyword{models}
\keyword{smooth}
\keyword{htest}

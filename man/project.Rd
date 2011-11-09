\name{project}
\alias{project}
\alias{project.ssanova}
\alias{project.ssanova9}
\alias{project.gssanova}
\alias{project.ssden}
\alias{project.ssden1}
\alias{project.sscden}
\alias{project.sscden1}
\alias{project.sshzd}
\alias{project.sscox}
\alias{project.sshzd1}
\alias{project.ssllrm}
\title{Projecting Smoothing Spline ANOVA Fits for Model Diagnostics}
\description{
    Calculate Kullback-Leibler projection of smoothing spline ANOVA fits
    for model diagnostics.
}
\usage{
project(object, ...)
\method{project}{ssanova}(object, include, ...)
\method{project}{ssanova9}(object, include, ...)
\method{project}{gssanova}(object, include, ...)
\method{project}{ssden}(object, include, mesh=FALSE, ...)
\method{project}{ssden1}(object, include, ...)
\method{project}{sscden}(object, include, ...)
\method{project}{sscden1}(object, include, ...)
\method{project}{sshzd}(object, include, mesh=FALSE, ...)
\method{project}{sscox}(object, include, ...)
\method{project}{sshzd1}(object, include, ...)
\method{project}{ssllrm}(object, include, ...)
}
\arguments{
    \item{object}{Object of class \code{"ssanova"}, \code{"gssanova"},
        \code{"ssden"}, \code{"ssden1"}, \code{"sscden"},
	\code{"sscden1"}, \code{"sshzd"}, \code{"sshzd1"}, or
	\code{"ssllrm"}.}
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
\note{
    \code{project.ssden1}, \code{project.sscden1}, and
    \code{project.sshzd1} calculates square error projections.
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

    For mixed-effect models from \code{ssanova} and \code{gssanova},
    the estimated random effects are treated as offset.
}
\author{Chong Gu, \email{chong@stat.purdue.edu}}
\references{
    Gu, C. (2004), Model diagnostics for smoothing spline ANOVA models.
    \emph{The Canadian Journal of Statistics}, \bold{32}, 347--358.
}
\seealso{
    Fitting functions \code{\link{ssanova}}, \code{\link{gssanova}},
    \code{\link{ssden}}, \code{\link{sshzd}}, and \code{\link{sshzd1}}.
}
\keyword{models}
\keyword{smooth}
\keyword{htest}

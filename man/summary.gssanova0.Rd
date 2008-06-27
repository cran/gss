\name{summary.gssanova0}
\alias{summary.gssanova0}
\title{Assessing Smoothing Spline ANOVA Fits with Non-Gaussian Responses}
\description{
    Calculate various summaries of smoothing spline ANOVA fits with
    non-Gaussian responses.
}
\usage{
\method{summary}{gssanova0}(object, diagnostics=FALSE, ...)
}
\arguments{
    \item{object}{Object of class \code{"gssanova"}.}
    \item{diagnostics}{Flag indicating if diagnostics are required.}
    \item{...}{Ignored.}
}
\value{
    \code{summary.gssanova0} returns a list object of \code{\link{class}}
    \code{"summary.gssanova0"} consisting of the following components.
    The entries \code{pi}, \code{kappa}, \code{cosines}, and
    \code{roughness} are only calculated if \code{diagnostics=TRUE}.
    \item{call}{Fitting call.}
    \item{family}{Error distribution.}
    \item{method}{Method for smoothing parameter selection.}
    \item{dispersion}{Assumed or estimated dispersion parameter.}
    \item{iter}{Number of performance-oriented iterations performed.}
    \item{fitted}{Fitted values on the link scale.}
    \item{residuals}{Working residuals on the link scale.}
    \item{rss}{Residual sum of squares.}
    \item{dev.resid}{Deviance residuals.}
    \item{deviance}{Deviance of the fit.}
    \item{dev.null}{Deviance of the null model.}
    \item{alpha}{Estimated size for \code{family="nbinomial"} with one
        column responses.  Estimated inverse scale of log life time for
	\code{family="nbinomial"}, \code{"lognorm"}, or
	\code{"loglogis"}.}
    \item{penalty}{Roughness penalty associated with the fit.}
    \item{pi}{"Percentage decomposition" of "explained variance" into
	model terms.}
    \item{kappa}{Concurvity diagnostics for model terms.  Virtually the
	square roots of variance inflation factors of a retrospective
	linear model.}
    \item{cosines}{Cosine diagnostics for practical significance of
	model terms.}
    \item{roughness}{Percentage decomposition of the roughness penalty
      \code{penalty} into model terms.}
}
\details{
    Similar to the iterated weighted least squares fitting of
    \code{\link{glm}}, penalized likelihood regression fit can be calculated
    through iterated penalized weighted least squares.

    The diagnostics are based on the "pseudo" Gaussian response model
    behind the weighted least squares problem at convergence.
}
\author{Chong Gu, \email{chong@stat.purdue.edu}}
\references{
    Gu, C. (1992), Diagnostics for nonparametric regression models with
    additive terms.  \emph{Journal of the American Statistical
	Association}, \bold{87}, 1051--1058.
}
\seealso{
    Fitting function \code{\link{gssanova0}} and methods
    \code{\link{predict.ssanova0}}, \code{\link{fitted.gssanova}}.
}
\keyword{models}
\keyword{regression}
\keyword{smooth}

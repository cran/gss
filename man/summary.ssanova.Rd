\name{summary.ssanova}
\alias{summary.ssanova}
\title{Assessing Smoothing Spline ANOVA Fits}
\description{
    Calculate various summaries of smoothing spline ANOVA fits.
}
\usage{
summary.ssanova(obj, diagnostics=FALSE)
}
\arguments{
    \item{obj}{Object of class \code{"ssanova"}.}
    \item{diagnostics}{Flag indicating if diagnostics are required.}
}
\value{
    \code{summary.ssanova} returns a list object of \code{\link{class}}
    \code{"summary.ssanova"} consisting of the following components.
    The entries \code{pi}, \code{kappa}, \code{cosines}, and
    \code{roughness} are only calculated if \code{diagnostics=TRUE}; see
    the reference below for details concerning the diagnostics.
    \item{call}{Fitting call.}
    \item{method}{Method for smoothing parameter selection.}
    \item{fitted}{Fitted values.}
    \item{residuals}{Residuals.}
    \item{sigma}{Assumed or estimated error standard deviation.}
    \item{r.squared}{Fraction of "explained variance" by the fitted model.}
    \item{rss}{Residual sum of squares.}
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
\author{Chong Gu, \email{chong@stat.purdue.edu}}
\references{
    Gu, C. (1992), Diagnostics for nonparametric regression models with
    additive terms.  \emph{Journal of the American Statistical
	Association}, \bold{87}, 1051--1058.
}
\seealso{
    Fitting function \code{\link{ssanova}} and methods
    \code{\link{predict.ssanova}}, \code{\link{fitted.ssanova}}.
}
\keyword{models}
\keyword{regression}
\keyword{smooth}

\name{summary.ssanova}
\title{Assessing Smoothing Spline ANOVA Fits}
\usage{
summary[.ssanova](obj, diagnostics=FALSE)
}
\arguments{
    \item{obj}{an object of class \code{"ssanova"}.}
    \item{diagnostics}{a logical flag.}
}
\description{
    \code{summary.ssanova} calculates various summaries of smoothing
    spline ANOVA fits.
}
\value{
    \code{summary.ssanova} returns a list object of \code{\link{class}}
    \code{"summary.ssanova"} consisting of the following components.
    The entries \code{kappa}, \code{cosines}, and \code{roughness} are
    only calculated for \code{diagnostics=TRUE}.
    \item{call}{the fitting call.}
    \item{method}{the smoothing parameter selection method.}
    \item{fitted}{the fitted values.}
    \item{residuals}{the residuals.}
    \item{sigma}{the square root of the estimated variance of the random
	error.}
    \item{r.squared}{the "fraction of variance explained by the model".}
    \item{rss}{the residual sum of squares.}
    \item{penalty}{the penalty associated with the fit.}
    \item{kappa}{the concurvity diagnostics for model terms.  These are
	virtually the variance inflation factors of a retrospective
	linear model.}
    \item{cosines}{the cosine diagnostics for practical significance of
	the model terms.}
    \item{roughness}{the roughness of individual model terms as
	percentages of the overall roughness, which is proportional to
	\code{penalty}.}
}
\author{Chong Gu, \email{chong@stat.purdue.edu}}
\seealso{
    The model fitting function \code{\link{ssanova}} and the predicting
    function \code{\link{predict.ssanova}}.
}
\keyword{ssanova}
\keyword{models}
\keyword{regression}
\keyword{smoothing}
\keyword{spline}

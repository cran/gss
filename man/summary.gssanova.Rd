\name{summary.gssanova}
\title{Assessing Smoothing Spline ANOVA Fits with Non Gaussian Responses}
\usage{
summary[.gssanova](obj, diagnostics=FALSE)
}
\arguments{
    \item{obj}{an object of class \code{"gssanova"}.}
    \item{diagnostics}{a logical flag.}
}
\description{
    \code{summary.gssanova} calculates various summaries of smoothing
    spline ANOVA fits with non Gaussian responses.
}
\value{
    \code{summary.gssanova} returns a list object of \code{\link{class}}
    \code{"summary.gssanova"} consisting of the following components.
    The entries \code{kappa}, \code{cosines}, and \code{roughness} are
    only calculated for \code{diagnostics=TRUE}.
    \item{call}{the fitting call.}
    \item{family}{the error distribution.}
    \item{method}{the smoothing parameter selection method.}
    \item{dispersion}{the assumed or estimated dispersion parameter.}
    \item{iter}{the number of performance-oriented iterations performed.}
    \item{fitted}{the fitted values on the scale of link.}
    \item{residuals}{the working residuals.}
    \item{rss}{the residual sum of squares.}
    \item{dev.resid}{the deviance residuals.}
    \item{deviance}{the deviance of the fit.}
    \item{dev.null}{the deviance of the null model.}
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
    The model fitting function \code{\link{gssanova}} and the predicting
    function \code{\link{predict.ssanova}}.
}
\keyword{gssanova}
\keyword{models}
\keyword{regression}
\keyword{smoothing}
\keyword{spline}

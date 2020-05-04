\name{predict9.gssanova}
\alias{predict9.gssanova}
\alias{predict9}
\title{Predicting from Smoothing Spline ANOVA Fits with Non-Gaussian Responses}
\description{
    Evaluate smoothing spline ANOVA fits with non-Gaussian responses at arbitrary
    points, with results on the response scale.
}
\usage{
\method{predict9}{gssanova}(object, newdata, ci=FALSE, level=.95, nu=NULL, ...)
}
\arguments{
    \item{object}{Object of class inheriting from \code{"gssanova"}.}
    \item{newdata}{Data frame or model frame in which to predict.}
    \item{ci}{Flag indicating if Bayesian confidence intervals are required.
              Ignored for \code{family="polr"}.} 
    \item{level}{Confidence level.  Ignored when \code{ci=FALSE}.}
    \item{nu}{Sizes for \code{"nbinomial"} fits with known sizes.  Ignored
              otherwise.}
    \item{...}{Ignored.}
}
\value{
    For \code{ci=FALSE}, \code{predict9.gssanova} returns a vector of the
    evaluated fit,

    For \code{ci=TRUE}, \code{predict9.gssanova} returns a list of three
    elements.
    \item{fit}{Vector of evaluated fit on response scale.}
    \item{lcl}{Vector of lower confidence limit on response scale.}
    \item{ucl}{Vector of upper confidence limit on response scale.}

    For \code{family="polr"}, \code{predict9.gssanova} returns a matrix of
    probabilities with each row adding up to 1.
}
\note{
    For mixed-effect models through \code{\link{gssanova}} or
    \code{\link{gssanova1}}, the Z matrix is set to 0 if not supplied.
    To supply the Z matrix, add an element \code{random=I(...)} in
    \code{newdata}, where the as-is function \code{I(...)} preserves the
    integrity of the Z matrix in data frame.

    Unlike on the link scale, partial sums make no sense on the response
    scale, so all terms are forced in here.
}
\seealso{
    Fitting functions \code{\link{gssanova}}, \code{\link{gssanova1}} and
    methods \code{\link{predict.ssanova}}, \code{\link{summary.gssanova}},
    \code{\link{project.gssanova}}, \code{\link{fitted.gssanova}}.
}
\references{
    Gu, C. (2013), \emph{Smoothing Spline ANOVA Models (2nd Ed)}.  New
    York: Springer-Verlag.
}
\keyword{models}
\keyword{regression}
\keyword{smooth}

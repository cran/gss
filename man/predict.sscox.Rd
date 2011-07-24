\name{predict.sscox}
\alias{predict.sscox}
\title{Evaluating Smoothing Spline ANOVA Estimate of Relative Risk}
\description{
    Evaluate terms in a smoothing spline ANOVA estimate of relative risk
    at arbitrary points.  Standard errors of the terms can be requested
    for use in constructing Bayesian confidence intervals.
}
\usage{
\method{predict}{sscox}(object, newdata, se.fit=FALSE,
                        include=c(object$terms$labels,object$lab.p), ...)
}
\arguments{
    \item{object}{Object of class \code{"sscox"}.}
    \item{newdata}{Data frame or model frame in which to predict.}
    \item{se.fit}{Flag indicating if standard errors are required.}
    \item{include}{List of model terms to be included in the
	prediction.}
    \item{...}{Ignored.}
}
\value{
    For \code{se.fit=FALSE}, \code{predict.sscox} returns a vector of
    the evaluated relative risk.

    For \code{se.fit=TRUE}, \code{predict.sscox} returns a list
    consisting of the following components.
    \item{fit}{Vector of evaluated relative risk.}
    \item{se.fit}{Vector of standard errors for log relative risk.}
}
\note{
    For mixed-effect models through \code{\link{sscox}}, the Z matrix is
    set to 0 if not supplied.  To supply the Z matrix, add a component
    \code{random=I(...)} in \code{newdata}, where the as-is function
    \code{I(...)} preserves the integrity of the Z matrix in data
    frame.
}
\seealso{
    Fitting functions \code{\link{sscox}} and method
    \code{\link{project.sscox}}.
}
\author{Chong Gu, \email{chong@stat.purdue.edu}}
\keyword{models}
\keyword{regression}
\keyword{smooth}
\keyword{survival}

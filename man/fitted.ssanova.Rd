\name{fitted.ssanova}
\alias{fitted.ssanova}
\alias{residuals.ssanova}
\alias{fitted.gssanova}
\alias{residuals.gssanova}
\title{Fitted Values and Residuals from Smoothing Spline ANOVA Fits}
\description{
    Methods for extracting fitted values and residuals from smoothing
    spline ANOVA fits.
}
\usage{
fitted.ssanova(object, ...)
residuals.ssanova(object, ...)

fitted.gssanova(object, ...)
residuals.gssanova(object, type="working", ...)
}
\arguments{
    \item{object}{Object of class \code{"ssanova"} or \code{"gssanova"}.}
    \item{type}{Type of residuals desired, with two alternatives
      \code{"working"} (default) or \code{"deviance"}.}
    \item{...}{Ignored.}
}
\details{
    The fitted values for \code{"gssanova"} objects are on the link
    scale, so are the \code{"working"} residuals.
}
\keyword{models}
\keyword{regression}
\keyword{smooth}

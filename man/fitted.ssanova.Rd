\name{fitted.ssanova}
\alias{residuals.ssanova}
\alias{fitted.gssanova}
\alias{residuals.gssanova}
\title{Fitted Values and Residuals from Smoothing Spline ANOVA Fits}
\description{
    Methods for extracting fitted values and residuals from smoothing
    spline ANOVA fits.
}
\usage{
fitted[.ssanova](obj)
residuals[.ssanova](obj)

fitted[.gssanova](obj)
residuals[.gssanova](obj,type="working")
}
\arguments{
    \item{obj}{Object of class \code{"ssanova"} or \code{"gssanova"}.}
    \item{type}{Type of residuals desired, with two alternatives
	\code{"working"} (default) or \code{"deviance"}.}
}
\details{
    The fitted values for \code{"gssanova"} objects are on the link
    scale, the same scale the \code{"working"} residuals are on.
}
\keyword{models}
\keyword{regression}
\keyword{smooth}

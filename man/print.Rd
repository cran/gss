\name{print}
\alias{print.ssanova}
\alias{print.ssanova0}
\alias{print.gssanova}
\alias{print.ssden}
\alias{print.sscden}
\alias{print.sshzd}
\alias{print.sscox}
\alias{print.ssllrm}
\alias{print.summary.ssanova}
\alias{print.summary.gssanova}
\alias{print.summary.gssanova0}
\title{Print Functions for Smoothing Spline ANOVA Models}
\description{
    Print functions for Smoothing Spline ANOVA models.
}
\usage{
\method{print}{ssanova}(x, ...)
\method{print}{ssanova0}(x, ...)
\method{print}{gssanova}(x, ...)
\method{print}{ssden}(x, ...)
\method{print}{sscden}(x, ...)
\method{print}{sshzd}(x, ...)
\method{print}{sscox}(x, ...)
\method{print}{ssllrm}(x, ...)
\method{print}{summary.ssanova}(x, digits=6, ...)
\method{print}{summary.gssanova}(x, digits=6, ...)
\method{print}{summary.gssanova0}(x, digits=6, ...)
}
\arguments{
    \item{x}{Object of class \code{ssanova}, \code{summary.ssanova},
        \code{summary.gssanova}, or \code{ssden}.}
    \item{digits}{Number of significant digits to be printed in
	values.}
    \item{...}{Ignored.}
}
\seealso{
    \code{\link{ssanova}}, \code{\link{ssanova0}},
    \code{\link{gssanova}}, \code{\link{gssanova0}},
    \code{\link{ssden}}, \code{\link{ssllrm}}, \code{\link{sshzd}},
    \code{\link{summary.ssanova}}, \code{\link{summary.gssanova}},
    \code{\link{summary.gssanova0}}.
}
\keyword{smooth}
\keyword{models}

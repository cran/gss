\name{family.gssanova}
\alias{mkdata.binomial}
\alias{dev.resid.binomial}
\alias{dev.null.binomial}
\alias{mkdata.poisson}
\alias{dev.resid.poisson}
\alias{dev.null.poisson}
\alias{mkdata.Gamma}
\alias{dev.resid.Gamma}
\alias{dev.null.Gamma}
\alias{mkdata.inverse.gaussian}
\alias{dev.resid.inverse.gaussian}
\alias{dev.null.inverse.gaussian}
\alias{mkdata.nbinomial}
\alias{dev.resid.nbinomial}
\alias{dev.null.nbinomial}
\alias{mkdata.weibull}
\alias{dev.resid.weibull}
\alias{dev.null.weibull}
\alias{mkdata.lognorm}
\alias{dev.resid.lognorm}
\alias{dev.null.lognorm}
\alias{mkdata.loglogis}
\alias{dev.resid.loglogis}
\alias{dev.null.loglogis}
\title{Utility Functions for Error Families}
\description{
    Utility functions for fitting Smoothing Spline ANOVA models with
    non Gaussian responses.
}
\usage{
mkdata.[family](y, eta, wt, offset)
dev.resid.[family](y, eta, wt)
dev.null.[family](y, wt, offset)
}
\details{
    \code{family} appearing in the function names is to be substituted
    by \code{binomial}, \code{poisson}, \code{Gamma},
    \code{inverse.gaussian}, \code{nbinomial}, \code{weibull},
    \code{lognorm}, or \code{loglogis}.
    
    Only one link is implemented for each \code{family}.  It is the
    logit link for \code{binomial}, and the log link for
    \code{poisson}, \code{Gamma}, and \code{inverse.gaussian}.
    For \code{nbinomial}, the working parameter is the logit of the
    probability \eqn{p}; see \code{\link{NegBinomial}}.  For
    \code{weibull}, \code{lognorm}, and \code{loglogis}, it is the
    location parameter for the log life time.

    \code{mkdata.x} create the pseudo data to be used in iterated
    penalized least squares fitting.  \code{dev.resid.x} calculate the
    deviance residuals.  \code{dev.null.x} calculate the deviance of the
    constant null model.

    For \code{weibull}, \code{lognorm}, and \code{loglogis}, an extra
    argument \code{alpha} appears in the argument lists.  The argument
    \code{alpha} also appears in \code{mkdata.nbinomial}.
}
\arguments{
    \item{y}{Model response.}
    \item{eta}{Fitted values on link scale.}
    \item{wt}{Model weights.}
    \item{offset}{Model offset.}
}
\note{
    For \code{binomial}, the response can be specified either as two
    columns of counts or as a column of sample proportion plus a column
    of weights, as in \code{\link{glm}}.

    For \code{nbinomial}, the response may be specified as two columns
    with the second being the known sizes, or simply as a single column
    with the common unknown size to be estimated through the maximum
    likelihood method.
}
\seealso{
    \code{\link{gssanova}}.
}
\keyword{smooth}

\name{predict.ssllrm}
\alias{predict.ssllrm}
\title{Evaluating Log-Linear Regression Model Fits}
\description{
    Evaluate conditional density in a log-linear regression model fit at
    arbitrary x, or contrast of log conditional density possibly with
    standard errors for constructing Bayesian confidence intervals.
}
\usage{
\method{predict}{ssllrm}(object, x, y=object$qd.pt, odds=NULL, se.odds=FALSE, ...)
}
\arguments{
    \item{object}{Object of class \code{"ssllrm"}.}
    \item{x}{Data frame of x values.}
    \item{y}{Data frame of y values; y-variables must be factors.}
    \item{odds}{Optional coefficients of contrast.}
    \item{se.odds}{Flag indicating if standard errors are required.
      Ignored when \code{odds=NULL}.}
    \item{...}{Ignored.}
}
\value{
    For \code{odds=NULL}, \code{predict.ssanova} returns a vector/matrix
    of the estimated \code{f(y|x)}.

    When \code{odds} is given, it should match \code{y} in length and
    the coefficients must add to zero; \code{predict.ssanova} then
    returns a vector of estimated "odds ratios" if \code{se.odds=FALSE}
    or a list consisting of the following components if
    \code{se.odds=TRUE}.
    \item{fit}{Vector of evaluated fit.}
    \item{se.fit}{Vector of standard errors.}
}
\seealso{
    Fitting function \code{\link{ssllrm}}.
}
\author{Chong Gu, \email{chong@stat.purdue.edu}}
\keyword{models}
\keyword{regression}
\keyword{smooth}

\name{mkrk.nominal}
\alias{mkrk.nominal}
\alias{mkrk.ordinal}
\title{
    Crafting Building Blocks for Discrete Splines
}
\description{
    Craft numerical functions to be used by \code{mkterm} to assemble
    model terms involving factors.
}
\usage{
mkrk.nominal(levels)
mkrk.ordinal(levels)
}
\arguments{
    \item{levels}{Levels of the factor.}
}
\details{
    For a nominal factor with levels \eqn{1,2,\dots,k}, the level means
    \eqn{f(i)} will be shrunk towards each other through a penalty
    proportional to
    \deqn{(f(1)-f(.))^2+\dots+(f(k)-f(.))^2}
    where \eqn{f(.)=(f(1)+\dots+f(k))/k}.

    For a ordinal factor with levels \eqn{1<2<\dots<k}, the level means
    \eqn{f(i)} will be shrunk towards each other through a penalty
    proportional to
    \deqn{(f(1)-f(2))^2+\dots+(f(k-1)-f(k))^2}
}
\value{
    A list of two elements.
    \item{fun}{Function definition.}
    \item{env}{Portable local constants derived from the arguments.}
}
\note{
    \code{mkrk.x} create a bivariate function
    \code{fun(x,y,env,outer=FALSE)}, where \code{x}, \code{y} are real
    arguments and local constants can be passed in through \code{env}.
}
\seealso{
    \code{\link{mkterm}}, \code{\link{mkrk.cubic}}, and
    \code{\link{mkrk.tp}}.
}
\keyword{internal}

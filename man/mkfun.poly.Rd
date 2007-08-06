\name{mkfun.poly}
\alias{mkfun.poly}
\alias{mkrk.cubic}
\alias{mkphi.cubic}
\alias{mkrk.cubic.per}
\alias{mkrk.linear}
\alias{mkrk.linear.per}
\title{
    Crafting Building Blocks for Polynomial Splines
}
\description{
    Craft numerical functions to be used by \code{\link{mkterm}} to
    assemble model terms.
}
\usage{
mkrk.cubic(range)
mkphi.cubic(range)
mkrk.cubic.per(range)
mkrk.linear(range)
mkrk.linear.per(range)
}
\arguments{
    \item{range}{Numerical vector whose minimum and maximum specify the
	range on which the function to be crafted is defined.}
}
\value{
    A list of two components.
    \item{fun}{Function definition.}
    \item{env}{Portable local constants derived from the argument.}
}
\note{
    \code{mkrk.x} create a bivariate function
    \code{fun(x,y,env,outer=FALSE)}, where \code{x}, \code{y} are real
    arguments and local constants can be passed in through \code{env}.

    \code{mkphi.cubic} creates a univariate function
    \code{fun(x,nu,env)}.
}
\seealso{
    \code{\link{mkterm}}, \code{\link{mkfun.tp}}, and
    \code{\link{mkrk.nominal}}.
}
\details{
    \code{mkrk.cubic}, \code{mkphi.cubic}, and \code{mkrk.linear}
    implement the polynomial spline construction in Gu (2002,
    Sec. 2.3.3) for \eqn{m=2,1}.

    
    \code{mkrk.cubic.per} and \code{mkrk.linear.per} implement the
    periodic polynomial spline construction in Gu (2002, Sec. 4.2.1) for
    \eqn{m=2,1}.
}
\author{Chong Gu, \email{chong@stat.purdue.edu}}
\references{
    Gu, C. (2002), \emph{Smoothing Spline ANOVA Models}.  New York:
    Springer-Verlag.
}
\keyword{internal}

\name{mkfun.poly}
\alias{mkrk.linear}
\alias{mkrk.cubic}
\alias{mkphi.cubic}
\title{
    Crafting Building Blocks for Polynomial Splines
}
\description{
    Craft numerical functions to be used by \code{\link{mkterm.linear}},
    \code{\link{mkterm.cubic}} to assemble model terms.
}
\usage{
mkrk.linear(range)
mkrk.cubic(range)
mkphi.cubic(range)
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
    \code{fun(x,nu=1,env)}.
}
\seealso{
    \code{\link{mkfun.tp}}, \code{\link{mkrk.factor}}.
}
\keyword{smooth}

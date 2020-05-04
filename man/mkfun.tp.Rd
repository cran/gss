\name{mkfun.tp}
\alias{mkfun.tp}
\alias{mkrk.tp}
\alias{mkrk.tp.p}
\alias{mkphi.tp}
\alias{mkphi.tp.p}
\alias{mkrk.sphere}
\title{
    Crafting Building Blocks for Thin-Plate and Spherical Splines
}
\description{
    Craft numerical functions to be used by \code{\link{mkterm}} to
    assemble model terms.
}
\usage{
mkrk.tp(dm, order, mesh, weight)
mkphi.tp(dm, order, mesh, weight)
mkrk.tp.p(dm, order)
mkphi.tp.p(dm, order)

mkrk.sphere(order)
}
\arguments{
    \item{dm}{Dimension of the variable \eqn{d}.}
    \item{order}{Order of the differential operator \eqn{m}.}
    \item{mesh}{Normalizing mesh.}
    \item{weight}{Normalizing weights.}
}
\details{
    \code{mkrk.tp}, \code{mkphi.tp}, \code{mkrk.tp.p}, and
    \code{mkphi.tp.p} implement the construction in Gu (2002,
    Sec. 4.4).  Thin-plate splines are defined for \eqn{2m>d}.

    \code{mkrk.tp.p} generates the pseudo kernel, and \code{mkphi.tp.p}
    generates the \eqn{(m+d-1)!/d!/(m-1)!} lower order polynomials with
    total order less than \eqn{m}.

    \code{mkphi.tp} generates normalized lower order polynomials
    orthonormal w.r.t. a norm specified by \code{mesh} and
    \code{weight}, and \code{mkrk.tp} conditions the pseudo kernel to
    generate the reproducing kernel orthogonal to the lower order
    polynomials w.r.t. the norm.

    \code{mkrk.sphere} implements the reproducing kernel construction of
    Wahba (1981) for \eqn{m=2,3,4}.
}
\value{
    A list of two elements.
    \item{fun}{Function definition.}
    \item{env}{Portable local constants derived from the arguments.}
}
\note{
    \code{mkrk.tp} and \code{mkrk.sphere} create a bivariate function
    \code{fun(x,y,env,outer=FALSE)}, where \code{x}, \code{y} are real
    arguments and local constants can be passed in through \code{env}.

    \code{mkphi.tp} creates a collection of univariate functions
    \code{fun(x,nu,env)}, where \code{x} is the argument and \code{nu}
    is the index.
}
\references{
    Gu, C. (2013), \emph{Smoothing Spline ANOVA Models (2nd Ed)}.  New
    York: Springer-Verlag.

    Wahba, G. (1981), Spline interpolation and smoothing on the sphere.
    \emph{SIAM Journal on Scientific and Statistical Computing},
    \bold{2}, 5--16.
}
\seealso{
    \code{\link{mkterm}}, \code{\link{mkfun.poly}}, and
    \code{\link{mkrk.nominal}}.
}
\keyword{internal}

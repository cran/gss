\name{mkfun.tp}
\alias{mkfun.tp}
\alias{mkrk.tp}
\alias{mkrk.tp.p}
\alias{mkphi.tp}
\alias{mkphi.tp.p}
\title{
    Crafting Building Blocks for Thin-Plate Splines
}
\description{
    Craft numerical functions to be used by \code{\link{mkterm.tp}} to
    assemble model terms.
}
\usage{
mkrk.tp(dm,order,mesh,weight)
mkphi.tp(dm,order,mesh,weight)

mkrk.tp.p(dm,order)
mkphi.tp.p(dm,order)
}
\arguments{
    \item{dm}{Dimension of the variable \eqn{d}.}
    \item{order}{Order of the differential operator \eqn{m}.}
    \item{mesh}{Normalizing mesh.}
    \item{weight}{Normalizing weights.}
}
\details{
    These are not to be called by the user.

    Thin-plate splines are defined for \eqn{2m>d}.

    \code{mkrk.tp.p} generates the pseudo kernel, and \code{mkphi.tp.p}
    generates the \eqn{(m+d-1)!/d!/(m-1)!} lower order polynomials with
    total order less than \eqn{m}.

    \code{mkphi.tp} generates normalized lower order polynomials
    orthonormal w.r.t. a norm specified by \code{mesh} and \code{weight}
    as described in the reference, and \code{mkrk.tp} conditions the
    pseudo kernel to generate the reproducing kernel orthogonal to the
    lower order polynomials w.r.t. the norm.
}
\value{
    A list of two components.
    \item{fun}{Function definition.}
    \item{env}{Portable local constants derived from the arguments.}
}
\note{
    \code{mkrk.tp} creates a bivariate function
    \code{fun(x,y,env,outer=FALSE)}, where \code{x}, \code{y} are real
    arguments and local constants can be passed in through \code{env}.

    \code{mkphi.tp} creates a collection of univariate functions
    \code{fun(x,nu,env)}, where \code{x} is the argument and \code{nu}
    is the index.
}
\references{
    Gu, C. and Wahba, G. (1993), Semiparametric analysis of variance
    with tensor product thin plate splines.  \emph{Journal of the Royal
	Statistical Society Ser. B}, \bold{55}, 353--368.
}
\seealso{
    \code{\link{mkfun.poly}}, \code{\link{mkrk.factor}}.
}
\keyword{internal}

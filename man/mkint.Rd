\name{mkint}
\alias{mkint}
\alias{mkint2}
\title{
    Generating Integrals of Basis Terms
}
\description{
    Generate integrals of basis terms for use in the ssden1 suite.
}
\usage{
mkint(mf, type, id.basis, quad, term, rho, rho.int)
mkint2(mf, type, id.basis, quad, term)
}
\arguments{
    \item{mf}{Model frame of the model formula.}
    \item{type}{List specifying the type of spline for each variable.}
    \item{id.basis}{Index of observations to be used as "knots."}
    \item{quad}{Quadratures on marginal domains, weighted by rho.}
    \item{term}{Model terms generated by \code{mkterm}.}
    \item{rho}{Marginal log(rho) on quadrature points.}
    \item{rho.int}{Marginal integrals of log(rho).}
}
\details{
    \code{mkint} calculates the first moments of basis functions with
    respect to the indepent joint density rho; \code{mkint2} calculates
    the second moments for use in \code{project.ssden1}.
}
\value{
    \code{mkint} returns a list of five elements.
    \item{s}{First moments of phi's.}
    \item{r}{First moments of rk's.}
    \item{s.rho}{Cross moments of phi*log(rho).}
    \item{r.rho}{Cross moments of rk*log(rho).}
    \item{var.type}{Types for variables.}

    \code{mkint2} returns a list of three elements.
    \item{ss}{Second moments of phi's.}
    \item{sr}{Cross moments of phi's and rk's.}
    \item{rr}{Second moments of rk.}
}
\keyword{internal}

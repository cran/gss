\name{smolyak}
\alias{smolyak.quad}
\alias{ccsmolyak.quad}
\alias{smolyak.size}
\alias{ccsmolyak.size}
\title{Generating Smolyak Cubature}
\description{
    Generate Smolyak cubatures using C routines modified from
    \code{smolyak.c} and \code{ccsmolyak.c} found in Knut Petras'
    SMOLPACK.
}
\usage{
smolyak.quad(d,k)
ccsmolyak.quad(d,k)

smolyak.size(d,k)
ccsmolyak.size(d,k)
}
\arguments{
    \item{d}{Dimension of unit cube.}
    \item{k}{Depth of algorithm.}
}
\value{
    \code{smolyak.quad} and \code{ccsmolyak.quad} return a list object
    with the following components.
    \item{pt}{Quadrature nodes in rows of matrix.}
    \item{wt}{Quadrature weights.}

    \code{smolyak.size} and \code{ccsmolyak.size} return a integer.
}
\details{
    \code{smolyak.quad} and \code{smolyak.size} concern the delayed
    Smolyak cubature.  \code{ccsmolyak.quad} and \code{ccsmolyak.size}
    concern the Clenshaw-Curtis Smolyak cubature.
}
\keyword{math}

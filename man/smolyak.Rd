\name{smolyak}
\alias{smolyak.quad}
\alias{smolyak.size}
\title{Generating Smolyak Cubature}
\description{
    Generate delayed Smolyak cubatures using C routines modified from
    \code{smolyak.c} found in Knut Petras' SMOLPACK.
}
\usage{
smolyak.quad(d, k)

smolyak.size(d, k)
}
\arguments{
    \item{d}{Dimension of unit cube.}
    \item{k}{Depth of algorithm.}
}
\value{
    \code{smolyak.quad} returns a list object with the following
    components.
    \item{pt}{Quadrature nodes in rows of matrix.}
    \item{wt}{Quadrature weights.}

    \code{smolyak.size} returns an integer.
}
\keyword{math}

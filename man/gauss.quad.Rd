\name{gauss.quad}
\alias{gauss.quad}
\title{Generating Gauss-Legendre Quadrature}
\description{
    Generate Gauss-Legendre quadratures using the FORTRAN routine
    \code{gaussq.f} found on NETLIB.
}
\usage{
gauss.quad(size, interval)
}
\arguments{
    \item{size}{Size of quadrature.}
    \item{interval}{Interval to be covered.}
}
\value{
    \code{gauss.quad} returns a list object with the following components.
    \item{pt}{Quadrature nodes.}
    \item{wt}{Quadrature weights.}
}
\keyword{math}

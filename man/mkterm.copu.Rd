\name{mkterm.copu}
\alias{mkterm.copu}
\title{
    Assembling Model Terms for Copula Density Estimation
}
\description{
    Assemble numerical functions for calculating model terms in copula
    density estimation.
}
\usage{
mkterm.copu(dm, order, symmetry, exclude)
}
\arguments{
    \item{dim}{Dimension of the domain.}
    \item{order}{Highest order of interactions allowed in log density.}
    \item{symmetry}{Flag indicating whether to enforce symmetry, or
        invariance under coordinate permutation.}
    \item{exclude}{Pair(s) of marginals whose interactions to be
        excluded in log density.}
}
\author{Chong Gu, \email{chong@stat.purdue.edu}}
\keyword{internal}

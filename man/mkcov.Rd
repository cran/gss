\name{mkcov}
\alias{mkcov}
\alias{mkcov.arma}
\alias{mkcov.long}
\alias{mkcov.known}
\title{
    Generating Covariance for Correlated Data
}
\description{
    Generate entries of covariance functions for correlated data.
}
\usage{
mkcov.arma(p, q, n)
mkcov.long(id)
mkcov.known(w)
}
\arguments{
    \item{p}{Order of AR terms.}
    \item{q}{Order of MA terms.}
    \item{n}{Dimension of covariance matrix.}

    \item{id}{Factor of subject ID.}

    \item{w}{Covariance matrix; only the upper triangular part is used.}
}
\details{
    \code{mkcov.arma} generates covariance functions for ARMA(p,q)
    model.

    \code{mkcov.long} generates covariance functions for longitudinal
    data.

    \code{mkcov.known} allows one to use a known covariance matrix in
    \code{ssanova9}.
}
\value{
    A list of three components.
    \item{fun}{Covariance matrix to be evaluated through
    \code{fun(gamma,env)} or \code{fun(env)}.}
    \item{env}{Constants in covariance function.}
    \item{init}{Initial values for correlation parameters.}
}
\note{
    One may pass \code{list(fun=...,env=...,init=...)} directly to the
    argument \code{cov} in calls to \code{\link{ssanova9}}, or make use
    of the \code{mkcov.x} functions through
    \code{cov=list("arma",c(p,q))}, \code{cov=list("long",id)}, or
    \code{cov=list("known",w)}.
}
\author{Chong Gu, \email{chong@stat.purdue.edu}}
\keyword{internal}

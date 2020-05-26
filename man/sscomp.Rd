\name{sscomp}
\alias{sscomp}
\alias{sscomp2}
\title{Composition Estimation}
\description{
    Estimate composition using multinomial counts.
}
\usage{
sscomp(x,wt=rep(1,length(x)),alpha=1.4)

sscomp2(x,alpha=1.4)
}
\arguments{
    \item{x}{Numerical vector or matrix of multinomial counts.}
    \item{wt}{Numerical vector of integration weights.}
    \item{alpha}{Parameter defining cross-validation score for smoothing
      parameter selection.}
}
\details{
    \code{sscomp} takes a vector \code{x} to estimate composition using
    density estimation on a nominal discrete domain; zero counts must be
    included in \code{x} to specify the domain.  \code{wt} mimicking the
    shape of the unknown density could improve performance.

    \code{sscomp2} takes a matrix \code{x}, collapses columns to
    estimate a density using \code{sscomp}, then using that as \code{wt}
    in further \code{sscomp} calls to estimate composition for each
    column.
}
\references{
    Gu, C. (2020), Composition estimation via shrinkage.
    \emph{manuscript}.
}

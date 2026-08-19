\name{retroCI}
\alias{retroCI}
\title{Calculating Retrospective Confidence Intervals using Clones of
  \code{ssanova} and \code{gssanova} Objects}
\description{
    Evaluate clones of \code{ssanova} or \code{gssanova} objects at
    arbitrary points, then obtain quantiles of the clone evaluations for
    retrospective confidence intervals.
}
\usage{
retroCI(object, newdata, limits=c(0.025,0.975),
 include=c(object$fit$terms$labels,object$fit$lab.p))
}
\arguments{
    \item{object}{Object of class \code{"clone"}.}
    \item{newdata}{Data frame or model frame in which to predict.}
    \item{limits}{To be passed to \code{\link{quantile}} for quantiles
      of clone evaluations.}
    \item{include}{List of model terms to be included in the
	prediction.  The \code{offset} term, if present, is to be
	specified by \code{"offset"}.}
}
\value{
    A matrix of \code{length(limits)} columns of clone evaluation
    quantiles on the \code{newdata} points.
}
\note{
    The function is essentially a batch version of
    \code{\link{predict.ssanova}} followed by a call to
    \code{\link{quantile}}, inheriting much of the syntax and the
    internal code from \code{predict.ssanova}.
}
\references{
    Gu, C. (2026), Retrospective statistical inference.  Manuscript.
}
\keyword{models}
\keyword{regression}
\keyword{smooth}

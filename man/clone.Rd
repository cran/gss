\name{clone}
\alias{clone}
\alias{clone.ssanova}
%\alias{project.ssanova9}
\alias{clone.gssanova}
%\alias{project.ssden}
%\alias{project.ssden1}
%\alias{project.sscden}
%\alias{project.sscden1}
%\alias{project.sshzd}
%\alias{project.sscox}
%\alias{project.sshzd1}
%\alias{project.ssllrm}
\title{Cloning Smoothing Spline ANOVA Fits for Retrospective Inference}
\description{
    Simulate replicates from smoothing spline ANOVA fits then calculate
    clones of the fits.
}
\usage{
clone(object, ...)
\method{clone}{ssanova}(object, nrep=1000, type=1, seed=NULL, eta=NULL,
fix=FALSE, ...)
%\method{project}{ssanova9}(object, include, ...)
\method{clone}{gssanova}(object, nrep=1000, type=1, seed=NULL, eta=NULL,
fix=FALSE, ...)
%\method{project}{ssden}(object, include, mesh=FALSE, ...)
%\method{project}{ssden1}(object, include, drop1=FALSE, ...)
%\method{project}{sscden}(object, include, ...)
%\method{project}{sscden1}(object, include, ...)
%\method{project}{sshzd}(object, include, mesh=FALSE, ...)
%\method{project}{sscox}(object, include, ...)
%\method{project}{sshzd1}(object, include, ...)
%\method{project}{ssllrm}(object, include, ...)
}
\arguments{
    \item{object}{Object of class \code{"ssanova"} or \code{"gssanova"}.}
%        \code{"ssden"}, \code{"ssden1"}, \code{"sscden"},
%	\code{"sscden1"}, \code{"sshzd"}, \code{"sshzd1"}, or
%	\code{"ssllrm"}.}
    \item{...}{Additional arguments.  Ignored in \code{clone.x}.}
    \item{nrep}{Number of clones to be generated.}
    \item{type}{type=1 tries to enforce effective constraints, type=2
      searches for optimal performer, both via smoothing parameter
      selection.}
    \item{seed}{Seed to be used in random number generation; if NULL
      seed will be set internally to ensure reproducible results.}
    \item{eta}{True function values to generate replicates from; if NULL
      the fitted value of \code{object} will be used.}
    \item{fix}{Logical flag indicating whether to fix theta's with
      multiple smoothing parameters.}
%    \item{mesh}{Flag indicating whether to return evaluations of the
%        projection.}
%    \item{drop1}{If TRUE, calculate \code{p<-length(include)} projections
%	with \code{include[-i]}, \code{i=1,...,p}.}
}
\value{
    The function returns a list object of class \code{clone}, to be use
    as input to \code{\link{retroCI}} to product retrospective
    confidence intervals.

    The input object is in element \code{fit}
}
\references{
    Gu, C. (2026), Retrospective statistical inference.  Manuscript.
}
\seealso{
    Fitting functions \code{\link{ssanova}}, \code{\link{gssanova}}.
%    \code{\link{ssden}}, \code{\link{sshzd}}, and \code{\link{sshzd1}}.
}
\keyword{models}
\keyword{smooth}
%\keyword{htest}

\name{project1}
\alias{project1}
\alias{project1.ssanova}
%\alias{project.ssanova9}
\alias{project1.gssanova}
%\alias{project.ssden}
%\alias{project.ssden1}
%\alias{project.sscden}
%\alias{project.sscden1}
%\alias{project.sshzd}
%\alias{project.sscox}
%\alias{project.sshzd1}
%\alias{project.ssllrm}
\title{Projecting Smoothing Spline ANOVA Fits for Retrospective Testing}
\description{
    Calculate Kullback-Leibler projection of smoothing spline ANOVA fits
    and calibrate using clone distribution under the null.
}
\usage{
project1(object, ...)
\method{project1}{ssanova}(object, include, nrep=200, type=1, ...)
%\method{project}{ssanova9}(object, include, ...)
\method{project1}{gssanova}(object, include, nrep=200, type=1, ...)
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
    \item{...}{Additional arguments.  Ignored in \code{project1.x}.}
    \item{include}{List of model terms to be included in the reduced
        model space.  The \code{partial} and \code{offset} terms, if
	present, are to be specified by \code{"partial"} and
	\code{"offset"}, respectively.}
    \item{nrep}{Number of clones to be simulated under the null .}
    \item{type}{Clone type to use; see \code{\link{clone}}.}
%    \item{mesh}{Flag indicating whether to return evaluations of the
%        projection.}
%    \item{drop1}{If TRUE, calculate \code{p<-length(include)} projections
%	with \code{include[-i]}, \code{i=1,...,p}.}
}
\value{
    The functions return a list consisting of the following elements.
    \item{ratio}{KL(fit0,fit1)/KL(fit0,null); the smaller the value, the
        more feasible the reduced model is.}
    \item{kl}{KL(fit0,fit1).}
%    For regression fits, the list also contains the following element.
    \item{check}{KL(fit0,fit1)/KL(fit0,null)+KL(fit1,null)/KL(fit0,null);
        a value closer to 1 is preferred.}
    \item{pvalue}{Retrospective p-value of ratio, calibrated against
        clone distribution under the null.}

%    For density and hazard fits, the list may contain the following
%    optional element.
%    \item{mesh}{The evaluations of the projection.}
}
%\note{
%    \code{project.ssden1}, \code{project.sscden1}, and
%    \code{project.sshzd1} calculates square error projections.
%}
\details{
    The entropy KL(fit0,null) can be decomposed as the sum of
    KL(fit0,fit1) and KL(fit1,null), where fit0 is the fit to be
    projected, fit1 is the projection in the reduced model space, and
    null is the constant fit.  The ratio KL(fit0,fit1)/KL(fit0,null)
    serves as a diagnostic of the feasibility of the reduced model.

    The projection itself is exactly as in \code{\link{project}}, but
    with a retrospective p-value attached, calibrating against the
    clone distribution under the null using \code{nrep} simulated
    clones.
}
\references{
    Gu, C. (2004), Model diagnostics for smoothing spline ANOVA models.
    \emph{The Canadian Journal of Statistics}, \bold{32}, 347--358.

    Gu, C. (2026), Retrospective statistical inference.  Manuscript.
}
\seealso{
    Fitting functions \code{\link{ssanova}}, \code{\link{gssanova}}, and
    cloning function \code{\link{clone}}.
%    \code{\link{ssden}}, \code{\link{sshzd}}, and \code{\link{sshzd1}}.
}
\keyword{models}
\keyword{smooth}
\keyword{htest}

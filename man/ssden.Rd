\name{ssden}
\alias{ssden}
\title{Estimating Probability Density Using Smoothing Splines}
\description{
    Estimate probability densities using smoothing spline ANOVA models
    with cubic spline, linear spline, or thin-plate spline marginals for
    numerical variables.  The symbolic model specification via
    \code{formula} follows the same rules as in \code{\link{lm}}, but
    with the response missing.
}
\usage{
ssden(formula, type="cubic", data=list(), alpha=1.4, weights=NULL,
      subset, na.action=na.omit, id.basis=NULL, nbasis=NULL, seed=NULL,
      domain=NULL, quadrature=NULL, ext=.05, order=2, prec=1e-7,
      maxiter=30)
}
\arguments{
    \item{formula}{Symbolic description of the model to be fit.}
    \item{type}{Type of numerical marginals to be used.  Supported are
	\code{type="cubic"} for cubic spline marginals,
	\code{type="linear"} for linear spline marginals, and
	\code{type="tp"} for thin-plate spline marginals.}
    \item{data}{Optional data frame containing the variables in the
        model.}
    \item{alpha}{Parameter defining cross-validation score for smoothing
        parameter selection.}
    \item{weights}{Optional vector of bin-counts for histogram data.}
    \item{subset}{Optional vector specifying a subset of observations
	to be used in the fitting process.}
    \item{na.action}{Function which indicates what should happen when
        the data contain NAs.}
    \item{id.basis}{Index of observations to be used as "knots."}
    \item{nbasis}{Number of "knots" to be used.  Ignored when
        \code{id.basis} is specified.}
    \item{seed}{Seed to be used for the random generation of "knots."
        Ignored when \code{id.basis} is specified.}
    \item{domain}{Data frame specifying marginal support of density.}
    \item{quadrature}{Quadrature for calculating integral.  Mandatory
        for \code{type="tp"}.}
    \item{ext}{For cubic spline and linear spline marginals, this option
	specifies how far to extend the domain beyond the minimum and
	the maximum as a percentage of the range.  The default
	\code{ext=.05} specifies marginal domains of lengths 110 percent
	of their respective ranges.  Evaluation outside of the domain
	will result in an error.  Ignored if \code{type="tp"} or
	\code{domain} are specified.}
    \item{order}{For thin-plate spline marginals, this option specifies
	the order of the marginal penalties.  Ignored if
	\code{type="cubic"} or \code{type="linear"} are specified.}
    \item{prec}{Precision requirement for internal iterations.}
    \item{maxiter}{Maximum number of iterations allowed for
	internal iterations.}
}
\details{
    The model specification via \code{formula} is for the log density.
    For example, \code{~x1*x2} prescribes a model of the form
    \deqn{
	log f(x1,x2) = g_{1}(x1) + g_{2}(x2) + g_{12}(x1,x2) + C
    }
    with the terms denoted by \code{"x1"}, \code{"x2"}, and
    \code{"x1:x2"}; the constant is determined by the fact that a
    density integrates to one.

    The selective term elimination may characterize (conditional)
    independence structures between variables.  For example,
    \code{~x1*x2+x1*x3} yields the conditional independence of x2 and x3
    given x1.  Currently, up to four variables are supported.

    Parallel to those in a \code{\link{ssanova}} object, the model terms
    are sums of unpenalized \emph{fixed effects} and the penalized
    \emph{random effects}.  Attached to every random effect there is a
    smoothing parameter, and the model complexity is largely determined
    by the number of smoothing parameters.

    The selection of smoothing parameters is through a cross-validation
    mechanism described in the references, with a parameter
    \code{alpha}; \code{alpha=1} is "unbiased" for the minimization of
    Kullback-Leibler loss but may yield severe undersmoothing, whereas
    larger \code{alpha} yields smoother estimates.

    A subset of the observations are selected as "knots."  Unless
    specified via \code{id.basis} or \code{nbasis}, the subset size is
    determined by \eqn{max(30,10n^(2/9))}, which is appropriate for
    \code{type="cubic"} but not necessarily for \code{type="linear"} or
    \code{type="tp"}.
}
\note{
    For \code{type="cubic"} and \code{type="linear"}, the quadrature
    will be generated if not provided by the user.  The default
    quadrature in 1-D is the 200-point Gauss-Legendre formula on the
    domain.  The default quadratures on 2-D, 3-D, and 4-D cubes are
    selected delayed Smolyak cubatures with 449, 2527, and 13697 points,
    on properly scaled product domains.  See \code{\link{gauss.quad}}
    and \code{\link{smolyak.quad}}.
}
\value{
    \code{ssden} returns a list object of \code{\link{class} "ssden"}.

    \code{\link{dssden}} and \code{\link{cdssden}} can be used to
    evaluate the estimated joint density and conditional density.
    \code{pssden}, \code{qssden}, \code{cpssden}, and \code{cqssden} can
    be used to evaluate (conditional) cdf and quantiles.
}
\seealso{
    \code{\link{dssden}} and \code{\link{cdssden}}.
}
\author{Chong Gu, \email{chong@stat.purdue.edu}}
\references{
    Gu, C. and Wang, J. (2002), \emph{Penalized Likelihood Density
    Estimation: Direct Cross-Validation and Scalable Approximation}.
    Available at \url{http://stat.purdue.edu/~chong/manu.html}.
  
    Gu, C. (2002), \emph{Smoothing Spline ANOVA Models}.  New York:
    Springer-Verlag.
}
\examples{
## 1-D estimate: Buffalo snowfall
data(buffalo)
buff.fit <- ssden(~buffalo,domain=data.frame(buffalo=c(0,150)))
plot(xx<-seq(0,150,len=101),dssden(buff.fit,xx),type="l")
plot(xx,pssden(buff.fit,xx),type="l")
plot(qq<-seq(0,1,len=51),qssden(buff.fit,qq),type="l")
## Clean up
\dontrun{rm(buffalo,buff.fit,xx,qq)
dev.off()}

## 2-D with triangular domain: AIDS incubation
data(aids)
## rectangular quadrature
quad.pt <- expand.grid(incu=((1:40)-.5)/40*100,infe=((1:40)-.5)/40*100)
quad.pt <- quad.pt[quad.pt$incu<=quad.pt$infe,]
quad.wt <- rep(1,nrow(quad.pt))
quad.wt[quad.pt$incu==quad.pt$infe] <- .5
quad.wt <- quad.wt/sum(quad.wt)*5e3
## additive model (pre-truncation independence)
aids.fit <- ssden(~incu+infe,data=aids,subset=age>=60,
                  domain=data.frame(incu=c(0,100),infe=c(0,100)),
                  quad=list(pt=quad.pt,wt=quad.wt))
## conditional (marginal) density and quantiles of infe
jk <- cdssden(aids.fit,xx<-seq(0,100,len=51),data.frame(incu=50))
plot(xx,jk$pdf,type="l")
cqssden(aids.fit,c(.05,.25,.5,.75,.95),data.frame(incu=50),jk$int)
## Clean up
\dontrun{rm(aids,quad.pt,quad.wt,aids.fit,jk,xx)
dev.off()}
}
\keyword{smooth}
\keyword{models}
\keyword{distribution}

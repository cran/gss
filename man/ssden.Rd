\name{ssden}
\alias{ssden}
\alias{ssden1}
\title{Estimating Probability Density Using Smoothing Splines}
\description{
    Estimate probability densities using smoothing spline ANOVA models.
    The symbolic model specification via \code{formula} follows the same
    rules as in \code{\link{lm}}, but with the response missing.
}
\usage{
ssden(formula, type=NULL, data=list(), alpha=1.4, weights=NULL,
      subset, na.action=na.omit, id.basis=NULL, nbasis=NULL, seed=NULL,
      domain=as.list(NULL), quad=NULL, qdsz.depth=NULL, bias=NULL,
      prec=1e-7, maxiter=30, skip.iter=FALSE)

ssden1(formula, type=NULL, data=list(), alpha=1.4, weights=NULL,
       subset, na.action=na.omit, id.basis=NULL, nbasis=NULL, seed=NULL,
       domain=as.list(NULL), quad=NULL, prec=1e-7, maxiter=30)
}
\arguments{
    \item{formula}{Symbolic description of the model to be fit.}
    \item{type}{List specifying the type of spline for each variable.
        See \code{\link{mkterm}} for details.}
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
    \item{quad}{Quadrature for calculating integral.  Mandatory if
        variables other than factors or numerical vectors are involved.}
    \item{qdsz.depth}{Depth to be used in \code{\link{smolyak.quad}} for
        the generation of quadrature.}
    \item{bias}{Input for sampling bias.}
    \item{prec}{Precision requirement for internal iterations.}
    \item{maxiter}{Maximum number of iterations allowed for
        internal iterations.}
    \item{skip.iter}{Flag indicating whether to use initial values of
        theta and skip theta iteration.  See \code{\link{ssanova}} for
	notes on skipping theta iteration.}
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
    given x1.

    Parallel to those in a \code{\link{ssanova}} object, the model terms
    are sums of unpenalized and penalized terms.  Attached to every
    penalized term there is a smoothing parameter, and the model
    complexity is largely determined by the number of smoothing
    parameters.

    The selection of smoothing parameters is through a cross-validation
    mechanism described in the references, with a parameter
    \code{alpha}; \code{alpha=1} is "unbiased" for the minimization of
    Kullback-Leibler loss but may yield severe undersmoothing, whereas
    larger \code{alpha} yields smoother estimates.

    A subset of the observations are selected as "knots."  Unless
    specified via \code{id.basis} or \code{nbasis}, the number of
    "knots" \eqn{q} is determined by \eqn{max(30,10n^{2/9})}, which is
    appropriate for the default cubic splines for numerical vectors.
}
\note{
    In \code{ssden}, default quadrature will be constructed for
    numerical vectors on a hyper cube, then outer product with factor
    levels will be taken if factors are involved.  The sides of the
    hyper cube are specified by \code{domain}; for \code{domain$x}
    missing, the default is
    \code{c(min(x),max(x))+c(-1,1)*(max(x)-mimn(x))*.05}.  In 1-D, the
    quadrature is the 200-point Gauss-Legendre formula returned from
    \code{\link{gauss.quad}}.  In multi-D, delayed Smolyak cubatures
    from \code{\link{smolyak.quad}} are used on cubes with the marginals
    properly transformed; see Gu and Wang (2003) for the marginal
    transformations.

    For reasonable execution time in higher dimensions, set
    \code{skip.iter=TRUE} in call to \code{ssden}.

    If you get an error message from \code{ssden} stating \code{"Newton
    iteration diverges"}, try to use a larger \code{qdsz.depth} which
    will execute slower, or switch to \code{ssden1}.  The default values
    of \code{qdsz.depth} for dimensions 4, 5, 6+ are 12, 11, 10.

    \code{ssden1} does not involve multi-D quadrature but does not
    perform as well as \code{ssden}.  It can be used in very high
    dimensions where \code{ssden} is infeasible.

    The results may vary from run to run.  For consistency, specify
    \code{id.basis} or set \code{seed}.
}
\value{
    \code{ssden} returns a list object of class \code{"ssden"}.
    \code{ssden1} returns a list object of class
    \code{c("ssden1","ssden")}.

    \code{\link{dssden}} and \code{\link{cdssden}} can be used to
    evaluate the estimated joint density and conditional density;
    \code{\link{pssden}}, \code{\link{qssden}}, \code{\link{cpssden}},
    and \code{\link{cqssden}} can be used to evaluate (conditional) cdf
    and quantiles.

    The method \code{\link{project.ssden}} can be used to calculate the
    Kullback-Leibler projection of \code{"ssden"} objects for model
    selection; \code{\link{project.ssden1}} can be used to calculate the
    square error projection of \code{"ssden1"} objects.
}
\author{Chong Gu, \email{chong@stat.purdue.edu}}
\references{
    Gu, C. (2013), \emph{Smoothing Spline ANOVA Models (2nd Ed)}.  New
    York: Springer-Verlag.

    Gu, C. and Wang, J. (2003), Penalized likelihood density
    estimation: Direct cross-validation and scalable approximation.
    \emph{Statistica Sinica}, \bold{13}, 811--826.
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
## conditional (marginal) density of infe
jk <- cdssden(aids.fit,xx<-seq(0,100,len=51),data.frame(incu=50))
plot(xx,jk$pdf,type="l")
## conditional (marginal) quantiles of infe (TIME-CONSUMING)
\dontrun{
cqssden(aids.fit,c(.05,.25,.5,.75,.95),data.frame(incu=50))
}
## Clean up
\dontrun{rm(aids,quad.pt,quad.wt,aids.fit,jk,xx)
dev.off()}

## One factor plus one vector
data(gastric)
gastric$trt
fit <- ssden(~futime*trt,data=gastric)
## conditional density
cdssden(fit,c("1","2"),cond=data.frame(futime=150))
## conditional quantiles
cqssden(fit,c(.05,.25,.5,.75,.95),data.frame(trt="1"))
## Clean up
\dontrun{rm(gastric,fit)}

## Sampling bias
## (X,T) is truncated to T<X<1 for T~U(0,1), so X is length-biased
rbias <- function(n) {
  t <- runif(n)
  x <- rnorm(n,.5,.15)
  ok <- (x>t)&(x<1)
  while(m<-sum(!ok)) {
    t[!ok] <- runif(m)
    x[!ok] <- rnorm(m,.5,.15)
    ok <- (x>t)&(x<1)
  }
  cbind(x,t)
}
xt <- rbias(100)
x <- xt[,1]; t <- xt[,2]
## length-biased
bias1 <- list(t=1,wt=1,fun=function(t,x){x[,]})
fit1 <- ssden(~x,domain=list(x=c(0,1)),bias=bias1)
plot(xx<-seq(0,1,len=101),dssden(fit1,xx),type="l")
## truncated
bias2 <- list(t=t,wt=rep(1/100,100),fun=function(t,x){x[,]>t})
fit2 <- ssden(~x,domain=list(x=c(0,1)),bias=bias2)
plot(xx,dssden(fit2,xx),type="l")
## Clean up
\dontrun{rm(rbias,xt,x,t,bias1,fit1,bias2,fit2)}
}
\keyword{smooth}
\keyword{models}
\keyword{distribution}

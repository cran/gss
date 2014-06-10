\name{ssanova9}
\alias{ssanova9}
\alias{para.arma}
\title{Fitting Smoothing Spline ANOVA Models with Correlated Data}
\description{
    Fit smoothing spline ANOVA models with correlated Gaussian data.
    The symbolic model specification via \code{formula} follows the same
    rules as in \code{\link{lm}}.
}
\usage{
ssanova9(formula, type=NULL, data=list(), subset, offset,
         na.action=na.omit, partial=NULL, method="v", alpha=1.4,
         varht=1, id.basis=NULL, nbasis=NULL, seed=NULL, cov,
         skip.iter=FALSE)

para.arma(fit)
}
\arguments{
    \item{formula}{Symbolic description of the model to be fit.}
    \item{type}{List specifying the type of spline for each variable.
        See \code{\link{mkterm}} for details.}
    \item{data}{Optional data frame containing the variables in the
	model.}
    \item{subset}{Optional vector specifying a subset of observations
	to be used in the fitting process.}
    \item{offset}{Optional offset term with known parameter 1.}
    \item{na.action}{Function which indicates what should happen when
	the data contain NAs.}
    \item{partial}{Optional symbolic description of parametric terms in
        partial spline models.}
    \item{method}{Method for smoothing parameter selection.  Supported
	are \code{method="v"} for V, \code{method="m"} for M, and
	\code{method="u"} for U; see the reference for definitions of U,
	V, and M.}
    \item{alpha}{Parameter modifying V or U; larger absolute values
        yield smoother fits.  Ignored when \code{method="m"} are
        specified.}
    \item{varht}{External variance estimate needed for
	\code{method="u"}.  Ignored when \code{method="v"} or
	\code{method="m"} are specified.}
    \item{id.basis}{Index designating selected "knots".}
    \item{nbasis}{Number of "knots" to be selected.  Ignored when
	\code{id.basis} is supplied.}
    \item{seed}{Seed to be used for the random generation of "knots".
	Ignored when \code{id.basis} is supplied.}
    \item{cov}{Input for covariance functions.  See \code{\link{mkcov}}
        for details.}
    \item{skip.iter}{Flag indicating whether to use initial values of
        theta and skip theta iteration.  See notes on skipping theta
	iteration.}

    \item{fit}{\code{ssanova9} fit with ARMA error.}
}
\details{
    The model specification via \code{formula} is intuitive.  For
    example, \code{y~x1*x2} yields a model of the form
    \deqn{
	y = C + f_{1}(x1) + f_{2}(x2) + f_{12}(x1,x2) + e
    }
    with the terms denoted by \code{"1"}, \code{"x1"}, \code{"x2"}, and
    \code{"x1:x2"}.

    The model terms are sums of unpenalized and penalized
    terms. Attached to every penalized term there is a smoothing
    parameter, and the model complexity is largely determined by the
    number of smoothing parameters.

    A subset of the observations are selected as "knots."  Unless
    specified via \code{id.basis} or \code{nbasis}, the number of
    "knots" \eqn{q} is determined by \eqn{max(30,10n^{2/9})}, which is
    appropriate for the default cubic splines for numerical vectors.

    Using \eqn{q} "knots," \code{ssanova} calculates an approximate
    solution to the penalized least squares problem using algorithms of
    the order \eqn{O(nq^{2})}, which for \eqn{q<<n} scale better than
    the \eqn{O(n^{3})} algorithms of \code{\link{ssanova0}}.  For the
    exact solution, one may set \eqn{q=n} in \code{ssanova}, but
    \code{\link{ssanova0}} would be much faster.
}
\section{Skipping Theta Iteration}{
    For the selection of multiple smoothing parameters,
    \code{\link{nlm}} is used to minimize the selection criterion such
    as the GCV score.  When the number of smoothing parameters is large,
    the process can be time-consuming due to the great amount of
    function evaluations involved.

    The starting values for the \code{nlm} iteration are obtained using
    Algorith 3.2 in Gu and Wahba (1991).  These starting values usually
    yield good estimates themselves, leaving the subsequent quasi-Newton
    iteration to pick up the "last 10\%" performance with extra effort
    many times of the initial one.  Thus, it is often a good idea to
    skip the iteration by specifying \code{skip.iter=TRUE}, especially
    in high-dimensions and/or with multi-way interactions.

    \code{skip.iter=TRUE} could be made the default in future releases.
}
\note{
    The results may vary from run to run. For consistency, specify
    \code{id.basis} or set \code{seed}.
}
\value{
    \code{ssanova9} returns a list object of class
    \code{c("ssanova9","ssanova")}.

    The method \code{\link{summary.ssanova9}} can be used to obtain
    summaries of the fits.  The method \code{\link{predict.ssanova}} can
    be used to evaluate the fits at arbitrary points along with standard
    errors.  The method \code{\link{project.ssanova9}} can be used to
    calculate the Kullback-Leibler projection for model selection.  The
    methods \code{\link{residuals.ssanova}} and
    \code{\link{fitted.ssanova}} extract the respective traits from the
    fits.

    \code{para.arma} returns the fitted ARMA coefficients for
    \code{cov=list("arma",c(p,q))} in the call to \code{ssanova9}.
}
\author{Chong Gu, \email{chong@stat.purdue.edu}}
\references{
    Han, C. and Gu, C. (2008), Optimal smoothing with correlated data,
    \emph{Sankhya}, \bold{70-A}, 38--72.

    Gu, C. (2013), \emph{Smoothing Spline ANOVA Models (2nd Ed)}.  New
    York: Springer-Verlag.

    Chong Gu (2014), Smoothing Spline ANOVA Models: R Package gss.
    \emph{Journal of Statistical Software}, 58(5), 1-25. URL
    http://www.jstatsoft.org/v58/i05/.
}
\examples{
x <- runif(100); y <- 5 + 3*sin(2*pi*x) + rnorm(x)
## independent fit
fit <- ssanova9(y~x,cov=list("known",diag(1,100)))
## AR(1) fit
fit <- ssanova9(y~x,cov=list("arma",c(1,0)))
para.arma(fit)
## MA(1) fit
e <- rnorm(101); e <- e[-1]-.5*e[-101]
x <- runif(100); y <- 5 + 3*sin(2*pi*x) + e
fit <- ssanova9(y~x,cov=list("arma",c(0,1)))
para.arma(fit)
## Clean up
\dontrun{rm(x,y,e,fit)}
}
\keyword{smooth}
\keyword{models}
\keyword{regression}

\name{predict.ssanova}
\title{Predicting from Smoothing Spline ANOVA Fits}
\usage{
predict[.ssanova](obj, newdata, se.fit=FALSE, include=obj$terms$labels)
}
\arguments{
    \item{obj}{an object of class \code{"ssanova"}.}
    \item{newdata}{a data frame or model frame.}
    \item{se.fit}{a logical flag.}
    \item{include}{a list of model terms to be included in the
	prediction.  The \code{partial} and \code{offset} terms, if
	present, are specified by \code{"partial"} and \code{"offset"},
	respectively.}
}
\description{
    \code{predict.ssanova} evaluates a smoothing spline ANOVA fit at
    arbitrary data points as specified in \code{newdata}.  The model
    terms to be included in the evaluation are specified via
    \code{include}, which defaults to all terms in the model.

    With the flag \code{se.fit=TRUE}, the standard errors of the fit are
    also calculated, which are to be used to construct Bayesian
    confidence intervals of the fit.
}
\value{
    For \code{se.fit=FALSE}, \code{predict.ssanova} returns a vector of
    the evaluated fit.

    For \code{se.fit=TRUE}, \code{predict.ssanova} returns a list
    consisting of the following components.
    \item{fit}{a vector of the evaluated fit.}
    \item{se.fit}{a vector of the corresponding standard errors.}
}
\seealso{
    The model fitting function \code{\link{ssanova}} and the summarizing
    function \code{\link{summary.ssanova}}.
}
\examples{
## Fit a model with thin-plate marginals, where geog is 2-D
data(lake.acid)
fit <- ssanova(ph~log(cal)*geog,"tp",lake.acid)

## Obtain estimates and standard errors on a grid
new <- data.frame(cal=1,geog=I(matrix(0,1,2)))
new <- model.frame(~log(cal)+geog,new)
predict(fit,new,se=TRUE)

## Evaluate the geog main effect
predict(fit,new,se=TRUE,inc="geog")

## Evaluate the sum of the geog main effect and the interaction
predict(fit,new,se=TRUE,inc=c("geog","log(cal):geog"))

## Evaluate the geog main effect on a grid
grid <- seq(-.04,.04,len=21)
new <- model.frame(~geog,list(geog=cbind(rep(grid,21),rep(grid,rep(21,21)))))
est <- predict(fit,new,se=TRUE,inc="geog")

## Plot the fit and standard error
par(pty="s")
contour(grid,grid,matrix(est$fit,21,21),col=1)
contour(grid,grid,matrix(est$se,21,21),add=TRUE,col=2)
}
\author{Chong Gu, \email{chong@stat.purdue.edu}}
\keyword{ssanova}
\keyword{models}
\keyword{regression}
\keyword{smoothing}
\keyword{spline}

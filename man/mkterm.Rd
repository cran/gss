\name{mkterm}
\alias{mkterm}
\title{
    Assembling Model Terms for Smoothing Spline ANOVA Models
}
\description{
    Assemble numerical functions for calculating model terms in a
    smoothing spline ANOVA model.
}
\usage{
mkterm(mf, type)
}
\arguments{
    \item{mf}{Model frame of the model formula.}
    \item{type}{List specifying the type of spline for each variable.}
}
\section{Background}{
    Tensor-product splines are constructed based on the model formula
    and the marginal reproducing kernels, as described in Gu (2002,
    Sec. 2.4).  The marginal variables can be factors, numerical
    vectors, and numerical matrices, as specified in the details
    section.

    One-way ANOVA decompositions are built in the supported marginal
    constructions, in which one has the constant, a "nonparametric
    contrast," and possibly also a "parametric contrast."  To the
    "nonparametric contrast" there corresponds a reproducing kernal
    \code{rk}, and to a "parametric contrast" there corresponds a set
    of null space basis \code{phi}.  The reproducing kernels and null
    space basis on the product domain can be constructed from the
    marginal \code{rk} and \code{phi} in a systematic manner.

    The marginal one-way ANOVA structures induce a multi-way ANOVA
    structure on the product domain, with model terms consisting of
    unpenalized "parametric contrasts" and/or penalized "nonparametric
    contrasts."  One only needs to construct \code{rk}'s and
    \code{phi}'s associated with the model terms implied by the model
    formula.
}
\details{
    For a \bold{factor} \code{x}, \code{type$x} is ignored;
    \code{\link{mkrk.ordinal}} is used if \code{is.ordered(x)==TRUE} and
    \code{\link{mkrk.nominal}} is used otherwise.  Factors with 3 or
    more levels are penalized.

    For a \bold{numerical vector} \code{x}, \code{type$x} is of the form
    \code{type.x} for \code{type.x}=\code{"cubic"}, \code{"linear"},
    or of the form \code{list(type.x, range)} for
    \code{type.x}=\code{"per"}, \code{"cubic.per"}, \code{"linear.per"},
    \code{"cubic"}, \code{"linear"}; \code{"per"} is short for
    \code{"cubic.per"}.  See \code{\link{mkfun.poly}} for the functions
    used.  For \code{type.x} missing, the \bold{default} is \code{"cubic"}.
    For \code{range} missing with \code{type.x}=\code{"cubic"},
    \code{"linear"}, the \bold{default} is
    \code{c(min(x),max(x))+c(-1,1)*(max(x)-mimn(x))*.05}.

    For a \bold{numerical matrix} \code{x}, \code{type$x} is of the form
    \code{type.x} or \code{list(type.x, order)} for
    \code{type.x}=\code{"tp"}, \code{"sphere"}, or of the form
    \code{list("tp",list(order=order,mesh=mesh,weight=weight))}.  See
    \code{\link{mkfun.tp}} for the functions used.  For \code{type.x}
    missing, the \bold{default} is \code{"tp"}.  For \code{order}
    missing, the \bold{default} is \code{2}.  For \code{mesh} and
    \code{weight} missing with \code{type.x}=\code{"tp"} and
    \code{order} given, the \bold{defaults} are \code{mesh}=\code{x} and
    \code{weight}=\code{1}.

    For a \bold{numerical vector} or \bold{numerical matrix} \code{x},
    one may also use \code{type$x} of the form
    \code{list("custom",list(nphi=nphi,mkphi=mkphi,mkrk=mkrk,env=env))};
    \code{nphi} is the null space dimension \emph{excluding the
      constant}, and \code{mkphi} is ignored if \code{nphi}=0.  See
    examples below.  This feature allows the use of other marginal
    constructions; one may modify \code{\link{mkphi.cubic}} or
    \code{\link{mkphi.tp.p}} for \code{mkphi} and modify
    \code{\link{mkrk.cubic}} or \code{\link{mkrk.sphere}} for
    \code{mkrk}.
}
\note{
    For a \bold{numerical vector} \code{x} in \code{\link{ssden}},
    the default \code{range} is \code{domain$x}.

    For a \bold{numerical matrix} \code{x} with
    \code{type.x}=\code{"sphere"}, it is assumed that
    \code{dim(x)[2]==2}, \code{x[,1]} between [-90,90] the latitude in
    degrees, and \code{x[,2]} between [-180,180] the longitude in
    degrees.

    For \bold{backward compatibility}, one may set \code{type="cubic"},
    \code{"linear"}, or \code{"tp"}, but then the default parameters can
    not be overridden; the type is simply duplicated for each variable.
}
\value{
    A list object with an element \code{labels} containing the labels
    of all model terms.  For each of the model terms, there is an
    element holding the numerical functions for calculating the
    unpenalized and penalized parts within the term.
}
\references{
    Gu, C. (2013), \emph{Smoothing Spline ANOVA Models (2nd Ed)}.  New
    York: Springer-Verlag.
}
\examples{
## cubic marginals
x1 <- rnorm(100); x2 <- rnorm(100); y <- 3+5*sin(x1-2*x2)+rnorm(x1)
fit <- ssanova0(y~x1*x2)
## the same fit
fit1 <- ssanova0(y~x1*x2,type=list(x1="cubic"))
## the same fit one more time
par <- list(nphi=1,mkphi=mkphi.cubic,mkrk=mkrk.cubic,
            env=c(min(x2),max(x2))+c(-1,1)*(max(x2)-min(x2))*.05)
fit2 <- ssanova0(y~x1*x2,type=list(x2=list("custom",par)))
## Clean up
\dontrun{rm(x1,x2,y,fit,fit1,par,fit2)}

## cubic and thin-plate marginals
x1 <- rnorm(100); x2 <- matrix(rnorm(200),100,2)
y <- 3+5*sin(x1-2*x2[,1]*x2[,2])+rnorm(x1)
fit <- ssanova0(y~x1*x2)
## the same fit
fit1 <- ssanova0(y~x1*x2,type=list(x2="tp"))
## the same fit one more time
mkphi.tp1 <- function(x) mkphi.tp(x$dm,x$ord,x$mesh,x$wt)
mkrk.tp1 <- function(x) mkrk.tp(x$dm,x$ord,x$mesh,x$wt)
env <- list(dm=2,ord=2,mesh=x2,wt=1)
par <- list(nphi=2,mkphi=mkphi.tp1,mkrk=mkrk.tp1,env=env)
fit2 <- ssanova0(y~x1*x2,type=list(x2=list("custom",par)))
## Clean up
\dontrun{rm(x1,x2,y,fit,fit1,mkphi.tp1,mkrk.tp1,env,par,fit2)}
}
\keyword{internal}

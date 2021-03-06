\name{mkran}
\alias{mkran}
\alias{mkran1}
\title{
    Generating Random Effects in Mixed-Effect Models
}
\description{
    Generate entries representing random effects in mixed-effect models.
}
\usage{
mkran(formula, data)
mkran1(ran1, ran2)
}
\arguments{
    \item{formula}{Symbolic description of the random effects.}
    \item{data}{Data frame containing the variables in the model.}
    \item{ran1}{Random effects in the form of the value of \code{mkran}}.
    \item{ran2}{Random effects in the form of the value of \code{mkran}}.
}
\details{
    \code{mkran} generates random effect terms from simple grouping
    variables, for use in nonparametric mixed-effect models as described
    in Gu and Ma (2005a, b).  The syntax of the formula resembles that
    of similar utilities for linear and nonlinear mixed-effect models,
    as described in Pinheiro and Bates (2000).

    Currently, \code{mkran} takes only two kinds of basic formulas,
    \code{~1|grp2} or \code{~grp1|grp2}.  Both \code{grp1} and
    \code{grp2} should be factors, and for the second formula, the
    levels of \code{grp2} should be nested under those of \code{grp1}.

    The Z matrix is determined by \code{grp2}.  When observations are
    ordered according to the levels of \code{grp2}, the Z matrix is
    block diagonal of 1 vectors.

    The Sigma matrix is diagonal.  For \code{~1|grp2}, it has one tuning
    parameter.  For \code{~grp1|grp2}, the number of parameters equals
    the number of levels of \code{grp1}, with each parameter shared by
    the \code{grp2} levels nested under the same \code{grp1} level.

    \code{mkran1} adds together two independent random effects, and can
    be used recursively to add more than two terms.  The arguments are
    of the form of the value of \code{mkran} or \code{mkran1}, which may
    or may not be created by \code{mkran} or \code{mkran1}.

    Multiple terms of random effects can also be specified via the likes
    of \code{mkran(~1|grp1+1|grp2,data)}, which is equivalent to
    \code{mkran1(mkran(~1|grp1,data),mkran(~1|grp2,data))}.
}
\value{
    A list of three elements.
    \item{z}{Z matrix.}
    \item{sigma}{Sigma matrix to be evaluated through
        \code{sigma$fun(para,sigma$env)}.}
    \item{init}{Initial parameter values.}
}
\note{
    One may pass a formula or a list to the argument \code{random} in
    calls to \code{\link{ssanova}} or\code{\link{gssanova}} to fit
    nonparametric mixed-effect models.  A formula will be converted to a
    list using \code{mkran}.  A list should be of the same form as the
    value of \code{mkran}.
}
\references{
    Gu, C. and Ma, P. (2005), Optimal smoothing in nonparametric
    mixed-effect models.  \emph{The Annals of Statistics}, \bold{33},
    1357--1379.

    Gu, C. and Ma, P. (2005), Generalized nonparametric mixed-effect
    models: computation and smoothing parameter selection.
    \emph{Journal of Computational and Graphical Statistics}, \bold{14},
    485--504.

    Pinheiro and Bates (2000), \emph{Mixed-Effects Models in S and
    S-PLUS}.  New York: Springer-Verlag.
}
\examples{
## Toy data
test <- data.frame(grp=as.factor(rep(1:2,c(2,3))))
## First formula
ran.test <- mkran(~1|grp,test)
ran.test$z
ran.test$sigma$fun(2,ran.test$sigma$env) # diag(10^(-2),2)
## Second formula
ran.test <- mkran(~grp|grp,test)
ran.test$z
ran.test$sigma$fun(c(1,2),ran.test$sigma$env) # diag(10^(-1),10^(-2))
## Clean up
\dontrun{rm(test,ran.test)}
}
\keyword{internal}

\name{bacteriuria}
\alias{bacteriuria}
\title{Treatment of Bacteriuria}
\description{
    Bacteriuria patients were randomly assigned to two treatment
    groups.  Weekly binary indicator of bacteriuria was recorded for
    every patient over 4 to 16 weeks.  A total of 72 patients were
    represented in the data, with 36 each in the two treatment groups.
}
\usage{data(bacteriuria)}
\format{
    A data frame containing 820 observations on the following variables.
    \tabular{ll}{
	\code{id} \tab Identification of patients, a factor.\cr
	\code{trt} \tab Treatments 1 or 2, a factor.\cr
	\code{time} \tab Weeks after randomization.\cr
	\code{infect} \tab Binary indicator of bacteriuria (bacteria in
	urine).
    }
}
\source{
    Joe, H. (1997), \emph{Multivariate Models and Dependence Concepts}.
    London: Chapman and Hall.
}
\references{
    Gu and Ma (2003b), \emph{Generalized Nonparametric Mixed-Effect
    Models: Computation and Smoothing Parameter Selection}.  Available
    at \url{http://stat.purdue.edu/~chong/manu.html}.
}
\keyword{datasets}

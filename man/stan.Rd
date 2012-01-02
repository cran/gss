\name{stan}
\alias{stan}
\title{Stanford Heart Transplant Data}
\description{
    Survival of patients from the Stanford heart transplant program.
}
\usage{data(stan)}
\format{
    A data frame containing 184 observations on the following variables.
    \tabular{ll}{
	\code{time} \tab Follow-up time after transplant, in days.\cr
	\code{status} \tab Censoring status.\cr
	\code{age} \tab Age at transplant.\cr
	\code{futime} \tab Square root of \code{time}.
    }
}
\source{
  Miller, R. G. and Halpern, J. (1982), Regression with censored data.
  \emph{Biometrika}, \bold{69}, 521--531.
}
\keyword{datasets}

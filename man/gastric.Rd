\name{gastric}
\alias{gastric}
\title{Gastric Cancer Data}
\description{
    Survival of gastric cancer patients under chemotherapy and
    chemotherapy-radiotherapy combination.
}
\usage{data(gastric)}
\format{
    A data frame containing 90 observations on the following variables.
    \tabular{ll}{
	\code{futime} \tab Follow-up time, in days.\cr
	\code{status} \tab Censoring status.\cr
	\code{trt} \tab Factor indicating the treatments: 1 --
	chemothrapy, 2 -- combination.
    }
}
\source{
  Moreau, T., O'Quigley, J., and Mesbah, M. (1985), A global
  goodness-of-fit statistic for the proportional hazards model.
  \emph{Applied Statistics}, \bold{34}, 212-218.
}
\keyword{datasets}

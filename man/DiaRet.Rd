\name{DiaRet}
\alias{DiaRet}
\title{Diabetic Retinopathy}
\description{
    Time to blindness of 197 diabetic retinopathy patients who received
    a laser treatment in one eye.
}
\usage{data(DiaRet)}
\format{
    A data frame containing 197 observations on the following variables.
    \tabular{ll}{
        \code{id} \tab Patient ID.\cr
	\code{time1} \tab Follow-up time of left eye.\cr
	\code{time2} \tab Follow-up time of right eye.\cr
	\code{status1} \tab Censoring indicator of left eye.\cr
	\code{status2} \tab Censoring indicator of right eye.\cr
	\code{trt1} \tab Treatment indicator of left eye.\cr
	\code{trt2} \tab Treatment indicator of right eye.\cr
	\code{type} \tab Type of diabetes.\cr
	\code{age} \tab Age of patient at diagnosis.\cr
	\code{time.t} \tab Follow-up time of treated eye.\cr
	\code{time.u} \tab Follow-up time of untreated eye.\cr
	\code{status.t} \tab Censoring indicator of treated eye.\cr
	\code{status.u} \tab Censoring indicator of untreated eye.
    }
}
\source{
    This is reformatted from the data frame \code{diabetes} in the R
    package \code{timereg} by Thomas H. Scheike.
}
\references{
    Huster, W.J., Brookmeyer, R., and Self, S.G. (1989), Modelling
    paired survival data with covariates.  \emph{Biometrics}, \bold{45},
    145--56.
}
\keyword{datasets}

\name{NO2}
\alias{NO2}
\title{Air Pollution and Road Traffic}
\description{
    A subset of 500 hourly observations collected by the Norwegian
    Public Roads Administration at Alnabru in Oslo, Norway, between
    October 2001 and August 2003.
}
\usage{data(NO2)}
\format{
    A data frame containing 500 observations on the following variables.
    \tabular{ll}{
	\code{no2} \tab Concentration of NO2, on log scale.\cr
	\code{cars} \tab Traffic volume of the hour, on log scale.\cr
	\code{temp} \tab Temperature 2 meters above ground, in Celsius.\cr
	\code{wind} \tab wind speed, meters/second.\cr
	\code{temp2} \tab Temperature difference between 25 and 2 meters
	above ground, in Celsius.\cr
	\code{wind2} \tab Wind direction, in degrees between 0 and 360.
    }
}
\source{
    Statlib Datasets Archive at \code{http://lib.stat.cmu.edu/datasets},
    contributed by Magne Aldrin.
}
\keyword{datasets}

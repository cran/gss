\name{ozone}
\alias{ozone}
\title{Ozone Concentration in Los Angeles Basin}
\description{
    Daily measurements of ozone concentration and eight meteorological
    quantities in the Los Angeles basin for 330 days of 1976.
}
\usage{data(ozone)}
\format{
    A data frame containing 330 observations on the following variables.
    \tabular{ll}{
	\code{upo3} \tab Upland ozone concentration, in ppm.\cr
	\code{vdht} \tab Vandenberg 500 millibar height, in meters.\cr
	\code{wdsp} \tab Wind speed, in miles per hour.\cr
	\code{hmdt} \tab Humidity.\cr
	\code{sbtp} \tab Sandburg Air Base temperature, in Celsius.\cr
	\code{ibht} \tab Inversion base height, in foot.\cr
	\code{dgpg} \tab Dagget pressure gradient, in mmhg.\cr
	\code{ibtp} \tab Inversion base temperature, in Fahrenheit.\cr
	\code{vsty} \tab Visibility, in miles.\cr
	\code{day} \tab Calendar day, between 1 and 366.
    }
}
\source{
    Unknown.
}
\references{
    Breiman, L. and Friedman, J. H. (1985).  Estimating optimal
    transformations for multiple regression and correlation.
    \emph{Journal of the American Statistical Association}, \bold{80},
    580--598.

    Hastie, T. and Tibshirani, R. (1990).  \emph{Generalized Additive
    Models}. Chapman and Hall.
}
\keyword{datasets}

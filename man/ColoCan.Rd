\name{ColoCan}
\alias{ColoCan}
\title{Colorectal Cancer Mortality Rate in Indiana Counties}
\description{
    County-wise death counts of colorectal cancer patients in Indiana
    during years 2000 through 2004.
}
\usage{data(ColoCan)}
\format{
    A data frame containing 184 observations on the following variables.
    \tabular{ll}{
	\code{event} \tab Death counts.\cr
	\code{pop} \tab Population from Census 2000.\cr
	\code{sex} \tab Gender of population.\cr
	\code{wrt} \tab Proportion of Whites.\cr
	\code{brt} \tab Proportion of Blacks.\cr
	\code{ort} \tab Proportion of other minorities.\cr
	\code{lat} \tab Latitude.\cr
	\code{lon} \tab Longitude.\cr
	\code{geog} \tab Geographic location, derived from \code{lat}
	and \code{lon}.\cr
	\code{scrn} \tab Colorectal cancer screening rate.\cr
	\code{name} \tab County name.
    }
}
\details{
    \code{geog} was generated from \code{lat} and \code{lon} using the
    code given in the example section.
}
\source{
    Dr. Tonglin Zhang.
}
\references{
    Zhang, T. and Lin, G. (2009), Cluster detection based on spatial
    associations and iterated residuals in generalized linear mixed
    models.  \emph{Biometrics}, \bold{65}, 353--360.
}
\examples{
## Converting latitude and longitude to x-y coordinates
## The 49th county is Marion, where Indianapolis is located.
\dontrun{ltln2xy <- function(latlon,latlon0) {
  lat <- latlon[,1]*pi/180; lon <- latlon[,2]*pi/180
  lt0 <- latlon0[1]*pi/180; ln0 <- latlon0[2]*pi/180
  x <- cos(lt0)*sin(lon-ln0); y <- sin(lat-lt0)
  cbind(x,y)
}
data(ColoCan)
latlon <- as.matrix(ColoCan[,c("lat","lon")])
ltln2xy(latlon,latlon[49,])
## Clean up
rm(ltln2xy,ColoCan,latlon)}
}
\keyword{datasets}

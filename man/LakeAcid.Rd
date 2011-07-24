\name{LakeAcidity}
\alias{LakeAcidity}
\title{Water Acidity in Lakes}
\description{
    Data extracted from the Eastern Lake Survey of 1984 conducted by the
    United States Environmental Protection Agency, concerning 112 lakes
    in the Blue Ridge.
}
\usage{data(LakeAcidity)}
\format{
    A data frame containing 112 observations on the following variables.
    \tabular{ll}{
	\code{ph} \tab Surface ph.\cr
	\code{cal} \tab Calcium concentration.\cr
	\code{lat} \tab Latitude.\cr
	\code{lon} \tab Longitude.\cr
	\code{geog} \tab Geographic location, derived from \code{lat}
	and \code{lon}
    }
}
\details{
    \code{geog} was generated from \code{lat} and \code{lon} using the
    code given in the Example section.
}
\source{
    Douglas, A. and Delampady, M. (1990), \emph{Eastern Lake Survey --
	Phase I: Documentation for the Data Base and the Derived Data
	sets.} Tech Report 160 (SIMS), Dept. Statistics, University of
        British Columbia.
}
\references{
    Gu, C. and Wahba, G. (1993), Semiparametric analysis of variance
    with tensor product thin plate splines.  \emph{Journal of the Royal
	Statistical Society Ser. B}, \bold{55}, 353--368.
}
\examples{
## Converting latitude and longitude to x-y coordinates
\dontrun{ltln2xy <- function(latlon,latlon0) {
  lat <- latlon[,1]*pi/180; lon <- latlon[,2]*pi/180
  lt0 <- latlon0[1]*pi/180; ln0 <- latlon0[2]*pi/180
  x <- cos(lt0)*sin(lon-ln0); y <- sin(lat-lt0)
  cbind(x,y)
}
data(LakeAcidity)
latlon <- as.matrix(LakeAcidity[,c("lat","lon")])
m.lat <- (min(latlon[,1])+max(latlon[,1]))/2
m.lon <- (min(latlon[,2])+max(latlon[,2]))/2
ltln2xy(latlon,c(m.lat,m.lon))
## Clean up
rm(ltln2xy,LakeAcidity,latlon,m.lat,m.lon)}
}
\keyword{datasets}

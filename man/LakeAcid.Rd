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
    A list containing 112 observations on the following variables.
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
\dontrun{convert <- function(lat, lon) {
    lat <- lat/180*pi
    lon <- lon/180*pi
    m.lat <- (max(lat)+min(lat))/2
    m.lon <- (max(lon)+min(lon))/2
    x <- cos(m.lat)*sin(m.lon-lon)
    y <- sin(lat-m.lat)
    cbind(x,y)
}
convert(LakeAcidity$lat,LakeAcidity$lon)}
}
\keyword{datasets}

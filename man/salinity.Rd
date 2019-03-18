\name{salinity}
\alias{salinity}
\docType{data}
\title{Salinity tolerance}
\description{
    72-hour acute salinity tolerance (LC50 values) of riverine macro-invertebrates. 
}
\usage{data(acidity)}
\format{
  A data frame with 108 observations on the following two variables:
  \describe{
    \item{\code{left}}{A numeric vector.}
    \item{\code{right}}{A numeric vector.}
  }
}
\source{
\code{fitdistrplus} R-package
}
\references{
Kefford, B.J., Nugegoda, D., Metzeling, L., Fields, E. 2006. Validating species sensitivity distributions using salinity tolerance of riverine macroinvertebrates in the southern Murray-darling Basin (Vitoria, Australia). Canadian Journal of Fisheries and Aquatic Science, 63, 1865-1877.}
\examples{
data(salinity)
hist(salinity$left)
}
\keyword{datasets}

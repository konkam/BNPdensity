\name{galaxy}
\alias{galaxy}
\docType{data}
\title{Galaxy Data Set}
\description{
Velocities of 82 galaxies diverging from our own galaxy.
}
\usage{data(galaxy)}
\format{
  A data frame with 82 observations on the following variable:
  \describe{
    \item{\code{velocity}}{A numeric vector.}
    }
}
\source{
\url{https://people.maths.bris.ac.uk/~mapjg/}
}
\references{
Roeder, K. (1990) "Density estimation with confidence sets exemplified by 
superclusters and voids in the galaxies". Journal of the American Statitical 
Association. 85, 617-624.
}
\examples{
data(galaxy)
hist(galaxy)
}
\keyword{datasets}

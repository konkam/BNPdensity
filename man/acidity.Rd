\name{acidity}
\alias{acidity}
\docType{data}
\title{Acidity Index Dataset}
\description{
    Concerns an acidity index measured in a sample of 155 lakes in north-central 
    Wisconsin. 
}
\usage{data(acidity)}
\format{
  A data frame with 155 observations on the following variable:
  \describe{
    \item{\code{acidity}}{A numeric vector.}
  }
}
\source{
\url{https://people.maths.bris.ac.uk/~mapjg/}
}
\references{
Crawford, S. L., DeGroot, M. H., Kadane, J. B. and Small, M. J. (1992). 
Modeling lake chemistry distributions: approximate Bayesian methods for estimating 
a finite mixture model. Technometrics, 34, 441-453.}
\examples{
data(acidity)
hist(acidity)
}
\keyword{datasets}

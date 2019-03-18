\name{enzyme}
\alias{enzyme}
\docType{data}
\title{Enzyme Dataset}
\description{
Concerns the distribution of enzymatic activity in the blood, 
for an enzyme involved in the metabolism of carcinogenetic substances, 
among a group of 245 unrelated individuals. 
}
\usage{data(enzyme)}
\format{
  A data frame with 244 observations on the following variable:
  \describe{
    \item{\code{enzyme}}{A numeric vector.}
  }
}
\source{
\url{https://people.maths.bris.ac.uk/~mapjg/}
}
\references{
Bechtel, Y. C., Bonaiti-Pellie, C., Poisson, N., Magnette, J. and Bechtel, P.R. (1993). 
A population and family study of N-acetyltransferase using caffeine urinary 
metabilites. Clin. Pharm. Therp., 54, 134-141.
}
\examples{
data(enzyme)
hist(enzyme)
}
\keyword{datasets}

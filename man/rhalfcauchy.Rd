\name{rhalfcauchy}
\alias{rhalfcauchy}
\title{
Random number generator half Cauchy}
\description{
Computes a random number.}
\usage{
rhalfcauchy(n, location = 0, scale = 1)
}
\details{
For internal use
}
\examples{
## The function is currently defined as
function (n, location = 0, scale = 1) 
{
    abs(rcauchy(n, location, scale))
  }
}
\keyword{internal}

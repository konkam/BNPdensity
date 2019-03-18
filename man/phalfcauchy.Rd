\name{phalfcauchy}
\alias{phalfcauchy}
\title{
Distribution function half Cauchy
}
\description{
Computes the cdf.
}
\usage{
phalfcauchy(q, location = 0, scale = 1)
}
\details{
For internal use}
\examples{
## The function is currently defined as
function (q, location = 0, scale = 1) 
{
    ifelse(x < 0, 0, 1) * (pcauchy(q, location, scale) - pcauchy(0, 
        location, scale))/(1 - pcauchy(0, location, scale))
  }
}
\keyword{internal}

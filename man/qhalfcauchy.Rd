\name{qhalfcauchy}
\alias{qhalfcauchy}
\title{
Quantile function half Cauchy
}
\description{
Computes the quantiles.
}
\usage{
qhalfcauchy(p, location = 0, scale = 1)
}
\details{
For internal use}
\examples{
## The function is currently defined as
function (p, location = 0, scale = 1) 
{
    qcauchy(p * (1 - pcauchy(0, location, scale)) + pcauchy(0, 
        location, scale), location, scale)
  }
}
\keyword{internal}

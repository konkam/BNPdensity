\name{rhalfnorm}
\alias{rhalfnorm}
\title{
Random number generator half Normal}
\description{
Computes a random number.
}
\usage{
rhalfnorm(n, mean = 0, sd = 1)
}
\details{
For internal use
}
\examples{
## The function is currently defined as
function (n, mean = 0, sd = 1) 
{
    abs(rnorm(n, mean, sd))
  }
}
\keyword{internal}

\name{rhalft}
\alias{rhalft}
\title{
Random number generator half Student-t}
\description{
Generates a random number.
}
\usage{
rhalft(n, df = 1, mean = 0, sd = 1)
}
\details{
For internal use
}
\examples{
## The function is currently defined as
function (n, df = 1, mean = 0, sd = 1) 
{
    abs(rt_(n, df, mean, sd))
  }
}
\keyword{internal}

\name{rt_}
\alias{rt_}
\title{
Random number generator non-standard Student-t}
\description{
Computes a random number.
}
\usage{
rt_(n, df, mean, sd)
}
\details{
For internal use
}
\examples{
## The function is currently defined as
function (n, df, mean, sd) 
{
    mean + sd * rt(n, df, ncp = 0)
  }
}
\keyword{internal}

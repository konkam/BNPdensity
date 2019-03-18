\name{dt_}
\alias{dt_}
\title{
Non-standard student-t densit
}
\description{
Computes the density.
}
\usage{
dt_(x, df, mean, sd)
}
\details{
For internal use
}
\examples{
## The function is currently defined as
function (x, df, mean, sd) 
{
    dt((x - mean)/sd, df, ncp = 0)/sd
  }
}
\keyword{internal}

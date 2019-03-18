\name{dhalfnorm}
\alias{dhalfnorm}
\title{
Density half normal
}
\description{
Computes the density.
}
\usage{
dhalfnorm(x, mean = 0, sd = 1)
}
\details{
For internal use
}
\examples{
## The function is currently defined as
function (x, mean = 0, sd = 1) 
{
    ifelse(x < 0, 0, 1) * dnorm(x, mean, sd)/(1 - pnorm(0, mean, 
        sd))
  }
}
\keyword{internal}

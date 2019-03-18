\name{qhalfnorm}
\alias{qhalfnorm}
\title{
Quantile function half Normal
}
\description{
Computes the quantiles.
}
\usage{
qhalfnorm(p, mean = 0, sd = 1)
}
\details{
For internal use
}
\examples{
## The function is currently defined as
function (p, mean = 0, sd = 1) 
{
    qnorm(p * (1 - pnorm(0, mean, sd)) + pnorm(0, mean, sd), 
        mean, sd)
  }
}
\keyword{internal}

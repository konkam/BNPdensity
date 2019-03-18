\name{qhalft}
\alias{qhalft}
\title{
Quantile function half Student-t
}
\description{
Computes the quantiles.
}
\usage{
qhalft(p, df = 1, mean = 0, sd = 1)
}
\details{
For internal use
}
\examples{
## The function is currently defined as
function (p, df = 1, mean = 0, sd = 1) 
{
    qt_(p * (1 - pt_(0, df, mean, sd)) + pt_(0, df, mean, sd), 
        df, mean, sd)
  }
}
\keyword{internal}

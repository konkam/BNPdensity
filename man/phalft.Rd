\name{phalft}
\alias{phalft}
\title{
Distribution function half Student-t
}
\description{
Computes the cumulative distribution function.
}
\usage{
phalft(q, df = 1, mean = 0, sd = 1)
}
\details{
For internal use
}
\examples{
## The function is currently defined as
function (q, df = 1, mean = 0, sd = 1) 
{
    ifelse(x < 0, 0, 1) * (pt_(q, df, mean, sd) - pt_(0, df, 
        mean, sd))/(1 - pt_(0, df, mean, sd))
  }
}
\keyword{internal}

\name{dhalft}
\alias{dhalft}
\title{
Density half Student-t
}
\description{
Computes the density. 
}
\usage{
dhalft(x, df = 1, mean = 0, sd = 1)
}
\details{
For internal use
}
\examples{
## The function is currently defined as
function (x, df = 1, mean = 0, sd = 1) 
{
    ifelse(x < 0, 0, 1) * dt_(x, df, mean, sd)/(1 - pt_(0, df, 
        mean, sd))
  }
}
\keyword{internal}
